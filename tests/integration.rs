mod common;

use common::{assert_geac_success, count_bam_reads, parquet_count, parquet_query_i32,
             parquet_query_str, write_bam, write_reference};
use tempfile::TempDir;

/// Sanity check: verify the BAM helper writes the expected number of reads.
#[test]
fn bam_helper_writes_expected_reads() {
    let dir = TempDir::new().unwrap();
    let _fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(), "s1.bam", "sample1", 200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    assert_eq!(count_bam_reads(&bam), 15, "expected 5 alt + 10 ref reads");
}

// ── geac collect ──────────────────────────────────────────────────────────────

/// Basic collect: 5 T alts + 10 A refs at chr1:50 → exactly 1 Parquet row
/// with the expected counts.
#[test]
fn collect_writes_valid_parquet() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(), "s1.bam", "sample1", 200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let out = dir.path().join("s1.parquet");

    assert_geac_success(&[
        "collect",
        "--input",    bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",   out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline", "raw",
    ]);

    assert!(out.exists(), "output Parquet not created");
    assert_eq!(parquet_count(&out), 1, "expected exactly 1 alt locus");
    assert_eq!(parquet_query_i32(&out, "alt_count",   "pos = 50"), 5);
    assert_eq!(parquet_query_i32(&out, "total_depth", "pos = 50"), 15);
    assert_eq!(parquet_query_str(&out, "alt_allele",  "pos = 50"), "T");
    assert_eq!(parquet_query_str(&out, "sample_id",   "pos = 50"), "sample1");
}

/// `--sample-id` flag overrides the SM tag read from the BAM header.
#[test]
fn collect_sample_id_override() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(), "s1.bam", "original_id", 200,
        vec![(50, b'T', 3, 10)],
        20,
    );
    let out = dir.path().join("s1.parquet");

    assert_geac_success(&[
        "collect",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
        "--sample-id", "overridden",
    ]);

    assert_eq!(parquet_query_str(&out, "sample_id", "pos = 50"), "overridden");
}

/// `--region` restricts output to only loci within the requested interval.
/// BAM has alts at pos 50 (0-based) and pos 150. With --region chr1:101-200
/// (1-based, covering 0-based 100–199) only pos 150 should appear.
#[test]
fn collect_region_restricts_output() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 300);
    let bam = write_bam(
        dir.path(), "s1.bam", "sample1", 300,
        vec![(50, b'T', 5, 10), (150, b'G', 3, 10)],
        20,
    );
    let out = dir.path().join("s1.parquet");

    assert_geac_success(&[
        "collect",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
        "--region",    "chr1:101-200",
    ]);

    assert_eq!(parquet_count(&out), 1, "region should restrict to 1 locus");
    assert_eq!(parquet_query_i32(&out, "alt_count", "pos = 150"), 3);
}

// ── geac merge ────────────────────────────────────────────────────────────────

/// merge combines two per-sample Parquets into a DuckDB with both samples present.
#[test]
fn merge_produces_duckdb_with_both_samples() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    for (name, sample) in [("s1.bam", "sampleA"), ("s2.bam", "sampleB")] {
        let bam = write_bam(
            dir.path(), name, sample, 200,
            vec![(50, b'T', 3, 10)],
            20,
        );
        let pq = dir.path().join(name.replace(".bam", ".parquet"));
        assert_geac_success(&[
            "collect",
            "--input",     bam.to_str().unwrap(),
            "--reference", fa.to_str().unwrap(),
            "--output",    pq.to_str().unwrap(),
            "--read-type", "raw",
            "--pipeline",  "raw",
        ]);
    }

    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge",
        "--output", db.to_str().unwrap(),
        dir.path().join("s1.parquet").to_str().unwrap(),
        dir.path().join("s2.parquet").to_str().unwrap(),
    ]);

    assert!(db.exists(), "DuckDB file not created");

    let conn = duckdb::Connection::open(&db).unwrap();
    let n_samples: i64 = conn
        .query_row("SELECT COUNT(DISTINCT sample_id) FROM alt_bases", [], |r| r.get(0))
        .unwrap();
    assert_eq!(n_samples, 2, "expected 2 distinct samples in merged DuckDB");
}

// ── geac qc ───────────────────────────────────────────────────────────────────

/// qc exits successfully and produces a TSV when --output is given.
#[test]
fn qc_exits_successfully() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(), "s1.bam", "sample1", 200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let pq = dir.path().join("s1.parquet");
    assert_geac_success(&[
        "collect",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    pq.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
    ]);

    let tsv = dir.path().join("qc.tsv");
    assert_geac_success(&[
        "qc",
        "--output", tsv.to_str().unwrap(),
        pq.to_str().unwrap(),
    ]);

    assert!(tsv.exists(), "QC TSV not created");
    let content = std::fs::read_to_string(&tsv).unwrap();
    assert!(content.contains("sample1"), "QC TSV should contain sample1");
}

// ── geac cohort ───────────────────────────────────────────────────────────────

/// cohort identifies a locus shared across two samples.
#[test]
fn cohort_finds_recurrent_locus() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    for (name, sample) in [("s1.bam", "sampleA"), ("s2.bam", "sampleB")] {
        let bam = write_bam(
            dir.path(), name, sample, 200,
            vec![(50, b'T', 3, 10)],
            20,
        );
        let pq = dir.path().join(name.replace(".bam", ".parquet"));
        assert_geac_success(&[
            "collect",
            "--input",     bam.to_str().unwrap(),
            "--reference", fa.to_str().unwrap(),
            "--output",    pq.to_str().unwrap(),
            "--read-type", "raw",
            "--pipeline",  "raw",
        ]);
    }

    let tsv = dir.path().join("cohort.tsv");
    assert_geac_success(&[
        "cohort",
        "--output",      tsv.to_str().unwrap(),
        "--min-samples", "2",
        dir.path().join("s1.parquet").to_str().unwrap(),
        dir.path().join("s2.parquet").to_str().unwrap(),
    ]);

    assert!(tsv.exists(), "cohort TSV not created");
    let content = std::fs::read_to_string(&tsv).unwrap();
    assert!(content.contains("chr1"), "cohort TSV should contain chr1");
}
