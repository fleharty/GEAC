mod common;

use common::{
    assert_geac_success, count_bam_reads, duckdb_count, duckdb_table_exists,
    parquet_count, parquet_query_f32, parquet_query_i32, parquet_query_i64, parquet_query_str,
    write_bam, write_bed, write_coverage_bam, write_paired_bam, write_reference,
    write_reference_with_base, CovRead,
};
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

// ── geac annotate-normal ───────────────────────────────────────────────────────

/// Basic annotate-normal: normal has 10 ref reads at the tumor alt position.
/// Expect an anchor row (normal_alt_allele IS NULL) with normal_depth=10, normal_alt_count=0.
#[test]
fn annotate_normal_produces_output() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    let tumor_bam = write_bam(dir.path(), "tumor.bam", "tumor1", 200,
        vec![(50, b'T', 5, 10)], 20);
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect", "--input", tumor_bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output", tumor_pq.to_str().unwrap(),
        "--read-type", "raw", "--pipeline", "raw",
    ]);

    // Normal: 10 ref (A) reads at pos 50, no alts.
    let normal_bam = write_bam(dir.path(), "normal.bam", "normal1", 200,
        vec![(50, b'T', 0, 10)], 20);
    let out_pq = dir.path().join("tumor.normal_evidence.parquet");

    assert_geac_success(&[
        "annotate-normal",
        "--tumor-parquet", tumor_pq.to_str().unwrap(),
        "--normal-bam",    normal_bam.to_str().unwrap(),
        "--reference",     fa.to_str().unwrap(),
        "--output",        out_pq.to_str().unwrap(),
    ]);

    assert!(out_pq.exists(), "normal_evidence Parquet not created");
    // Anchor row must exist (normal_alt_allele IS NULL).
    assert!(
        parquet_count(&out_pq) >= 1,
        "expected at least one row in normal_evidence Parquet"
    );
    assert_eq!(
        parquet_query_i32(&out_pq, "normal_depth",     "pos = 50 AND normal_alt_allele IS NULL"),
        10,
        "normal_depth should be 10"
    );
    assert_eq!(
        parquet_query_i32(&out_pq, "normal_alt_count", "pos = 50 AND normal_alt_allele IS NULL"),
        0,
        "normal_alt_count should be 0 for all-ref normal"
    );
}

/// When the normal also carries the tumor alt allele at the same position,
/// annotate-normal emits a per-allele row with normal_alt_count > 0.
#[test]
fn annotate_normal_detects_normal_alt() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    let tumor_bam = write_bam(dir.path(), "tumor.bam", "tumor1", 200,
        vec![(50, b'T', 5, 7)], 20);
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect", "--input", tumor_bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output", tumor_pq.to_str().unwrap(),
        "--read-type", "raw", "--pipeline", "raw",
    ]);

    // Normal: 3 T alts + 7 ref reads at pos 50.
    let normal_bam = write_bam(dir.path(), "normal.bam", "normal1", 200,
        vec![(50, b'T', 3, 7)], 20);
    let out_pq = dir.path().join("tumor.normal_evidence.parquet");

    assert_geac_success(&[
        "annotate-normal",
        "--tumor-parquet", tumor_pq.to_str().unwrap(),
        "--normal-bam",    normal_bam.to_str().unwrap(),
        "--reference",     fa.to_str().unwrap(),
        "--output",        out_pq.to_str().unwrap(),
    ]);

    // Per-allele row for T must have normal_alt_count = 3.
    assert_eq!(
        parquet_query_i32(&out_pq, "normal_alt_count",
            "pos = 50 AND normal_alt_allele = 'T'"),
        3,
        "normal_alt_count should be 3 for the T allele row"
    );
}

/// When the normal has zero coverage at a tumor locus, annotate-normal emits
/// an anchor row with normal_depth=0.
#[test]
fn annotate_normal_zero_coverage() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    let tumor_bam = write_bam(dir.path(), "tumor.bam", "tumor1", 200,
        vec![(50, b'T', 5, 10)], 20);
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect", "--input", tumor_bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output", tumor_pq.to_str().unwrap(),
        "--read-type", "raw", "--pipeline", "raw",
    ]);

    // Normal: no reads at all.
    let normal_bam = write_bam(dir.path(), "normal.bam", "normal1", 200, vec![], 20);
    let out_pq = dir.path().join("tumor.normal_evidence.parquet");

    assert_geac_success(&[
        "annotate-normal",
        "--tumor-parquet", tumor_pq.to_str().unwrap(),
        "--normal-bam",    normal_bam.to_str().unwrap(),
        "--reference",     fa.to_str().unwrap(),
        "--output",        out_pq.to_str().unwrap(),
    ]);

    assert_eq!(
        parquet_query_i32(&out_pq, "normal_depth",
            "pos = 50 AND normal_alt_allele IS NULL"),
        0,
        "zero-coverage position should have normal_depth = 0"
    );
}

// ── geac annotate-pon ─────────────────────────────────────────────────────────

/// Helper: build a PoN DuckDB from synthetic normal samples.
fn build_pon_db(
    dir: &std::path::Path,
    fa: &std::path::Path,
    normals: &[(&str, &str, Vec<common::Locus>)],
    db_name: &str,
) -> std::path::PathBuf {
    let mut parquet_paths: Vec<std::path::PathBuf> = Vec::new();
    for (bam_name, sample, loci) in normals {
        let bam = write_bam(dir, bam_name, sample, 200, loci.clone(), 20);
        let pq = dir.join(format!("{sample}.parquet"));
        assert_geac_success(&[
            "collect", "--input", bam.to_str().unwrap(),
            "--reference", fa.to_str().unwrap(),
            "--output", pq.to_str().unwrap(),
            "--read-type", "raw", "--pipeline", "raw",
        ]);
        parquet_paths.push(pq);
    }
    let db = dir.join(db_name);
    let mut args = vec!["merge", "--output", db.to_str().unwrap()];
    let path_strs: Vec<String> = parquet_paths.iter()
        .map(|p| p.to_str().unwrap().to_owned())
        .collect();
    for s in &path_strs { args.push(s.as_str()); }
    assert_geac_success(&args);
    db
}

/// When the tumor alt allele is absent from all PoN samples, n_pon_samples = 0.
#[test]
fn annotate_pon_locus_absent_from_pon() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // Tumor: T alt at pos 50.
    let tumor_bam = write_bam(dir.path(), "tumor.bam", "tumor1", 200,
        vec![(50, b'T', 5, 10)], 20);
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect", "--input", tumor_bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output", tumor_pq.to_str().unwrap(),
        "--read-type", "raw", "--pipeline", "raw",
    ]);

    // PoN: 2 normals with G alt at pos 100 only — no overlap with tumor locus.
    let pon_db = build_pon_db(dir.path(), &fa, &[
        ("pon1.bam", "pon1", vec![(100, b'G', 2, 10)]),
        ("pon2.bam", "pon2", vec![(100, b'G', 2, 10)]),
    ], "pon.duckdb");

    let out_pq = dir.path().join("tumor.pon_evidence.parquet");
    assert_geac_success(&[
        "annotate-pon",
        "--tumor-parquet", tumor_pq.to_str().unwrap(),
        "--pon-db",        pon_db.to_str().unwrap(),
        "--output",        out_pq.to_str().unwrap(),
    ]);

    assert!(out_pq.exists(), "pon_evidence Parquet not created");
    assert_eq!(
        parquet_query_i64(&out_pq, "n_pon_samples",
            "pos = 50 AND tumor_alt_allele = 'T'"),
        0,
        "tumor locus absent from PoN should have n_pon_samples = 0"
    );
}

/// When 2 of 3 PoN samples carry the same alt at the tumor locus, n_pon_samples=2
/// and pon_total_samples=3.
#[test]
fn annotate_pon_locus_present_in_pon() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // Tumor: T alt at pos 50.
    let tumor_bam = write_bam(dir.path(), "tumor.bam", "tumor1", 200,
        vec![(50, b'T', 5, 10)], 20);
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect", "--input", tumor_bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output", tumor_pq.to_str().unwrap(),
        "--read-type", "raw", "--pipeline", "raw",
    ]);

    // PoN: 3 normals. 2 have T at pos 50; all 3 have G at pos 100 (ensures
    // pon_total_samples = 3, since COUNT(DISTINCT sample_id) in alt_bases = 3).
    let pon_db = build_pon_db(dir.path(), &fa, &[
        ("pon1.bam", "pon1", vec![(50, b'T', 3, 10), (100, b'G', 2, 10)]),
        ("pon2.bam", "pon2", vec![(50, b'T', 2, 10), (100, b'G', 2, 10)]),
        ("pon3.bam", "pon3", vec![(100, b'G', 1, 10)]),
    ], "pon.duckdb");

    let out_pq = dir.path().join("tumor.pon_evidence.parquet");
    assert_geac_success(&[
        "annotate-pon",
        "--tumor-parquet", tumor_pq.to_str().unwrap(),
        "--pon-db",        pon_db.to_str().unwrap(),
        "--output",        out_pq.to_str().unwrap(),
    ]);

    assert_eq!(
        parquet_query_i64(&out_pq, "n_pon_samples",
            "pos = 50 AND tumor_alt_allele = 'T'"),
        2,
        "2 of 3 PoN samples carry T at pos 50"
    );
    assert_eq!(
        parquet_query_i64(&out_pq, "pon_total_samples",
            "pos = 50 AND tumor_alt_allele = 'T'"),
        3,
        "PoN has 3 total samples"
    );
}

// ── geac merge suffix routing ─────────────────────────────────────────────────

/// .reads.parquet files (from --reads-output) are routed to the alt_reads table.
#[test]
fn merge_routes_reads_parquet_to_alt_reads_table() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(dir.path(), "s1.bam", "sample1", 200,
        vec![(50, b'T', 5, 10)], 20);

    // --reads-output produces s1.locus.parquet + s1.reads.parquet.
    let out_stem = dir.path().join("s1.parquet");
    assert_geac_success(&[
        "collect",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out_stem.to_str().unwrap(),
        "--read-type", "raw", "--pipeline", "raw",
        "--reads-output",
    ]);

    let locus_pq = dir.path().join("s1.locus.parquet");
    let reads_pq = dir.path().join("s1.reads.parquet");
    assert!(locus_pq.exists(), "locus parquet not created");
    assert!(reads_pq.exists(), "reads parquet not created");

    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge", "--output", db.to_str().unwrap(),
        locus_pq.to_str().unwrap(),
        reads_pq.to_str().unwrap(),
    ]);

    assert!(duckdb_table_exists(&db, "alt_reads"), "alt_reads table not created");
    assert!(duckdb_count(&db, "alt_reads") > 0, "alt_reads table is empty");
}

/// .normal_evidence.parquet files are routed to the normal_evidence table.
#[test]
fn merge_routes_normal_evidence_parquet() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    let tumor_bam = write_bam(dir.path(), "tumor.bam", "tumor1", 200,
        vec![(50, b'T', 5, 10)], 20);
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect", "--input", tumor_bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output", tumor_pq.to_str().unwrap(),
        "--read-type", "raw", "--pipeline", "raw",
    ]);

    let normal_bam = write_bam(dir.path(), "normal.bam", "normal1", 200,
        vec![(50, b'T', 0, 10)], 20);
    let ne_pq = dir.path().join("tumor.normal_evidence.parquet");
    assert_geac_success(&[
        "annotate-normal",
        "--tumor-parquet", tumor_pq.to_str().unwrap(),
        "--normal-bam",    normal_bam.to_str().unwrap(),
        "--reference",     fa.to_str().unwrap(),
        "--output",        ne_pq.to_str().unwrap(),
    ]);

    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge", "--output", db.to_str().unwrap(),
        tumor_pq.to_str().unwrap(),
        ne_pq.to_str().unwrap(),
    ]);

    assert!(duckdb_table_exists(&db, "normal_evidence"), "normal_evidence table not created");
    assert!(duckdb_count(&db, "normal_evidence") > 0, "normal_evidence table is empty");
}

/// .pon_evidence.parquet files are routed to the pon_evidence table.
#[test]
fn merge_routes_pon_evidence_parquet() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    let tumor_bam = write_bam(dir.path(), "tumor.bam", "tumor1", 200,
        vec![(50, b'T', 5, 10)], 20);
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect", "--input", tumor_bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output", tumor_pq.to_str().unwrap(),
        "--read-type", "raw", "--pipeline", "raw",
    ]);

    let pon_db = build_pon_db(dir.path(), &fa, &[
        ("pon1.bam", "pon1", vec![(100, b'G', 2, 10)]),
    ], "pon.duckdb");

    let pe_pq = dir.path().join("tumor.pon_evidence.parquet");
    assert_geac_success(&[
        "annotate-pon",
        "--tumor-parquet", tumor_pq.to_str().unwrap(),
        "--pon-db",        pon_db.to_str().unwrap(),
        "--output",        pe_pq.to_str().unwrap(),
    ]);

    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge", "--output", db.to_str().unwrap(),
        tumor_pq.to_str().unwrap(),
        pe_pq.to_str().unwrap(),
    ]);

    assert!(duckdb_table_exists(&db, "pon_evidence"), "pon_evidence table not created");
    assert!(duckdb_count(&db, "pon_evidence") > 0, "pon_evidence table is empty");
}

// ── geac coverage ─────────────────────────────────────────────────────────────

/// Basic smoke test: 15 regular reads at pos 50 → total_depth = 15, frac_dup = 0.
#[test]
fn coverage_basic_depth() {
    let dir = TempDir::new().unwrap();
    let fa  = write_reference(dir.path(), 200);
    // 15 regular reads all starting at pos 40, covering positions 40–59 (read_len=20)
    let reads: Vec<CovRead> = (0..15).map(|_| CovRead::regular(40)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
    ]);

    assert!(out.exists(), "coverage Parquet not created");
    // Pos 40 is the start of reads; all 15 reads cover it.
    assert_eq!(parquet_query_i32(&out, "total_depth",    "pos = 40"), 15);
    assert_eq!(parquet_query_i32(&out, "raw_read_depth", "pos = 40"), 15);
    // All reads are regular (no dups), so frac_dup = 0.
    let frac_dup = parquet_query_f32(&out, "frac_dup", "pos = 40");
    assert!((frac_dup - 0.0).abs() < 1e-6, "frac_dup should be 0.0, got {frac_dup}");
}

/// Duplicate reads are counted in raw_read_depth but excluded from total_depth.
/// frac_dup = n_dup / raw_read_depth.
#[test]
fn coverage_frac_dup_excludes_duplicates() {
    let dir = TempDir::new().unwrap();
    let fa  = write_reference(dir.path(), 200);

    // 10 regular + 5 duplicate reads at pos 40
    let mut reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    reads.extend((0..5).map(|_| CovRead::duplicate(40)));
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
    ]);

    assert_eq!(parquet_query_i32(&out, "raw_read_depth", "pos = 40"), 15,
        "raw_read_depth should count all reads including dups");
    assert_eq!(parquet_query_i32(&out, "total_depth",    "pos = 40"), 10,
        "total_depth should exclude duplicate reads");

    let frac_dup = parquet_query_f32(&out, "frac_dup", "pos = 40");
    // 5 dups / 15 total = 0.3333...
    assert!((frac_dup - 5.0_f32 / 15.0).abs() < 1e-4,
        "frac_dup should be 5/15 ≈ {:.4}, got {frac_dup:.4}", 5.0_f32 / 15.0);
}

/// MAPQ=0 reads are tracked in frac_mapq0 and excluded from total_depth when
/// --min-map-qual is set above 0.
#[test]
fn coverage_mapq0_tracked_and_excluded() {
    let dir = TempDir::new().unwrap();
    let fa  = write_reference(dir.path(), 200);

    // 10 regular (mapq=60) + 5 mapq=0 reads at pos 40
    let mut reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    reads.extend((0..5).map(|_| CovRead::mapq0(40)));
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",       bam.to_str().unwrap(),
        "--reference",   fa.to_str().unwrap(),
        "--output",      out.to_str().unwrap(),
        "--read-type",   "raw",
        "--pipeline",    "raw",
        "--min-map-qual", "20",   // mapq=0 reads fail the filter → excluded from total_depth
    ]);

    // raw_read_depth counts non-secondary, non-supplementary reads (both regular and mapq0)
    assert_eq!(parquet_query_i32(&out, "raw_read_depth", "pos = 40"), 15);
    // total_depth only counts reads passing --min-map-qual (10 reads with mapq=60)
    assert_eq!(parquet_query_i32(&out, "total_depth",    "pos = 40"), 10,
        "mapq=0 reads should be excluded from total_depth with --min-map-qual 20");

    // frac_mapq0 = 5 mapq0 reads / 15 non-dup reads
    let frac_mapq0 = parquet_query_f32(&out, "frac_mapq0", "pos = 40");
    assert!((frac_mapq0 - 5.0_f32 / 15.0).abs() < 1e-4,
        "frac_mapq0 should be 5/15 ≈ {:.4}, got {frac_mapq0:.4}", 5.0_f32 / 15.0);
}

/// GC content is computed from the reference sequence, not from reads.
/// With an all-C reference, every position should have gc_content = 1.0.
#[test]
fn coverage_gc_content_computed_from_reference() {
    let dir = TempDir::new().unwrap();
    // All-C reference → GC fraction = 1.0
    let fa  = write_reference_with_base(dir.path(), 200, b'C');
    let reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
        "--gc-window", "10",
    ]);

    let gc = parquet_query_f32(&out, "gc_content", "pos = 40");
    assert!((gc - 1.0).abs() < 1e-6, "gc_content should be 1.0 for all-C reference, got {gc}");
}

/// With an all-A reference, gc_content should be 0.0 everywhere.
#[test]
fn coverage_gc_content_zero_for_all_a_reference() {
    let dir = TempDir::new().unwrap();
    let fa  = write_reference(dir.path(), 200);
    let reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
    ]);

    let gc = parquet_query_f32(&out, "gc_content", "pos = 40");
    assert!((gc - 0.0).abs() < 1e-6, "gc_content should be 0.0 for all-A reference, got {gc}");
}

/// When --targets is provided, zero-depth positions within the target interval
/// are emitted with total_depth = 0 (dropout capture).
#[test]
fn coverage_targets_emits_zero_depth_positions() {
    let dir = TempDir::new().unwrap();
    let fa  = write_reference(dir.path(), 200);

    // Reads at pos 40, read_len=20 → cover positions 40–59 (depth > 0)
    let reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);

    // Target: positions 10–59 (50 bases). Positions 10–39 have no reads → zero depth.
    let bed = write_bed(dir.path(), "targets.bed", &[(10, 60)]);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
        "--targets",   bed.to_str().unwrap(),
    ]);

    // All 50 target positions should appear in the output
    assert_eq!(parquet_count(&out), 50, "expected one row per target position");

    // Positions with reads have depth > 0
    assert!(parquet_query_i32(&out, "total_depth", "pos = 40") > 0,
        "pos 40 should have reads");

    // Zero-depth positions are recorded with total_depth = 0
    assert_eq!(parquet_query_i32(&out, "total_depth", "pos = 10"), 0,
        "pos 10 (no reads) should have total_depth = 0");
    assert_eq!(parquet_query_i32(&out, "total_depth", "pos = 39"), 0,
        "pos 39 (no reads) should have total_depth = 0");
}

/// Without --targets, only positions covered by reads are emitted.
/// Positions outside the covered region do not appear.
#[test]
fn coverage_no_targets_skips_zero_depth() {
    let dir = TempDir::new().unwrap();
    let fa  = write_reference(dir.path(), 200);

    // Reads at pos 40, read_len=20 → cover positions 40–59 only
    let reads: Vec<CovRead> = (0..5).map(|_| CovRead::regular(40)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
    ]);

    // Only covered positions are emitted (positions 40–59 = 20 rows)
    let n = parquet_count(&out);
    assert_eq!(n, 20, "expected 20 covered positions, got {n}");
}

/// insert_size is collected from properly-paired R1 reads; mean_insert_size reflects
/// the TLEN set on those reads and n_insert_size_obs counts them.
#[test]
fn coverage_insert_size_from_paired_reads() {
    let dir = TempDir::new().unwrap();
    let fa  = write_reference(dir.path(), 500);

    // 8 proper-pair R1 reads at pos 40 with TLEN = 300
    let reads: Vec<CovRead> = (0..8).map(|_| CovRead::r1_paired(40, 300)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 500, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
    ]);

    assert_eq!(parquet_query_i32(&out, "n_insert_size_obs", "pos = 40"), 8,
        "expected 8 insert size observations");
    let mean_ins = parquet_query_f32(&out, "mean_insert_size", "pos = 40");
    assert!((mean_ins - 300.0).abs() < 1.0,
        "mean_insert_size should be 300.0, got {mean_ins}");
}

/// .coverage.parquet files are routed to the coverage table by geac merge.
#[test]
fn merge_routes_coverage_parquet_to_coverage_table() {
    let dir = TempDir::new().unwrap();
    let fa  = write_reference(dir.path(), 200);

    // Build a regular locus parquet (required by merge; coverage-only merge not supported)
    let locus_bam = write_bam(dir.path(), "locus.bam", "sample1", 200,
        vec![(50, b'T', 5, 10)], 20);
    let locus_pq  = dir.path().join("locus.parquet");
    assert_geac_success(&[
        "collect",
        "--input",     locus_bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    locus_pq.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
    ]);

    // Build a coverage parquet
    let reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    let cov_bam = write_coverage_bam(dir.path(), "cov.bam", "sample1", 200, reads, 20);
    let cov_pq  = dir.path().join("sample1.coverage.parquet");
    assert_geac_success(&[
        "coverage",
        "--input",     cov_bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    cov_pq.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
    ]);

    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge", "--output", db.to_str().unwrap(),
        locus_pq.to_str().unwrap(),
        cov_pq.to_str().unwrap(),
    ]);

    assert!(duckdb_table_exists(&db, "coverage"), "coverage table not created");
    assert!(duckdb_count(&db, "coverage") > 0,    "coverage table is empty");
}

// ── fwd/rev alt count read-level semantics ────────────────────────────────────

/// For overlapping pairs where both R1 (forward) and R2 (reverse) carry the alt,
/// fwd_alt_count and rev_alt_count each equal the number of alt pairs.
/// This verifies that fwd + rev = 2 * alt_count for standard overlapping pairs,
/// distinguishing read-level counts from the fragment-level alt_count.
#[test]
fn collect_fwd_rev_alt_counts_are_read_level() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    // 2 pairs where both R1+R2 carry alt, 8 pairs where both carry ref
    let bam = write_paired_bam(
        dir.path(), "paired.bam", "sample1", 200,
        vec![(50, b'T', 2, 8, 0)],
        20,
    );
    let out = dir.path().join("paired.parquet");

    assert_geac_success(&[
        "collect",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
    ]);

    assert_eq!(parquet_query_i32(&out, "alt_count",     "pos = 50"), 2,
        "fragment-level alt count: 2 pairs");
    assert_eq!(parquet_query_i32(&out, "fwd_alt_count", "pos = 50"), 2,
        "R1 forward reads carrying alt: 2");
    assert_eq!(parquet_query_i32(&out, "rev_alt_count", "pos = 50"), 2,
        "R2 reverse reads carrying alt: 2");
}

/// For pairs where only R2 carries the alt (R2-only artefact pattern),
/// fwd_alt_count must be 0 and rev_alt_count must equal the number of such pairs.
/// Before the fix, rev_alt_count was always 0 because strand was attributed using
/// R1's orientation for the whole fragment instead of each read's own is_reverse flag.
#[test]
fn collect_r2_only_artefact_visible_in_rev_alt_count() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    // 2 R2-only alt pairs (R1=ref, R2=alt) + 8 both-ref pairs
    let bam = write_paired_bam(
        dir.path(), "r2only.bam", "sample1", 200,
        vec![(50, b'T', 0, 8, 2)],
        20,
    );
    let out = dir.path().join("r2only.parquet");

    assert_geac_success(&[
        "collect",
        "--input",     bam.to_str().unwrap(),
        "--reference", fa.to_str().unwrap(),
        "--output",    out.to_str().unwrap(),
        "--read-type", "raw",
        "--pipeline",  "raw",
    ]);

    assert_eq!(parquet_query_i32(&out, "alt_count",     "pos = 50"), 2,
        "fragment-level alt count: 2 R2-only alt pairs");
    assert_eq!(parquet_query_i32(&out, "fwd_alt_count", "pos = 50"), 0,
        "R1 forward reads do not carry alt: fwd_alt_count = 0");
    assert_eq!(parquet_query_i32(&out, "rev_alt_count", "pos = 50"), 2,
        "R2 reverse reads carry alt: rev_alt_count = 2");
}
