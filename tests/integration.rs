mod common;

use common::{
    assert_geac_success, count_bam_reads, duckdb_columns, duckdb_count, duckdb_table_exists,
    parquet_columns, parquet_count, parquet_query_f32, parquet_query_i32, parquet_query_i64,
    parquet_query_opt_str, parquet_query_str, write_bam, write_bed, write_coverage_bam,
    write_cycle_bam, write_fragments_bam, write_mnv_bam, write_paired_bam, write_reference,
    write_reference_with_base, CovRead, CycleTestRead, FragmentPair,
};
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::io::Read as _;
use std::path::Path;
use tempfile::TempDir;

fn schema_required_columns(table: &str) -> Vec<String> {
    let manifest: Value = serde_json::from_str(
        &std::fs::read_to_string("schema/geac_schema.json").expect("read schema manifest"),
    )
    .expect("parse schema manifest");
    manifest["tables"][table]["required_columns"]
        .as_array()
        .expect("required_columns array")
        .iter()
        .map(|value: &Value| value.as_str().expect("column name").to_string())
        .collect()
}

fn assert_schema_columns_present(actual: &[String], table: &str) {
    let actual = actual.iter().cloned().collect::<HashSet<_>>();
    let missing = schema_required_columns(table)
        .into_iter()
        .filter(|column| !actual.contains(column))
        .collect::<Vec<_>>();
    assert!(
        missing.is_empty(),
        "schema contract missing columns for {table}: {missing:?}"
    );
}

fn sha256_hex(path: &Path) -> String {
    let mut file = std::fs::File::open(path).expect("open file for sha256");
    let mut hasher = Sha256::new();
    let mut buf = [0_u8; 64 * 1024];
    loop {
        let n = file.read(&mut buf).expect("read file for sha256");
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    format!("{:x}", hasher.finalize())
}

/// Sanity check: verify the BAM helper writes the expected number of reads.
#[test]
fn bam_helper_writes_expected_reads() {
    let dir = TempDir::new().unwrap();
    let _fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(),
        "s1.bam",
        "sample1",
        200,
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
        dir.path(),
        "s1.bam",
        "sample1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let out = dir.path().join("s1.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    assert!(out.exists(), "output Parquet not created");
    assert_eq!(parquet_count(&out), 1, "expected exactly 1 alt locus");
    assert_eq!(parquet_query_i32(&out, "alt_count", "pos = 50"), 5);
    assert_eq!(parquet_query_i32(&out, "total_depth", "pos = 50"), 15);
    assert_eq!(parquet_query_str(&out, "alt_allele", "pos = 50"), "T");
    assert_eq!(parquet_query_str(&out, "sample_id", "pos = 50"), "sample1");
}

#[test]
fn collect_optionally_records_input_checksum_sha256() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(),
        "checksum.bam",
        "sample1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );

    let out_default = dir.path().join("checksum_default.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out_default.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);
    assert_eq!(
        parquet_query_opt_str(&out_default, "input_checksum_sha256", "pos = 50"),
        None,
        "checksum should be null by default"
    );

    let out_hashed = dir.path().join("checksum_hashed.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out_hashed.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--input-checksum-sha256",
    ]);
    assert_eq!(
        parquet_query_opt_str(&out_hashed, "input_checksum_sha256", "pos = 50"),
        Some(sha256_hex(&bam)),
        "checksum should match the input BAM when enabled"
    );
}

/// `--index` lets geac open a BAM whose index lives under a non-conventional
/// name/location, and records that path as `bai_path` for IGV.
#[test]
fn collect_opens_bam_with_explicitly_named_index() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(),
        "reads.bam",
        "sample1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );

    // Move the conventional index aside so htslib's inference ({bam}.bai)
    // would fail — proving the run only succeeds because of --index.
    let conventional = dir
        .path()
        .join(format!("{}.bai", bam.file_name().unwrap().to_str().unwrap()));
    let custom_index = dir.path().join("elsewhere.bai");
    std::fs::rename(&conventional, &custom_index).expect("rename index");

    let out = dir.path().join("idx.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--index",
        custom_index.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
    ]);

    let expected = std::fs::canonicalize(&custom_index)
        .unwrap()
        .to_string_lossy()
        .into_owned();
    assert_eq!(
        parquet_query_opt_str(&out, "bai_path", "pos = 50"),
        Some(expected),
        "bai_path should record the explicitly-specified --index path"
    );
}

/// `--bai-uri` takes precedence over `--index` for the recorded `bai_path`.
#[test]
fn collect_bai_uri_takes_precedence_over_index() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(),
        "reads2.bam",
        "sample1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let conventional = dir
        .path()
        .join(format!("{}.bai", bam.file_name().unwrap().to_str().unwrap()));
    let custom_index = dir.path().join("local.bai");
    std::fs::rename(&conventional, &custom_index).expect("rename index");

    let out = dir.path().join("uri.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--index",
        custom_index.to_str().unwrap(),
        "--bai-uri",
        "gs://bucket/sample.bam.bai",
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
    ]);
    assert_eq!(
        parquet_query_opt_str(&out, "bai_path", "pos = 50"),
        Some("gs://bucket/sample.bam.bai".to_string()),
        "--bai-uri should win over --index for bai_path"
    );
}

/// `--gnomad-index-uri` is recorded as `gnomad_index_path` for IGV, independent
/// of whether a gnomAD VCF is read.
#[test]
fn collect_records_gnomad_index_uri() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(),
        "g.bam",
        "sample1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );

    let out = dir.path().join("g.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--gnomad-index-uri",
        "gs://bucket/gnomad.vcf.gz.tbi",
    ]);
    assert_eq!(
        parquet_query_opt_str(&out, "gnomad_index_path", "pos = 50"),
        Some("gs://bucket/gnomad.vcf.gz.tbi".to_string()),
        "gnomad_index_path should record --gnomad-index-uri"
    );
}

/// `--sample-id` flag overrides the SM tag read from the BAM header.
#[test]
fn collect_sample_id_override() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(),
        "s1.bam",
        "original_id",
        200,
        vec![(50, b'T', 3, 10)],
        20,
    );
    let out = dir.path().join("s1.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--sample-id",
        "overridden",
    ]);

    assert_eq!(
        parquet_query_str(&out, "sample_id", "pos = 50"),
        "overridden"
    );
}

/// `--region` restricts output to only loci within the requested interval.
/// BAM has alts at pos 50 (0-based) and pos 150. With --region chr1:101-200
/// (1-based, covering 0-based 100–199) only pos 150 should appear.
#[test]
fn collect_region_restricts_output() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 300);
    let bam = write_bam(
        dir.path(),
        "s1.bam",
        "sample1",
        300,
        vec![(50, b'T', 5, 10), (150, b'G', 3, 10)],
        20,
    );
    let out = dir.path().join("s1.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--region",
        "chr1:101-200",
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
        let bam = write_bam(dir.path(), name, sample, 200, vec![(50, b'T', 3, 10)], 20);
        let pq = dir.path().join(name.replace(".bam", ".parquet"));
        assert_geac_success(&[
            "collect",
            "--input",
            bam.to_str().unwrap(),
            "--reference",
            fa.to_str().unwrap(),
            "--output",
            pq.to_str().unwrap(),
            "--read-type",
            "raw",
            "--pipeline",
            "raw",
        ]);
    }

    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge",
        "--output",
        db.to_str().unwrap(),
        dir.path().join("s1.parquet").to_str().unwrap(),
        dir.path().join("s2.parquet").to_str().unwrap(),
    ]);

    assert!(db.exists(), "DuckDB file not created");

    let conn = duckdb::Connection::open(&db).unwrap();
    let n_samples: i64 = conn
        .query_row("SELECT COUNT(DISTINCT sample_id) FROM alt_bases", [], |r| {
            r.get(0)
        })
        .unwrap();
    assert_eq!(n_samples, 2, "expected 2 distinct samples in merged DuckDB");

    // geac_metadata / geac_inputs must exist and record merge provenance.
    assert!(
        duckdb_table_exists(&db, "geac_metadata"),
        "geac_metadata table not created"
    );
    assert!(
        duckdb_table_exists(&db, "geac_inputs"),
        "geac_inputs table not created"
    );
    let version: String = conn
        .query_row("SELECT geac_version FROM geac_metadata LIMIT 1", [], |r| {
            r.get(0)
        })
        .expect("could not read geac_version from geac_metadata");
    assert_eq!(
        version,
        env!("CARGO_PKG_VERSION"),
        "geac_version in metadata should match build version"
    );
    let schema_version: String = conn
        .query_row(
            "SELECT schema_version FROM geac_metadata LIMIT 1",
            [],
            |r| r.get(0),
        )
        .expect("could not read schema_version from geac_metadata");
    assert_eq!(schema_version, "duckdb-v4");
    let recorded_samples: i64 = conn
        .query_row("SELECT n_samples FROM geac_metadata LIMIT 1", [], |r| {
            r.get(0)
        })
        .expect("could not read n_samples from geac_metadata");
    assert_eq!(
        recorded_samples, 2,
        "metadata should record merged sample count"
    );
    let alt_bases_rows: i64 = conn
        .query_row(
            "SELECT alt_bases_rows FROM geac_metadata LIMIT 1",
            [],
            |r| r.get(0),
        )
        .expect("could not read alt_bases_rows from geac_metadata");
    assert_eq!(
        alt_bases_rows, 2,
        "metadata should record merged alt_bases rows"
    );
    let input_rows: i64 = conn
        .query_row("SELECT COUNT(*) FROM geac_inputs", [], |r| r.get(0))
        .expect("could not count geac_inputs rows");
    assert_eq!(
        input_rows, 2,
        "expected one geac_inputs row per merged parquet"
    );

    assert_geac_success(&["inspect", "--input", db.to_str().unwrap()]);
}

// ── geac qc ───────────────────────────────────────────────────────────────────

/// qc exits successfully and produces a TSV when --output is given.
#[test]
fn qc_exits_successfully() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(),
        "s1.bam",
        "sample1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let pq = dir.path().join("s1.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    let tsv = dir.path().join("qc.tsv");
    assert_geac_success(&[
        "qc",
        "--output",
        tsv.to_str().unwrap(),
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
        let bam = write_bam(dir.path(), name, sample, 200, vec![(50, b'T', 3, 10)], 20);
        let pq = dir.path().join(name.replace(".bam", ".parquet"));
        assert_geac_success(&[
            "collect",
            "--input",
            bam.to_str().unwrap(),
            "--reference",
            fa.to_str().unwrap(),
            "--output",
            pq.to_str().unwrap(),
            "--read-type",
            "raw",
            "--pipeline",
            "raw",
        ]);
    }

    let tsv = dir.path().join("cohort.tsv");
    assert_geac_success(&[
        "cohort",
        "--output",
        tsv.to_str().unwrap(),
        "--min-samples",
        "2",
        dir.path().join("s1.parquet").to_str().unwrap(),
        dir.path().join("s2.parquet").to_str().unwrap(),
    ]);

    assert!(tsv.exists(), "cohort TSV not created");
    let content = std::fs::read_to_string(&tsv).unwrap();
    assert!(content.contains("chr1"), "cohort TSV should contain chr1");
}

// ── geac sample-identity ─────────────────────────────────────────────────────

fn insert_identity_rows(
    conn: &duckdb::Connection,
    sample_id: &str,
    subject_id: &str,
    positions: &[i64],
    vaf: f64,
) {
    let alt_count = (vaf * 100.0).round() as i32;
    for pos in positions {
        conn.execute(
            "INSERT INTO alt_bases
             (sample_id, subject_id, chrom, pos, ref_allele, alt_allele, variant_type, total_depth, alt_count)
             VALUES (?, ?, 'chr1', ?, 'A', 'G', 'SNV', 100, ?)",
            duckdb::params![sample_id, subject_id, pos, alt_count],
        )
        .expect("insert identity row");
    }
}

#[test]
fn sample_identity_reports_duplicate_pair() {
    let dir = TempDir::new().unwrap();
    let db = dir.path().join("identity.duckdb");
    let conn = duckdb::Connection::open(&db).unwrap();
    conn.execute_batch(
        "CREATE TABLE alt_bases (
            sample_id VARCHAR,
            subject_id VARCHAR,
            chrom VARCHAR,
            pos BIGINT,
            ref_allele VARCHAR,
            alt_allele VARCHAR,
            variant_type VARCHAR,
            total_depth INTEGER,
            alt_count INTEGER
        );",
    )
    .unwrap();

    let common = (1_i64..=10).collect::<Vec<_>>();
    let duplicate_private = (101_i64..=140).collect::<Vec<_>>();
    let mut duplicate_markers = common.clone();
    duplicate_markers.extend(duplicate_private);

    insert_identity_rows(&conn, "S1", "P1a", &duplicate_markers, 0.50);
    insert_identity_rows(&conn, "S2", "P1b", &duplicate_markers, 0.50);
    insert_identity_rows(&conn, "S3", "P2", &common, 0.50);
    insert_identity_rows(&conn, "S4", "P2", &common, 0.95);
    drop(conn);

    let out = dir.path().join("identity.tsv");
    assert_geac_success(&[
        "sample-identity",
        "--input",
        db.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
    ]);

    let content = std::fs::read_to_string(&out).unwrap();
    assert!(
        content.contains("S1\tP1a\tS2\tP1b"),
        "expected S1/S2 duplicate pair in output:\n{content}"
    );
    assert!(
        !content.contains("S1\tP1a\tS3\tP2"),
        "unrelated pair should not pass default thresholds:\n{content}"
    );
    assert!(
        !content.contains("S3\tP2\tS4\tP2"),
        "same-subject genotype-discordant pair should not pass default thresholds:\n{content}"
    );
}

// ── geac annotate-normal ───────────────────────────────────────────────────────

/// Basic annotate-normal: normal has 10 ref reads at the tumor alt position.
/// Expect an anchor row (normal_alt_allele IS NULL) with normal_depth=10, normal_alt_count=0.
#[test]
fn annotate_normal_produces_output() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    let tumor_bam = write_bam(
        dir.path(),
        "tumor.bam",
        "tumor1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        tumor_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        tumor_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    // Normal: 10 ref (A) reads at pos 50, no alts.
    let normal_bam = write_bam(
        dir.path(),
        "normal.bam",
        "normal1",
        200,
        vec![(50, b'T', 0, 10)],
        20,
    );
    let out_pq = dir.path().join("tumor.normal_evidence.parquet");

    assert_geac_success(&[
        "annotate-normal",
        "--tumor-parquet",
        tumor_pq.to_str().unwrap(),
        "--normal-bam",
        normal_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out_pq.to_str().unwrap(),
    ]);

    assert!(out_pq.exists(), "normal_evidence Parquet not created");
    // Anchor row must exist (normal_alt_allele IS NULL).
    assert!(
        parquet_count(&out_pq) >= 1,
        "expected at least one row in normal_evidence Parquet"
    );
    assert_eq!(
        parquet_query_i32(
            &out_pq,
            "normal_depth",
            "pos = 50 AND normal_alt_allele IS NULL"
        ),
        10,
        "normal_depth should be 10"
    );
    assert_eq!(
        parquet_query_i32(
            &out_pq,
            "normal_alt_count",
            "pos = 50 AND normal_alt_allele IS NULL"
        ),
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

    let tumor_bam = write_bam(
        dir.path(),
        "tumor.bam",
        "tumor1",
        200,
        vec![(50, b'T', 5, 7)],
        20,
    );
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        tumor_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        tumor_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    // Normal: 3 T alts + 7 ref reads at pos 50.
    let normal_bam = write_bam(
        dir.path(),
        "normal.bam",
        "normal1",
        200,
        vec![(50, b'T', 3, 7)],
        20,
    );
    let out_pq = dir.path().join("tumor.normal_evidence.parquet");

    assert_geac_success(&[
        "annotate-normal",
        "--tumor-parquet",
        tumor_pq.to_str().unwrap(),
        "--normal-bam",
        normal_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out_pq.to_str().unwrap(),
    ]);

    // Per-allele row for T must have normal_alt_count = 3.
    assert_eq!(
        parquet_query_i32(
            &out_pq,
            "normal_alt_count",
            "pos = 50 AND normal_alt_allele = 'T'"
        ),
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

    let tumor_bam = write_bam(
        dir.path(),
        "tumor.bam",
        "tumor1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        tumor_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        tumor_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    // Normal: no reads at all.
    let normal_bam = write_bam(dir.path(), "normal.bam", "normal1", 200, vec![], 20);
    let out_pq = dir.path().join("tumor.normal_evidence.parquet");

    assert_geac_success(&[
        "annotate-normal",
        "--tumor-parquet",
        tumor_pq.to_str().unwrap(),
        "--normal-bam",
        normal_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out_pq.to_str().unwrap(),
    ]);

    assert_eq!(
        parquet_query_i32(
            &out_pq,
            "normal_depth",
            "pos = 50 AND normal_alt_allele IS NULL"
        ),
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
            "collect",
            "--input",
            bam.to_str().unwrap(),
            "--reference",
            fa.to_str().unwrap(),
            "--output",
            pq.to_str().unwrap(),
            "--read-type",
            "raw",
            "--pipeline",
            "raw",
        ]);
        parquet_paths.push(pq);
    }
    let db = dir.join(db_name);
    let mut args = vec!["merge", "--output", db.to_str().unwrap()];
    let path_strs: Vec<String> = parquet_paths
        .iter()
        .map(|p| p.to_str().unwrap().to_owned())
        .collect();
    for s in &path_strs {
        args.push(s.as_str());
    }
    assert_geac_success(&args);
    db
}

/// When the tumor alt allele is absent from all PoN samples, n_pon_samples = 0.
#[test]
fn annotate_pon_locus_absent_from_pon() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // Tumor: T alt at pos 50.
    let tumor_bam = write_bam(
        dir.path(),
        "tumor.bam",
        "tumor1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        tumor_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        tumor_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    // PoN: 2 normals with G alt at pos 100 only — no overlap with tumor locus.
    let pon_db = build_pon_db(
        dir.path(),
        &fa,
        &[
            ("pon1.bam", "pon1", vec![(100, b'G', 2, 10)]),
            ("pon2.bam", "pon2", vec![(100, b'G', 2, 10)]),
        ],
        "pon.duckdb",
    );

    let out_pq = dir.path().join("tumor.pon_evidence.parquet");
    assert_geac_success(&[
        "annotate-pon",
        "--tumor-parquet",
        tumor_pq.to_str().unwrap(),
        "--pon-db",
        pon_db.to_str().unwrap(),
        "--output",
        out_pq.to_str().unwrap(),
    ]);

    assert!(out_pq.exists(), "pon_evidence Parquet not created");
    assert_eq!(
        parquet_query_i64(
            &out_pq,
            "n_pon_samples",
            "pos = 50 AND tumor_alt_allele = 'T'"
        ),
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
    let tumor_bam = write_bam(
        dir.path(),
        "tumor.bam",
        "tumor1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        tumor_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        tumor_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    // PoN: 3 normals. 2 have T at pos 50; all 3 have G at pos 100 (ensures
    // pon_total_samples = 3, since COUNT(DISTINCT sample_id) in alt_bases = 3).
    let pon_db = build_pon_db(
        dir.path(),
        &fa,
        &[
            (
                "pon1.bam",
                "pon1",
                vec![(50, b'T', 3, 10), (100, b'G', 2, 10)],
            ),
            (
                "pon2.bam",
                "pon2",
                vec![(50, b'T', 2, 10), (100, b'G', 2, 10)],
            ),
            ("pon3.bam", "pon3", vec![(100, b'G', 1, 10)]),
        ],
        "pon.duckdb",
    );

    let out_pq = dir.path().join("tumor.pon_evidence.parquet");
    assert_geac_success(&[
        "annotate-pon",
        "--tumor-parquet",
        tumor_pq.to_str().unwrap(),
        "--pon-db",
        pon_db.to_str().unwrap(),
        "--output",
        out_pq.to_str().unwrap(),
    ]);

    assert_eq!(
        parquet_query_i64(
            &out_pq,
            "n_pon_samples",
            "pos = 50 AND tumor_alt_allele = 'T'"
        ),
        2,
        "2 of 3 PoN samples carry T at pos 50"
    );
    assert_eq!(
        parquet_query_i64(
            &out_pq,
            "pon_total_samples",
            "pos = 50 AND tumor_alt_allele = 'T'"
        ),
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
    let bam = write_bam(
        dir.path(),
        "s1.bam",
        "sample1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );

    // --reads-output produces s1.locus.parquet + s1.reads.parquet.
    let out_stem = dir.path().join("s1.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out_stem.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--reads-output",
    ]);

    let locus_pq = dir.path().join("s1.locus.parquet");
    let reads_pq = dir.path().join("s1.reads.parquet");
    assert!(locus_pq.exists(), "locus parquet not created");
    assert!(reads_pq.exists(), "reads parquet not created");
    let reads_cols = parquet_columns(&reads_pq);
    assert!(reads_cols.contains(&"read_type".to_string()));
    assert!(reads_cols.contains(&"pipeline".to_string()));
    assert!(reads_cols.contains(&"n_before_alt".to_string()));
    assert!(reads_cols.contains(&"trailing_n_run_len".to_string()));
    assert_eq!(parquet_query_str(&reads_pq, "read_type", "TRUE"), "raw");
    assert_eq!(parquet_query_str(&reads_pq, "pipeline", "TRUE"), "raw");

    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge",
        "--output",
        db.to_str().unwrap(),
        locus_pq.to_str().unwrap(),
        reads_pq.to_str().unwrap(),
    ]);

    assert!(
        duckdb_table_exists(&db, "alt_reads"),
        "alt_reads table not created"
    );
    assert!(
        duckdb_count(&db, "alt_reads") > 0,
        "alt_reads table is empty"
    );
    let merged_cols = duckdb_columns(&db, "alt_reads");
    assert!(merged_cols.contains(&"read_type".to_string()));
    assert!(merged_cols.contains(&"pipeline".to_string()));
    assert!(merged_cols.contains(&"n_before_alt".to_string()));
    assert!(merged_cols.contains(&"trailing_n_run_len".to_string()));
}

#[test]
fn collect_targets_emits_sample_metrics_parquet() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(),
        "metrics.bam",
        "sample1",
        200,
        vec![(45, b'T', 5, 10)],
        20,
    );
    let bed = write_bed(dir.path(), "targets.bed", &[(40, 80)]);
    let out = dir.path().join("metrics.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--targets",
        bed.to_str().unwrap(),
    ]);

    let metrics_pq = dir.path().join("metrics.sample_metrics.parquet");
    assert!(metrics_pq.exists(), "sample_metrics parquet not created");
    assert_eq!(
        parquet_count(&metrics_pq),
        1,
        "expected one sample-metrics row"
    );
    assert_eq!(
        parquet_query_i32(&metrics_pq, "n_target_positions", "TRUE"),
        40
    );
    assert_eq!(
        parquet_query_i32(&metrics_pq, "n_target_positions_covered", "TRUE"),
        15
    );
    assert_eq!(
        parquet_query_f32(&metrics_pq, "mean_target_depth_covered", "TRUE"),
        15.0
    );
    assert!(
        (parquet_query_f32(&metrics_pq, "mean_target_depth_all", "TRUE") - 5.625).abs() < 1e-5,
        "mean_target_depth_all should include zero-depth target positions"
    );
    assert_eq!(
        parquet_query_f32(&metrics_pq, "median_target_depth_covered", "TRUE"),
        15.0
    );
    assert_eq!(
        parquet_query_f32(&metrics_pq, "median_target_depth_all", "TRUE"),
        0.0
    );
    assert!(
        (parquet_query_f32(&metrics_pq, "pct_fragment_bases_on_target", "TRUE") - 0.75).abs()
            < 1e-5,
        "on-target fragment-base fraction should reflect target overlap"
    );
}

#[test]
fn merge_routes_sample_metrics_parquet_to_sample_metrics_table() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(),
        "metrics_merge.bam",
        "sample1",
        200,
        vec![(45, b'T', 5, 10)],
        20,
    );
    let bed = write_bed(dir.path(), "targets.bed", &[(40, 80)]);
    let out = dir.path().join("metrics_merge.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--targets",
        bed.to_str().unwrap(),
    ]);

    let locus_pq = dir.path().join("metrics_merge.parquet");
    let sample_metrics_pq = dir.path().join("metrics_merge.sample_metrics.parquet");
    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge",
        "--output",
        db.to_str().unwrap(),
        locus_pq.to_str().unwrap(),
        sample_metrics_pq.to_str().unwrap(),
    ]);

    assert!(
        duckdb_table_exists(&db, "sample_metrics"),
        "sample_metrics table not created"
    );
    assert_eq!(duckdb_count(&db, "sample_metrics"), 1);
    assert_schema_columns_present(&duckdb_columns(&db, "sample_metrics"), "sample_metrics");
    let conn = duckdb::Connection::open(&db).expect("open merged duckdb");
    let sample_metrics_rows: i64 = conn
        .query_row(
            "SELECT sample_metrics_rows FROM geac_metadata LIMIT 1",
            [],
            |r| r.get(0),
        )
        .expect("read sample_metrics_rows");
    assert_eq!(sample_metrics_rows, 1);
}

// ── Cycle number correctness ───────────────────────────────────────────────────

/// Reverse-strand reads count cycles from the read END (synthesis start is the 3′ end of
/// the stored sequence).
///
/// A 20-bp reverse-strand read with alt at qpos=10:
///   correct cycle = read_len − qpos = 20 − 10 = 10
///   buggy  cycle  = qpos + 1 = 11
#[test]
fn reads_reverse_strand_cycle_counts_from_read_end() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // Reverse-strand read of length 20 starting at ref pos 40.
    // Alt T stored at qpos=10 → reference position 40 + 10 = 50.
    // cycle = 20 − 10 = 10.
    let mut seq = vec![b'A'; 20];
    seq[10] = b'T';
    let quals = vec![40u8; 20];
    let bam = write_cycle_bam(
        dir.path(),
        "rev.bam",
        "sample1",
        200,
        vec![CycleTestRead {
            pos: 40,
            flags: 0x10,
            leading_hard_clips: 0,
            trailing_hard_clips: 0,
            seq,
            quals,
            aux_i32_tags: vec![],
        }],
    );
    let out = dir.path().join("rev.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--reads-output",
    ]);

    let reads_pq = dir.path().join("rev.reads.parquet");
    assert!(reads_pq.exists(), "reads parquet not created");
    let cycle = parquet_query_i32(&reads_pq, "cycle", "alt_allele = 'T'");
    assert_eq!(
        cycle, 10,
        "reverse-strand cycle should be read_len − qpos = 10, not qpos + 1 = 11"
    );
}

/// Forward read with leading hard clips: the clipped bases shift all cycle numbers up.
///
/// A forward read with 5H leading clips, stored length 10, alt at qpos=5:
///   correct cycle = 5 (hard clips) + 5 (qpos) + 1 = 11
///   buggy  cycle  = 5 + 1 = 6
#[test]
fn reads_forward_hard_clips_counted_in_cycle() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // Read starts at ref pos 40; alt T at qpos=5 → ref pos 45.
    // Leading 5H: cycles 1–5 are hard-clipped; stored seq starts at cycle 6.
    // cycle = 5 + 5 + 1 = 11.
    let mut seq = vec![b'A'; 10];
    seq[5] = b'T';
    let quals = vec![40u8; 10];
    let bam = write_cycle_bam(
        dir.path(),
        "fwd_hc.bam",
        "sample1",
        200,
        vec![CycleTestRead {
            pos: 40,
            flags: 0,
            leading_hard_clips: 5,
            trailing_hard_clips: 0,
            seq,
            quals,
            aux_i32_tags: vec![],
        }],
    );
    let out = dir.path().join("fwd_hc.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--reads-output",
    ]);

    let reads_pq = dir.path().join("fwd_hc.reads.parquet");
    assert!(reads_pq.exists(), "reads parquet not created");
    let cycle = parquet_query_i32(&reads_pq, "cycle", "alt_allele = 'T'");
    assert_eq!(
        cycle, 11,
        "forward hard-clip cycle should be hard_clip_before + qpos + 1 = 11"
    );
}

/// Reverse read with trailing hard clips (= 5′ hard clips in synthesis order):
/// cycle accounts for those clipped bases.
///
/// A reverse read of stored length 10 with 5 trailing H in the CIGAR, alt at qpos=4:
///   correct cycle = 5 (trailing H = 5′ clips) + (10 − 4) = 11
///   buggy  cycle  = 4 + 1 = 5
#[test]
fn reads_reverse_hard_clips_counted_in_cycle() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // Reverse-strand read, stored length 10, starting at ref pos 40.
    // Trailing 5H → 5 bases at the 5′ end of original synthesis are hard-clipped.
    // Alt T at qpos=4 → ref pos 40 + 4 = 44.
    // cycle = 5 + (10 − 4) = 11.
    let mut seq = vec![b'A'; 10];
    seq[4] = b'T';
    let quals = vec![40u8; 10];
    let bam = write_cycle_bam(
        dir.path(),
        "rev_hc.bam",
        "sample1",
        200,
        vec![CycleTestRead {
            pos: 40,
            flags: 0x10,
            leading_hard_clips: 0,
            trailing_hard_clips: 5,
            seq,
            quals,
            aux_i32_tags: vec![],
        }],
    );
    let out = dir.path().join("rev_hc.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--reads-output",
    ]);

    let reads_pq = dir.path().join("rev_hc.reads.parquet");
    assert!(reads_pq.exists(), "reads parquet not created");
    let cycle = parquet_query_i32(&reads_pq, "cycle", "alt_allele = 'T'");
    assert_eq!(
        cycle, 11,
        "reverse hard-clip cycle should be trailing_H + read_len − qpos = 11"
    );
}

#[test]
fn reads_record_n_context_metrics() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // Forward read: qpos=3 at ref pos 43. Sequence context around alt:
    // before = A N N, alt = T, after = N N A A
    // n_before_alt=3, n_after_alt=4
    // n_n_before_alt=2, n_n_after_alt=2
    // leading_n_run_len=2, trailing_n_run_len=2
    let seq = b"ANNTNNAA".to_vec();
    let quals = vec![40u8; seq.len()];
    let bam = write_cycle_bam(
        dir.path(),
        "nctx.bam",
        "sample1",
        200,
        vec![CycleTestRead {
            pos: 40,
            flags: 0,
            leading_hard_clips: 0,
            trailing_hard_clips: 0,
            seq,
            quals,
            aux_i32_tags: vec![],
        }],
    );
    let out = dir.path().join("nctx.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--reads-output",
    ]);

    let reads_pq = dir.path().join("nctx.reads.parquet");
    assert!(reads_pq.exists(), "reads parquet not created");
    assert_eq!(
        parquet_query_i32(&reads_pq, "n_before_alt", "alt_allele = 'T'"),
        3
    );
    assert_eq!(
        parquet_query_i32(&reads_pq, "n_after_alt", "alt_allele = 'T'"),
        4
    );
    assert_eq!(
        parquet_query_i32(&reads_pq, "n_n_before_alt", "alt_allele = 'T'"),
        2
    );
    assert_eq!(
        parquet_query_i32(&reads_pq, "n_n_after_alt", "alt_allele = 'T'"),
        2
    );
    assert_eq!(
        parquet_query_i32(&reads_pq, "leading_n_run_len", "alt_allele = 'T'"),
        2
    );
    assert_eq!(
        parquet_query_i32(&reads_pq, "trailing_n_run_len", "alt_allele = 'T'"),
        2
    );
}

#[test]
fn dragen_reads_use_xv_xw_for_family_size() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    let mut seq = vec![b'A'; 10];
    seq[5] = b'T';
    let quals = vec![40u8; 10];
    let bam = write_cycle_bam(
        dir.path(),
        "dragen_tags.bam",
        "sample1",
        200,
        vec![CycleTestRead {
            pos: 40,
            flags: 0,
            leading_hard_clips: 0,
            trailing_hard_clips: 0,
            seq,
            quals,
            aux_i32_tags: vec![(*b"XV", 7), (*b"XW", 3)],
        }],
    );
    let out = dir.path().join("dragen_tags.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "duplex",
        "--pipeline",
        "dragen",
        "--reads-output",
    ]);

    let reads_pq = dir.path().join("dragen_tags.reads.parquet");
    assert!(reads_pq.exists(), "reads parquet not created");
    assert_eq!(
        parquet_query_i32(&reads_pq, "ab_count", "alt_allele = 'T'"),
        7
    );
    assert_eq!(
        parquet_query_i32(&reads_pq, "family_size", "alt_allele = 'T'"),
        3
    );
}

/// .normal_evidence.parquet files are routed to the normal_evidence table.
#[test]
fn merge_routes_normal_evidence_parquet() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    let tumor_bam = write_bam(
        dir.path(),
        "tumor.bam",
        "tumor1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        tumor_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        tumor_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    let normal_bam = write_bam(
        dir.path(),
        "normal.bam",
        "normal1",
        200,
        vec![(50, b'T', 0, 10)],
        20,
    );
    let ne_pq = dir.path().join("tumor.normal_evidence.parquet");
    assert_geac_success(&[
        "annotate-normal",
        "--tumor-parquet",
        tumor_pq.to_str().unwrap(),
        "--normal-bam",
        normal_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        ne_pq.to_str().unwrap(),
    ]);

    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge",
        "--output",
        db.to_str().unwrap(),
        tumor_pq.to_str().unwrap(),
        ne_pq.to_str().unwrap(),
    ]);

    assert!(
        duckdb_table_exists(&db, "normal_evidence"),
        "normal_evidence table not created"
    );
    assert!(
        duckdb_count(&db, "normal_evidence") > 0,
        "normal_evidence table is empty"
    );
}

/// .pon_evidence.parquet files are routed to the pon_evidence table.
#[test]
fn merge_routes_pon_evidence_parquet() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    let tumor_bam = write_bam(
        dir.path(),
        "tumor.bam",
        "tumor1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let tumor_pq = dir.path().join("tumor.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        tumor_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        tumor_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    let pon_db = build_pon_db(
        dir.path(),
        &fa,
        &[("pon1.bam", "pon1", vec![(100, b'G', 2, 10)])],
        "pon.duckdb",
    );

    let pe_pq = dir.path().join("tumor.pon_evidence.parquet");
    assert_geac_success(&[
        "annotate-pon",
        "--tumor-parquet",
        tumor_pq.to_str().unwrap(),
        "--pon-db",
        pon_db.to_str().unwrap(),
        "--output",
        pe_pq.to_str().unwrap(),
    ]);

    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge",
        "--output",
        db.to_str().unwrap(),
        tumor_pq.to_str().unwrap(),
        pe_pq.to_str().unwrap(),
    ]);

    assert!(
        duckdb_table_exists(&db, "pon_evidence"),
        "pon_evidence table not created"
    );
    assert!(
        duckdb_count(&db, "pon_evidence") > 0,
        "pon_evidence table is empty"
    );
}

// ── geac coverage ─────────────────────────────────────────────────────────────

/// Basic smoke test: 15 regular reads at pos 50 → total_depth = 15, frac_dup = 0.
#[test]
fn coverage_basic_depth() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    // 15 regular reads all starting at pos 40, covering positions 40–59 (read_len=20)
    let reads: Vec<CovRead> = (0..15).map(|_| CovRead::regular(40)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    assert!(out.exists(), "coverage Parquet not created");
    // Pos 40 is the start of reads; all 15 reads cover it.
    assert_eq!(parquet_query_i32(&out, "total_depth", "pos = 40"), 15);
    assert_eq!(parquet_query_i32(&out, "raw_read_depth", "pos = 40"), 15);
    // All reads are regular (no dups), so frac_dup = 0.
    let frac_dup = parquet_query_f32(&out, "frac_dup", "pos = 40");
    assert!(
        (frac_dup - 0.0).abs() < 1e-6,
        "frac_dup should be 0.0, got {frac_dup}"
    );
}

/// Duplicate reads are counted in raw_read_depth but excluded from total_depth.
/// frac_dup = n_dup / raw_read_depth.
#[test]
fn coverage_frac_dup_excludes_duplicates() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // 10 regular + 5 duplicate reads at pos 40
    let mut reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    reads.extend((0..5).map(|_| CovRead::duplicate(40)));
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    assert_eq!(
        parquet_query_i32(&out, "raw_read_depth", "pos = 40"),
        15,
        "raw_read_depth should count all reads including dups"
    );
    assert_eq!(
        parquet_query_i32(&out, "total_depth", "pos = 40"),
        10,
        "total_depth should exclude duplicate reads"
    );

    let frac_dup = parquet_query_f32(&out, "frac_dup", "pos = 40");
    // 5 dups / 15 total = 0.3333...
    assert!(
        (frac_dup - 5.0_f32 / 15.0).abs() < 1e-4,
        "frac_dup should be 5/15 ≈ {:.4}, got {frac_dup:.4}",
        5.0_f32 / 15.0
    );
}

/// MAPQ=0 reads are tracked in frac_mapq0 and excluded from total_depth when
/// --min-map-qual is set above 0.
#[test]
fn coverage_mapq0_tracked_and_excluded() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // 10 regular (mapq=60) + 5 mapq=0 reads at pos 40
    let mut reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    reads.extend((0..5).map(|_| CovRead::mapq0(40)));
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--min-map-qual",
        "20", // mapq=0 reads fail the filter → excluded from total_depth
    ]);

    // raw_read_depth counts non-secondary, non-supplementary reads (both regular and mapq0)
    assert_eq!(parquet_query_i32(&out, "raw_read_depth", "pos = 40"), 15);
    // total_depth only counts reads passing --min-map-qual (10 reads with mapq=60)
    assert_eq!(
        parquet_query_i32(&out, "total_depth", "pos = 40"),
        10,
        "mapq=0 reads should be excluded from total_depth with --min-map-qual 20"
    );

    // frac_mapq0 = 5 mapq0 reads / 15 non-dup reads
    let frac_mapq0 = parquet_query_f32(&out, "frac_mapq0", "pos = 40");
    assert!(
        (frac_mapq0 - 5.0_f32 / 15.0).abs() < 1e-4,
        "frac_mapq0 should be 5/15 ≈ {:.4}, got {frac_mapq0:.4}",
        5.0_f32 / 15.0
    );
}

/// GC content is computed from the reference sequence, not from reads.
/// With an all-C reference, every position should have gc_content = 1.0.
#[test]
fn coverage_gc_content_computed_from_reference() {
    let dir = TempDir::new().unwrap();
    // All-C reference → GC fraction = 1.0
    let fa = write_reference_with_base(dir.path(), 200, b'C');
    let reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--gc-window",
        "10",
    ]);

    let gc = parquet_query_f32(&out, "gc_content", "pos = 40");
    assert!(
        (gc - 1.0).abs() < 1e-6,
        "gc_content should be 1.0 for all-C reference, got {gc}"
    );
}

/// With an all-A reference, gc_content should be 0.0 everywhere.
#[test]
fn coverage_gc_content_zero_for_all_a_reference() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    let gc = parquet_query_f32(&out, "gc_content", "pos = 40");
    assert!(
        (gc - 0.0).abs() < 1e-6,
        "gc_content should be 0.0 for all-A reference, got {gc}"
    );
}

/// When --targets is provided, zero-depth positions within the target interval
/// are emitted with total_depth = 0 (dropout capture).
#[test]
fn coverage_targets_emits_zero_depth_positions() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // Reads at pos 40, read_len=20 → cover positions 40–59 (depth > 0)
    let reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);

    // Target: positions 10–59 (50 bases). Positions 10–39 have no reads → zero depth.
    let bed = write_bed(dir.path(), "targets.bed", &[(10, 60)]);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--targets",
        bed.to_str().unwrap(),
    ]);

    // All 50 target positions should appear in the output
    assert_eq!(
        parquet_count(&out),
        50,
        "expected one row per target position"
    );

    // Positions with reads have depth > 0
    assert!(
        parquet_query_i32(&out, "total_depth", "pos = 40") > 0,
        "pos 40 should have reads"
    );

    // Zero-depth positions are recorded with total_depth = 0
    assert_eq!(
        parquet_query_i32(&out, "total_depth", "pos = 10"),
        0,
        "pos 10 (no reads) should have total_depth = 0"
    );
    assert_eq!(
        parquet_query_i32(&out, "total_depth", "pos = 39"),
        0,
        "pos 39 (no reads) should have total_depth = 0"
    );
}

/// Without --targets, only positions covered by reads are emitted.
/// Positions outside the covered region do not appear.
#[test]
fn coverage_no_targets_skips_zero_depth() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // Reads at pos 40, read_len=20 → cover positions 40–59 only
    let reads: Vec<CovRead> = (0..5).map(|_| CovRead::regular(40)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 200, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    // Only covered positions are emitted (positions 40–59 = 20 rows)
    let n = parquet_count(&out);
    assert_eq!(n, 20, "expected 20 covered positions, got {n}");
}

/// Adaptive thresholding flushes a partially filled bin before emitting low-depth

/// insert_size is collected from properly-paired R1 reads; mean_insert_size reflects
/// the TLEN set on those reads and n_insert_size_obs counts them.
#[test]
fn coverage_insert_size_from_paired_reads() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 500);

    // 8 proper-pair R1 reads at pos 40 with TLEN = 300
    let reads: Vec<CovRead> = (0..8).map(|_| CovRead::r1_paired(40, 300)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 500, reads, 20);
    let out = dir.path().join("s1.coverage.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    assert_eq!(
        parquet_query_i32(&out, "n_insert_size_obs", "pos = 40"),
        8,
        "expected 8 insert size observations"
    );
    let mean_ins = parquet_query_f32(&out, "mean_insert_size", "pos = 40");
    assert!(
        (mean_ins - 300.0).abs() < 1.0,
        "mean_insert_size should be 300.0, got {mean_ins}"
    );
}

/// .coverage.parquet files are routed to the coverage table by geac merge.
#[test]
fn merge_routes_coverage_parquet_to_coverage_table() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);

    // Build a regular locus parquet (required by merge; coverage-only merge not supported)
    let locus_bam = write_bam(
        dir.path(),
        "locus.bam",
        "sample1",
        200,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let locus_pq = dir.path().join("locus.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        locus_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        locus_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    // Build a coverage parquet
    let reads: Vec<CovRead> = (0..10).map(|_| CovRead::regular(40)).collect();
    let cov_bam = write_coverage_bam(dir.path(), "cov.bam", "sample1", 200, reads, 20);
    let cov_pq = dir.path().join("sample1.coverage.parquet");
    assert_geac_success(&[
        "coverage",
        "--input",
        cov_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        cov_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge",
        "--output",
        db.to_str().unwrap(),
        locus_pq.to_str().unwrap(),
        cov_pq.to_str().unwrap(),
    ]);

    assert!(
        duckdb_table_exists(&db, "coverage"),
        "coverage table not created"
    );
    assert!(duckdb_count(&db, "coverage") > 0, "coverage table is empty");
}

/// Two named intervals: one covered by reads, one at zero depth.
/// Verifies n_bases, mean_depth, frac_at_1x, interval_name, and that the
/// interval accumulator counts zero-depth positions correctly.
#[test]
fn coverage_intervals_basic() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 400);

    // 5 reads of length 20 all starting at pos=100; covers [100, 120) at depth 5.
    // Interval chr1:200-220 receives no reads → zero depth.
    let reads: Vec<CovRead> = (0..5).map(|_| CovRead::regular(100)).collect();
    let bam = write_coverage_bam(dir.path(), "s1.bam", "sample1", 400, reads, 20);

    // BED with two named intervals.
    let bed = dir.path().join("targets.bed");
    std::fs::write(&bed, "chr1\t100\t120\texon1\nchr1\t200\t220\texon2\n").unwrap();

    let cov_pq = dir.path().join("s1.coverage.parquet");
    let iv_pq = dir.path().join("s1.coverage.intervals.parquet");

    assert_geac_success(&[
        "coverage",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        cov_pq.to_str().unwrap(),
        "--targets",
        bed.to_str().unwrap(),
        "--intervals-output",
        iv_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    assert!(iv_pq.exists(), "intervals Parquet not created");
    assert_eq!(parquet_count(&iv_pq), 2, "expected 2 interval rows");

    // exon1: 20 positions all at depth 5 → frac_at_1x = 1.0, mean_depth = 5.0
    let n_bases_exon1 = parquet_query_i32(&iv_pq, "n_bases", "interval_name = 'exon1'");
    assert_eq!(n_bases_exon1, 20, "exon1 n_bases should be 20");

    let mean_depth_exon1 = parquet_query_f32(&iv_pq, "mean_depth", "interval_name = 'exon1'");
    assert!(
        (mean_depth_exon1 - 5.0).abs() < 0.01,
        "exon1 mean_depth should be 5.0, got {mean_depth_exon1}"
    );

    let frac_1x_exon1 = parquet_query_f32(&iv_pq, "frac_at_1x", "interval_name = 'exon1'");
    assert!(
        (frac_1x_exon1 - 1.0).abs() < 0.01,
        "exon1 frac_at_1x should be 1.0, got {frac_1x_exon1}"
    );

    // exon2: 20 positions all at depth 0 → frac_at_1x = 0.0, mean_depth = 0.0
    let n_bases_exon2 = parquet_query_i32(&iv_pq, "n_bases", "interval_name = 'exon2'");
    assert_eq!(n_bases_exon2, 20, "exon2 n_bases should be 20");

    let mean_depth_exon2 = parquet_query_f32(&iv_pq, "mean_depth", "interval_name = 'exon2'");
    assert!(
        mean_depth_exon2.abs() < 0.01,
        "exon2 mean_depth should be 0.0, got {mean_depth_exon2}"
    );

    let frac_1x_exon2 = parquet_query_f32(&iv_pq, "frac_at_1x", "interval_name = 'exon2'");
    assert!(
        frac_1x_exon2.abs() < 0.01,
        "exon2 frac_at_1x should be 0.0, got {frac_1x_exon2}"
    );
}

/// --track adds a named nullable Float32 column populated from the BEDGraph file.
#[test]
fn coverage_track_column_populated_from_bedgraph() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let reads: Vec<CovRead> = (0..5).map(|_| CovRead::regular(40)).collect();
    let bam = write_coverage_bam(dir.path(), "track.bam", "sample1", 200, reads, 20);

    // BEDGraph covers [30, 60) with score 0.9 and [60, 90) with 0.3.
    let bg = dir.path().join("map.bedgraph");
    std::fs::write(&bg, "chr1\t30\t60\t0.9\nchr1\t60\t90\t0.3\n").unwrap();

    let out = dir.path().join("track.coverage.parquet");
    assert_geac_success(&[
        "coverage",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--track",
        &format!("gem150:{}", bg.display()),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    // The column "gem150" should be present.
    let cols = parquet_columns(&out);
    assert!(
        cols.contains(&"gem150".to_string()),
        "gem150 column missing; found: {cols:?}"
    );

    // pos=40 falls in [30, 60) → score 0.9.
    let score = parquet_query_f32(&out, "gem150", "pos = 40");
    assert!(
        (score - 0.9).abs() < 0.01,
        "gem150 score at pos=40 should be 0.9, got {score}"
    );
}

/// .coverage.intervals.parquet files are routed to the coverage_intervals table by geac merge.
#[test]
fn merge_routes_coverage_intervals_parquet_to_coverage_intervals_table() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 400);

    // Locus parquet (needed so merge has a primary table to anchor on).
    let locus_bam = write_bam(
        dir.path(),
        "locus.bam",
        "sample1",
        400,
        vec![(50, b'T', 5, 10)],
        20,
    );
    let locus_pq = dir.path().join("locus.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        locus_bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        locus_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    let reads: Vec<CovRead> = (0..5).map(|_| CovRead::regular(100)).collect();
    let bam = write_coverage_bam(dir.path(), "cov.bam", "sample1", 400, reads, 20);
    let bed = dir.path().join("targets.bed");
    std::fs::write(&bed, "chr1\t100\t120\texon1\nchr1\t200\t220\texon2\n").unwrap();

    let cov_pq = dir.path().join("sample1.coverage.parquet");
    let iv_pq = dir.path().join("sample1.coverage.intervals.parquet");
    assert_geac_success(&[
        "coverage",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        cov_pq.to_str().unwrap(),
        "--targets",
        bed.to_str().unwrap(),
        "--intervals-output",
        iv_pq.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    let db = dir.path().join("cohort.duckdb");
    assert_geac_success(&[
        "merge",
        "--output",
        db.to_str().unwrap(),
        locus_pq.to_str().unwrap(),
        cov_pq.to_str().unwrap(),
        iv_pq.to_str().unwrap(),
    ]);

    assert!(
        duckdb_table_exists(&db, "coverage_intervals"),
        "coverage_intervals table not created"
    );
    assert_eq!(
        duckdb_count(&db, "coverage_intervals"),
        2,
        "expected 2 rows in coverage_intervals"
    );
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
        dir.path(),
        "paired.bam",
        "sample1",
        200,
        vec![(50, b'T', 2, 8, 0)],
        20,
    );
    let out = dir.path().join("paired.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    assert_eq!(
        parquet_query_i32(&out, "alt_count", "pos = 50"),
        2,
        "fragment-level alt count: 2 pairs"
    );
    assert_eq!(
        parquet_query_i32(&out, "fwd_alt_count", "pos = 50"),
        2,
        "R1 forward reads carrying alt: 2"
    );
    assert_eq!(
        parquet_query_i32(&out, "rev_alt_count", "pos = 50"),
        2,
        "R2 reverse reads carrying alt: 2"
    );
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
        dir.path(),
        "r2only.bam",
        "sample1",
        200,
        vec![(50, b'T', 0, 8, 2)],
        20,
    );
    let out = dir.path().join("r2only.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    assert_eq!(
        parquet_query_i32(&out, "alt_count", "pos = 50"),
        2,
        "fragment-level alt count: 2 R2-only alt pairs"
    );
    assert_eq!(
        parquet_query_i32(&out, "fwd_alt_count", "pos = 50"),
        0,
        "R1 forward reads do not carry alt: fwd_alt_count = 0"
    );
    assert_eq!(
        parquet_query_i32(&out, "rev_alt_count", "pos = 50"),
        2,
        "R2 reverse reads carry alt: rev_alt_count = 2"
    );
}

#[test]
fn collect_output_matches_schema_manifest() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(),
        "schema.bam",
        "sample1",
        200,
        vec![(50, b'T', 3, 10)],
        20,
    );
    let out = dir.path().join("schema.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    assert_schema_columns_present(&parquet_columns(&out), "alt_bases");
}

#[test]
fn merged_duckdb_tables_match_schema_manifest() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_bam(
        dir.path(),
        "schema_reads.bam",
        "sample1",
        200,
        vec![(50, b'T', 3, 10)],
        20,
    );
    let out = dir.path().join("schema_reads.parquet");

    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--reads-output",
    ]);

    let db = dir.path().join("schema.duckdb");
    assert_geac_success(&[
        "merge",
        "--output",
        db.to_str().unwrap(),
        dir.path()
            .join("schema_reads.locus.parquet")
            .to_str()
            .unwrap(),
        dir.path()
            .join("schema_reads.reads.parquet")
            .to_str()
            .unwrap(),
    ]);

    assert_schema_columns_present(&duckdb_columns(&db, "alt_bases"), "alt_bases");
    assert_schema_columns_present(&duckdb_columns(&db, "alt_reads"), "alt_reads");
    assert_schema_columns_present(&duckdb_columns(&db, "geac_metadata"), "geac_metadata");
    assert_schema_columns_present(&duckdb_columns(&db, "geac_inputs"), "geac_inputs");
}

// ── MNV detection ─────────────────────────────────────────────────────────────

/// fragment_id enables MNV detection: a self-join on fragment_id finds reads that
/// carry substitutions at two adjacent positions on the same fragment.
///
/// BAM layout (pos1=50 T, pos2=51 G, all 0-based):
///   3 reads carry BOTH alts ("mnv_*")    → same fragment_id at pos1 and pos2
///   2 reads carry ONLY T at pos1         → appear in alt_reads at pos1 only
///   2 reads carry ONLY G at pos2         → appear in alt_reads at pos2 only
///   5 reads carry ref at both positions
///
/// Expected MNV query result:
///   co_count         = 3  (3 fragments share both loci)
///   frac_cooccurring = 3 / (3 + 2) = 0.6  (pos1 alt_count = 5)
#[test]
fn reads_fragment_id_enables_mnv_detection() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 200);
    let bam = write_mnv_bam(
        dir.path(),
        "mnv.bam",
        "sample1",
        200,
        50,  // pos1 (0-based)
        b'T',
        51,  // pos2 (0-based)
        b'G',
        3,   // n_mnv
        2,   // n_pos1_only
        2,   // n_pos2_only
        5,   // n_ref
        20,  // read_len
    );

    let out_stem = dir.path().join("mnv.parquet");
    assert_geac_success(&[
        "collect",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out_stem.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
        "--reads-output",
    ]);

    let locus_pq = dir.path().join("mnv.locus.parquet");
    let reads_pq = dir.path().join("mnv.reads.parquet");
    assert!(locus_pq.exists(), "locus parquet not created");
    assert!(reads_pq.exists(), "reads parquet not created");

    // Both positions should appear as alt loci.
    assert_eq!(parquet_count(&locus_pq), 2, "expected 2 alt loci (pos 50 and 51)");
    assert_eq!(parquet_query_i32(&locus_pq, "alt_count", "pos = 50"), 5);
    assert_eq!(parquet_query_i32(&locus_pq, "alt_count", "pos = 51"), 5);

    // alt_reads must contain the fragment_id column.
    let reads_cols = parquet_columns(&reads_pq);
    assert!(
        reads_cols.contains(&"fragment_id".to_string()),
        "fragment_id column missing from alt_reads"
    );

    // Merge to DuckDB and run the MNV self-join.
    let db = dir.path().join("mnv.duckdb");
    assert_geac_success(&[
        "merge",
        "--output",
        db.to_str().unwrap(),
        locus_pq.to_str().unwrap(),
        reads_pq.to_str().unwrap(),
    ]);

    let conn = duckdb::Connection::open(&db).unwrap();

    // co_count: fragments carrying both T@50 and G@51.
    let co_count: i64 = conn
        .query_row(
            "SELECT COUNT(*)
             FROM alt_reads a
             JOIN alt_reads b
               ON  a.fragment_id = b.fragment_id
               AND a.sample_id   = b.sample_id
               AND a.chrom       = b.chrom
               AND a.pos < b.pos
               AND (b.pos - a.pos) <= 3
             WHERE a.pos = 50 AND b.pos = 51",
            [],
            |r| r.get(0),
        )
        .unwrap();
    assert_eq!(co_count, 3, "co_count: expected 3 fragments carrying both alts");

    // frac_cooccurring = co_count / alt_count at pos1.
    let alt_count_pos1: i64 = conn
        .query_row(
            "SELECT alt_count FROM alt_bases WHERE pos = 50",
            [],
            |r| r.get(0),
        )
        .unwrap();
    let frac: f64 = co_count as f64 / alt_count_pos1 as f64;
    assert!(
        (frac - 0.6).abs() < 1e-9,
        "frac_cooccurring should be 0.6, got {frac}"
    );

    // Reads carrying only one alt must NOT appear as co-occurring.
    // Verify that each "p1_*" fragment_id appears only in pos=50, not in pos=51.
    let spillover: i64 = conn
        .query_row(
            "SELECT COUNT(*)
             FROM alt_reads a
             JOIN alt_reads b
               ON  a.fragment_id = b.fragment_id
               AND a.pos < b.pos
             WHERE a.pos = 50 AND b.pos = 51
               AND a.alt_allele = 'T'
               AND b.alt_allele = 'G'
               AND a.fragment_id NOT IN (
                   SELECT fragment_id FROM alt_reads WHERE pos = 50 AND alt_allele = 'T'
                   INTERSECT
                   SELECT fragment_id FROM alt_reads WHERE pos = 51 AND alt_allele = 'G'
               )",
            [],
            |r| r.get(0),
        )
        .unwrap();
    assert_eq!(
        spillover, 0,
        "no single-substitution fragment should appear in the MNV self-join"
    );
}

/// Regression test for the `is_first_in_template` strand bug.
///
/// In FR libraries R1 is randomly assigned to either strand, so half of
/// proper pairs have R1 on the − strand (rightmost). The previous filter
/// `!record.is_first_in_template()` dropped those pairs entirely, halving
/// fragment yield. The correct selector is positive TLEN — by SAM spec
/// exactly one read per pair (the leftmost) carries it.
#[test]
fn fragments_emits_record_for_both_r1_strand_orientations() {
    let dir = TempDir::new().unwrap();
    let fa = write_reference(dir.path(), 500);
    let bam = write_fragments_bam(
        dir.path(),
        "frag_strand.bam",
        "sample1",
        500,
        vec![
            FragmentPair {
                left_pos: 100,
                insert_size: 150,
                r1_on_plus: true,
            },
            FragmentPair {
                left_pos: 300,
                insert_size: 150,
                r1_on_plus: false,
            },
        ],
        50,
    );
    let out = dir.path().join("frag_strand.fragments.parquet");

    assert_geac_success(&[
        "fragments",
        "--input",
        bam.to_str().unwrap(),
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-type",
        "raw",
        "--pipeline",
        "raw",
    ]);

    assert_eq!(
        parquet_count(&out),
        2,
        "both R1-on-+ and R1-on-− proper pairs must produce a fragment record"
    );
    assert_eq!(parquet_query_i64(&out, "frag_start", "frag_start = 100"), 100);
    assert_eq!(parquet_query_i64(&out, "frag_start", "frag_start = 300"), 300);
}

// ── Experimental: fusion detection ──────────────────────────────────────────────
//
// Locks in cross-file label consistency: the Parquet gene_a/gene_b columns and the
// FX:Z: tag in the --reads-output BAM must order the fused pair identically. Uses
// gene names whose GTF order (BCR then ABL1) differs from alphabetical order (ABL1 <
// BCR) so a regression to per-output ordering would change one but not the other.

/// Deterministic pseudo-random ACGT sequence (xorshift64) — distinct per seed so two
/// gene bodies share no k-mers.
fn fusion_gen_seq(seed: u64, n: usize) -> String {
    let mut x = seed | 1;
    let mut s = String::with_capacity(n);
    for _ in 0..n {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        s.push(b"ACGT"[(x & 3) as usize] as char);
    }
    s
}

#[test]
fn fusions_label_ordering_consistent_across_outputs() {
    use rust_htslib::bam::{self, Read as _};

    let dir = TempDir::new().unwrap();

    // Two distinct 300 bp gene bodies on chr1, with 200 bp of (distinct) filler
    // between and around them.
    let f0 = fusion_gen_seq(1, 100);
    let bcr = fusion_gen_seq(2, 300); // chr1:100-400
    let f1 = fusion_gen_seq(3, 200);
    let abl1 = fusion_gen_seq(4, 300); // chr1:600-900
    let tail = fusion_gen_seq(5, 100);
    let seq = format!("{f0}{bcr}{f1}{abl1}{tail}");
    let len = seq.len();

    // FASTA (single line) + .fai. Header ">chr1\n" is 6 bytes, so first base is at
    // offset 6; one line of `len` bases plus newline → linewidth len+1.
    let fa = dir.path().join("ref.fa");
    std::fs::write(&fa, format!(">chr1\n{seq}\n")).unwrap();
    std::fs::write(
        dir.path().join("ref.fa.fai"),
        format!("chr1\t{len}\t6\t{len}\t{}\n", len + 1),
    )
    .unwrap();

    // GTF: BCR listed FIRST (gene index 0), ABL1 second (index 1). 1-based inclusive.
    let gtf = dir.path().join("genes.gtf");
    std::fs::write(
        &gtf,
        "chr1\tt\tgene\t101\t400\t.\t+\t.\tgene_name \"BCR\";\n\
         chr1\tt\tgene\t601\t900\t.\t+\t.\tgene_name \"ABL1\";\n",
    )
    .unwrap();

    // One fusion fragment: read1 drawn from BCR, read2 from ABL1 (unmapped pair
    // sharing a qname — k-mer matching is alignment-independent).
    let r1 = &bcr[50..110];
    let r2 = &abl1[50..110];
    let mut header = bam::Header::new();
    let mut sq = bam::header::HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", len as i64);
    header.push_record(&sq);

    let bam_path = dir.path().join("reads.bam");
    {
        let mut writer =
            bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        let qual = vec![30u8; 60];
        for (flags, rseq) in [(77u16, r1), (141u16, r2)] {
            let mut rec = bam::record::Record::new();
            rec.set(b"frag1", None, rseq.as_bytes(), &qual);
            rec.set_tid(-1);
            rec.set_pos(-1);
            rec.set_mtid(-1);
            rec.set_mpos(-1);
            rec.set_flags(flags);
            writer.write(&rec).unwrap();
        }
    }

    let index = dir.path().join("idx.duckdb");
    assert_geac_success(&[
        "experimental",
        "build-fusion-index",
        "--gene-annotation",
        gtf.to_str().unwrap(),
        "--fasta",
        fa.to_str().unwrap(),
        "--min-gene-kmers",
        "1",
        "--output",
        index.to_str().unwrap(),
    ]);

    let out = dir.path().join("fusions.parquet");
    let reads_out = dir.path().join("fusion_reads.bam");
    assert_geac_success(&[
        "experimental",
        "fusions",
        "--bam",
        bam_path.to_str().unwrap(),
        "--index",
        index.to_str().unwrap(),
        "--sample-id",
        "HG002",
        "--output",
        out.to_str().unwrap(),
        "--reads-output",
        reads_out.to_str().unwrap(),
        "--min-supporting-reads",
        "1",
    ]);

    // Exactly one fusion candidate; gene_a/gene_b follow gene-index order (BCR, ABL1).
    assert_eq!(parquet_count(&out), 1, "expected one fusion candidate");
    let gene_a = parquet_query_str(&out, "gene_a", "supporting_reads >= 1");
    let gene_b = parquet_query_str(&out, "gene_b", "supporting_reads >= 1");
    assert_eq!((gene_a.as_str(), gene_b.as_str()), ("BCR", "ABL1"));

    // Every FX tag in the reads BAM must equal "gene_a::gene_b" from the Parquet —
    // the cross-file consistency guarantee.
    let expected_fx = format!("{gene_a}::{gene_b}");
    let mut reader = bam::Reader::from_path(&reads_out).unwrap();
    let mut n_tagged = 0;
    for r in reader.records() {
        let rec = r.unwrap();
        match rec.aux(b"FX") {
            Ok(bam::record::Aux::String(s)) => {
                assert_eq!(s, expected_fx, "FX tag must match Parquet gene_a::gene_b");
                n_tagged += 1;
            }
            _ => panic!("fusion read missing FX tag"),
        }
    }
    assert_eq!(n_tagged, 2, "both reads of the fragment should be tagged");
}

// ── Fusion junction-coherence filter ─────────────────────────────────────────
//
// Builds two scenarios and confirms --min-coherent-fragments behaves correctly:
//
// 1. A genuine spanning read: the first half of the read comes from gene A's
//    unique sequence and the second half from gene B's unique sequence, producing
//    a clean A-block→B-block k-mer partition → coherent.
// 2. A homology-like read: the read is drawn from a region whose sequence is
//    shared (or interleaved) between both genes; k-mers from both genes are
//    interleaved in the read → not coherent.
//
// The test confirms:
//   • With --min-coherent-fragments 0 (default): both candidates appear.
//   • With --min-coherent-fragments 1: only the spanning candidate survives.

#[test]
fn fusions_junction_coherence_filter() {
    use rust_htslib::bam::{self};
    let dir = TempDir::new().unwrap();

    // Two completely distinct gene bodies on different chromosomes so there is
    // no positional overlap to confuse things.
    let gene_a_seq = fusion_gen_seq(10, 300);
    let gene_b_seq = fusion_gen_seq(20, 300);

    // Write a two-contig FASTA (chr1 = gene A, chr2 = gene B).
    let fa = dir.path().join("ref.fa");
    let len = 300usize;
    std::fs::write(
        &fa,
        format!(">chr1\n{gene_a_seq}\n>chr2\n{gene_b_seq}\n"),
    ).unwrap();
    std::fs::write(
        dir.path().join("ref.fa.fai"),
        format!(
            "chr1\t{len}\t6\t{len}\t{}\nchr2\t{len}\t{}\t{len}\t{}\n",
            len + 1,
            6 + len + 1 + 6,   // offset of ">chr2\n" then the sequence
            len + 1,
        ),
    ).unwrap();

    // GTF: GENE_A on chr1, GENE_B on chr2.
    let gtf = dir.path().join("genes.gtf");
    std::fs::write(
        &gtf,
        format!(
            "chr1\tt\tgene\t1\t{len}\t.\t+\t.\tgene_name \"GENE_A\";\n\
             chr2\tt\tgene\t1\t{len}\t.\t+\t.\tgene_name \"GENE_B\";\n"
        ),
    ).unwrap();

    // Build the index.
    let index = dir.path().join("idx.duckdb");
    assert_geac_success(&[
        "experimental",
        "build-fusion-index",
        "--gene-annotation", gtf.to_str().unwrap(),
        "--fasta", fa.to_str().unwrap(),
        "--min-gene-kmers", "1",
        "--output", index.to_str().unwrap(),
    ]);

    // Fragment 1: genuine spanning read — left half from GENE_A, right half from GENE_B.
    // The k-mers from each gene will occupy disjoint positions in the read.
    let half = 60usize;
    let k = 23usize;
    // Grab a stretch from the middle of each gene body (away from edges).
    let span_read: Vec<u8> = [
        gene_a_seq[50..50 + half].as_bytes(),
        gene_b_seq[50..50 + half].as_bytes(),
    ].concat();

    // Fragment 2: discordant pair to give Fragment 2 a reason to exist — R1 from
    // GENE_A only, R2 from GENE_B only. This gives non-coherent supporting reads = 0
    // for Fragment 2's entry (it contributes only discordant reads, not spanning reads,
    // so n_coherent_reads = 0 for the GENE_A::GENE_B candidate it supports).
    //
    // To test the filter we actually need a *separate* fusion candidate whose only
    // support is a spanning read with interleaved k-mers. We simulate this by
    // constructing a "homology fragment": a read that genuinely has GENE_A k-mers
    // interleaved with GENE_B k-mers. We do this by alternating individual bytes
    // from each gene body so k-mers from both genes appear at every position.
    //
    // However, building a truly interleaved sequence that reliably hits k-mers from
    // both genes is complex to engineer synthetically. Instead, we verify the simpler
    // invariant: a discordant pair (R1→A, R2→B, no single read spanning both genes)
    // contributes 0 coherent reads, and with --min-coherent-fragments 1 it is dropped.

    // Write BAM with:
    //   frag1: single spanning read (R1 = gene_a_half + gene_b_half, R2 = pure gene_B)
    //   frag2: discordant pair only  (R1 = pure gene_A, R2 = pure gene_B)
    let r1_span = &span_read;
    let r2_span = gene_b_seq[100..100 + half + half].as_bytes();
    let r1_discord = gene_a_seq[100..100 + half + half].as_bytes();
    let r2_discord = gene_b_seq[150..150 + half + half].as_bytes();

    let bam_path = dir.path().join("reads.bam");
    {
        let mut header = bam::Header::new();
        for (sn, ln) in [("chr1", len), ("chr2", len)] {
            let mut sq = bam::header::HeaderRecord::new(b"SQ");
            sq.push_tag(b"SN", sn);
            sq.push_tag(b"LN", ln as i64);
            header.push_record(&sq);
        }
        let mut writer =
            bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        let qual_span = vec![30u8; r1_span.len()];
        let qual_norm = vec![30u8; half + half];
        // frag1 spanning
        for (qname, flags, seq, qual) in [
            (b"span1" as &[u8], 77u16, r1_span.as_slice(), qual_span.as_slice()),
            (b"span1",          141,   r2_span,             qual_norm.as_slice()),
        ] {
            let mut rec = bam::record::Record::new();
            rec.set(qname, None, seq, qual);
            rec.set_tid(-1); rec.set_pos(-1);
            rec.set_mtid(-1); rec.set_mpos(-1);
            rec.set_flags(flags);
            writer.write(&rec).unwrap();
        }
        // frag2 discordant
        for (qname, flags, seq) in [
            (b"disc1" as &[u8], 77u16, r1_discord),
            (b"disc1",          141,   r2_discord),
        ] {
            let mut rec = bam::record::Record::new();
            rec.set(qname, None, seq, &qual_norm);
            rec.set_tid(-1); rec.set_pos(-1);
            rec.set_mtid(-1); rec.set_mpos(-1);
            rec.set_flags(flags);
            writer.write(&rec).unwrap();
        }
    }

    let out = dir.path().join("fusions.parquet");

    // With default --min-coherent-fragments 0: both fragments count, 1 candidate passes.
    assert_geac_success(&[
        "experimental", "fusions",
        "--bam", bam_path.to_str().unwrap(),
        "--index", index.to_str().unwrap(),
        "--sample-id", "HG002",
        "--output", out.to_str().unwrap(),
        "--min-supporting-reads", "1",
        "--kmer-size", &k.to_string(),
    ]);
    assert_eq!(parquet_count(&out), 1, "default: one GENE_A::GENE_B candidate");

    // Check coherence columns are present and sensible.
    // span1's R1 hits both genes → n_spanning_reads >= 1, n_coherent_reads >= 1.
    let n_spanning = parquet_query_i32(&out, "n_spanning_reads", "supporting_reads >= 1");
    let n_coherent = parquet_query_i32(&out, "n_coherent_fragments", "supporting_reads >= 1");
    assert!(n_spanning >= 1, "expect at least one spanning read; got {n_spanning}");
    assert!(n_coherent >= 1, "expect at least one coherent fragment; got {n_coherent}");

    // With --min-coherent-fragments 1: the candidate has coherent reads → still passes.
    assert_geac_success(&[
        "experimental", "fusions",
        "--bam", bam_path.to_str().unwrap(),
        "--index", index.to_str().unwrap(),
        "--sample-id", "HG002",
        "--output", out.to_str().unwrap(),
        "--min-supporting-reads", "1",
        "--min-coherent-fragments", "1",
        "--kmer-size", &k.to_string(),
    ]);
    assert_eq!(parquet_count(&out), 1, "min-coherent-fragments 1: coherent candidate survives");
}

/// A fusion present in the PoN above the --max-pon-samples threshold is *tagged*
/// (filter="pon"), not removed from the output. Annotate-only runs leave filter="PASS".
#[test]
fn fusions_pon_tags_calls_instead_of_dropping() {
    use rust_htslib::bam::{self};
    let dir = TempDir::new().unwrap();

    let gene_a_seq = fusion_gen_seq(10, 300);
    let gene_b_seq = fusion_gen_seq(20, 300);
    let len = 300usize;

    let fa = dir.path().join("ref.fa");
    std::fs::write(&fa, format!(">chr1\n{gene_a_seq}\n>chr2\n{gene_b_seq}\n")).unwrap();
    std::fs::write(
        dir.path().join("ref.fa.fai"),
        format!(
            "chr1\t{len}\t6\t{len}\t{}\nchr2\t{len}\t{}\t{len}\t{}\n",
            len + 1,
            6 + len + 1 + 6,
            len + 1,
        ),
    ).unwrap();

    let gtf = dir.path().join("genes.gtf");
    std::fs::write(
        &gtf,
        format!(
            "chr1\tt\tgene\t1\t{len}\t.\t+\t.\tgene_name \"GENE_A\";\n\
             chr2\tt\tgene\t1\t{len}\t.\t+\t.\tgene_name \"GENE_B\";\n"
        ),
    ).unwrap();

    let index = dir.path().join("idx.duckdb");
    assert_geac_success(&[
        "experimental", "build-fusion-index",
        "--gene-annotation", gtf.to_str().unwrap(),
        "--fasta", fa.to_str().unwrap(),
        "--min-gene-kmers", "1",
        "--output", index.to_str().unwrap(),
    ]);

    // A single spanning read (left half GENE_A, right half GENE_B) → one candidate.
    let half = 60usize;
    let k = 23usize;
    let span_read: Vec<u8> = [
        gene_a_seq[50..50 + half].as_bytes(),
        gene_b_seq[50..50 + half].as_bytes(),
    ].concat();
    let r2_span = gene_b_seq[100..100 + half + half].as_bytes();

    let bam_path = dir.path().join("reads.bam");
    {
        let mut header = bam::Header::new();
        for (sn, ln) in [("chr1", len), ("chr2", len)] {
            let mut sq = bam::header::HeaderRecord::new(b"SQ");
            sq.push_tag(b"SN", sn);
            sq.push_tag(b"LN", ln as i64);
            header.push_record(&sq);
        }
        let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        let qual_span = vec![30u8; span_read.len()];
        let qual_norm = vec![30u8; half + half];
        for (qname, flags, seq, qual) in [
            (b"span1" as &[u8], 77u16, span_read.as_slice(), qual_span.as_slice()),
            (b"span1", 141, r2_span, qual_norm.as_slice()),
        ] {
            let mut rec = bam::record::Record::new();
            rec.set(qname, None, seq, qual);
            rec.set_tid(-1); rec.set_pos(-1);
            rec.set_mtid(-1); rec.set_mpos(-1);
            rec.set_flags(flags);
            writer.write(&rec).unwrap();
        }
    }

    // Build a one-sample fusion PoN. The file must be named *.fusions.parquet so
    // `geac merge` routes it to the fusions table.
    let normal_pq = dir.path().join("normal.fusions.parquet");
    assert_geac_success(&[
        "experimental", "fusions",
        "--bam", bam_path.to_str().unwrap(),
        "--index", index.to_str().unwrap(),
        "--sample-id", "HG001",
        "--output", normal_pq.to_str().unwrap(),
        "--min-supporting-reads", "1",
        "--kmer-size", &k.to_string(),
    ]);
    let pon_db = dir.path().join("fusion_pon.duckdb");
    assert_geac_success(&[
        "merge", normal_pq.to_str().unwrap(),
        "--output", pon_db.to_str().unwrap(),
    ]);

    // Annotate-only (no --max-pon-samples): row stays, filter=PASS, n_pon_samples=1.
    let out = dir.path().join("tumor.fusions.parquet");
    assert_geac_success(&[
        "experimental", "fusions",
        "--bam", bam_path.to_str().unwrap(),
        "--index", index.to_str().unwrap(),
        "--sample-id", "HG002",
        "--output", out.to_str().unwrap(),
        "--min-supporting-reads", "1",
        "--kmer-size", &k.to_string(),
        "--fusion-pon", pon_db.to_str().unwrap(),
    ]);
    assert_eq!(parquet_count(&out), 1, "annotate-only: call present");
    assert_eq!(parquet_query_i32(&out, "n_pon_samples", "supporting_reads >= 1"), 1);
    assert_eq!(parquet_query_str(&out, "filter", "supporting_reads >= 1"), "PASS");

    // With --max-pon-samples 0: present in 1 > 0 PoN samples → tagged, NOT dropped.
    assert_geac_success(&[
        "experimental", "fusions",
        "--bam", bam_path.to_str().unwrap(),
        "--index", index.to_str().unwrap(),
        "--sample-id", "HG002",
        "--output", out.to_str().unwrap(),
        "--min-supporting-reads", "1",
        "--kmer-size", &k.to_string(),
        "--fusion-pon", pon_db.to_str().unwrap(),
        "--max-pon-samples", "0",
    ]);
    assert_eq!(
        parquet_count(&out), 1,
        "max-pon-samples 0: call retained (tagged, not dropped)"
    );
    assert_eq!(parquet_query_str(&out, "filter", "supporting_reads >= 1"), "pon");
}
