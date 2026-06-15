use std::collections::{BTreeMap, BTreeSet};

use anyhow::{bail, Context, Result};
use duckdb::Connection;

use crate::cli::InspectArgs;

const REQUIRED_ALT_BASE_COLUMNS: &[&str] = &[
    "sample_id",
    "chrom",
    "pos",
    "ref_allele",
    "alt_allele",
    "variant_type",
    "total_depth",
    "alt_count",
    "ref_count",
    "fwd_depth",
    "rev_depth",
    "fwd_alt_count",
    "rev_alt_count",
    "fwd_ref_count",
    "rev_ref_count",
    "overlap_depth",
    "overlap_alt_agree",
    "overlap_alt_disagree",
    "overlap_ref_agree",
    "read_type",
    "pipeline",
];

const OPTIONAL_TABLES: &[&str] = &[
    "alt_reads",
    "normal_evidence",
    "pon_evidence",
    "coverage",
    "coverage_intervals",
    "sample_metrics",
    "fragments",
    "fusions",
];

const RESOURCE_COLUMNS: &[&str] = &[
    "bam_path",
    "bai_path",
    "variants_path",
    "gnomad_path",
    "gnomad_index_path",
    "targets_path",
];

#[derive(Debug, Default)]
struct InspectReport {
    errors: Vec<String>,
    warnings: Vec<String>,
    info: Vec<String>,
}

impl InspectReport {
    fn error(&mut self, msg: impl Into<String>) {
        self.errors.push(msg.into());
    }

    fn warn(&mut self, msg: impl Into<String>) {
        self.warnings.push(msg.into());
    }

    fn note(&mut self, msg: impl Into<String>) {
        self.info.push(msg.into());
    }
}

pub fn run_inspect(args: &InspectArgs) -> Result<()> {
    let conn = Connection::open(&args.input)
        .with_context(|| format!("failed to open DuckDB: {}", args.input.display()))?;

    let mut report = InspectReport::default();
    report.note(format!("input: {}", args.input.display()));

    let tables = table_names(&conn)?;
    if tables.is_empty() {
        report.error("database contains no user tables");
    }

    inspect_metadata(&conn, &tables, &mut report)?;
    inspect_tables(&conn, &tables, &mut report)?;
    inspect_alt_bases(&conn, &tables, &mut report)?;
    inspect_samples(&conn, &tables, &mut report)?;

    print_report(&report);

    if !report.errors.is_empty() {
        bail!("inspect found {} error(s)", report.errors.len());
    }
    if args.strict && !report.warnings.is_empty() {
        bail!(
            "inspect found {} warning(s) in --strict mode",
            report.warnings.len()
        );
    }
    Ok(())
}

fn table_names(conn: &Connection) -> Result<BTreeSet<String>> {
    let mut stmt = conn
        .prepare(
            "SELECT table_name FROM information_schema.tables \
             WHERE table_schema = 'main' ORDER BY table_name",
        )
        .context("failed to inspect DuckDB tables")?;
    let rows = stmt
        .query_map([], |row| row.get::<_, String>(0))
        .context("failed to query DuckDB tables")?;
    Ok(rows.filter_map(|r| r.ok()).collect())
}

fn table_columns(conn: &Connection, table: &str) -> Result<BTreeSet<String>> {
    let mut stmt = conn
        .prepare(&format!("DESCRIBE SELECT * FROM {table} LIMIT 0"))
        .with_context(|| format!("failed to describe table {table}"))?;
    let rows = stmt
        .query_map([], |row| row.get::<_, String>(0))
        .with_context(|| format!("failed to query columns for {table}"))?;
    Ok(rows.filter_map(|r| r.ok()).collect())
}

fn table_count(conn: &Connection, table: &str) -> Result<i64> {
    conn.query_row(&format!("SELECT COUNT(*) FROM {table}"), [], |row| {
        row.get(0)
    })
    .with_context(|| format!("failed to count rows in {table}"))
}

fn inspect_metadata(
    conn: &Connection,
    tables: &BTreeSet<String>,
    report: &mut InspectReport,
) -> Result<()> {
    if !tables.contains("geac_metadata") {
        report.warn("missing geac_metadata provenance table");
    } else {
        let cols = table_columns(conn, "geac_metadata")?;
        for required in ["geac_version", "schema_version", "created_at_utc"] {
            if !cols.contains(required) {
                report.warn(format!("geac_metadata missing column {required}"));
            }
        }
        let n = table_count(conn, "geac_metadata")?;
        if n != 1 {
            report.warn(format!("geac_metadata has {n} rows; expected 1"));
        }
    }

    if !tables.contains("geac_inputs") {
        report.warn("missing geac_inputs provenance table");
    } else {
        let n = table_count(conn, "geac_inputs")?;
        report.note(format!("geac_inputs rows: {n}"));
    }
    Ok(())
}

fn inspect_tables(
    conn: &Connection,
    tables: &BTreeSet<String>,
    report: &mut InspectReport,
) -> Result<()> {
    if !tables.contains("alt_bases") {
        report.error("missing required alt_bases table");
        return Ok(());
    }

    let mut present_optional = Vec::new();
    for table in OPTIONAL_TABLES {
        if tables.contains(*table) {
            let n = table_count(conn, table)?;
            present_optional.push(format!("{table} ({n})"));
        }
    }
    if present_optional.is_empty() {
        report.note("optional feature tables: none");
    } else {
        report.note(format!(
            "optional feature tables: {}",
            present_optional.join(", ")
        ));
    }
    Ok(())
}

fn inspect_alt_bases(
    conn: &Connection,
    tables: &BTreeSet<String>,
    report: &mut InspectReport,
) -> Result<()> {
    if !tables.contains("alt_bases") {
        return Ok(());
    }

    let cols = table_columns(conn, "alt_bases")?;
    let missing: Vec<&str> = REQUIRED_ALT_BASE_COLUMNS
        .iter()
        .copied()
        .filter(|c| !cols.contains(*c))
        .collect();
    if missing.is_empty() {
        report.note("alt_bases required columns: ok");
    } else {
        report.error(format!(
            "alt_bases missing required columns: {}",
            missing.join(", ")
        ));
    }

    let n_rows = table_count(conn, "alt_bases")?;
    let n_samples: i64 = conn
        .query_row(
            "SELECT COUNT(DISTINCT sample_id) FROM alt_bases",
            [],
            |row| row.get(0),
        )
        .context("failed to count distinct alt_bases samples")?;
    report.note(format!("alt_bases rows: {n_rows}; samples: {n_samples}"));

    if n_rows == 0 {
        report.warn("alt_bases is empty");
    }

    for col in ["sample_id", "chrom", "alt_allele"] {
        if cols.contains(col) {
            let n_null: i64 = conn
                .query_row(
                    &format!("SELECT COUNT(*) FROM alt_bases WHERE {col} IS NULL"),
                    [],
                    |row| row.get(0),
                )
                .with_context(|| format!("failed to count null {col} values"))?;
            if n_null > 0 {
                report.warn(format!("alt_bases.{col} has {n_null} NULL values"));
            }
        }
    }

    Ok(())
}

fn inspect_samples(
    conn: &Connection,
    tables: &BTreeSet<String>,
    report: &mut InspectReport,
) -> Result<()> {
    if !tables.contains("samples") {
        report.warn("missing samples summary table; rerun geac merge with a current version");
        return Ok(());
    }

    let cols = table_columns(conn, "samples")?;
    let n = table_count(conn, "samples")?;
    report.note(format!("samples rows: {n}"));

    inspect_sample_metadata_conflicts(conn, &cols, report)?;
    inspect_resource_paths(conn, &cols, report)?;
    Ok(())
}

fn inspect_sample_metadata_conflicts(
    conn: &Connection,
    cols: &BTreeSet<String>,
    report: &mut InspectReport,
) -> Result<()> {
    let mut available = Vec::new();
    for col in ["subject_id", "sample_type", "read_type", "pipeline"] {
        if cols.contains(col) {
            available.push(col);
        }
    }
    if available.is_empty() {
        return Ok(());
    }

    let select_parts: Vec<String> = available
        .iter()
        .map(|col| format!("COUNT(DISTINCT {col}) FILTER (WHERE {col} IS NOT NULL) AS n_{col}"))
        .collect();
    let sql = format!(
        "SELECT sample_id, {} FROM samples GROUP BY sample_id",
        select_parts.join(", ")
    );
    let mut stmt = conn
        .prepare(&sql)
        .context("failed to prepare sample metadata conflict query")?;
    let mut rows = stmt
        .query([])
        .context("failed to query sample metadata conflicts")?;

    let mut conflicts: BTreeMap<&str, i64> = BTreeMap::new();
    while let Some(row) = rows.next().context("failed to read conflict row")? {
        for (idx, col) in available.iter().enumerate() {
            let n_distinct: i64 = row.get(idx + 1)?;
            if n_distinct > 1 {
                *conflicts.entry(*col).or_default() += 1;
            }
        }
    }

    for (col, n_samples) in conflicts {
        report.warn(format!(
            "{n_samples} sample_id value(s) have conflicting {col} values"
        ));
    }
    Ok(())
}

fn inspect_resource_paths(
    conn: &Connection,
    cols: &BTreeSet<String>,
    report: &mut InspectReport,
) -> Result<()> {
    for col in RESOURCE_COLUMNS {
        if !cols.contains(*col) {
            report.warn(format!("samples missing resource column {col}"));
            continue;
        }

        let n_non_null: i64 = conn
            .query_row(
                &format!("SELECT COUNT(*) FROM samples WHERE {col} IS NOT NULL AND {col} <> ''"),
                [],
                |row| row.get(0),
            )
            .with_context(|| format!("failed to count samples.{col}"))?;
        if *col == "bam_path" && n_non_null == 0 {
            report.warn("no embedded BAM/CRAM paths found; IGV requires an external manifest");
        }
        if n_non_null == 0 {
            continue;
        }

        let n_localized: i64 = conn
            .query_row(
                &format!(
                    "SELECT COUNT(*) FROM samples \
                     WHERE {col} IS NOT NULL AND {col} <> '' \
                       AND {col} NOT LIKE 'gs://%' \
                       AND {col} NOT LIKE 's3://%' \
                       AND {col} NOT LIKE 'http://%' \
                       AND {col} NOT LIKE 'https://%'"
                ),
                [],
                |row| row.get(0),
            )
            .with_context(|| format!("failed to classify samples.{col} paths"))?;
        if n_localized > 0 {
            report.warn(format!(
                "samples.{col} has {n_localized} non-URI path(s); cloud/IGV use may need explicit URI inputs or a manifest"
            ));
        }
    }
    Ok(())
}

fn print_report(report: &InspectReport) {
    println!("GEAC inspect");
    println!("============");

    if !report.info.is_empty() {
        println!("\nInfo:");
        for msg in &report.info {
            println!("  - {msg}");
        }
    }

    if !report.warnings.is_empty() {
        println!("\nWarnings:");
        for msg in &report.warnings {
            println!("  - {msg}");
        }
    }

    if !report.errors.is_empty() {
        println!("\nErrors:");
        for msg in &report.errors {
            println!("  - {msg}");
        }
    }

    println!(
        "\nSummary: {} error(s), {} warning(s)",
        report.errors.len(),
        report.warnings.len()
    );
}
