use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;

use anyhow::{bail, Context, Result};
use duckdb::Connection;
use tracing::info;

use crate::cli::MergeArgs;

const DUCKDB_SCHEMA_VERSION: &str = "duckdb-v2";

struct TableSpec {
    table: &'static str,
    suffix: Option<&'static str>,
    index_sql: Option<&'static str>,
    rebuild_samples_summary: bool,
}

struct InputProvenance {
    input_path: String,
    input_kind: &'static str,
    source_kind: &'static str,
    file_size_bytes: Option<i64>,
    modified_at_epoch_seconds: Option<f64>,
    sample_count: Option<i64>,
    row_count: Option<i64>,
}

const TABLE_SPECS: &[TableSpec] = &[
    TableSpec {
        table: "alt_bases",
        suffix: None,
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_chrom_pos ON alt_bases (chrom, pos);
             CREATE INDEX IF NOT EXISTS idx_sample_id ON alt_bases (sample_id);",
        ),
        rebuild_samples_summary: true,
    },
    TableSpec {
        table: "alt_reads",
        suffix: Some(".reads.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_reads_locus \
             ON alt_reads (sample_id, chrom, pos, alt_allele);",
        ),
        rebuild_samples_summary: false,
    },
    TableSpec {
        table: "normal_evidence",
        suffix: Some(".normal_evidence.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_ne_locus \
             ON normal_evidence (tumor_sample_id, chrom, pos, tumor_alt_allele);",
        ),
        rebuild_samples_summary: false,
    },
    TableSpec {
        table: "pon_evidence",
        suffix: Some(".pon_evidence.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_pe_locus \
             ON pon_evidence (tumor_sample_id, chrom, pos, tumor_alt_allele);",
        ),
        rebuild_samples_summary: false,
    },
    TableSpec {
        table: "coverage",
        suffix: Some(".coverage.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_coverage_locus \
             ON coverage (sample_id, chrom, pos);",
        ),
        rebuild_samples_summary: false,
    },
    TableSpec {
        table: "coverage_intervals",
        suffix: Some(".coverage.intervals.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_coverage_intervals_locus \
             ON coverage_intervals (sample_id, chrom, start, \"end\");",
        ),
        rebuild_samples_summary: false,
    },
];

fn escape_path(path: &Path) -> String {
    escape_sql_literal(&path.display().to_string())
}

fn escape_sql_literal(value: &str) -> String {
    value.replace('\'', "''")
}

fn shell_quote_arg(arg: &str) -> String {
    if !arg.is_empty()
        && arg
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || "-_./:=,".contains(c))
    {
        arg.to_string()
    } else {
        format!("'{}'", arg.replace('\'', "'\"'\"'"))
    }
}

fn command_line_string() -> String {
    std::env::args_os()
        .map(|arg| shell_quote_arg(&arg.to_string_lossy()))
        .collect::<Vec<_>>()
        .join(" ")
}

fn file_size_bytes(path: &Path) -> Option<i64> {
    let len = fs::metadata(path).ok()?.len();
    i64::try_from(len).ok()
}

fn modified_at_epoch_seconds(path: &Path) -> Option<f64> {
    let modified = fs::metadata(path).ok()?.modified().ok()?;
    let duration = modified.duration_since(UNIX_EPOCH).ok()?;
    Some(duration.as_secs_f64())
}

fn sample_id_column(table: &str) -> Option<&'static str> {
    match table {
        "alt_bases" | "alt_reads" | "coverage" | "coverage_intervals" => Some("sample_id"),
        "normal_evidence" | "pon_evidence" => Some("tumor_sample_id"),
        _ => None,
    }
}

fn dst_table_row_count(conn: &Connection, table: &str) -> Result<i64> {
    if !dst_table_exists(conn, table)? {
        return Ok(0);
    }

    conn.query_row(&format!("SELECT COUNT(*) FROM {table}"), [], |row| {
        row.get(0)
    })
    .with_context(|| format!("failed to count rows in {table}"))
}

fn merged_sample_count(conn: &Connection) -> Result<i64> {
    if dst_table_exists(conn, "samples")? {
        return conn
            .query_row("SELECT COUNT(DISTINCT sample_id) FROM samples", [], |row| {
                row.get(0)
            })
            .context("failed to count samples from samples table");
    }

    if dst_table_exists(conn, "alt_bases")? {
        return conn
            .query_row(
                "SELECT COUNT(DISTINCT sample_id) FROM alt_bases",
                [],
                |row| row.get(0),
            )
            .context("failed to count samples from alt_bases");
    }

    if dst_table_exists(conn, "coverage")? {
        return conn
            .query_row(
                "SELECT COUNT(DISTINCT sample_id) FROM coverage",
                [],
                |row| row.get(0),
            )
            .context("failed to count samples from coverage");
    }

    Ok(0)
}

/// Returns true if `table` exists in the attached database aliased as `alias`.
fn src_table_exists(conn: &Connection, alias: &str, table: &str) -> Result<bool> {
    let count: i64 = conn
        .query_row(
            &format!(
                "SELECT COUNT(*) FROM duckdb_tables() \
                 WHERE database_name = '{alias}' AND schema_name = 'main' \
                 AND table_name = '{table}'"
            ),
            [],
            |row| row.get(0),
        )
        .with_context(|| format!("failed to query table list in attached database '{alias}'"))?;
    Ok(count > 0)
}

/// Returns true if `table` exists in the output (current) database.
fn dst_table_exists(conn: &Connection, table: &str) -> Result<bool> {
    let count: i64 = conn
        .query_row(
            &format!(
                "SELECT COUNT(*) FROM duckdb_tables() \
                 WHERE database_name = current_database() AND schema_name = 'main' \
                 AND table_name = '{table}'"
            ),
            [],
            |row| row.get(0),
        )
        .with_context(|| format!("failed to query output table list for '{table}'"))?;
    Ok(count > 0)
}

/// Copy a table from an attached source database into the output database.
/// Creates the table if it doesn't exist yet; otherwise inserts rows.
fn copy_table(conn: &Connection, alias: &str, table: &str) -> Result<()> {
    if dst_table_exists(conn, table)? {
        conn.execute_batch(&format!(
            "INSERT INTO {table} BY NAME SELECT * FROM {alias}.{table};"
        ))
        .with_context(|| format!("failed to insert '{table}' from attached '{alias}'"))?;
    } else {
        conn.execute_batch(&format!(
            "CREATE TABLE {table} AS SELECT * FROM {alias}.{table};"
        ))
        .with_context(|| format!("failed to create '{table}' from attached '{alias}'"))?;
    }
    Ok(())
}

fn parquet_spec_for_path(path: &Path) -> &'static TableSpec {
    let file_name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default();
    TABLE_SPECS
        .iter()
        .find(|spec| {
            spec.suffix
                .is_some_and(|suffix| file_name.ends_with(suffix))
        })
        .unwrap_or(&TABLE_SPECS[0])
}

fn classify_inputs(inputs: &[PathBuf]) -> (Vec<&PathBuf>, BTreeMap<&'static str, Vec<&PathBuf>>) {
    let duckdb_inputs = inputs
        .iter()
        .filter(|path| path.extension().and_then(|e| e.to_str()) == Some("duckdb"))
        .collect::<Vec<_>>();

    let mut parquet_groups: BTreeMap<&'static str, Vec<&PathBuf>> = TABLE_SPECS
        .iter()
        .map(|spec| (spec.table, Vec::new()))
        .collect();

    for input in inputs
        .iter()
        .filter(|path| path.extension().and_then(|e| e.to_str()) != Some("duckdb"))
    {
        let spec = parquet_spec_for_path(input);
        parquet_groups.get_mut(spec.table).unwrap().push(input);
    }

    (duckdb_inputs, parquet_groups)
}

fn collect_parquet_input_provenance(
    conn: &Connection,
    spec: &TableSpec,
    path: &Path,
) -> Result<InputProvenance> {
    let escaped = escape_path(path);
    let row_count = conn
        .query_row(
            &format!("SELECT COUNT(*) FROM read_parquet('{escaped}')"),
            [],
            |row| row.get(0),
        )
        .with_context(|| format!("failed to count rows in {}", path.display()))?;

    let sample_count = if let Some(column) = sample_id_column(spec.table) {
        Some(
            conn.query_row(
                &format!("SELECT COUNT(DISTINCT {column}) FROM read_parquet('{escaped}')"),
                [],
                |row| row.get(0),
            )
            .with_context(|| {
                format!(
                    "failed to count distinct samples in {} from column {column}",
                    path.display()
                )
            })?,
        )
    } else {
        None
    };

    Ok(InputProvenance {
        input_path: path.display().to_string(),
        input_kind: spec.table,
        source_kind: "parquet",
        file_size_bytes: file_size_bytes(path),
        modified_at_epoch_seconds: modified_at_epoch_seconds(path),
        sample_count,
        row_count: Some(row_count),
    })
}

fn merge_parquet_group(
    conn: &Connection,
    spec: &TableSpec,
    inputs: &[&PathBuf],
    input_provenance: &mut Vec<InputProvenance>,
) -> Result<()> {
    if inputs.is_empty() {
        return Ok(());
    }

    info!(
        table = spec.table,
        n_files = inputs.len(),
        "reading Parquet files"
    );

    let first_escaped = escape_path(inputs[0]);
    conn.execute_batch(&format!(
        "CREATE TABLE {} AS SELECT * FROM read_parquet('{first_escaped}');",
        spec.table
    ))
    .with_context(|| format!("failed to create {} table from first file", spec.table))?;

    for (idx, path) in inputs[1..].iter().enumerate() {
        let escaped = escape_path(path);
        info!(
            table = spec.table,
            file = %path.display(),
            idx = idx + 1,
            total = inputs.len(),
            "inserting Parquet file"
        );
        conn.execute_batch(&format!(
            "INSERT INTO {} BY NAME SELECT * FROM read_parquet('{escaped}');",
            spec.table
        ))
        .with_context(|| format!("failed to insert {} from {}", spec.table, path.display()))?;
    }

    let n_rows: i64 = conn
        .query_row(&format!("SELECT COUNT(*) FROM {}", spec.table), [], |row| {
            row.get(0)
        })
        .with_context(|| format!("failed to count rows in {}", spec.table))?;
    info!(table = spec.table, n_rows, "table created from Parquet");

    for path in inputs {
        input_provenance.push(collect_parquet_input_provenance(conn, spec, path)?);
    }

    Ok(())
}

fn merge_duckdb_inputs(
    conn: &Connection,
    duckdb_inputs: &[&PathBuf],
    input_provenance: &mut Vec<InputProvenance>,
) -> Result<()> {
    for (i, db_path) in duckdb_inputs.iter().enumerate() {
        let alias = format!("_src{i}");
        let escaped = escape_path(db_path);
        info!(file = %db_path.display(), "attaching source DuckDB");

        conn.execute_batch(&format!("ATTACH '{escaped}' AS {alias} (READ_ONLY);"))
            .with_context(|| format!("failed to attach {}", db_path.display()))?;

        let mut total_rows = 0_i64;
        for spec in TABLE_SPECS {
            if src_table_exists(conn, &alias, spec.table)? {
                let n_src: i64 = conn
                    .query_row(
                        &format!("SELECT COUNT(*) FROM {alias}.{}", spec.table),
                        [],
                        |row| row.get(0),
                    )
                    .with_context(|| format!("failed to count rows in {alias}.{}", spec.table))?;

                info!(
                    table = spec.table,
                    n_src,
                    file = %db_path.display(),
                    "copying table from source DuckDB"
                );
                total_rows += n_src;
                copy_table(conn, &alias, spec.table)?;
            }
        }

        let sample_count = if src_table_exists(conn, &alias, "samples")? {
            Some(
                conn.query_row(
                    &format!("SELECT COUNT(DISTINCT sample_id) FROM {alias}.samples"),
                    [],
                    |row| row.get(0),
                )
                .with_context(|| {
                    format!(
                        "failed to count samples in attached database {}",
                        db_path.display()
                    )
                })?,
            )
        } else if src_table_exists(conn, &alias, "alt_bases")? {
            Some(
                conn.query_row(
                    &format!("SELECT COUNT(DISTINCT sample_id) FROM {alias}.alt_bases"),
                    [],
                    |row| row.get(0),
                )
                .with_context(|| {
                    format!(
                        "failed to count samples in attached database {}",
                        db_path.display()
                    )
                })?,
            )
        } else if src_table_exists(conn, &alias, "coverage")? {
            Some(
                conn.query_row(
                    &format!("SELECT COUNT(DISTINCT sample_id) FROM {alias}.coverage"),
                    [],
                    |row| row.get(0),
                )
                .with_context(|| {
                    format!(
                        "failed to count samples in attached database {}",
                        db_path.display()
                    )
                })?,
            )
        } else {
            None
        };

        input_provenance.push(InputProvenance {
            input_path: db_path.display().to_string(),
            input_kind: "cohort",
            source_kind: "duckdb",
            file_size_bytes: file_size_bytes(db_path),
            modified_at_epoch_seconds: modified_at_epoch_seconds(db_path),
            sample_count,
            row_count: Some(total_rows),
        });

        conn.execute_batch(&format!("DETACH {alias};"))
            .with_context(|| format!("failed to detach {}", db_path.display()))?;
    }

    Ok(())
}

fn rebuild_samples_summary(conn: &Connection) -> Result<()> {
    info!("(re)building samples summary table...");
    conn.execute_batch(
        "DROP TABLE IF EXISTS samples;
         CREATE TABLE samples AS
         SELECT
             sample_id,
             batch,
             COUNT(*)                                   AS n_alt_loci,
             SUM(alt_count)                             AS total_alt_reads,
             COUNT(DISTINCT chrom || ':' || pos)        AS n_positions,
             MIN(pos)                                   AS min_pos,
             MAX(pos)                                   AS max_pos,
             ANY_VALUE(read_type)                       AS read_type,
             ANY_VALUE(pipeline)                        AS pipeline
         FROM alt_bases
         GROUP BY sample_id, batch
         ORDER BY sample_id, batch;",
    )
    .context("failed to build samples summary table")?;

    let n_samples: i64 = conn
        .query_row("SELECT COUNT(*) FROM samples", [], |row| row.get(0))
        .context("failed to count samples")?;
    info!(n_samples, "samples summary table built");
    Ok(())
}

fn build_indices_and_derived_tables(conn: &Connection) -> Result<()> {
    for spec in TABLE_SPECS {
        if dst_table_exists(conn, spec.table)? {
            if let Some(index_sql) = spec.index_sql {
                info!(table = spec.table, "creating indices");
                conn.execute_batch(index_sql)
                    .with_context(|| format!("failed to create indices on {}", spec.table))?;
            }

            if spec.rebuild_samples_summary {
                rebuild_samples_summary(conn)?;
            }
        }
    }

    Ok(())
}

fn write_input_provenance(conn: &Connection, inputs: &[InputProvenance]) -> Result<()> {
    info!("writing geac_inputs table...");
    conn.execute_batch(
        "DROP TABLE IF EXISTS geac_inputs;
         CREATE TABLE geac_inputs (
             input_path VARCHAR,
             input_kind VARCHAR,
             source_kind VARCHAR,
             file_size_bytes BIGINT,
             modified_at TIMESTAMPTZ,
             checksum_sha256 VARCHAR,
             sample_count BIGINT,
             row_count BIGINT
         );",
    )
    .context("failed to create geac_inputs table")?;

    for input in inputs {
        let modified_at_sql = input
            .modified_at_epoch_seconds
            .map(|seconds| format!("to_timestamp({seconds})"))
            .unwrap_or_else(|| "NULL".to_string());
        let file_size_sql = input
            .file_size_bytes
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NULL".to_string());
        let sample_count_sql = input
            .sample_count
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NULL".to_string());
        let row_count_sql = input
            .row_count
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NULL".to_string());

        conn.execute_batch(&format!(
            "INSERT INTO geac_inputs (
                 input_path, input_kind, source_kind, file_size_bytes,
                 modified_at, checksum_sha256, sample_count, row_count
             ) VALUES (
                 '{}', '{}', '{}', {file_size_sql},
                 {modified_at_sql}, NULL, {sample_count_sql}, {row_count_sql}
             );",
            escape_sql_literal(&input.input_path),
            input.input_kind,
            input.source_kind,
        ))
        .with_context(|| format!("failed to record provenance for {}", input.input_path))?;
    }

    Ok(())
}

fn write_metadata(
    conn: &Connection,
    args: &MergeArgs,
    parquet_groups: &BTreeMap<&'static str, Vec<&PathBuf>>,
    duckdb_inputs: &[&PathBuf],
) -> Result<()> {
    info!("writing geac_metadata table...");

    let n_samples = merged_sample_count(conn)?;
    let alt_bases_rows = dst_table_row_count(conn, "alt_bases")?;
    let alt_reads_rows = dst_table_row_count(conn, "alt_reads")?;
    let normal_evidence_rows = dst_table_row_count(conn, "normal_evidence")?;
    let pon_evidence_rows = dst_table_row_count(conn, "pon_evidence")?;
    let coverage_rows = dst_table_row_count(conn, "coverage")?;
    let samples_rows = dst_table_row_count(conn, "samples")?;

    conn.execute_batch(&format!(
        "DROP TABLE IF EXISTS geac_metadata;
         CREATE TABLE geac_metadata (
             schema_version VARCHAR,
             geac_version VARCHAR,
             created_at TIMESTAMPTZ,
             command_line VARCHAR,
             output_path VARCHAR,
             platform_os VARCHAR,
             platform_arch VARCHAR,
             platform_family VARCHAR,
             n_alt_bases_inputs BIGINT,
             n_alt_reads_inputs BIGINT,
             n_normal_evidence_inputs BIGINT,
             n_pon_evidence_inputs BIGINT,
             n_coverage_inputs BIGINT,
             n_duckdb_inputs BIGINT,
             n_samples BIGINT,
             alt_bases_rows BIGINT,
             alt_reads_rows BIGINT,
             normal_evidence_rows BIGINT,
             pon_evidence_rows BIGINT,
             coverage_rows BIGINT,
             samples_rows BIGINT
         );
         INSERT INTO geac_metadata VALUES (
             '{}',
             '{}',
             current_timestamp,
             '{}',
             '{}',
             '{}',
             '{}',
             '{}',
             {},
             {},
             {},
             {},
             {},
             {},
             {},
             {},
             {},
             {},
             {},
             {},
             {}
         );",
        DUCKDB_SCHEMA_VERSION,
        env!("CARGO_PKG_VERSION"),
        escape_sql_literal(&command_line_string()),
        escape_sql_literal(&args.output.display().to_string()),
        std::env::consts::OS,
        std::env::consts::ARCH,
        std::env::consts::FAMILY,
        parquet_groups.get("alt_bases").map_or(0, Vec::len),
        parquet_groups.get("alt_reads").map_or(0, Vec::len),
        parquet_groups.get("normal_evidence").map_or(0, Vec::len),
        parquet_groups.get("pon_evidence").map_or(0, Vec::len),
        parquet_groups.get("coverage").map_or(0, Vec::len),
        duckdb_inputs.len(),
        n_samples,
        alt_bases_rows,
        alt_reads_rows,
        normal_evidence_rows,
        pon_evidence_rows,
        coverage_rows,
        samples_rows,
    ))
    .context("failed to write geac_metadata table")?;
    Ok(())
}

pub fn merge(args: &MergeArgs) -> Result<()> {
    if args.output.exists() {
        bail!(
            "output file already exists: {}. Remove it or choose a different path.",
            args.output.display()
        );
    }

    if args.inputs.is_empty() {
        bail!("no input files provided");
    }

    for path in &args.inputs {
        if !path.exists() {
            bail!("input file not found: {}", path.display());
        }
    }

    let (duckdb_inputs, parquet_groups) = classify_inputs(&args.inputs);

    info!(
        n_locus_files = parquet_groups.get("alt_bases").map_or(0, Vec::len),
        n_reads_files = parquet_groups.get("alt_reads").map_or(0, Vec::len),
        n_normal_files = parquet_groups.get("normal_evidence").map_or(0, Vec::len),
        n_pon_files = parquet_groups.get("pon_evidence").map_or(0, Vec::len),
        n_coverage_files = parquet_groups.get("coverage").map_or(0, Vec::len),
        n_duckdb_files = duckdb_inputs.len(),
        output = %args.output.display(),
        "merging files into DuckDB"
    );

    let conn = Connection::open(&args.output).with_context(|| {
        format!(
            "failed to create DuckDB database: {}",
            args.output.display()
        )
    })?;

    let mut input_provenance = Vec::new();

    for spec in TABLE_SPECS {
        let inputs = parquet_groups
            .get(spec.table)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        merge_parquet_group(&conn, spec, inputs, &mut input_provenance)?;
    }

    merge_duckdb_inputs(&conn, &duckdb_inputs, &mut input_provenance)?;
    build_indices_and_derived_tables(&conn)?;
    write_input_provenance(&conn, &input_provenance)?;
    write_metadata(&conn, args, &parquet_groups, &duckdb_inputs)?;

    info!(output = %args.output.display(), "merge complete");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classify_inputs_routes_by_suffix() {
        let inputs = vec![
            PathBuf::from("sample.parquet"),
            PathBuf::from("sample.reads.parquet"),
            PathBuf::from("sample.normal_evidence.parquet"),
            PathBuf::from("sample.pon_evidence.parquet"),
            PathBuf::from("sample.coverage.parquet"),
            PathBuf::from("cohort.duckdb"),
        ];

        let (duckdb_inputs, parquet_groups) = classify_inputs(&inputs);

        let duckdb_names = duckdb_inputs
            .iter()
            .map(|path| path.display().to_string())
            .collect::<Vec<_>>();
        assert_eq!(duckdb_names, vec!["cohort.duckdb"]);

        let group_names = |table: &str| {
            parquet_groups[table]
                .iter()
                .map(|path| path.display().to_string())
                .collect::<Vec<_>>()
        };
        assert_eq!(group_names("alt_bases"), vec!["sample.parquet"]);
        assert_eq!(group_names("alt_reads"), vec!["sample.reads.parquet"]);
        assert_eq!(
            group_names("normal_evidence"),
            vec!["sample.normal_evidence.parquet"]
        );
        assert_eq!(
            group_names("pon_evidence"),
            vec!["sample.pon_evidence.parquet"]
        );
        assert_eq!(group_names("coverage"), vec!["sample.coverage.parquet"]);
    }

    #[test]
    fn table_registry_builds_indices_and_samples_table() {
        let conn = Connection::open_in_memory().expect("in-memory duckdb");
        conn.execute_batch(
            "CREATE TABLE alt_bases (
                 sample_id VARCHAR,
                 batch VARCHAR,
                 chrom VARCHAR,
                 pos BIGINT,
                 alt_count INTEGER,
                 read_type VARCHAR,
                 pipeline VARCHAR
             );
             INSERT INTO alt_bases VALUES ('s1', 'b1', 'chr1', 10, 2, 'raw', 'raw');
             CREATE TABLE alt_reads (
                 sample_id VARCHAR,
                 chrom VARCHAR,
                 pos BIGINT,
                 alt_allele VARCHAR
             );
             INSERT INTO alt_reads VALUES ('s1', 'chr1', 10, 'T');",
        )
        .expect("seed tables");

        build_indices_and_derived_tables(&conn).expect("build derived tables");

        assert!(dst_table_exists(&conn, "samples").expect("samples exists"));

        let n_samples: i64 = conn
            .query_row("SELECT COUNT(*) FROM samples", [], |row| row.get(0))
            .expect("count samples");
        assert_eq!(n_samples, 1);

        let n_indices: i64 = conn
            .query_row(
                "SELECT COUNT(*) FROM duckdb_indexes() WHERE table_name IN ('alt_bases', 'alt_reads')",
                [],
                |row| row.get(0),
            )
            .expect("count indices");
        assert!(
            n_indices >= 2,
            "expected registry-driven indices to be created"
        );
    }
}
