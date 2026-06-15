use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;

use anyhow::{bail, Context, Result};
use duckdb::{Connection, OptionalExt};
use tracing::{info, warn};

use crate::cli::MergeArgs;

const DUCKDB_SCHEMA_VERSION: &str = "duckdb-v4";

struct TableSpec {
    table: &'static str,
    suffix: Option<&'static str>,
    index_sql: Option<&'static str>,
    rebuild_samples_summary: bool,
    /// When true, files are registered as an external DuckDB VIEW over read_parquet()
    /// rather than ingested into a table. Use for large outputs where copying into
    /// DuckDB would be prohibitively expensive (e.g. WGS-scale fragment records).
    view_only: bool,
}

struct InputProvenance {
    input_path: String,
    input_kind: &'static str,
    source_kind: &'static str,
    file_size_bytes: Option<i64>,
    modified_at_epoch_seconds: Option<f64>,
    checksum_sha256: Option<String>,
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
        view_only: false,
    },
    TableSpec {
        table: "alt_reads",
        suffix: Some(".reads.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_reads_locus \
             ON alt_reads (sample_id, chrom, pos, alt_allele);",
        ),
        rebuild_samples_summary: false,
        view_only: false,
    },
    TableSpec {
        table: "normal_evidence",
        suffix: Some(".normal_evidence.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_ne_locus \
             ON normal_evidence (tumor_sample_id, chrom, pos, tumor_alt_allele);",
        ),
        rebuild_samples_summary: false,
        view_only: false,
    },
    TableSpec {
        table: "pon_evidence",
        suffix: Some(".pon_evidence.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_pe_locus \
             ON pon_evidence (tumor_sample_id, chrom, pos, tumor_alt_allele);",
        ),
        rebuild_samples_summary: false,
        view_only: false,
    },
    TableSpec {
        table: "coverage",
        suffix: Some(".coverage.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_coverage_locus \
             ON coverage (sample_id, chrom, pos);",
        ),
        rebuild_samples_summary: false,
        view_only: false,
    },
    TableSpec {
        table: "coverage_intervals",
        suffix: Some(".coverage.intervals.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_coverage_intervals_locus \
             ON coverage_intervals (sample_id, chrom, start, \"end\");",
        ),
        rebuild_samples_summary: false,
        view_only: false,
    },
    TableSpec {
        table: "sample_metrics",
        suffix: Some(".sample_metrics.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_sample_metrics_sample \
             ON sample_metrics (sample_id);",
        ),
        rebuild_samples_summary: false,
        view_only: false,
    },
    TableSpec {
        table: "fragments",
        suffix: Some(".fragments.parquet"),
        index_sql: Some(
            "CREATE INDEX IF NOT EXISTS idx_fragments_locus \
             ON fragments (sample_id, chrom, midpoint);",
        ),
        rebuild_samples_summary: false,
        view_only: false,
    },
    TableSpec {
        table: "fusions",
        suffix: Some(".fusions.parquet"),
        index_sql: Some("CREATE INDEX IF NOT EXISTS idx_fusions_pair ON fusions (gene_a, gene_b);"),
        rebuild_samples_summary: false,
        view_only: false,
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
        "alt_bases" | "alt_reads" | "coverage" | "coverage_intervals" | "sample_metrics"
        | "fragments" => Some("sample_id"),
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

/// Read the `schema_version` recorded in an attached source database's
/// `geac_metadata` table. Returns `None` when the input predates schema-version
/// stamping (no `geac_metadata` table, or a NULL/absent value) so callers can
/// distinguish "incompatible" from "unverifiable".
fn attached_schema_version(conn: &Connection, alias: &str) -> Result<Option<String>> {
    if !src_table_exists(conn, alias, "geac_metadata")? {
        return Ok(None);
    }
    let version: Option<String> = conn
        .query_row(
            &format!("SELECT schema_version FROM {alias}.geac_metadata LIMIT 1"),
            [],
            |row| row.get(0),
        )
        .optional()
        .with_context(|| format!("failed to read schema_version from attached database '{alias}'"))?
        .flatten();
    Ok(version)
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

/// Column names present in a Parquet file (via DuckDB's `DESCRIBE`).
fn parquet_column_names(
    conn: &Connection,
    path: &Path,
) -> Result<std::collections::HashSet<String>> {
    let escaped = escape_path(path);
    let mut stmt = conn
        .prepare(&format!("DESCRIBE SELECT * FROM read_parquet('{escaped}')"))
        .with_context(|| format!("failed to read Parquet schema from {}", path.display()))?;
    let mut rows = stmt.query([])?;
    let mut cols = std::collections::HashSet::new();
    while let Some(row) = rows.next()? {
        // DESCRIBE columns: column_name, column_type, null, key, default, extra.
        cols.insert(row.get::<_, String>(0)?);
    }
    Ok(cols)
}

/// Guard the `alt_bases` catch-all: any Parquet whose name matches no known GEAC
/// suffix lands here, so a mis-named fusion/coverage/etc. Parquet would otherwise
/// be loaded as a locus table and fail later with a cryptic binder error on a
/// missing provenance column. Detect the mismatch up front and tell the user how
/// to fix it (rename to the right suffix).
fn validate_locus_inputs(conn: &Connection, inputs: &[&PathBuf]) -> Result<()> {
    for path in inputs {
        let cols = parquet_column_names(conn, path)?;

        // A genuine `geac collect` locus Parquet is keyed on chrom + pos. If those
        // are present we accept the file regardless of its other columns.
        if cols.contains("chrom") && cols.contains("pos") {
            continue;
        }

        // Identify the most likely intended type from its signature columns so the
        // error names the exact suffix to use.
        let (kind, suffix) = if cols.contains("gene_a") && cols.contains("gene_b") {
            ("fusion-caller", ".fusions.parquet")
        } else if cols.contains("insert_size") || cols.contains("midpoint") {
            ("fragments", ".fragments.parquet")
        } else if cols.contains("mean_coverage") || cols.contains("coverage") {
            ("coverage", ".coverage.parquet")
        } else {
            bail!(
                "{} does not look like a `geac collect` locus Parquet (missing chrom/pos \
                 columns) and its name matches no known GEAC output suffix. Rename it to the \
                 correct suffix (e.g. .fusions.parquet, .coverage.parquet, .fragments.parquet) \
                 so `geac merge` routes it to the right table.",
                path.display()
            );
        };

        bail!(
            "{path} looks like a {kind} Parquet but was routed to the locus table because its \
             name does not end in '{suffix}'. Rename it to '<sample>{suffix}' (or regenerate it \
             with that --output) so `geac merge` builds the correct table.",
            path = path.display(),
        );
    }
    Ok(())
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

    // The checksum column is only present in alt_bases Parquets.
    let checksum_sha256: Option<String> = if spec.table == "alt_bases" {
        conn.query_row(
            &format!(
                "SELECT input_checksum_sha256 \
                 FROM read_parquet('{escaped}') \
                 WHERE input_checksum_sha256 IS NOT NULL \
                 LIMIT 1"
            ),
            [],
            |row| row.get(0),
        )
        .optional()
        .with_context(|| format!("failed to read checksum from {}", path.display()))?
        .flatten()
    } else {
        None
    };

    Ok(InputProvenance {
        input_path: path.display().to_string(),
        input_kind: spec.table,
        source_kind: "parquet",
        file_size_bytes: file_size_bytes(path),
        modified_at_epoch_seconds: modified_at_epoch_seconds(path),
        checksum_sha256,
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

    if spec.view_only {
        return merge_parquet_view(conn, spec, inputs, input_provenance);
    }

    // The locus table is the catch-all for any Parquet whose name matches no known
    // suffix. Validate those inputs really are locus Parquets before building the
    // table, so a mis-named fusion/coverage Parquet fails with an actionable message.
    if spec.table == "alt_bases" {
        validate_locus_inputs(conn, inputs)?;
    }

    info!(
        table = spec.table,
        n_files = inputs.len(),
        "reading Parquet files"
    );

    // Build the table schema as the union of every input file's columns so
    // older Parquets (missing newly-added columns) and newer ones can coexist
    // — INSERT … BY NAME would otherwise reject sources with extra columns.
    let path_list = inputs
        .iter()
        .map(|p| format!("'{}'", escape_path(p)))
        .collect::<Vec<_>>()
        .join(",");
    conn.execute_batch(&format!(
        "CREATE TABLE {} AS SELECT * FROM read_parquet([{path_list}], union_by_name=true) WHERE FALSE;",
        spec.table
    ))
    .with_context(|| format!("failed to create {} table schema from inputs", spec.table))?;

    for (idx, path) in inputs.iter().enumerate() {
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

fn merge_parquet_view(
    conn: &Connection,
    spec: &TableSpec,
    inputs: &[&PathBuf],
    input_provenance: &mut Vec<InputProvenance>,
) -> Result<()> {
    info!(
        table = spec.table,
        n_files = inputs.len(),
        "registering Parquet files as external view"
    );

    let path_list = inputs
        .iter()
        .map(|p| format!("'{}'", escape_path(p)))
        .collect::<Vec<_>>()
        .join(", ");

    conn.execute_batch(&format!(
        "CREATE VIEW {} AS SELECT * FROM read_parquet([{path_list}]);",
        spec.table
    ))
    .with_context(|| format!("failed to create view {} over Parquet files", spec.table))?;

    info!(table = spec.table, "external view created");

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

        // Refuse to merge inputs written by an incompatible DuckDB schema version.
        // Copying mismatched schemas via `INSERT ... BY NAME` silently drops or
        // NULL-fills columns instead of failing, so validate up front.
        match attached_schema_version(conn, &alias)? {
            Some(version) if version != DUCKDB_SCHEMA_VERSION => {
                bail!(
                    "DuckDB input {} has schema_version '{}', but this GEAC build merges \
                     '{}'. Merging mismatched schema versions can silently drop columns. \
                     Regenerate the input with a matching GEAC version, or upgrade GEAC.",
                    db_path.display(),
                    version,
                    DUCKDB_SCHEMA_VERSION,
                );
            }
            Some(_) => {}
            None => {
                warn!(
                    file = %db_path.display(),
                    expected = DUCKDB_SCHEMA_VERSION,
                    "attached DuckDB has no schema_version in geac_metadata; \
                     cannot verify schema compatibility before merge"
                );
            }
        }

        let mut total_rows = 0_i64;
        for spec in TABLE_SPECS {
            if spec.view_only {
                continue;
            }
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
            checksum_sha256: None,
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

    // Discover which columns actually exist in alt_bases (older Parquet files may
    // predate the resource-path columns; using union_by_name=true during merge
    // means they may be absent).
    let alt_bases_cols: std::collections::HashSet<String> = conn
        .prepare("DESCRIBE SELECT * FROM alt_bases LIMIT 0")
        .and_then(|mut stmt| {
            stmt.query_map([], |row| row.get::<_, String>(0))
                .map(|rows| rows.filter_map(|r| r.ok()).collect::<Vec<String>>())
        })
        .unwrap_or_default()
        .into_iter()
        .collect();

    let resource_cols = [
        "bam_path",
        "bai_path",
        "variants_path",
        "gnomad_path",
        "gnomad_index_path",
        "targets_path",
    ];
    let resource_selects: String = resource_cols
        .iter()
        .map(|col| {
            if alt_bases_cols.contains(*col) {
                format!("ANY_VALUE({col})                        AS {col}")
            } else {
                format!("NULL::VARCHAR                           AS {col}")
            }
        })
        .collect::<Vec<_>>()
        .join(",\n             ");

    let sql = format!(
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
             ANY_VALUE(pipeline)                        AS pipeline,
             {resource_selects}
         FROM alt_bases
         GROUP BY sample_id, batch
         ORDER BY sample_id, batch;"
    );

    conn.execute_batch(&sql)
        .context("failed to build samples summary table")?;

    let n_samples: i64 = conn
        .query_row("SELECT COUNT(*) FROM samples", [], |row| row.get(0))
        .context("failed to count samples")?;
    info!(n_samples, "samples summary table built");
    Ok(())
}

fn build_indices_and_derived_tables(conn: &Connection) -> Result<()> {
    for spec in TABLE_SPECS {
        if spec.view_only {
            continue;
        }
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
        let checksum_sql = input
            .checksum_sha256
            .as_deref()
            .map(|h| format!("'{}'", escape_sql_literal(h)))
            .unwrap_or_else(|| "NULL".to_string());

        conn.execute_batch(&format!(
            "INSERT INTO geac_inputs (
                 input_path, input_kind, source_kind, file_size_bytes,
                 modified_at, checksum_sha256, sample_count, row_count
             ) VALUES (
                 '{}', '{}', '{}', {file_size_sql},
                 {modified_at_sql}, {checksum_sql}, {sample_count_sql}, {row_count_sql}
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
    let sample_metrics_rows = dst_table_row_count(conn, "sample_metrics")?;
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
             n_sample_metrics_inputs BIGINT,
             n_duckdb_inputs BIGINT,
             n_samples BIGINT,
             alt_bases_rows BIGINT,
             alt_reads_rows BIGINT,
             normal_evidence_rows BIGINT,
             pon_evidence_rows BIGINT,
             coverage_rows BIGINT,
             sample_metrics_rows BIGINT,
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
        parquet_groups.get("sample_metrics").map_or(0, Vec::len),
        duckdb_inputs.len(),
        n_samples,
        alt_bases_rows,
        alt_reads_rows,
        normal_evidence_rows,
        pon_evidence_rows,
        coverage_rows,
        sample_metrics_rows,
        samples_rows,
    ))
    .context("failed to write geac_metadata table")?;
    Ok(())
}

fn alt_bases_column_names(conn: &Connection) -> Result<std::collections::HashSet<String>> {
    conn.prepare("DESCRIBE SELECT * FROM alt_bases LIMIT 0")
        .and_then(|mut stmt| {
            stmt.query_map([], |row| row.get::<_, String>(0))
                .map(|rows| {
                    rows.filter_map(|r| r.ok())
                        .collect::<std::collections::HashSet<String>>()
                })
        })
        .context("failed to describe alt_bases columns")
}

fn write_on_target_tsv(conn: &Connection, path: &std::path::Path) -> Result<()> {
    if !dst_table_exists(conn, "alt_bases")? {
        info!("alt_bases table not found; skipping on-target TSV");
        return Ok(());
    }

    let cols = alt_bases_column_names(conn)?;

    if !cols.contains("on_target") {
        info!("on_target column absent from alt_bases (no --targets BED was used); skipping on-target TSV");
        return Ok(());
    }

    // Columns in explorer display order: always-present first, then conditional.
    // Each entry is (sql_expression, output_alias, optional_guard_column).
    // When optional_guard_column is Some(c), the column is included only if c exists.
    let column_specs: &[(&str, &str, Option<&str>)] = &[
        ("sample_id", "sample_id", None),
        ("chrom", "chrom", None),
        ("pos + 1", "pos", None),
        ("ref_allele", "ref_allele", Some("ref_allele")),
        ("alt_allele", "alt_allele", None),
        ("variant_type", "variant_type", Some("variant_type")),
        ("ROUND(alt_count * 1.0 / total_depth, 4)", "vaf", None),
        ("alt_count", "alt_count", None),
        ("ref_count", "ref_count", Some("ref_count")),
        ("total_depth", "total_depth", None),
        ("fwd_alt_count", "fwd_alt_count", Some("fwd_alt_count")),
        ("rev_alt_count", "rev_alt_count", Some("rev_alt_count")),
        (
            "overlap_alt_agree",
            "overlap_alt_agree",
            Some("overlap_alt_agree"),
        ),
        (
            "overlap_alt_disagree",
            "overlap_alt_disagree",
            Some("overlap_alt_disagree"),
        ),
        ("variant_called", "variant_called", Some("variant_called")),
        ("variant_filter", "variant_filter", Some("variant_filter")),
        ("gene", "gene", Some("gene")),
        ("gnomad_af", "gnomad_af", Some("gnomad_af")),
    ];

    let select_parts: Vec<String> = column_specs
        .iter()
        .filter(|(_, _, guard)| guard.is_none_or(|c| cols.contains(c)))
        .map(|(expr, alias, _)| {
            if *expr == *alias {
                expr.to_string()
            } else {
                format!("{expr} AS {alias}")
            }
        })
        .collect();

    let select_sql = select_parts.join(",\n            ");
    let escaped_path = escape_path(path);

    conn.execute_batch(&format!(
        "COPY (
            SELECT
                {select_sql}
            FROM alt_bases
            WHERE on_target = true
            ORDER BY chrom, pos, alt_allele, sample_id
        ) TO '{escaped_path}' (FORMAT CSV, HEADER true, DELIMITER '\t');"
    ))
    .with_context(|| format!("failed to write on-target TSV to {}", path.display()))?;

    let n_rows: i64 = conn
        .query_row(
            "SELECT COUNT(*) FROM alt_bases WHERE on_target = true",
            [],
            |row| row.get(0),
        )
        .context("failed to count on-target rows")?;
    info!(n_rows, path = %path.display(), "on-target TSV written");

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

    if let Some(tsv_path) = &args.on_target_tsv {
        write_on_target_tsv(&conn, tsv_path)?;
    }

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

        // fragments routing test
        let frag_inputs = vec![PathBuf::from("sample.fragments.parquet")];
        let (_, frag_groups) = classify_inputs(&frag_inputs);
        assert_eq!(
            frag_groups["fragments"]
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>(),
            vec!["sample.fragments.parquet"]
        );

        // fusions routing test — must not fall through to the alt_bases catch-all
        let fusion_inputs = vec![PathBuf::from("sample.fusions.parquet")];
        let (_, fusion_groups) = classify_inputs(&fusion_inputs);
        assert_eq!(
            fusion_groups["fusions"]
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>(),
            vec!["sample.fusions.parquet"]
        );
        assert!(fusion_groups["alt_bases"].is_empty());
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

    /// Write a single-row Parquet with the given `SELECT` body to `path`.
    ///
    /// Serialized via a static mutex because DuckDB's parquet extension is
    /// installed to a shared filesystem path (~/.duckdb/extensions/…) on first
    /// use. Without serialization, parallel tests race to install the extension
    /// and DuckDB v1.2.2 fails when two threads hit the install step concurrently.
    fn write_parquet(conn: &Connection, path: &Path, select_body: &str) {
        static LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
        let _guard = LOCK.lock().unwrap();
        let escaped = escape_path(path);
        conn.execute_batch(&format!(
            "COPY (SELECT {select_body}) TO '{escaped}' (FORMAT PARQUET);"
        ))
        .expect("write parquet");
    }

    #[test]
    fn validate_locus_inputs_rejects_misnamed_fusion_parquet() {
        let conn = Connection::open_in_memory().expect("in-memory duckdb");
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("Control_1_simplex.parquet");
        write_parquet(
            &conn,
            &path,
            "'EWSR1' AS gene_a, 'FLI1' AS gene_b, 5 AS supporting_reads",
        );

        let err = validate_locus_inputs(&conn, &[&path])
            .expect_err("misnamed fusion parquet must be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains(".fusions.parquet"),
            "message should name the suffix: {msg}"
        );
        assert!(
            msg.contains("fusion-caller"),
            "message should identify the type: {msg}"
        );
    }

    #[test]
    fn validate_locus_inputs_accepts_real_locus_parquet() {
        let conn = Connection::open_in_memory().expect("in-memory duckdb");
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("sample.parquet");
        write_parquet(&conn, &path, "'chr1' AS chrom, 100 AS pos, 3 AS alt_count");

        validate_locus_inputs(&conn, &[&path]).expect("genuine locus parquet must pass");
    }

    #[test]
    fn validate_locus_inputs_rejects_unknown_parquet() {
        let conn = Connection::open_in_memory().expect("in-memory duckdb");
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("mystery.parquet");
        write_parquet(&conn, &path, "1 AS foo, 2 AS bar");

        let err = validate_locus_inputs(&conn, &[&path])
            .expect_err("unrecognized parquet must be rejected");
        assert!(
            err.to_string().contains("missing chrom/pos"),
            "message should explain why: {err}"
        );
    }

    /// Build a DuckDB file at `path` containing only the supplied
    /// `geac_metadata`/table SQL, used to exercise the schema-version guard.
    fn write_source_duckdb(path: &Path, setup_sql: &str) {
        let src = Connection::open(path).expect("open source duckdb");
        src.execute_batch(setup_sql)
            .expect("populate source duckdb");
    }

    #[test]
    fn attached_schema_version_reads_metadata() {
        let dir = tempfile::tempdir().expect("tempdir");
        let src_path = dir.path().join("src.duckdb");
        write_source_duckdb(
            &src_path,
            "CREATE TABLE geac_metadata (schema_version VARCHAR);
             INSERT INTO geac_metadata VALUES ('duckdb-v4');",
        );

        let conn = Connection::open_in_memory().expect("in-memory duckdb");
        let escaped = escape_path(&src_path);
        conn.execute_batch(&format!("ATTACH '{escaped}' AS _src0 (READ_ONLY);"))
            .expect("attach source");

        assert_eq!(
            attached_schema_version(&conn, "_src0").expect("read schema_version"),
            Some("duckdb-v4".to_string())
        );
    }

    #[test]
    fn attached_schema_version_none_when_metadata_absent() {
        let dir = tempfile::tempdir().expect("tempdir");
        let src_path = dir.path().join("legacy.duckdb");
        write_source_duckdb(&src_path, "CREATE TABLE alt_bases (chrom VARCHAR);");

        let conn = Connection::open_in_memory().expect("in-memory duckdb");
        let escaped = escape_path(&src_path);
        conn.execute_batch(&format!("ATTACH '{escaped}' AS _src0 (READ_ONLY);"))
            .expect("attach source");

        assert_eq!(
            attached_schema_version(&conn, "_src0").expect("read schema_version"),
            None,
            "legacy inputs without geac_metadata must read as unverifiable, not error"
        );
    }

    #[test]
    fn merge_duckdb_inputs_rejects_mismatched_schema_version() {
        let dir = tempfile::tempdir().expect("tempdir");
        let src_path = dir.path().join("old.duckdb");
        write_source_duckdb(
            &src_path,
            "CREATE TABLE geac_metadata (schema_version VARCHAR);
             INSERT INTO geac_metadata VALUES ('duckdb-v1');
             CREATE TABLE alt_bases (chrom VARCHAR);",
        );

        let conn = Connection::open_in_memory().expect("in-memory duckdb");
        let mut provenance = Vec::new();
        let inputs = vec![&src_path];
        let err = merge_duckdb_inputs(&conn, &inputs, &mut provenance)
            .expect_err("mismatched schema version must be rejected before copying");
        let msg = err.to_string();
        assert!(
            msg.contains("duckdb-v1"),
            "message should name the input version: {msg}"
        );
        assert!(
            msg.contains(DUCKDB_SCHEMA_VERSION),
            "message should name the expected version: {msg}"
        );
    }

    #[test]
    fn all_table_specs_present_in_schema_json() {
        let raw = std::fs::read_to_string("schema/geac_schema.json")
            .expect("read schema/geac_schema.json");
        let manifest: serde_json::Value =
            serde_json::from_str(&raw).expect("parse schema/geac_schema.json");
        let schema_tables: std::collections::HashSet<&str> = manifest["tables"]
            .as_object()
            .expect("schema tables must be an object")
            .keys()
            .map(String::as_str)
            .collect();

        let missing: Vec<&str> = TABLE_SPECS
            .iter()
            .map(|s| s.table)
            .filter(|t| !schema_tables.contains(t))
            .collect();

        assert!(
            missing.is_empty(),
            "TABLE_SPECS tables missing from schema/geac_schema.json: {missing:?}\n\
             Add them to the \"tables\" map in schema/geac_schema.json."
        );
    }
}
