use anyhow::{bail, Context, Result};
use duckdb::Connection;
use tracing::info;

use crate::cli::MergeArgs;

/// Returns true if `table` exists in the attached database aliased as `alias`.
fn src_table_exists(conn: &Connection, alias: &str, table: &str) -> Result<bool> {
    let count: i64 = conn
        .query_row(
            &format!(
                "SELECT COUNT(*) FROM {alias}.information_schema.tables \
                 WHERE table_name = '{table}' AND table_schema = 'main'"
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
                "SELECT COUNT(*) FROM information_schema.tables \
                 WHERE table_name = '{table}' AND table_schema = 'main'"
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
            "INSERT INTO {table} SELECT * FROM {alias}.{table};"
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

    // Validate all inputs exist before opening the database
    for path in &args.inputs {
        if !path.exists() {
            bail!("input file not found: {}", path.display());
        }
    }

    // Partition inputs into DuckDB files and Parquet files.
    let is_duckdb = |p: &&std::path::PathBuf| {
        p.extension().and_then(|e| e.to_str()) == Some("duckdb")
    };
    let duckdb_inputs: Vec<&std::path::PathBuf> = args.inputs.iter().filter(is_duckdb).collect();
    let parquet_inputs: Vec<&std::path::PathBuf> = args.inputs.iter().filter(|p| !is_duckdb(p)).collect();

    // Partition Parquet inputs by filename suffix.
    fn suffix(p: &&std::path::PathBuf, sfx: &str) -> bool {
        p.file_name().and_then(|n| n.to_str()).is_some_and(|n| n.ends_with(sfx))
    }
    let reads_inputs:    Vec<&std::path::PathBuf> = parquet_inputs.iter().copied().filter(|p| suffix(p, ".reads.parquet")).collect();
    let normal_inputs:   Vec<&std::path::PathBuf> = parquet_inputs.iter().copied().filter(|p| suffix(p, ".normal_evidence.parquet")).collect();
    let pon_inputs:      Vec<&std::path::PathBuf> = parquet_inputs.iter().copied().filter(|p| suffix(p, ".pon_evidence.parquet")).collect();
    let coverage_inputs: Vec<&std::path::PathBuf> = parquet_inputs.iter().copied().filter(|p| suffix(p, ".coverage.parquet")).collect();
    let locus_inputs:    Vec<&std::path::PathBuf> = parquet_inputs.iter().copied().filter(|p| {
        !suffix(p, ".reads.parquet")
            && !suffix(p, ".normal_evidence.parquet")
            && !suffix(p, ".pon_evidence.parquet")
            && !suffix(p, ".coverage.parquet")
    }).collect();

    info!(
        n_locus_files    = locus_inputs.len(),
        n_reads_files    = reads_inputs.len(),
        n_normal_files   = normal_inputs.len(),
        n_pon_files      = pon_inputs.len(),
        n_coverage_files = coverage_inputs.len(),
        n_duckdb_files   = duckdb_inputs.len(),
        output = %args.output.display(),
        "merging files into DuckDB"
    );

    let conn = Connection::open(&args.output)
        .with_context(|| format!("failed to create DuckDB database: {}", args.output.display()))?;

    // ── Phase 1: Parquet inputs ────────────────────────────────────────────────

    // alt_bases
    if !locus_inputs.is_empty() {
        info!("reading locus Parquet files and creating alt_bases table...");

        let first = locus_inputs[0];
        let first_escaped = first.display().to_string().replace('\'', "''");
        conn.execute_batch(&format!(
            "CREATE TABLE alt_bases AS SELECT * FROM read_parquet('{first_escaped}');"
        ))
        .context("failed to create alt_bases table from first file")?;

        for (idx, path) in locus_inputs[1..].iter().enumerate() {
            let escaped = path.display().to_string().replace('\'', "''");
            info!(file = %path.display(), idx = idx + 1, total = locus_inputs.len(), "inserting locus file");
            conn.execute_batch(&format!(
                "INSERT INTO alt_bases SELECT * FROM read_parquet('{escaped}');"
            ))
            .with_context(|| format!("failed to insert {}", path.display()))?;
        }

        let n_rows: i64 = conn
            .query_row("SELECT COUNT(*) FROM alt_bases", [], |row| row.get(0))
            .context("failed to count rows")?;
        info!(n_rows, "alt_bases table created from Parquet");
    }

    // alt_reads
    if !reads_inputs.is_empty() {
        info!("reading reads Parquet files and creating alt_reads table...");

        let first = reads_inputs[0];
        let first_escaped = first.display().to_string().replace('\'', "''");
        conn.execute_batch(&format!(
            "CREATE TABLE alt_reads AS SELECT * FROM read_parquet('{first_escaped}');"
        ))
        .context("failed to create alt_reads table from first reads file")?;

        for (idx, path) in reads_inputs[1..].iter().enumerate() {
            let escaped = path.display().to_string().replace('\'', "''");
            info!(file = %path.display(), idx = idx + 1, total = reads_inputs.len(), "inserting reads file");
            conn.execute_batch(&format!(
                "INSERT INTO alt_reads SELECT * FROM read_parquet('{escaped}');"
            ))
            .with_context(|| format!("failed to insert reads file {}", path.display()))?;
        }

        let n: i64 = conn
            .query_row("SELECT COUNT(*) FROM alt_reads", [], |row| row.get(0))
            .context("failed to count alt_reads rows")?;
        info!(n, "alt_reads table created from Parquet");
    }

    // normal_evidence
    if !normal_inputs.is_empty() {
        info!("reading normal evidence Parquet files...");

        let first = normal_inputs[0];
        let first_escaped = first.display().to_string().replace('\'', "''");
        conn.execute_batch(&format!(
            "CREATE TABLE normal_evidence AS SELECT * FROM read_parquet('{first_escaped}');"
        ))
        .context("failed to create normal_evidence table")?;

        for (idx, path) in normal_inputs[1..].iter().enumerate() {
            let escaped = path.display().to_string().replace('\'', "''");
            info!(file = %path.display(), idx = idx + 1, total = normal_inputs.len(), "inserting normal evidence file");
            conn.execute_batch(&format!(
                "INSERT INTO normal_evidence SELECT * FROM read_parquet('{escaped}');"
            ))
            .with_context(|| format!("failed to insert normal evidence file {}", path.display()))?;
        }

        let n: i64 = conn
            .query_row("SELECT COUNT(*) FROM normal_evidence", [], |row| row.get(0))
            .context("failed to count normal_evidence rows")?;
        info!(n, "normal_evidence table created from Parquet");
    }

    // pon_evidence
    if !pon_inputs.is_empty() {
        info!("reading PoN evidence Parquet files...");

        let first = pon_inputs[0];
        let first_escaped = first.display().to_string().replace('\'', "''");
        conn.execute_batch(&format!(
            "CREATE TABLE pon_evidence AS SELECT * FROM read_parquet('{first_escaped}');"
        ))
        .context("failed to create pon_evidence table")?;

        for (idx, path) in pon_inputs[1..].iter().enumerate() {
            let escaped = path.display().to_string().replace('\'', "''");
            info!(file = %path.display(), idx = idx + 1, total = pon_inputs.len(), "inserting PoN evidence file");
            conn.execute_batch(&format!(
                "INSERT INTO pon_evidence SELECT * FROM read_parquet('{escaped}');"
            ))
            .with_context(|| format!("failed to insert PoN evidence file {}", path.display()))?;
        }

        let n: i64 = conn
            .query_row("SELECT COUNT(*) FROM pon_evidence", [], |row| row.get(0))
            .context("failed to count pon_evidence rows")?;
        info!(n, "pon_evidence table created from Parquet");
    }

    // coverage
    if !coverage_inputs.is_empty() {
        info!("reading coverage Parquet files...");

        let first = coverage_inputs[0];
        let first_escaped = first.display().to_string().replace('\'', "''");
        conn.execute_batch(&format!(
            "CREATE TABLE coverage AS SELECT * FROM read_parquet('{first_escaped}');"
        ))
        .context("failed to create coverage table")?;

        for (idx, path) in coverage_inputs[1..].iter().enumerate() {
            let escaped = path.display().to_string().replace('\'', "''");
            info!(file = %path.display(), idx = idx + 1, total = coverage_inputs.len(), "inserting coverage file");
            conn.execute_batch(&format!(
                "INSERT INTO coverage SELECT * FROM read_parquet('{escaped}');"
            ))
            .with_context(|| format!("failed to insert coverage file {}", path.display()))?;
        }

        let n: i64 = conn
            .query_row("SELECT COUNT(*) FROM coverage", [], |row| row.get(0))
            .context("failed to count coverage rows")?;
        info!(n, "coverage table created from Parquet");
    }

    // ── Phase 2: DuckDB inputs ─────────────────────────────────────────────────
    // Known data tables to copy.  The derived `samples` table is intentionally
    // excluded — it will be rebuilt from the merged alt_bases at the end.
    const DATA_TABLES: &[&str] = &[
        "alt_bases",
        "alt_reads",
        "normal_evidence",
        "pon_evidence",
        "coverage",
    ];

    for (i, db_path) in duckdb_inputs.iter().enumerate() {
        let alias = format!("_src{i}");
        let escaped = db_path.display().to_string().replace('\'', "''");
        info!(file = %db_path.display(), "attaching source DuckDB");

        conn.execute_batch(&format!("ATTACH '{escaped}' AS {alias} (READ_ONLY);"))
            .with_context(|| format!("failed to attach {}", db_path.display()))?;

        for &table in DATA_TABLES {
            if src_table_exists(&conn, &alias, table)? {
                let n_src: i64 = conn
                    .query_row(
                        &format!("SELECT COUNT(*) FROM {alias}.{table}"),
                        [],
                        |row| row.get(0),
                    )
                    .with_context(|| format!("failed to count rows in {alias}.{table}"))?;

                info!(table, n_src, file = %db_path.display(), "copying table from source DuckDB");
                copy_table(&conn, &alias, table)?;
            }
        }

        conn.execute_batch(&format!("DETACH {alias};"))
            .with_context(|| format!("failed to detach {}", db_path.display()))?;
    }

    // ── Phase 3: Indices and derived tables ───────────────────────────────────
    if dst_table_exists(&conn, "alt_bases")? {
        info!("creating indices on alt_bases...");
        conn.execute_batch(
            "CREATE INDEX IF NOT EXISTS idx_chrom_pos ON alt_bases (chrom, pos);
             CREATE INDEX IF NOT EXISTS idx_sample_id ON alt_bases (sample_id);",
        )
        .context("failed to create alt_bases indices")?;

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
    }

    if dst_table_exists(&conn, "alt_reads")? {
        info!("creating index on alt_reads...");
        conn.execute_batch(
            "CREATE INDEX IF NOT EXISTS idx_reads_locus \
             ON alt_reads (sample_id, chrom, pos, alt_allele);",
        )
        .context("failed to create alt_reads index")?;
    }

    if dst_table_exists(&conn, "normal_evidence")? {
        info!("creating index on normal_evidence...");
        conn.execute_batch(
            "CREATE INDEX IF NOT EXISTS idx_ne_locus \
             ON normal_evidence (tumor_sample_id, chrom, pos, tumor_alt_allele);",
        )
        .context("failed to create normal_evidence index")?;
    }

    if dst_table_exists(&conn, "pon_evidence")? {
        info!("creating index on pon_evidence...");
        conn.execute_batch(
            "CREATE INDEX IF NOT EXISTS idx_pe_locus \
             ON pon_evidence (tumor_sample_id, chrom, pos, tumor_alt_allele);",
        )
        .context("failed to create pon_evidence index")?;
    }

    if dst_table_exists(&conn, "coverage")? {
        info!("creating index on coverage...");
        conn.execute_batch(
            "CREATE INDEX IF NOT EXISTS idx_coverage_locus \
             ON coverage (sample_id, chrom, pos);",
        )
        .context("failed to create coverage index")?;
    }

    info!(output = %args.output.display(), "merge complete");
    Ok(())
}
