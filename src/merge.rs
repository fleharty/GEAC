use anyhow::{bail, Context, Result};
use duckdb::Connection;
use tracing::info;

use crate::cli::MergeArgs;

pub fn merge(args: &MergeArgs) -> Result<()> {
    if args.output.exists() {
        bail!(
            "output file already exists: {}. Remove it or choose a different path.",
            args.output.display()
        );
    }

    if args.inputs.is_empty() {
        bail!("no input Parquet files provided");
    }

    // Validate all inputs exist before opening the database
    for path in &args.inputs {
        if !path.exists() {
            bail!("input file not found: {}", path.display());
        }
    }

    // Split inputs: files ending in ".reads.parquet" go to alt_reads; everything else to alt_bases.
    let (reads_inputs, locus_inputs): (Vec<_>, Vec<_>) = args.inputs.iter().partition(|p| {
        p.file_name()
            .and_then(|n| n.to_str())
            .is_some_and(|n| n.ends_with(".reads.parquet"))
    });

    if locus_inputs.is_empty() {
        bail!("no locus Parquet files provided (files ending in .reads.parquet are excluded from the locus table)");
    }

    info!(
        n_locus_files = locus_inputs.len(),
        n_reads_files = reads_inputs.len(),
        output = %args.output.display(),
        "merging Parquet files into DuckDB"
    );

    let conn = Connection::open(&args.output)
        .with_context(|| format!("failed to create DuckDB database: {}", args.output.display()))?;

    // ── alt_bases table ───────────────────────────────────────────────────────
    // Create the table from the first file to establish the schema, then
    // INSERT the remaining files one at a time.  This avoids the multi-file
    // schema-inference path in DuckDB which can segfault when given a large
    // array literal or when files have subtly different nullable/type metadata.
    info!("reading locus Parquet files and creating alt_bases table...");

    let first = locus_inputs[0];
    let first_escaped = first.display().to_string().replace('\'', "''");
    conn.execute_batch(&format!(
        "CREATE TABLE alt_bases AS SELECT * FROM read_parquet('{first_escaped}');"
    ))
    .context("failed to create alt_bases table from first file")?;

    for (idx, path) in locus_inputs[1..].iter().enumerate() {
        let escaped = path.display().to_string().replace('\'', "''");
        info!(
            file = %path.display(),
            idx = idx + 1,
            total = locus_inputs.len(),
            "inserting locus file"
        );
        conn.execute_batch(&format!(
            "INSERT INTO alt_bases SELECT * FROM read_parquet('{escaped}');"
        ))
        .with_context(|| format!("failed to insert {}", path.display()))?;
    }

    let n_rows: i64 = conn
        .query_row("SELECT COUNT(*) FROM alt_bases", [], |row| row.get(0))
        .context("failed to count rows")?;

    info!(n_rows, "alt_bases table created");

    // Create indices for common query patterns
    info!("creating indices...");
    conn.execute_batch(
        "CREATE INDEX idx_chrom_pos ON alt_bases (chrom, pos);
         CREATE INDEX idx_sample_id ON alt_bases (sample_id);",
    )
    .context("failed to create indices")?;

    // Create a samples summary table
    conn.execute_batch(
        "CREATE TABLE samples AS
         SELECT
             sample_id,
             COUNT(*)                                   AS n_alt_loci,
             SUM(alt_count)                             AS total_alt_reads,
             COUNT(DISTINCT chrom || ':' || pos)        AS n_positions,
             MIN(pos)                                   AS min_pos,
             MAX(pos)                                   AS max_pos,
             ANY_VALUE(read_type)                       AS read_type,
             ANY_VALUE(pipeline)                        AS pipeline
         FROM alt_bases
         GROUP BY sample_id
         ORDER BY sample_id;",
    )
    .context("failed to create samples summary table")?;

    let n_samples: i64 = conn
        .query_row("SELECT COUNT(*) FROM samples", [], |row| row.get(0))
        .context("failed to count samples")?;

    info!(n_samples, "samples summary table created");

    // ── alt_reads table (optional) ────────────────────────────────────────────
    if !reads_inputs.is_empty() {
        info!("reading reads Parquet files and creating alt_reads table...");

        let first_reads = reads_inputs[0];
        let first_reads_escaped = first_reads.display().to_string().replace('\'', "''");
        conn.execute_batch(&format!(
            "CREATE TABLE alt_reads AS SELECT * FROM read_parquet('{first_reads_escaped}');"
        ))
        .context("failed to create alt_reads table from first reads file")?;

        for (idx, path) in reads_inputs[1..].iter().enumerate() {
            let escaped = path.display().to_string().replace('\'', "''");
            info!(
                file = %path.display(),
                idx = idx + 1,
                total = reads_inputs.len(),
                "inserting reads file"
            );
            conn.execute_batch(&format!(
                "INSERT INTO alt_reads SELECT * FROM read_parquet('{escaped}');"
            ))
            .with_context(|| format!("failed to insert reads file {}", path.display()))?;
        }

        let n_read_rows: i64 = conn
            .query_row("SELECT COUNT(*) FROM alt_reads", [], |row| row.get(0))
            .context("failed to count alt_reads rows")?;

        info!(n_read_rows, "alt_reads table created");

        conn.execute_batch(
            "CREATE INDEX idx_reads_locus ON alt_reads (sample_id, chrom, pos, alt_allele);",
        )
        .context("failed to create alt_reads index")?;
    }

    info!(output = %args.output.display(), "merge complete");

    Ok(())
}
