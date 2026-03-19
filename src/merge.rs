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

    info!(
        n_files = args.inputs.len(),
        output = %args.output.display(),
        "merging Parquet files into DuckDB"
    );

    let conn = Connection::open(&args.output)
        .with_context(|| format!("failed to create DuckDB database: {}", args.output.display()))?;

    // Create the table from the first file to establish the schema, then
    // INSERT the remaining files one at a time.  This avoids the multi-file
    // schema-inference path in DuckDB which can segfault when given a large
    // array literal or when files have subtly different nullable/type metadata.
    info!("reading Parquet files and creating alt_bases table...");

    let first = &args.inputs[0];
    let first_escaped = first.display().to_string().replace('\'', "''");
    conn.execute_batch(&format!(
        "CREATE TABLE alt_bases AS SELECT * FROM read_parquet('{first_escaped}');"
    ))
    .context("failed to create alt_bases table from first file")?;

    for (idx, path) in args.inputs[1..].iter().enumerate() {
        let escaped = path.display().to_string().replace('\'', "''");
        info!(
            file = %path.display(),
            idx = idx + 1,
            total = args.inputs.len(),
            "inserting"
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
    info!(output = %args.output.display(), "merge complete");

    Ok(())
}
