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

    // Build a DuckDB array literal: ['file1.parquet', 'file2.parquet', ...]
    // Paths are escaped to handle spaces and special characters.
    let file_list: Vec<String> = args
        .inputs
        .iter()
        .map(|p| format!("'{}'", p.display().to_string().replace('\'', "''")))
        .collect();
    let file_array = file_list.join(", ");

    // Create the main table by reading all Parquet files in parallel.
    // DuckDB infers the schema from the Parquet files automatically.
    info!("reading Parquet files and creating alt_bases table...");
    conn.execute_batch(&format!(
        "CREATE TABLE alt_bases AS
         SELECT * FROM read_parquet([{file_array}]);"
    ))
    .context("failed to create alt_bases table")?;

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
