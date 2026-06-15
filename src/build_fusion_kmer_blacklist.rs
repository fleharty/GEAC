use anyhow::{Context, Result};
use duckdb::Connection;
use tracing::info;

use crate::cli::BuildFusionKmerBlacklistArgs;

pub fn build_fusion_kmer_blacklist(args: &BuildFusionKmerBlacklistArgs) -> Result<()> {
    if args.kmer_hits.is_empty() {
        anyhow::bail!("no --kmer-hits files provided");
    }
    for path in &args.kmer_hits {
        if !path.exists() {
            anyhow::bail!("kmer-hits file not found: {}", path.display());
        }
    }
    if args.output.exists() {
        anyhow::bail!(
            "output file already exists: {}. Remove it or choose a different path.",
            args.output.display()
        );
    }

    // Build a comma-separated quoted list of input paths for DuckDB's read_csv().
    let path_list = args
        .kmer_hits
        .iter()
        .map(|p| format!("'{}'", p.display().to_string().replace('\'', "''")))
        .collect::<Vec<_>>()
        .join(", ");

    let n_inputs = args.kmer_hits.len();
    let min_samples = args.min_pon_samples;
    let output = args.output.display().to_string();

    info!(
        n_inputs,
        min_pon_samples = min_samples,
        "aggregating k-mer hits from PoN samples..."
    );

    let conn = Connection::open_in_memory().context("failed to open in-memory DuckDB")?;

    // Read all TSVs, group by kmer_hash, count distinct samples, filter by threshold,
    // and write directly to Parquet.
    let sql = format!(
        "COPY (
            SELECT
                kmer_hash,
                COUNT(DISTINCT sample_id) AS n_pon_samples
            FROM read_csv([{path_list}], header=true, delim='\\t', union_by_name=true)
            GROUP BY kmer_hash
            HAVING COUNT(DISTINCT sample_id) >= {min_samples}
         ) TO '{output}' (FORMAT PARQUET);",
    );

    conn.execute_batch(&sql)
        .context("failed to aggregate k-mer hits and write blacklist Parquet")?;

    let n_blacklisted: i64 = conn
        .query_row(
            &format!("SELECT COUNT(*) FROM read_parquet('{output}')"),
            [],
            |r| r.get(0),
        )
        .unwrap_or(0);

    info!(
        n_blacklisted_kmers = n_blacklisted,
        output = %args.output.display(),
        "k-mer blacklist written"
    );
    Ok(())
}
