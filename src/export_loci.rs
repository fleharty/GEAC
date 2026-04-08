use std::io::Write;

use anyhow::{Context, Result};
use duckdb::Connection;
use tracing::info;

use crate::cli::ExportLociArgs;

/// Export a TSV of distinct loci from a cohort DuckDB or single-sample Parquet.
///
/// The output has two columns: `chrom` and `pos` (0-based), one row per distinct
/// locus that passes the VAF and sample-count filters.  The file is suitable as
/// input to `geac locus-depth --loci`.
pub fn export_loci(args: &ExportLociArgs) -> Result<()> {
    let db = Connection::open_in_memory().context("failed to open in-memory DuckDB")?;

    let input_str = args.input.display().to_string().replace('\'', "''");
    let ext = args
        .input
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("");

    // Accept both a DuckDB file and a Parquet file.
    match ext {
        "duckdb" => {
            db.execute_batch(&format!(
                "ATTACH '{input_str}' AS cohort (READ_ONLY); \
                 CREATE TABLE alt_bases AS SELECT * FROM cohort.alt_bases;"
            ))
            .context("failed to attach cohort DuckDB")?;
        }
        _ => {
            db.execute_batch(&format!(
                "CREATE TABLE alt_bases AS SELECT * FROM read_parquet('{input_str}');"
            ))
            .context("failed to load Parquet")?;
        }
    }

    // Build WHERE clause for VAF filter.
    let mut conditions = vec![
        "total_depth > 0".to_string(),
        format!(
            "alt_count * 1.0 / total_depth >= {}",
            args.min_vaf
        ),
    ];

    if let Some(max_vaf) = args.max_vaf {
        conditions.push(format!("alt_count * 1.0 / total_depth <= {}", max_vaf));
    }

    if !args.variant_types.is_empty() {
        let quoted: Vec<String> = args
            .variant_types
            .iter()
            .map(|t| format!("'{}'", t.replace('\'', "''")))
            .collect();
        conditions.push(format!("variant_type IN ({})", quoted.join(", ")));
    }

    let where_clause = conditions.join(" AND ");

    // If --min-samples > 1 we need a HAVING clause on the aggregated count.
    let sql = if args.min_samples > 1 {
        format!(
            "SELECT chrom, pos \
             FROM (
                 SELECT chrom, pos, COUNT(DISTINCT sample_id) AS n_samples \
                 FROM alt_bases \
                 WHERE {where_clause} \
                 GROUP BY chrom, pos
             ) \
             WHERE n_samples >= {} \
             ORDER BY chrom, pos",
            args.min_samples
        )
    } else {
        format!(
            "SELECT DISTINCT chrom, pos \
             FROM alt_bases \
             WHERE {where_clause} \
             ORDER BY chrom, pos"
        )
    };

    let mut stmt = db.prepare(&sql).context("failed to prepare loci query")?;

    let rows = stmt
        .query_map([], |row| {
            Ok((row.get::<_, String>(0)?, row.get::<_, i64>(1)?))
        })
        .context("failed to query loci")?;

    let out_file = std::fs::File::create(&args.output)
        .with_context(|| format!("failed to create output file: {}", args.output.display()))?;
    let mut writer = std::io::BufWriter::new(out_file);

    writeln!(writer, "chrom\tpos").context("failed to write header")?;

    let mut n_loci: u64 = 0;
    for row in rows {
        let (chrom, pos) = row.context("failed to read locus row")?;
        writeln!(writer, "{chrom}\t{pos}").context("failed to write locus row")?;
        n_loci += 1;
    }

    writer.flush().context("failed to flush output")?;

    info!(
        n_loci,
        output = %args.output.display(),
        "loci exported"
    );

    Ok(())
}
