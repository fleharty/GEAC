use anyhow::{bail, Context, Result};
use duckdb::Connection;
use tracing::info;

use crate::cli::AnnotatePonArgs;
use crate::record::PonEvidence;

/// Cross-annotate tumor alt-base loci against a Panel of Normals DuckDB.
///
/// Opens the PoN DuckDB (which must contain an `alt_bases` table produced by
/// `geac merge`) and joins the tumor locus Parquet against it using DuckDB's
/// inline `read_parquet()` function.  No BAM pileup is required.
///
/// For each tumor alt locus the query returns:
///   n_pon_samples     — PoN samples that carry this exact (chrom, pos, alt_allele)
///   pon_total_samples — total distinct samples in the PoN (denominator)
///   max_pon_vaf       — highest alt_count/total_depth seen in any PoN sample (NULL if absent)
///   mean_pon_vaf      — mean alt_count/total_depth across PoN samples (NULL if absent)
pub fn annotate_pon(args: &AnnotatePonArgs) -> Result<Vec<PonEvidence>> {
    if !args.tumor_parquet.exists() {
        bail!("tumor Parquet not found: {}", args.tumor_parquet.display());
    }
    if !args.pon_db.exists() {
        bail!("PoN DuckDB not found: {}", args.pon_db.display());
    }

    // Open the PoN DuckDB. The connection is read-only in practice — the SQL
    // below only SELECTs from it.
    let conn = Connection::open(&args.pon_db)
        .with_context(|| format!("failed to open PoN DuckDB: {}", args.pon_db.display()))?;

    // Verify the PoN database has an alt_bases table.
    conn.execute_batch("SELECT 1 FROM alt_bases LIMIT 0")
        .context("PoN DuckDB does not contain an alt_bases table — was it built with `geac merge`?")?;

    let tumor_path = args.tumor_parquet.display().to_string().replace('\'', "''");

    info!(
        tumor_parquet = %args.tumor_parquet.display(),
        pon_db        = %args.pon_db.display(),
        "running PoN annotation query"
    );

    let sql = format!(
        r#"
        WITH
        pon_total AS (
            SELECT COUNT(DISTINCT sample_id) AS n
            FROM alt_bases
        ),
        pon_agg AS (
            SELECT
                chrom,
                CAST(pos AS BIGINT)                               AS pos,
                alt_allele,
                COUNT(DISTINCT sample_id)                         AS n_pon_samples,
                MAX(alt_count * 1.0 / NULLIF(total_depth, 0))    AS max_pon_vaf,
                AVG(alt_count * 1.0 / NULLIF(total_depth, 0))    AS mean_pon_vaf
            FROM alt_bases
            GROUP BY chrom, pos, alt_allele
        ),
        tumor AS (
            SELECT DISTINCT
                sample_id           AS tumor_sample_id,
                chrom,
                CAST(pos AS BIGINT) AS pos,
                alt_allele          AS tumor_alt_allele
            FROM read_parquet('{tumor_path}')
        )
        SELECT
            t.tumor_sample_id,
            t.chrom,
            t.pos,
            t.tumor_alt_allele,
            COALESCE(p.n_pon_samples, 0)  AS n_pon_samples,
            pt.n                           AS pon_total_samples,
            p.max_pon_vaf,
            p.mean_pon_vaf
        FROM tumor t
        CROSS JOIN pon_total pt
        LEFT JOIN pon_agg p
               ON p.chrom     = t.chrom
              AND p.pos        = t.pos
              AND p.alt_allele = t.tumor_alt_allele
        ORDER BY t.chrom, t.pos, t.tumor_alt_allele
        "#
    );

    let mut stmt = conn.prepare(&sql).context("failed to prepare PoN annotation query")?;

    let rows = stmt
        .query_map([], |row| {
            Ok(PonEvidence {
                tumor_sample_id:   row.get(0)?,
                chrom:             row.get(1)?,
                pos:               row.get(2)?,
                tumor_alt_allele:  row.get(3)?,
                n_pon_samples:     row.get(4)?,
                pon_total_samples: row.get(5)?,
                max_pon_vaf:       row.get(6)?,
                mean_pon_vaf:      row.get(7)?,
            })
        })
        .context("failed to execute PoN annotation query")?;

    let mut results = Vec::new();
    for row in rows {
        results.push(row.context("failed to read PoN evidence row")?);
    }

    info!(n_records = results.len(), "PoN annotation complete");
    Ok(results)
}
