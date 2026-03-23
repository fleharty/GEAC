use anyhow::{Context, Result};
use duckdb::Connection;

use crate::cli::CohortArgs;

pub fn run_cohort(args: &CohortArgs) -> Result<()> {
    anyhow::ensure!(!args.inputs.is_empty(), "no input Parquet files provided");

    let con = Connection::open_in_memory().context("failed to open DuckDB")?;

    // ── Load all Parquet files into a view ─────────────────────────────────────
    let file_list: Vec<String> = args
        .inputs
        .iter()
        .map(|p| format!("'{}'", p.display().to_string().replace('\'', "''")))
        .collect();
    let source = if file_list.len() == 1 {
        format!("read_parquet({}, union_by_name=true)", file_list[0])
    } else {
        format!("read_parquet([{}], union_by_name=true)", file_list.join(", "))
    };

    let filter = if args.on_target_only {
        "WHERE on_target = true"
    } else {
        ""
    };

    con.execute_batch(&format!(
        "CREATE VIEW data AS SELECT * FROM {source} {filter}"
    ))
    .context("failed to create data view")?;

    // ── Detect optional annotation columns ────────────────────────────────────
    let schema_cols: Vec<String> = {
        let mut stmt = con.prepare("DESCRIBE SELECT * FROM data LIMIT 0")?;
        stmt.query_map([], |row| row.get::<_, String>(0))?
            .collect::<Result<Vec<_>, _>>()?
    };
    let has = |col: &str| schema_cols.iter().any(|c| c == col);

    // Extra annotation columns included in the file output only.
    let annot_select = [
        ("gene",           "ANY_VALUE(gene)           AS gene"),
        ("on_target",      "ANY_VALUE(on_target)      AS on_target"),
        ("homopolymer_len","ANY_VALUE(homopolymer_len) AS homopolymer_len"),
        ("str_len",        "ANY_VALUE(str_len)         AS str_len"),
        ("trinuc_context", "ANY_VALUE(trinuc_context)  AS trinuc_context"),
    ]
    .iter()
    .filter(|(col, _)| has(col))
    .map(|(_, expr)| format!(", {expr}"))
    .collect::<String>();

    // ── Shared CTE used by both queries ────────────────────────────────────────
    // The loci CTE computes per-locus aggregates.
    // The outer SELECT joins in the total sample count.
    let where_clause = format!(
        "l.n_samples >= {min_samples}
         AND l.n_samples * 1.0 / t.n_total >= {min_frac}",
        min_samples = args.min_samples,
        min_frac    = args.min_sample_fraction,
    );

    let cte = format!(
        "WITH totals AS (
             SELECT COUNT(DISTINCT sample_id) AS n_total FROM data
         ),
         loci AS (
             SELECT
                 chrom,
                 pos,
                 ref_allele,
                 alt_allele,
                 variant_type,
                 COUNT(DISTINCT sample_id)                              AS n_samples,
                 ROUND(AVG(alt_count * 1.0 / total_depth), 6)          AS mean_vaf,
                 ROUND(MAX(alt_count * 1.0 / total_depth), 6)          AS max_vaf,
                 ROUND(AVG(total_depth), 1)                            AS mean_depth,
                 ROUND(AVG(alt_count), 2)                              AS mean_alt_count,
                 ROUND(
                     AVG(fwd_alt_count * 1.0
                         / NULLIF(fwd_alt_count + rev_alt_count, 0)),
                     4
                 )                                                     AS mean_strand_balance,
                 ROUND(
                     MAX(fwd_alt_count * 1.0
                         / NULLIF(fwd_alt_count + rev_alt_count, 0)),
                     4
                 )                                                     AS max_strand_balance
                 {annot_select}
             FROM data
             GROUP BY chrom, pos, ref_allele, alt_allele, variant_type
         )"
    );

    // ── Top-N display query (explicit columns, stable indices) ─────────────────
    // Keep annotation columns out of this query so column indices are fixed.
    let display_sql = format!(
        "{cte}
         SELECT
             l.chrom,
             l.pos,
             l.ref_allele,
             l.alt_allele,
             l.variant_type,
             l.n_samples,
             t.n_total                                      AS total_samples,
             ROUND(l.n_samples * 1.0 / t.n_total, 4)       AS sample_fraction,
             l.mean_vaf,
             l.max_vaf,
             l.mean_depth,
             l.mean_strand_balance,
             l.max_strand_balance
         FROM loci l CROSS JOIN totals t
         WHERE {where_clause}
         ORDER BY sample_fraction DESC, chrom, pos
         LIMIT {top_n}",
        top_n = args.top_n,
    );

    struct Row {
        chrom:                String,
        pos:                  i64,
        ref_allele:           String,
        alt_allele:           String,
        variant_type:         String,
        n_samples:            i64,
        total_samples:        i64,
        sample_fraction:      f64,
        mean_vaf:             f64,
        max_vaf:              f64,
        mean_depth:           f64,
        mean_strand_balance:  f64,
        max_strand_balance:   f64,
    }

    let mut stmt = con.prepare(&display_sql).context("failed to prepare display query")?;
    let rows: Vec<Row> = stmt
        .query_map([], |row| {
            Ok(Row {
                chrom:               row.get(0)?,
                pos:                 row.get(1)?,
                ref_allele:          row.get(2)?,
                alt_allele:          row.get(3)?,
                variant_type:        row.get(4)?,
                n_samples:           row.get(5)?,
                total_samples:       row.get(6)?,
                sample_fraction:     row.get(7)?,
                mean_vaf:            row.get(8)?,
                max_vaf:             row.get(9)?,
                mean_depth:          row.get(10)?,
                mean_strand_balance: row.get(11)?,
                max_strand_balance:  row.get(12)?,
            })
        })?
        .collect::<Result<Vec<_>, _>>()
        .context("failed to fetch display rows")?;

    // ── Print stdout summary ───────────────────────────────────────────────────
    let total_samples: i64 = con
        .query_row(
            "SELECT COUNT(DISTINCT sample_id) FROM data",
            [],
            |r| r.get(0),
        )
        .context("failed to count total samples")?;

    println!(
        "Cohort: {} samples  |  min_samples={}  |  min_fraction={:.2}  |  showing top {}",
        total_samples, args.min_samples, args.min_sample_fraction, args.top_n
    );
    println!();
    println!(
        "{:<6}  {:>12}  {:<4}  {:<12}  {:<9}  {:>9}  {:>8}  {:>10}  {:>10}  {:>10}  {:>12}  {:>12}",
        "Chrom", "Pos", "Ref", "Alt", "Type",
        "Samples", "Frac", "Mean VAF", "Max VAF", "Mean Depth",
        "Mean Strand", "Max Strand",
    );
    println!("{}", "─".repeat(122));

    for r in &rows {
        println!(
            "{:<6}  {:>12}  {:<4}  {:<12}  {:<9}  {:>4}/{:<4}  {:>8.4}  {:>10.6}  {:>10.6}  {:>10.1}  {:>12.4}  {:>12.4}",
            r.chrom, r.pos, r.ref_allele, r.alt_allele, r.variant_type,
            r.n_samples, r.total_samples, r.sample_fraction,
            r.mean_vaf, r.max_vaf, r.mean_depth,
            r.mean_strand_balance, r.max_strand_balance,
        );
    }

    if rows.is_empty() {
        println!("(no loci meet the filter criteria)");
    }

    // ── Write output file (includes annotation columns) ────────────────────────
    let full_query = format!(
        "{cte}
         SELECT
             l.*,
             t.n_total                                AS total_samples,
             ROUND(l.n_samples * 1.0 / t.n_total, 4) AS sample_fraction
         FROM loci l CROSS JOIN totals t
         WHERE {where_clause}
         ORDER BY sample_fraction DESC, chrom, pos",
    );

    let output_str = args.output.display().to_string().replace('\'', "''");
    let fmt = match args.output.extension().and_then(|e| e.to_str()).unwrap_or("") {
        "parquet" => "FORMAT PARQUET".to_string(),
        _         => "FORMAT CSV, DELIMITER '\\t', HEADER true".to_string(),
    };

    con.execute_batch(&format!("COPY ({full_query}) TO '{output_str}' ({fmt})"))
        .with_context(|| format!("failed to write output: {}", args.output.display()))?;

    let n_loci: i64 = con
        .query_row(
            &format!("SELECT COUNT(*) FROM ({full_query})"),
            [],
            |r| r.get(0),
        )
        .context("failed to count output loci")?;

    println!();
    println!("{} loci written to {}", n_loci, args.output.display());

    Ok(())
}
