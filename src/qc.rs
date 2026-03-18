use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use duckdb::Connection;

use crate::cli::QcArgs;

const SBS6: [&str; 6] = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"];

struct OverallStats {
    sample_id:           String,
    n_snv:               i64,
    n_insertion:         i64,
    n_deletion:          i64,
    mean_depth:          f64,
    mean_vaf:            f64,
    strand_balance:      f64,
    overlap_concordance: Option<f64>,
}

struct SubstStats {
    n_loci:         i64,
    mean_vaf:       f64,
    strand_balance: f64,
}

pub fn run_qc(args: &QcArgs) -> Result<()> {
    let con = Connection::open_in_memory().context("failed to open DuckDB")?;

    // Build the parquet source expression.
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
    .context("failed to create data view — check that input files exist and are valid Parquet")?;

    // ── Overall stats per sample ───────────────────────────────────────────────
    let overall_sql = "
        SELECT
            sample_id,
            COUNT(*) FILTER (WHERE variant_type = 'SNV')        AS n_snv,
            COUNT(*) FILTER (WHERE variant_type = 'insertion')  AS n_insertion,
            COUNT(*) FILTER (WHERE variant_type = 'deletion')   AS n_deletion,
            ROUND(AVG(total_depth), 1)                          AS mean_depth,
            ROUND(AVG(alt_count * 1.0 / total_depth), 6)       AS mean_vaf,
            ROUND(
                AVG(fwd_alt_count * 1.0 / NULLIF(fwd_alt_count + rev_alt_count, 0)),
                4
            ) AS strand_balance,
            ROUND(
                SUM(overlap_alt_agree) * 1.0
                    / NULLIF(SUM(overlap_alt_agree + overlap_alt_disagree), 0),
                4
            ) AS overlap_concordance
        FROM data
        GROUP BY sample_id
        ORDER BY sample_id
    ";

    let mut stmt = con.prepare(overall_sql).context("failed to prepare overall query")?;
    let overall: Vec<OverallStats> = stmt
        .query_map([], |row| {
            Ok(OverallStats {
                sample_id:           row.get(0)?,
                n_snv:               row.get(1)?,
                n_insertion:         row.get(2)?,
                n_deletion:          row.get(3)?,
                mean_depth:          row.get(4)?,
                mean_vaf:            row.get(5)?,
                strand_balance:      row.get(6)?,
                overlap_concordance: row.get(7)?,
            })
        })?
        .collect::<Result<Vec<_>, _>>()
        .context("failed to fetch overall stats")?;

    if overall.is_empty() {
        anyhow::bail!("no records found — is the input empty or the on-target filter too strict?");
    }

    // ── Per-substitution stats per sample (SBS6, pyrimidine-normalised) ────────
    let subst_sql = "
        SELECT
            sample_id,
            CASE
                WHEN ref_allele IN ('C', 'T') THEN ref_allele || '>' || alt_allele
                WHEN ref_allele = 'G' AND alt_allele = 'T' THEN 'C>A'
                WHEN ref_allele = 'G' AND alt_allele = 'C' THEN 'C>G'
                WHEN ref_allele = 'G' AND alt_allele = 'A' THEN 'C>T'
                WHEN ref_allele = 'A' AND alt_allele = 'T' THEN 'T>A'
                WHEN ref_allele = 'A' AND alt_allele = 'G' THEN 'T>C'
                WHEN ref_allele = 'A' AND alt_allele = 'C' THEN 'T>G'
            END AS substitution,
            COUNT(*) AS n_loci,
            ROUND(AVG(alt_count * 1.0 / total_depth), 6) AS mean_vaf,
            ROUND(
                AVG(fwd_alt_count * 1.0 / NULLIF(fwd_alt_count + rev_alt_count, 0)),
                4
            ) AS strand_balance
        FROM data
        WHERE variant_type = 'SNV'
        GROUP BY sample_id, substitution
        HAVING substitution IS NOT NULL
        ORDER BY sample_id, substitution
    ";

    let mut stmt = con.prepare(subst_sql).context("failed to prepare substitution query")?;
    let subst_rows: Vec<(String, String, SubstStats)> = stmt
        .query_map([], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, String>(1)?,
                SubstStats {
                    n_loci:         row.get(2)?,
                    mean_vaf:       row.get(3)?,
                    strand_balance: row.get(4)?,
                },
            ))
        })?
        .collect::<Result<Vec<_>, _>>()
        .context("failed to fetch substitution stats")?;

    let mut subst_by_sample: HashMap<String, HashMap<String, SubstStats>> = HashMap::new();
    for (sample, subst, stats) in subst_rows {
        subst_by_sample
            .entry(sample)
            .or_default()
            .insert(subst, stats);
    }

    // ── Print human-readable report ────────────────────────────────────────────
    let sep = "═".repeat(64);
    let thin = "─".repeat(64);

    for s in &overall {
        println!("{sep}");
        println!("  GEAC QC  —  {}", s.sample_id);
        println!("{sep}");
        println!("  Overall");
        println!("    SNV loci            : {:>10}", s.n_snv);
        println!("    Insertion loci      : {:>10}", s.n_insertion);
        println!("    Deletion loci       : {:>10}", s.n_deletion);
        println!("    Mean depth          : {:>10.1}", s.mean_depth);
        println!("    Mean VAF            : {:>10.6}", s.mean_vaf);
        println!(
            "    Strand balance      : {:>10.4}  (0.5 = perfect)",
            s.strand_balance
        );
        match s.overlap_concordance {
            Some(c) => println!(
                "    Overlap concordance : {:>10.4}  (1.0 = perfect)",
                c
            ),
            None => println!(
                "    Overlap concordance :          —  (no overlapping fragment pairs)"
            ),
        }

        println!();
        println!("  SNV Substitution Profile (SBS6, pyrimidine-normalised)");
        println!("  {thin}");
        println!(
            "    {:<5}  {:>10}  {:>12}  {:>12}",
            "Subst", "Loci", "Mean VAF", "Strand bal."
        );
        println!("  {thin}");

        let empty = HashMap::new();
        let smap = subst_by_sample.get(&s.sample_id).unwrap_or(&empty);

        for sub in SBS6 {
            if let Some(ss) = smap.get(sub) {
                println!(
                    "    {:<5}  {:>10}  {:>12.6}  {:>12.4}",
                    sub, ss.n_loci, ss.mean_vaf, ss.strand_balance
                );
            } else {
                println!("    {:<5}  {:>10}  {:>12}  {:>12}", sub, 0, "—", "—");
            }
        }
        println!("{sep}");
        println!();
    }

    // ── Optional TSV output ────────────────────────────────────────────────────
    if let Some(out_path) = &args.output {
        write_tsv(&overall, &subst_by_sample, out_path)?;
        eprintln!("TSV written to {}", out_path.display());
    }

    Ok(())
}

fn write_tsv(
    overall: &[OverallStats],
    subst_by_sample: &HashMap<String, HashMap<String, SubstStats>>,
    path: &Path,
) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("failed to create TSV: {}", path.display()))?;

    // Header
    let mut header: Vec<String> = vec![
        "sample_id".into(),
        "n_snv".into(),
        "n_insertion".into(),
        "n_deletion".into(),
        "mean_depth".into(),
        "mean_vaf".into(),
        "strand_balance".into(),
        "overlap_concordance".into(),
    ];
    for sub in SBS6 {
        let prefix = sub.replace('>', "_");
        header.push(format!("{prefix}_n"));
        header.push(format!("{prefix}_mean_vaf"));
        header.push(format!("{prefix}_strand_balance"));
    }
    writeln!(f, "{}", header.join("\t"))?;

    let empty = HashMap::new();
    for s in overall {
        let oc = s
            .overlap_concordance
            .map_or("NA".to_string(), |v| format!("{v:.4}"));

        let mut row: Vec<String> = vec![
            s.sample_id.clone(),
            s.n_snv.to_string(),
            s.n_insertion.to_string(),
            s.n_deletion.to_string(),
            format!("{:.1}", s.mean_depth),
            format!("{:.6}", s.mean_vaf),
            format!("{:.4}", s.strand_balance),
            oc,
        ];

        let smap = subst_by_sample.get(&s.sample_id).unwrap_or(&empty);
        for sub in SBS6 {
            if let Some(ss) = smap.get(sub) {
                row.push(ss.n_loci.to_string());
                row.push(format!("{:.6}", ss.mean_vaf));
                row.push(format!("{:.4}", ss.strand_balance));
            } else {
                row.extend(["0".into(), "NA".into(), "NA".into()]);
            }
        }

        writeln!(f, "{}", row.join("\t"))?;
    }

    Ok(())
}
