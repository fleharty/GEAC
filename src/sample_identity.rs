use std::collections::HashMap;
use std::io::Write;

use anyhow::{bail, Context, Result};
use duckdb::Connection;
use rayon::prelude::*;
use tracing::info;

use crate::cli::SampleIdentityArgs;

#[derive(Debug)]
struct SampleFingerprint {
    sample_id: String,
    subject_id: Option<String>,
    present: Vec<u64>,
    hom_alt: Vec<u64>,
    n_nonref: u32,
}

#[derive(Debug)]
struct IdentityPair {
    sample_a: String,
    sample_b: String,
    subject_a: Option<String>,
    subject_b: Option<String>,
    n_nonref_a: u32,
    n_nonref_b: u32,
    n_shared: u32,
    n_shared_concordant: u32,
    jaccard: f64,
    concordance: f64,
    flag: String,
}

pub fn run_sample_identity(args: &SampleIdentityArgs) -> Result<()> {
    validate_args(args)?;

    let conn = Connection::open(&args.input)
        .with_context(|| format!("failed to open cohort DuckDB: {}", args.input.display()))?;
    ensure_alt_bases(&conn)?;

    let cols = table_columns(&conn, "alt_bases")?;
    let has_subject = cols.iter().any(|c| c == "subject_id");
    let use_gnomad = args.common_gnomad_only && cols.iter().any(|c| c == "gnomad_af");
    if args.common_gnomad_only && !use_gnomad {
        bail!("--common-gnomad-only requested, but alt_bases has no gnomad_af column");
    }

    let (samples, n_markers) = load_fingerprints(&conn, args, has_subject, use_gnomad)?;
    if samples.len() < 2 {
        bail!(
            "need at least two samples with usable fingerprint markers; found {}",
            samples.len()
        );
    }
    if args.max_samples > 0 && samples.len() > args.max_samples {
        bail!(
            "sample identity would compare {} samples, above --max-samples {}; raise --max-samples or prefilter the cohort",
            samples.len(),
            args.max_samples
        );
    }

    info!(
        n_samples = samples.len(),
        n_markers, "comparing sample identity fingerprints"
    );

    let mut pairs = compare_fingerprints(&samples, args);
    pairs.sort_by(|a, b| {
        b.jaccard
            .partial_cmp(&a.jaccard)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                b.concordance
                    .partial_cmp(&a.concordance)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| a.sample_a.cmp(&b.sample_a))
            .then_with(|| a.sample_b.cmp(&b.sample_b))
    });

    write_pairs_tsv(&args.output, &pairs)?;
    info!(
        n_pairs = pairs.len(),
        output = %args.output.display(),
        "sample identity output written"
    );
    Ok(())
}

fn validate_args(args: &SampleIdentityArgs) -> Result<()> {
    if args.min_depth < 1 {
        bail!("--min-depth must be >= 1");
    }
    if !(0.0..=1.0).contains(&args.min_vaf) {
        bail!("--min-vaf must be between 0 and 1");
    }
    if !(0.0..=1.0).contains(&args.hom_vaf) {
        bail!("--hom-vaf must be between 0 and 1");
    }
    if args.min_vaf >= args.hom_vaf {
        bail!("--min-vaf must be lower than --hom-vaf");
    }
    if args.min_recurrence < 1 {
        bail!("--min-recurrence must be >= 1");
    }
    if args.max_markers < 1 {
        bail!("--max-markers must be >= 1");
    }
    if !(0.0..=1.0).contains(&args.min_jaccard) {
        bail!("--min-jaccard must be between 0 and 1");
    }
    if !(0.0..=1.0).contains(&args.min_concordance) {
        bail!("--min-concordance must be between 0 and 1");
    }
    if !(0.0..=1.0).contains(&args.gnomad_min_af)
        || !(0.0..=1.0).contains(&args.gnomad_max_af)
        || args.gnomad_min_af > args.gnomad_max_af
    {
        bail!("gnomAD AF bounds must satisfy 0 <= --gnomad-min-af <= --gnomad-max-af <= 1");
    }
    Ok(())
}

fn ensure_alt_bases(conn: &Connection) -> Result<()> {
    let exists: i64 = conn
        .query_row(
            "SELECT COUNT(*) FROM information_schema.tables WHERE table_schema = 'main' AND table_name = 'alt_bases'",
            [],
            |row| row.get(0),
        )
        .context("failed to inspect DuckDB tables")?;
    if exists == 0 {
        bail!("input DuckDB does not contain an alt_bases table");
    }
    Ok(())
}

fn table_columns(conn: &Connection, table: &str) -> Result<Vec<String>> {
    let mut stmt = conn
        .prepare(&format!("DESCRIBE SELECT * FROM {table} LIMIT 0"))
        .with_context(|| format!("failed to describe {table} table"))?;
    stmt.query_map([], |row| row.get::<_, String>(0))?
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("failed to read {table} columns"))
}

fn load_fingerprints(
    conn: &Connection,
    args: &SampleIdentityArgs,
    has_subject: bool,
    use_gnomad: bool,
) -> Result<(Vec<SampleFingerprint>, usize)> {
    let subject_select = if has_subject {
        "ANY_VALUE(subject_id) AS subject_id"
    } else {
        "NULL::VARCHAR AS subject_id"
    };
    let gnomad_clause = if use_gnomad {
        format!(
            "AND gnomad_af BETWEEN {} AND {}",
            args.gnomad_min_af, args.gnomad_max_af
        )
    } else {
        String::new()
    };
    let sql = format!(
        "
        WITH geno AS (
            SELECT
                sample_id,
                {subject_select},
                chrom,
                pos,
                alt_allele,
                chrom || ':' || pos || ':' || alt_allele AS marker_id,
                CASE
                    WHEN alt_count * 1.0 / NULLIF(total_depth, 0) >= {hom_vaf}
                    THEN 2 ELSE 1
                END AS geno
            FROM alt_bases
            WHERE variant_type = 'SNV'
              AND total_depth >= {min_depth}
              AND alt_count * 1.0 / NULLIF(total_depth, 0) >= {min_vaf}
              {gnomad_clause}
            GROUP BY sample_id, chrom, pos, alt_allele, total_depth, alt_count
        ),
        marker_counts AS (
            SELECT chrom, pos, alt_allele, marker_id, COUNT(DISTINCT sample_id) AS n_samples
            FROM geno
            GROUP BY chrom, pos, alt_allele, marker_id
        ),
        markers AS (
            SELECT marker_id
            FROM marker_counts
            WHERE n_samples >= {min_recurrence}
            ORDER BY n_samples DESC, marker_id
            LIMIT {max_markers}
        )
        SELECT g.sample_id, g.subject_id, g.marker_id, MAX(g.geno) AS geno
        FROM geno g
        JOIN markers m USING (marker_id)
        GROUP BY g.sample_id, g.subject_id, g.marker_id
        ORDER BY g.sample_id, g.marker_id
        ",
        hom_vaf = args.hom_vaf,
        min_depth = args.min_depth,
        min_vaf = args.min_vaf,
        min_recurrence = args.min_recurrence,
        max_markers = args.max_markers,
    );

    let mut stmt = conn
        .prepare(&sql)
        .context("failed to prepare sample identity marker query")?;
    let rows = stmt
        .query_map([], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, Option<String>>(1)?,
                row.get::<_, String>(2)?,
                row.get::<_, i32>(3)?,
            ))
        })
        .context("failed to query sample identity markers")?;

    let mut marker_idx: HashMap<String, usize> = HashMap::new();
    let mut sample_idx: HashMap<String, usize> = HashMap::new();
    let mut pending: Vec<(usize, usize, i32)> = Vec::new();
    let mut sample_meta: Vec<(String, Option<String>)> = Vec::new();

    for row in rows {
        let (sample_id, subject_id, marker_id, geno) =
            row.context("failed to read sample identity marker row")?;
        let next_sample_idx = sample_meta.len();
        let sidx = *sample_idx.entry(sample_id.clone()).or_insert_with(|| {
            sample_meta.push((sample_id, subject_id));
            next_sample_idx
        });
        let next_marker_idx = marker_idx.len();
        let midx = *marker_idx.entry(marker_id).or_insert(next_marker_idx);
        pending.push((sidx, midx, geno));
    }

    let n_markers = marker_idx.len();
    if n_markers == 0 {
        bail!("no usable sample identity markers under the current thresholds");
    }
    let n_words = (n_markers + 63) / 64;
    let mut samples = sample_meta
        .into_iter()
        .map(|(sample_id, subject_id)| SampleFingerprint {
            sample_id,
            subject_id,
            present: vec![0; n_words],
            hom_alt: vec![0; n_words],
            n_nonref: 0,
        })
        .collect::<Vec<_>>();

    for (sidx, midx, geno) in pending {
        let word = midx / 64;
        let bit = 1_u64 << (midx % 64);
        let sample = &mut samples[sidx];
        if sample.present[word] & bit == 0 {
            sample.present[word] |= bit;
            sample.n_nonref += 1;
        }
        if geno >= 2 {
            sample.hom_alt[word] |= bit;
        }
    }

    samples.retain(|s| s.n_nonref > 0);
    samples.sort_by(|a, b| a.sample_id.cmp(&b.sample_id));
    Ok((samples, n_markers))
}

fn compare_fingerprints(
    samples: &[SampleFingerprint],
    args: &SampleIdentityArgs,
) -> Vec<IdentityPair> {
    (0..samples.len())
        .into_par_iter()
        .map(|i| {
            let mut local = Vec::new();
            for j in (i + 1)..samples.len() {
                if let Some(pair) = compare_pair(&samples[i], &samples[j], args) {
                    local.push(pair);
                }
            }
            local
        })
        .reduce(Vec::new, |mut acc, mut local| {
            acc.append(&mut local);
            acc
        })
}

fn compare_pair(
    a: &SampleFingerprint,
    b: &SampleFingerprint,
    args: &SampleIdentityArgs,
) -> Option<IdentityPair> {
    let mut n_shared = 0_u32;
    let mut n_shared_concordant = 0_u32;

    for idx in 0..a.present.len() {
        let shared = a.present[idx] & b.present[idx];
        if shared == 0 {
            continue;
        }
        n_shared += shared.count_ones();

        let a_hom = a.hom_alt[idx];
        let b_hom = b.hom_alt[idx];
        let hom_match = shared & a_hom & b_hom;
        let het_match = shared & !a_hom & !b_hom;
        n_shared_concordant += hom_match.count_ones() + het_match.count_ones();
    }

    if n_shared == 0 {
        return None;
    }
    let denom = a.n_nonref + b.n_nonref - n_shared;
    if denom == 0 {
        return None;
    }
    let jaccard = n_shared as f64 / denom as f64;
    let concordance = n_shared_concordant as f64 / n_shared as f64;
    let passes = jaccard >= args.min_jaccard && concordance >= args.min_concordance;
    if !args.all_pairs && !passes {
        return None;
    }

    let flag = if passes {
        match (&a.subject_id, &b.subject_id) {
            (Some(sa), Some(sb)) if sa == sb => "EXPECTED_MATCH",
            _ => "SAME_INDIVIDUAL",
        }
    } else {
        ""
    }
    .to_string();

    Some(IdentityPair {
        sample_a: a.sample_id.clone(),
        sample_b: b.sample_id.clone(),
        subject_a: a.subject_id.clone(),
        subject_b: b.subject_id.clone(),
        n_nonref_a: a.n_nonref,
        n_nonref_b: b.n_nonref,
        n_shared,
        n_shared_concordant,
        jaccard,
        concordance,
        flag,
    })
}

fn write_pairs_tsv(path: &std::path::Path, pairs: &[IdentityPair]) -> Result<()> {
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(path)
            .with_context(|| format!("failed to create output TSV: {}", path.display()))?,
    );
    writeln!(
        out,
        "sample_a\tsubject_a\tsample_b\tsubject_b\tn_nonref_a\tn_nonref_b\tn_shared\tn_shared_concordant\tjaccard\tconcordance\tflag"
    )?;
    for p in pairs {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}",
            p.sample_a,
            p.subject_a.as_deref().unwrap_or(""),
            p.sample_b,
            p.subject_b.as_deref().unwrap_or(""),
            p.n_nonref_a,
            p.n_nonref_b,
            p.n_shared,
            p.n_shared_concordant,
            p.jaccard,
            p.concordance,
            p.flag
        )?;
    }
    Ok(())
}
