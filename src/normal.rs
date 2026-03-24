use std::collections::HashMap;

use anyhow::{Context, Result};
use duckdb::Connection;
use rust_htslib::bam::{self, Read};
use tracing::info;

use crate::bam::{open_bam, read_group_sample_id};
use crate::cli::AnnotateNormalArgs;
use crate::record::NormalEvidence;

/// Cross-annotate tumor alt-base loci against a paired normal BAM/CRAM.
///
/// Reads the tumor locus Parquet to get the set of (chrom, pos, ref_allele, alt_allele)
/// tuples, then pileups the normal BAM at each position.  For every tumor locus an anchor
/// row (normal_alt_allele = None) is always emitted; additional rows are appended for each
/// non-ref allele observed in the normal at SNV positions.
pub fn annotate_normal(args: &AnnotateNormalArgs) -> Result<Vec<NormalEvidence>> {
    // ── Step 1: load tumor loci from Parquet ─────────────────────────────────
    let db = Connection::open_in_memory().context("failed to open in-memory DuckDB")?;

    let tumor_path = args.tumor_parquet.display().to_string().replace('\'', "''");
    db.execute_batch(&format!(
        "CREATE TABLE tumor AS SELECT * FROM read_parquet('{tumor_path}');"
    ))
    .context("failed to load tumor Parquet")?;

    let tumor_sample_id: String = db
        .query_row("SELECT ANY_VALUE(sample_id) FROM tumor", [], |row| row.get(0))
        .context("failed to read tumor sample_id from Parquet")?;

    info!(tumor_sample_id = %tumor_sample_id, "loaded tumor Parquet");

    // Query distinct (chrom, pos, ref_allele, alt_allele) tuples.
    // ref_allele is used to filter normal pileup bases to non-ref only.
    let mut stmt = db
        .prepare(
            "SELECT DISTINCT chrom, CAST(pos AS BIGINT), ref_allele, alt_allele \
             FROM tumor ORDER BY chrom, pos",
        )
        .context("failed to prepare tumor loci query")?;

    // targets: chrom -> pos -> (ref_allele, Vec<tumor_alt_allele>)
    let mut targets: HashMap<String, HashMap<i64, (String, Vec<String>)>> = HashMap::new();

    let rows = stmt
        .query_map([], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, i64>(1)?,
                row.get::<_, String>(2)?,
                row.get::<_, String>(3)?,
            ))
        })
        .context("failed to query tumor loci")?;

    let mut n_loci = 0usize;
    for row in rows {
        let (chrom, pos, ref_allele, alt_allele) = row.context("failed to read tumor locus row")?;
        let entry = targets
            .entry(chrom)
            .or_default()
            .entry(pos)
            .or_insert_with(|| (ref_allele, Vec::new()));
        entry.1.push(alt_allele);
        n_loci += 1;
    }

    info!(n_loci, n_chroms = targets.len(), "loaded tumor loci");

    // ── Step 2: open normal BAM ───────────────────────────────────────────────
    let mut bam = open_bam(&args.normal_bam, &args.reference)?;

    let normal_sample_id = match &args.normal_sample_id {
        Some(id) => id.clone(),
        None => read_group_sample_id(bam.header()).context(
            "--normal-sample-id was not provided and no SM tag found in normal BAM/CRAM header",
        )?,
    };

    info!(normal_sample_id = %normal_sample_id, "opened normal BAM");

    // ── Step 3: pileup the normal BAM at each tumor locus ────────────────────
    let mut results: Vec<NormalEvidence> = Vec::new();

    let mut sorted_chroms: Vec<String> = targets.keys().cloned().collect();
    sorted_chroms.sort();

    for chrom in &sorted_chroms {
        let pos_map = &targets[chrom];
        let mut sorted_positions: Vec<i64> = pos_map.keys().copied().collect();
        sorted_positions.sort();

        for &pos in &sorted_positions {
            let (ref_allele, tumor_alt_alleles) = &pos_map[&pos];

            // Fetch a 1-bp window around this position.
            // htslib region strings are 1-based inclusive.
            let region = format!("{}:{}-{}", chrom, pos + 1, pos + 1);
            let found = match bam.fetch(region.as_str()) {
                Ok(()) => {
                    pileup_position(
                        &mut bam,
                        chrom,
                        pos,
                        ref_allele,
                        tumor_alt_alleles,
                        &tumor_sample_id,
                        &normal_sample_id,
                        args.min_base_qual,
                        args.min_map_qual,
                        args.include_duplicates,
                        args.include_secondary,
                        args.include_supplementary,
                        &mut results,
                    )?
                }
                // If the chromosome is not in the normal BAM index (e.g. the normal
                // was aligned to a different reference), treat as zero-depth.
                Err(_) => false,
            };

            if !found {
                // No coverage — emit anchor rows only with depth = 0.
                for tumor_alt in tumor_alt_alleles {
                    results.push(NormalEvidence {
                        tumor_sample_id: tumor_sample_id.clone(),
                        chrom: chrom.clone(),
                        pos,
                        tumor_alt_allele: tumor_alt.clone(),
                        normal_sample_id: normal_sample_id.clone(),
                        normal_alt_allele: None,
                        normal_depth: 0,
                        normal_alt_count: 0,
                    });
                }
            }
        }
    }

    info!(n_records = results.len(), "normal annotation complete");
    Ok(results)
}

/// Pileup a single position in the normal BAM and append NormalEvidence records.
///
/// Returns `true` if the position was found in the pileup (has any coverage),
/// `false` otherwise.
#[allow(clippy::too_many_arguments)]
fn pileup_position(
    bam: &mut bam::IndexedReader,
    chrom: &str,
    pos: i64,
    ref_allele: &str,
    tumor_alt_alleles: &[String],
    tumor_sample_id: &str,
    normal_sample_id: &str,
    min_base_qual: u8,
    min_map_qual: u8,
    include_duplicates: bool,
    include_secondary: bool,
    include_supplementary: bool,
    results: &mut Vec<NormalEvidence>,
) -> Result<bool> {
    // SNV positions have a single-char ref allele.  We only emit per-allele rows
    // for SNV positions; indel positions get the anchor row only.
    let is_snv_pos = ref_allele.len() == 1;
    let ref_base: char = ref_allele.chars().next().unwrap_or('N');

    let mut found = false;

    for pileup_result in bam.pileup() {
        let pileup = pileup_result.context("error reading pileup from normal BAM")?;
        if pileup.pos() as i64 != pos {
            continue;
        }

        found = true;

        let (normal_depth, base_counts) = tally_position(
            &pileup,
            min_base_qual,
            min_map_qual,
            include_duplicates,
            include_secondary,
            include_supplementary,
        );

        for tumor_alt in tumor_alt_alleles {
            // Anchor row — always written.
            results.push(NormalEvidence {
                tumor_sample_id: tumor_sample_id.to_owned(),
                chrom: chrom.to_owned(),
                pos,
                tumor_alt_allele: tumor_alt.clone(),
                normal_sample_id: normal_sample_id.to_owned(),
                normal_alt_allele: None,
                normal_depth,
                normal_alt_count: 0,
            });

            // Non-ref allele rows — SNV positions only.
            if is_snv_pos {
                for (&base, &count) in &base_counts {
                    if base == ref_base || base == 'N' {
                        continue;
                    }
                    results.push(NormalEvidence {
                        tumor_sample_id: tumor_sample_id.to_owned(),
                        chrom: chrom.to_owned(),
                        pos,
                        tumor_alt_allele: tumor_alt.clone(),
                        normal_sample_id: normal_sample_id.to_owned(),
                        normal_alt_allele: Some(base.to_string()),
                        normal_depth,
                        normal_alt_count: count,
                    });
                }
            }
        }

        break; // Only one pileup position per fetch
    }

    Ok(found)
}

/// Tally observed bases at a single pileup column, de-duplicating by fragment (query name).
///
/// Returns `(total_depth, HashMap<base, count>)`.
/// `total_depth` counts fragments that have a non-N base passing all filters.
fn tally_position(
    pileup: &rust_htslib::bam::pileup::Pileup,
    min_base_qual: u8,
    min_map_qual: u8,
    include_duplicates: bool,
    include_secondary: bool,
    include_supplementary: bool,
) -> (i32, HashMap<char, i32>) {
    // Collect one base per fragment query name (first read seen wins).
    let mut by_qname: HashMap<Vec<u8>, char> = HashMap::new();

    for alignment in pileup.alignments() {
        if alignment.is_del() || alignment.is_refskip() {
            continue;
        }

        let record = alignment.record();

        if (!include_duplicates && record.is_duplicate())
            || (!include_secondary && record.is_secondary())
            || (!include_supplementary && record.is_supplementary())
        {
            continue;
        }

        if record.mapq() < min_map_qual {
            continue;
        }

        let qpos = match alignment.qpos() {
            Some(p) => p,
            None => continue,
        };

        let base_qual = record.qual()[qpos];
        if base_qual < min_base_qual {
            continue;
        }

        let base = record.seq()[qpos].to_ascii_uppercase() as char;
        by_qname.entry(record.qname().to_vec()).or_insert(base);
    }

    let mut base_counts: HashMap<char, i32> = HashMap::new();
    let mut total_depth: i32 = 0;

    for base in by_qname.values() {
        if *base != 'N' {
            total_depth += 1;
            *base_counts.entry(*base).or_insert(0) += 1;
        }
    }

    (total_depth, base_counts)
}
