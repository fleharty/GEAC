use std::io::BufRead;
use std::time::Instant;

use anyhow::{Context, Result};
use rust_htslib::bam::Read;
use tracing::info;

use crate::bam::{open_bam, read_group_sample_id};
use crate::cli::LociDepthArgs;
use crate::record::LocusDepthRecord;

/// Collect total depth at a fixed set of loci from a BAM/CRAM.
///
/// Reads a TSV produced by `geac export-loci` (columns: chrom, pos, 0-based),
/// pileups the BAM at each position, and returns one `LocusDepthRecord` per locus.
/// Loci with zero depth (e.g. missing chromosome in the BAM) are still emitted
/// with `total_depth = 0` so that the output set is exhaustive over the input loci.
pub fn collect_locus_depth(args: &LociDepthArgs) -> Result<Vec<LocusDepthRecord>> {
    // ── Step 1: load loci TSV ─────────────────────────────────────────────────
    let loci = load_loci(&args.loci)?;
    info!(n_loci = loci.len(), loci = %args.loci.display(), "loaded loci");

    // ── Step 2: open BAM ─────────────────────────────────────────────────────
    let mut bam = open_bam(&args.input, &args.reference)?;

    let sample_id = match &args.sample_id {
        Some(id) => id.clone(),
        None => read_group_sample_id(bam.header()).context(
            "--sample-id was not provided and no SM tag found in BAM/CRAM header @RG line",
        )?,
    };

    info!(sample_id = %sample_id, input = %args.input.display(), "opened BAM");

    // ── Step 3: pileup each locus ─────────────────────────────────────────────
    let mut results: Vec<LocusDepthRecord> = Vec::with_capacity(loci.len());
    let start = Instant::now();
    let report_interval = args.progress_interval;
    let mut last_report = Instant::now();

    for (i, (chrom, pos)) in loci.iter().enumerate() {
        // htslib region strings are 1-based inclusive; our positions are 0-based.
        let region = format!("{}:{}-{}", chrom, pos + 1, pos + 1);

        let (total_depth, fwd_depth, rev_depth) = match bam.fetch(region.as_str()) {
            Ok(()) => tally_position(
                &mut bam,
                *pos,
                args.min_base_qual,
                args.min_map_qual,
                args.include_duplicates,
                args.include_secondary,
                args.include_supplementary,
            )?,
            // Chromosome absent from index — emit zero depth.
            Err(_) => (0, 0, 0),
        };

        results.push(LocusDepthRecord {
            sample_id: sample_id.clone(),
            chrom: chrom.clone(),
            pos: *pos,
            total_depth,
            fwd_depth,
            rev_depth,
        });

        if report_interval > 0
            && last_report.elapsed().as_secs() >= report_interval
        {
            let elapsed = start.elapsed().as_secs_f64();
            eprintln!(
                "[{elapsed:>6.0}s] locus-depth: {i}/{} loci processed (current: {chrom}:{})",
                loci.len(),
                pos + 1
            );
            last_report = Instant::now();
        }
    }

    let elapsed = start.elapsed().as_secs_f64();
    eprintln!("[done] locus-depth: elapsed={elapsed:.1}s loci={}", loci.len());
    info!(n_records = results.len(), "locus-depth collection complete");
    Ok(results)
}

/// Parse a TSV produced by `geac export-loci`.
///
/// Expected header: `chrom\tpos` (0-based positions).
/// Lines starting with `#` are ignored.
fn load_loci(path: &std::path::Path) -> Result<Vec<(String, i64)>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("failed to open loci file: {}", path.display()))?;
    let reader = std::io::BufReader::new(file);

    let mut loci: Vec<(String, i64)> = Vec::new();
    let mut first = true;

    for (lineno, line_result) in reader.lines().enumerate() {
        let line = line_result
            .with_context(|| format!("error reading loci file at line {}", lineno + 1))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // Skip header row.
        if first {
            first = false;
            if line.starts_with("chrom") {
                continue;
            }
        }

        let mut fields = line.splitn(2, '\t');
        let chrom = fields
            .next()
            .with_context(|| format!("missing chrom field at line {}", lineno + 1))?
            .to_string();
        let pos_str = fields
            .next()
            .with_context(|| format!("missing pos field at line {}", lineno + 1))?;
        let pos: i64 = pos_str
            .parse()
            .with_context(|| format!("invalid pos '{}' at line {}", pos_str, lineno + 1))?;

        loci.push((chrom, pos));
    }

    Ok(loci)
}

/// Pileup a single position and tally total / forward / reverse depth.
///
/// Returns `(total_depth, fwd_depth, rev_depth)`.
/// De-duplicates by fragment query name (first read seen wins).
#[allow(clippy::too_many_arguments)]
fn tally_position(
    bam: &mut rust_htslib::bam::IndexedReader,
    pos: i64,
    min_base_qual: u8,
    min_map_qual: u8,
    include_duplicates: bool,
    include_secondary: bool,
    include_supplementary: bool,
) -> Result<(i32, i32, i32)> {
    use std::collections::HashMap;

    // qname -> is_reverse
    let mut by_qname: HashMap<Vec<u8>, bool> = HashMap::new();

    for pileup_result in bam.pileup() {
        let pileup = pileup_result.context("error reading pileup")?;
        if pileup.pos() as i64 != pos {
            continue;
        }

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

            let base = record.seq()[qpos].to_ascii_uppercase();
            if base == b'N' {
                continue;
            }

            by_qname
                .entry(record.qname().to_vec())
                .or_insert(record.is_reverse());
        }

        break; // One pileup column per fetch call
    }

    let mut total_depth: i32 = 0;
    let mut fwd_depth: i32 = 0;
    let mut rev_depth: i32 = 0;

    for is_rev in by_qname.values() {
        total_depth += 1;
        if *is_rev {
            rev_depth += 1;
        } else {
            fwd_depth += 1;
        }
    }

    Ok((total_depth, fwd_depth, rev_depth))
}
