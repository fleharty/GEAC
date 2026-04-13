use std::collections::HashSet;
use std::time::Instant;

use anyhow::{Context, Result};
use rust_htslib::bam::Read;
use tracing::info;

use crate::cli::CollectArgs;
use crate::gene_annotations::GeneAnnotations;
use crate::gnomad::GnomadIndex;
use crate::record::{RefBase, RefRead};
use crate::repeat::compute_repeat_metrics;
use crate::targets::TargetIntervals;

use super::compute_input_sha256;
use super::pileup::{tally_pileup, true_cycle};
use super::ref_utils::{open_bam, read_group_sample_id, RefCache};

/// Collect reference-site records for target positions where `alt_positions` has no entry.
///
/// For each target position not already covered by an alt record, the BAM is pileup-fetched
/// and a `RefBase` locus record is emitted (if depth > 0).  Per-read `RefRead` records are
/// always emitted alongside (equivalent to `--reads-output` for alt records).
///
/// The caller is responsible for ensuring `--targets` was provided; this function
/// panics if `target_intervals` is `None`.
pub fn collect_ref_bases(
    args: &CollectArgs,
    alt_positions: &HashSet<(String, i64)>,
    target_intervals: &TargetIntervals,
    gene_annots: Option<&GeneAnnotations>,
    mut gnomad: Option<&mut GnomadIndex>,
) -> Result<(Vec<RefBase>, Vec<RefRead>)> {
    let input_checksum_sha256: Option<String> = if args.input_checksum_sha256 {
        Some(compute_input_sha256(&args.input)?)
    } else {
        None
    };
    let mut bam = open_bam(&args.input, &args.reference)?;
    let mut ref_cache = RefCache::new(&args.reference)?;

    let sample_id = match &args.sample_id {
        Some(id) => id.clone(),
        None => read_group_sample_id(bam.header()).context(
            "--sample-id was not provided and no SM tag found in BAM/CRAM header @RG line",
        )?,
    };

    let bam_targets: Vec<(String, usize)> = {
        let header = bam.header();
        (0..header.target_count())
            .map(|tid| {
                let name = std::str::from_utf8(header.tid2name(tid))
                    .unwrap_or("unknown")
                    .to_string();
                let len = header.target_len(tid).unwrap_or(0) as usize;
                (name, len)
            })
            .collect()
    };

    let n_target_positions = target_intervals.total_bases();
    let n_alt_positions = alt_positions.len();
    let n_ref_positions = n_target_positions.saturating_sub(n_alt_positions);
    info!(
        n_target_positions,
        n_alt_positions, n_ref_positions, "collecting ref-site records"
    );

    let mut ref_bases: Vec<RefBase> = Vec::new();
    let mut ref_reads: Vec<RefRead> = Vec::new();

    let start = Instant::now();
    let report_interval = args.progress_interval;
    let mut last_report = Instant::now();
    let mut n_visited: u64 = 0;

    target_intervals.for_each_position(|chrom, pos| {
        if alt_positions.contains(&(chrom.to_string(), pos)) {
            return;
        }

        n_visited += 1;

        if report_interval > 0 && last_report.elapsed().as_secs() >= report_interval {
            let elapsed = start.elapsed().as_secs_f64();
            eprintln!(
                "[{elapsed:>6.0}s] ref-sites: {n_visited}/{n_ref_positions} non-alt positions processed (current: {chrom}:{})",
                pos + 1
            );
            last_report = Instant::now();
        }

        // htslib region strings are 1-based inclusive; our positions are 0-based.
        let region = format!("{}:{}-{}", chrom, pos + 1, pos + 1);
        if bam.fetch(region.as_str()).is_err() {
            // Chromosome absent from BAM index — skip silently.
            return;
        }

        // Look up the reference base.
        let tid = bam_targets.iter().position(|(n, _)| n == chrom);
        let tid = match tid {
            Some(t) => t,
            None => return,
        };
        let ref_base = match ref_cache.get(&bam_targets, tid, pos as usize) {
            Ok((_, b)) => b,
            Err(_) => return,
        };
        if ref_base == 'N' {
            return;
        }

        let mut found = false;
        for pileup_result in bam.pileup() {
            let pileup = match pileup_result {
                Ok(p) => p,
                Err(_) => break,
            };
            if pileup.pos() as i64 != pos {
                continue;
            }
            found = true;

            let result = tally_pileup(
                &pileup,
                args.pipeline,
                args.min_base_qual,
                args.min_map_qual,
                args.include_duplicates,
                args.include_secondary,
                args.include_supplementary,
                ref_base,
                false,
                true, // collect_all_reads — captures every non-N read
            );

            if result.total_depth == 0 {
                break;
            }

            let ref_tally = result.bases.get(&ref_base);
            let ref_count = ref_tally.map_or(0, |t| t.total);
            let fwd_ref_count = ref_tally.map_or(0, |t| t.fwd);
            let rev_ref_count = ref_tally.map_or(0, |t| t.rev);
            let overlap_ref_agree = ref_tally.map_or(0, |t| t.overlap_alt_agree);

            let on_target = Some(target_intervals.contains(chrom, pos));
            let gene = gene_annots.and_then(|g| g.get(chrom, pos)).map(|a| a.gene);
            let repeat =
                compute_repeat_metrics(ref_cache.current_seq(), pos as usize, args.repeat_window);

            let gnomad_af = if let Some(ref mut g) = gnomad {
                g.get(chrom, pos, &ref_base.to_string(), &ref_base.to_string())
                    .unwrap_or(None)
            } else {
                None
            };

            let homopolymer_len = if repeat.homopolymer_len > 0 { Some(repeat.homopolymer_len) } else { None };
            let str_period = if repeat.str_period > 0 { Some(repeat.str_period) } else { None };
            let str_len = if repeat.str_len > 0 { Some(repeat.str_len) } else { None };

            ref_bases.push(RefBase {
                sample_id: sample_id.clone(),
                chrom: chrom.to_string(),
                pos,
                ref_allele: ref_base.to_string(),
                total_depth: result.total_depth,
                fwd_depth: result.fwd_depth,
                rev_depth: result.rev_depth,
                ref_count,
                fwd_ref_count,
                rev_ref_count,
                overlap_depth: result.overlap_depth,
                overlap_ref_agree,
                read_type: args.read_type,
                pipeline: args.pipeline,
                batch: args.batch.clone(),
                label1: args.label1.clone(),
                label2: args.label2.clone(),
                label3: args.label3.clone(),
                on_target,
                gene,
                homopolymer_len,
                str_period,
                str_len,
                gnomad_af,
                input_checksum_sha256: input_checksum_sha256.clone(),
            });

            // Emit one RefRead per passing read (all bases, including ref).
            for details_vec in result.read_details.values() {
                for detail in details_vec {
                    ref_reads.push(RefRead {
                        sample_id: sample_id.clone(),
                        chrom: chrom.to_string(),
                        pos,
                        read_type: args.read_type,
                        pipeline: args.pipeline,
                        cycle: true_cycle(
                            detail.qpos,
                            detail.read_len,
                            detail.is_reverse,
                            detail.hard_clip_before,
                        ),
                        read_length: detail.read_len as i32,
                        is_read1: detail.is_first_in_pair,
                        ab_count: detail.ab_count,
                        ba_count: detail.ba_count,
                        family_size: detail.family_size,
                        base_qual: detail.base_qual as i32,
                        map_qual: detail.map_qual as i32,
                        insert_size: detail.insert_size,
                        n_before_alt: detail.n_before_alt as i32,
                        n_after_alt: detail.n_after_alt as i32,
                        n_n_before_alt: detail.n_n_before_alt as i32,
                        n_n_after_alt: detail.n_n_after_alt as i32,
                        leading_n_run_len: detail.leading_n_run_len as i32,
                        trailing_n_run_len: detail.trailing_n_run_len as i32,
                        input_checksum_sha256: input_checksum_sha256.clone(),
                    });
                }
            }

            break;
        }
        let _ = found; // zero-depth positions are simply skipped
    });

    let elapsed = start.elapsed().as_secs_f64();
    eprintln!("[done] ref-sites: elapsed={elapsed:.1}s");
    info!(
        n_ref_bases = ref_bases.len(),
        n_ref_reads = ref_reads.len(),
        "ref-site collection complete"
    );

    Ok((ref_bases, ref_reads))
}
