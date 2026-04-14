mod builders;
mod indel;
mod pileup;
mod ref_collect;
mod ref_utils;

use std::sync::atomic::Ordering;
use std::time::Instant;
use std::{fs::File, io::Read as _, path::Path};

use anyhow::{Context, Result};
use rust_htslib::bam::Read;
use sha2::{Digest, Sha256};

use crate::cli::CollectArgs;
use crate::gene_annotations::GeneAnnotations;
use crate::gnomad::GnomadIndex;
use crate::progress::ProgressReporter;
use crate::record::{AltBase, AltRead, SampleMetricsRecord};
use crate::repeat::compute_repeat_metrics;
use crate::targets::TargetIntervals;
use crate::vcf::VariantAnnotator;

use builders::LocusContext;
use indel::tally_indels;
use pileup::{locus_n_context_summary, tally_pileup, PileupResult};

pub use ref_collect::collect_ref_bases;
pub(crate) use ref_utils::RefCache;
pub(crate) use ref_utils::gc_frac;
pub use ref_utils::{open_bam, read_group_sample_id};

/// Process a BAM/CRAM file and return all alt base records (and optionally per-read detail records).
///
/// When `args.reads_output` is true, the second element of the returned tuple contains one
/// `AltRead` record per read (fragment) that supports an alt base. When false, the second
/// element is always empty.
pub fn collect_alt_bases(
    args: &CollectArgs,
    annotator: Option<&dyn VariantAnnotator>,
    target_intervals: Option<&TargetIntervals>,
    gene_annots: Option<&GeneAnnotations>,
    mut gnomad: Option<&mut GnomadIndex>,
) -> Result<(Vec<AltBase>, Vec<AltRead>, Option<SampleMetricsRecord>)> {
    let input_checksum_sha256 = if args.input_checksum_sha256 {
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

    let targets: Vec<(String, usize)> = {
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

    if let Some(region) = &args.region {
        bam.fetch(region.as_str()).with_context(|| {
            format!(
                "failed to fetch region '{region}': check that the region is valid and the BAM is indexed"
            )
        })?;
    }

    let start = Instant::now();
    let (reporter, progress) = ProgressReporter::start(args.progress_interval);
    let collect_reads = args.reads_output;
    let mut records: Vec<AltBase> = Vec::new();
    let mut read_records: Vec<AltRead> = Vec::new();
    let region_filter = args.region.as_deref().and_then(parse_region_spec);
    let mut sample_metrics_acc = target_intervals.map(|_| SampleMetricsAccumulator::default());

    for pileup in bam.pileup() {
        let pileup = pileup.context("error reading pileup")?;
        let tid = pileup.tid() as usize;
        let pos = pileup.pos() as i64;

        let (chrom, ref_base) = ref_cache.get(&targets, tid, pos as usize)?;
        if ref_base == 'N' {
            continue;
        }

        let PileupResult {
            bases,
            total_depth,
            fwd_depth,
            rev_depth,
            overlap_depth,
            read_details,
        } = tally_pileup(
            &pileup,
            args.pipeline,
            args.min_base_qual,
            args.min_map_qual,
            args.include_duplicates,
            args.include_secondary,
            args.include_supplementary,
            ref_base,
            collect_reads,
            false,
        );

        progress.positions_processed.fetch_add(1, Ordering::Relaxed);
        progress
            .reads_processed
            .fetch_add(total_depth as u64, Ordering::Relaxed);
        progress.update_locus(&chrom, pos);

        if total_depth == 0 {
            continue;
        }

        let on_target = target_intervals.map(|t| t.contains(&chrom, pos));
        if let Some(acc) = sample_metrics_acc.as_mut() {
            acc.observe_position(
                &chrom,
                pos,
                total_depth,
                on_target == Some(true),
                region_filter.as_ref(),
            );
        }
        let gene = gene_annots.and_then(|g| g.get(&chrom, pos)).map(|a| a.gene);
        let repeat =
            compute_repeat_metrics(ref_cache.current_seq(), pos as usize, args.repeat_window);
        let trinuc_context = {
            let seq = ref_cache.current_seq();
            let p = pos as usize;
            if p > 0 && p + 1 < seq.len() {
                Some(format!(
                    "{}{}{}",
                    seq[p - 1] as char,
                    seq[p] as char,
                    seq[p + 1] as char
                ))
            } else {
                None
            }
        };

        let ref_tally = bases.get(&ref_base);
        let locus = LocusContext::new(
            args,
            &sample_id,
            &chrom,
            pos,
            ref_base,
            total_depth,
            fwd_depth,
            rev_depth,
            overlap_depth,
            ref_tally,
            on_target,
            gene,
            &repeat,
            trinuc_context,
            input_checksum_sha256.clone(),
        );

        for (base, tally) in &bases {
            if *base == ref_base || *base == 'N' || tally.total == 0 {
                continue;
            }

            progress.alt_bases_found.fetch_add(1, Ordering::Relaxed);

            let alt_allele = base.to_string();
            let (variant_called, variant_filter) =
                vcf_annotation(annotator, &chrom, pos, &alt_allele);
            let gnomad_af = if let Some(ref mut g) = gnomad {
                g.get(&chrom, pos, &ref_base.to_string(), &alt_allele)?
            } else {
                None
            };

            let n_ctx_summary = if collect_reads {
                read_details
                    .get(base)
                    .map(|details| locus_n_context_summary(details))
            } else {
                None
            };

            locus.push_snv_record(
                &mut records,
                *base,
                tally,
                variant_called,
                variant_filter,
                gnomad_af,
                n_ctx_summary,
            );

            if collect_reads {
                if let Some(details) = read_details.get(base) {
                    let seq = ref_cache.current_seq();
                    for detail in details {
                        let frag_gc = detail.insert_size.and_then(|ins| {
                            gc_frac(seq, detail.frag_start as usize, ins as usize)
                        });
                        read_records.push(locus.build_alt_read(&alt_allele, detail, frag_gc));
                    }
                }
            }
        }

        let (indels, indel_read_details) = tally_indels(
            &pileup,
            args.pipeline,
            pos,
            ref_cache.current_seq(),
            args.min_map_qual,
            args.include_duplicates,
            args.include_secondary,
            args.include_supplementary,
            collect_reads,
        );

        for indel in indels.values() {
            if indel.total == 0 {
                continue;
            }

            progress.alt_bases_found.fetch_add(1, Ordering::Relaxed);

            let (variant_called, variant_filter) =
                vcf_annotation(annotator, &chrom, pos, &indel.alt_allele);
            let gnomad_af = if let Some(ref mut g) = gnomad {
                g.get(&chrom, pos, &indel.ref_allele, &indel.alt_allele)?
            } else {
                None
            };

            let n_ctx_summary = if collect_reads {
                indel_read_details
                    .get(&indel.alt_allele)
                    .map(|details| locus_n_context_summary(details))
            } else {
                None
            };

            locus.push_indel_record(
                &mut records,
                indel,
                variant_called,
                variant_filter,
                gnomad_af,
                n_ctx_summary,
            );
        }

        if collect_reads {
            for (alt_allele, details) in &indel_read_details {
                let seq = ref_cache.current_seq();
                for detail in details {
                    let frag_gc = detail.insert_size.and_then(|ins| {
                        gc_frac(seq, detail.frag_start as usize, ins as usize)
                    });
                    read_records.push(locus.build_alt_read(alt_allele, detail, frag_gc));
                }
            }
        }
    }

    reporter.finish(start);

    let sample_metrics = if let (Some(ti), Some(acc)) = (target_intervals, sample_metrics_acc) {
        Some(acc.build(
            &sample_id,
            args.batch.clone(),
            args.read_type,
            args.pipeline,
            input_checksum_sha256,
            ti,
            region_filter.as_ref(),
        ))
    } else {
        None
    };

    Ok((records, read_records, sample_metrics))
}

#[derive(Debug, Clone)]
struct RegionSpec {
    chrom: String,
    start: i64,
    end: Option<i64>,
}

impl RegionSpec {
    fn contains(&self, chrom: &str, pos: i64) -> bool {
        if chrom != self.chrom || pos < self.start {
            return false;
        }
        self.end.is_none_or(|end| pos < end)
    }
}

fn parse_region_spec(region: &str) -> Option<RegionSpec> {
    let (chrom, coords) = region.split_once(':').unwrap_or((region, ""));
    if coords.is_empty() {
        return Some(RegionSpec {
            chrom: chrom.to_string(),
            start: 0,
            end: None,
        });
    }
    let (start_s, end_s) = coords.split_once('-').unwrap_or((coords, ""));
    let start = start_s
        .replace(',', "")
        .parse::<i64>()
        .ok()?
        .saturating_sub(1);
    let end = if end_s.is_empty() {
        None
    } else {
        Some(end_s.replace(',', "").parse::<i64>().ok()?)
    };
    Some(RegionSpec {
        chrom: chrom.to_string(),
        start,
        end,
    })
}

#[derive(Default)]
struct SampleMetricsAccumulator {
    total_fragment_bases: i64,
    on_target_fragment_bases: i64,
    covered_target_depths: Vec<i32>,
}

impl SampleMetricsAccumulator {
    fn observe_position(
        &mut self,
        chrom: &str,
        pos: i64,
        total_depth: i32,
        on_target: bool,
        region: Option<&RegionSpec>,
    ) {
        if region.is_some_and(|r| !r.contains(chrom, pos)) {
            return;
        }
        self.total_fragment_bases += i64::from(total_depth);
        if on_target {
            self.on_target_fragment_bases += i64::from(total_depth);
            self.covered_target_depths.push(total_depth);
        }
    }

    fn build(
        mut self,
        sample_id: &str,
        batch: Option<String>,
        read_type: crate::record::ReadType,
        pipeline: crate::record::Pipeline,
        input_checksum_sha256: Option<String>,
        target_intervals: &TargetIntervals,
        region: Option<&RegionSpec>,
    ) -> SampleMetricsRecord {
        let mut n_target_positions = 0_i32;
        target_intervals.for_each_position(|chrom, pos| {
            if region.is_none_or(|r| r.contains(chrom, pos)) {
                n_target_positions += 1;
            }
        });

        self.covered_target_depths.sort_unstable();
        let n_target_positions_covered = self.covered_target_depths.len() as i32;
        let zero_target_positions =
            n_target_positions.saturating_sub(n_target_positions_covered) as usize;
        let sum_target_depths: i64 = self
            .covered_target_depths
            .iter()
            .map(|&d| i64::from(d))
            .sum();

        let mean_target_depth_covered = if n_target_positions_covered > 0 {
            Some(sum_target_depths as f32 / n_target_positions_covered as f32)
        } else {
            None
        };
        let mean_target_depth_all = if n_target_positions > 0 {
            Some(sum_target_depths as f32 / n_target_positions as f32)
        } else {
            None
        };
        let median_target_depth_covered = median_sorted_i32(&self.covered_target_depths);
        let median_target_depth_all =
            median_sorted_with_leading_zeros(&self.covered_target_depths, zero_target_positions);
        let pct_fragment_bases_on_target = if self.total_fragment_bases > 0 {
            Some(self.on_target_fragment_bases as f32 / self.total_fragment_bases as f32)
        } else {
            None
        };

        SampleMetricsRecord {
            sample_id: sample_id.to_string(),
            batch,
            read_type,
            pipeline,
            input_checksum_sha256,
            n_target_positions,
            n_target_positions_covered,
            mean_target_depth_covered,
            mean_target_depth_all,
            median_target_depth_covered,
            median_target_depth_all,
            pct_fragment_bases_on_target,
        }
    }
}

fn median_sorted_i32(sorted: &[i32]) -> Option<f32> {
    if sorted.is_empty() {
        return None;
    }
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        Some(sorted[mid] as f32)
    } else {
        Some((sorted[mid - 1] as f32 + sorted[mid] as f32) / 2.0)
    }
}

fn median_sorted_with_leading_zeros(sorted_positive: &[i32], zero_count: usize) -> Option<f32> {
    let total = zero_count + sorted_positive.len();
    if total == 0 {
        return None;
    }

    fn value_at(sorted_positive: &[i32], zero_count: usize, idx: usize) -> f32 {
        if idx < zero_count {
            0.0
        } else {
            sorted_positive[idx - zero_count] as f32
        }
    }

    let mid = total / 2;
    if total % 2 == 1 {
        Some(value_at(sorted_positive, zero_count, mid))
    } else {
        Some(
            (value_at(sorted_positive, zero_count, mid - 1)
                + value_at(sorted_positive, zero_count, mid))
                / 2.0,
        )
    }
}

pub(super) fn compute_input_sha256(path: &Path) -> Result<String> {
    let mut file = File::open(path)
        .with_context(|| format!("failed to open input for SHA-256: {}", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = [0_u8; 64 * 1024];

    loop {
        let read_n = file
            .read(&mut buffer)
            .with_context(|| format!("failed to read input for SHA-256: {}", path.display()))?;
        if read_n == 0 {
            break;
        }
        hasher.update(&buffer[..read_n]);
    }

    Ok(format!("{:x}", hasher.finalize()))
}

fn vcf_annotation(
    annotator: Option<&dyn VariantAnnotator>,
    chrom: &str,
    pos: i64,
    alt_allele: &str,
) -> (Option<bool>, Option<String>) {
    match annotator {
        None => (None, None),
        Some(ann) => match ann.get(chrom, pos, alt_allele) {
            Some(a) => (Some(true), Some(a.filter.clone())),
            None => (Some(false), None),
        },
    }
}
