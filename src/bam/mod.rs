mod builders;
mod indel;
mod pileup;
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
use crate::record::{
    AltBase, AltRead, SampleMetricsPartialRecord, SampleMetricsRecord, VariantType,
};
use crate::region::{parse_region_input, RegionInput};
use crate::repeat::compute_repeat_metrics;
use crate::sample_metrics_math;
use crate::targets::TargetIntervals;
use crate::vcf::VariantAnnotator;

use builders::LocusContext;
use indel::tally_indels;
use pileup::{locus_n_context_summary, tally_pileup, PileupResult};

pub(crate) use ref_utils::gc_frac;
pub(crate) use ref_utils::resolve_max_pileup_depth;
pub(crate) use ref_utils::RefCache;
pub use ref_utils::{open_bam, read_group_sample_id};

/// Sample-metrics output of a collect run: the final record, the combinable
/// per-shard partial (when `--sample-metrics-partial` is set), or none (no targets).
pub enum SampleMetricsOutput {
    None,
    Final(SampleMetricsRecord),
    Partial(SampleMetricsPartialRecord),
}

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
) -> Result<(Vec<AltBase>, Vec<AltRead>, SampleMetricsOutput)> {
    let input_checksum_sha256 = if args.input_checksum_sha256 {
        Some(compute_input_sha256(&args.input)?)
    } else {
        None
    };

    let mut bam = open_bam(&args.input, &args.reference, args.index.as_deref())?;
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

    let exclude_tags = crate::bam::pileup::parse_exclude_tags(&args.exclude_tag)?;
    let family_size_scheme = args.family_size_scheme()?;

    let region_input = parse_region_input(args.region.as_deref())?;
    // Each fetch unit is fetched independently and its emitted pileup columns are
    // clipped to `clip`. For an interval list we fetch the *merged* (disjoint)
    // intervals one at a time and clip each to its own bounds, so boundary
    // over-emission from htslib (columns for reads spanning an interval edge)
    // never leaks into — or double-counts against — an adjacent interval. This is
    // what makes splitting one sample across interval shards exact.
    let fetch_units: Vec<FetchUnit> = match &region_input {
        None => vec![FetchUnit {
            fetch: None,
            clip: None,
        }],
        Some(RegionInput::Single(s)) => vec![FetchUnit {
            fetch: Some(s.clone()),
            clip: parse_region_spec(s),
        }],
        Some(RegionInput::Intervals(ivs)) => ivs
            .merged_intervals()
            .into_iter()
            .map(|(chrom, start, end)| FetchUnit {
                fetch: Some(format!("{}:{}-{}", chrom, start + 1, end)),
                clip: Some(RegionSpec {
                    chrom: chrom.to_string(),
                    start: i64::from(start),
                    end: Some(i64::from(end)),
                }),
            })
            .collect(),
    };
    // Region membership used when counting target positions for sample metrics:
    // the union of all requested intervals (or the single region).
    let metrics_region = match &region_input {
        None => MetricsRegion::All,
        Some(RegionInput::Single(s)) => parse_region_spec(s)
            .map(MetricsRegion::Single)
            .unwrap_or(MetricsRegion::All),
        Some(RegionInput::Intervals(t)) => MetricsRegion::Intervals(t),
    };

    let start = Instant::now();
    let (reporter, progress) = ProgressReporter::start(args.progress_interval);
    let collect_reads = args.reads_output;
    let mut records: Vec<AltBase> = Vec::new();
    let mut read_records: Vec<AltRead> = Vec::new();
    let mut sample_metrics_acc = target_intervals.map(|_| SampleMetricsAccumulator::default());

    for unit in &fetch_units {
        if let Some(r) = unit.fetch.as_deref() {
            bam.fetch(r).with_context(|| {
                format!(
                    "failed to fetch region '{r}': check that the region is valid and the BAM is indexed"
                )
            })?;
        }

        let mut plp = bam.pileup();
        plp.set_max_depth(resolve_max_pileup_depth(args.max_pileup_depth));
        for pileup in plp {
            let pileup = pileup.context("error reading pileup")?;
            let tid = pileup.tid() as usize;
            let pos = pileup.pos() as i64;

            let (chrom, ref_base) = ref_cache.get(&targets, tid, pos as usize)?;
            // Clip to the current fetch interval: htslib may emit columns just
            // outside it for reads spanning the edge. Skipping them keeps interval
            // shards a true partition (no boundary double-counting).
            if unit.clip.as_ref().is_some_and(|c| !c.contains(&chrom, pos)) {
                continue;
            }
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
                &family_size_scheme,
                args.min_base_qual,
                args.min_map_qual,
                args.include_duplicates,
                args.include_secondary,
                args.include_supplementary,
                &exclude_tags,
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

            // on_target for SNVs/insertions is evaluated at the anchor position.
            // Deletions use range overlap instead (see indel loop below).
            let on_target = target_intervals.map(|t| t.contains(&chrom, pos));
            if let Some(acc) = sample_metrics_acc.as_mut() {
                // Positions are already clipped to the fetch interval above, so the
                // accumulator observes only in-region columns.
                acc.observe_position(total_depth, on_target == Some(true));
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
                &family_size_scheme,
                pos,
                ref_cache.current_seq(),
                args.min_map_qual,
                args.include_duplicates,
                args.include_secondary,
                args.include_supplementary,
                &exclude_tags,
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
                    // gnomAD lookup needs the anchor base (single char), not the
                    // deleted bases that IndelCount.ref_allele happens to store for
                    // deletions. Pass ref_base — matches the SNV call above and the
                    // AltBase.ref_allele contract downstream.
                    g.get(&chrom, pos, &ref_base.to_string(), &indel.alt_allele)?
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

                // For deletions, check whether the deleted range overlaps any target
                // interval rather than testing only the anchor. alt_allele is "-<bases>",
                // so pos + alt_allele.len() is the exclusive end of the deletion span.
                let indel_on_target = if indel.variant_type == VariantType::Deletion {
                    target_intervals
                        .map(|t| t.overlaps(&chrom, pos, pos + indel.alt_allele.len() as i64))
                } else {
                    on_target
                };

                locus.push_indel_record(
                    &mut records,
                    indel,
                    indel_on_target,
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
                        let frag_gc = detail
                            .insert_size
                            .and_then(|ins| gc_frac(seq, detail.frag_start as usize, ins as usize));
                        read_records.push(locus.build_alt_read(alt_allele, detail, frag_gc));
                    }
                }
            }
        }
    }

    reporter.finish(start);

    let sample_metrics = match (target_intervals, sample_metrics_acc) {
        (Some(ti), Some(acc)) if args.sample_metrics_partial => {
            SampleMetricsOutput::Partial(acc.build_partial(
                &sample_id,
                args.subject_id.clone(),
                args.sample_type.clone(),
                args.batch.clone(),
                args.read_type,
                args.pipeline.clone(),
                input_checksum_sha256,
                ti,
                &metrics_region,
            ))
        }
        (Some(ti), Some(acc)) => SampleMetricsOutput::Final(acc.build(
            &sample_id,
            args.subject_id.clone(),
            args.sample_type.clone(),
            args.batch.clone(),
            args.read_type,
            args.pipeline.clone(),
            input_checksum_sha256,
            ti,
            &metrics_region,
        )),
        _ => SampleMetricsOutput::None,
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

/// One BAM fetch plus the clip applied to its emitted pileup columns.
struct FetchUnit {
    /// htslib region string to fetch, or `None` to scan the whole BAM.
    fetch: Option<String>,
    /// Restrict emitted positions to this interval, or `None` for no clip.
    clip: Option<RegionSpec>,
}

/// Region membership test for counting target positions in sample metrics.
/// This is the union of all requested intervals (or the single region), as
/// opposed to the per-fetch `RegionSpec` clip used while emitting loci.
enum MetricsRegion<'a> {
    All,
    Single(RegionSpec),
    Intervals(&'a TargetIntervals),
}

impl MetricsRegion<'_> {
    fn contains(&self, chrom: &str, pos: i64) -> bool {
        match self {
            MetricsRegion::All => true,
            MetricsRegion::Single(r) => r.contains(chrom, pos),
            MetricsRegion::Intervals(t) => t.contains(chrom, pos),
        }
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
    /// Observe one pileup column. Callers must pre-clip positions to the region;
    /// the accumulator does no region filtering itself.
    fn observe_position(&mut self, total_depth: i32, on_target: bool) {
        self.total_fragment_bases += i64::from(total_depth);
        if on_target {
            self.on_target_fragment_bases += i64::from(total_depth);
            self.covered_target_depths.push(total_depth);
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn build(
        self,
        sample_id: &str,
        subject_id: Option<String>,
        sample_type: Option<String>,
        batch: Option<String>,
        read_type: crate::record::ReadType,
        pipeline: Option<String>,
        input_checksum_sha256: Option<String>,
        target_intervals: &TargetIntervals,
        region: &MetricsRegion,
    ) -> SampleMetricsRecord {
        let n_target_positions = count_target_positions(target_intervals, region);
        let hist = sample_metrics_math::histogram_from_depths(&self.covered_target_depths);
        let stats = sample_metrics_math::depth_stats(
            &hist,
            n_target_positions,
            self.total_fragment_bases,
            self.on_target_fragment_bases,
        );

        SampleMetricsRecord {
            sample_id: sample_id.to_string(),
            subject_id,
            sample_type,
            batch,
            read_type,
            pipeline,
            input_checksum_sha256,
            n_target_positions,
            n_target_positions_covered: stats.n_target_positions_covered,
            mean_target_depth_covered: stats.mean_target_depth_covered,
            mean_target_depth_all: stats.mean_target_depth_all,
            median_target_depth_covered: stats.median_target_depth_covered,
            median_target_depth_all: stats.median_target_depth_all,
            pct_fragment_bases_on_target: stats.pct_fragment_bases_on_target,
        }
    }

    /// Emit combinable sufficient statistics for this shard instead of final
    /// metrics. Used by `geac collect --sample-metrics-partial`; merged later by
    /// `geac aggregate-metrics`.
    #[allow(clippy::too_many_arguments)]
    fn build_partial(
        self,
        sample_id: &str,
        subject_id: Option<String>,
        sample_type: Option<String>,
        batch: Option<String>,
        read_type: crate::record::ReadType,
        pipeline: Option<String>,
        input_checksum_sha256: Option<String>,
        target_intervals: &TargetIntervals,
        region: &MetricsRegion,
    ) -> SampleMetricsPartialRecord {
        let n_target_positions = count_target_positions(target_intervals, region);
        let hist = sample_metrics_math::histogram_from_depths(&self.covered_target_depths);
        let (hist_depth, hist_count): (Vec<i32>, Vec<i64>) = hist.into_iter().unzip();

        SampleMetricsPartialRecord {
            sample_id: sample_id.to_string(),
            subject_id,
            sample_type,
            batch,
            read_type,
            pipeline,
            input_checksum_sha256,
            n_target_positions,
            total_fragment_bases: self.total_fragment_bases,
            on_target_fragment_bases: self.on_target_fragment_bases,
            hist_depth,
            hist_count,
        }
    }
}

/// Count target positions that fall within the metrics region (the union of all
/// requested intervals, or the whole targets file when unrestricted).
fn count_target_positions(target_intervals: &TargetIntervals, region: &MetricsRegion) -> i64 {
    let mut n = 0_i64;
    target_intervals.for_each_position(|chrom, pos| {
        if region.contains(chrom, pos) {
            n += 1;
        }
    });
    n
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
