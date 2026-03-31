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
use crate::record::{AltBase, AltRead};
use crate::repeat::compute_repeat_metrics;
use crate::targets::TargetIntervals;
use crate::vcf::VariantAnnotator;

use builders::LocusContext;
use indel::tally_indels;
use pileup::{tally_pileup, PileupResult};

pub(crate) use ref_utils::RefCache;
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
) -> Result<(Vec<AltBase>, Vec<AltRead>)> {
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
            args.min_base_qual,
            args.min_map_qual,
            args.include_duplicates,
            args.include_secondary,
            args.include_supplementary,
            ref_base,
            collect_reads,
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
        let gene = gene_annots.and_then(|g| g.get(&chrom, pos).map(str::to_owned));
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

            locus.push_snv_record(
                &mut records,
                *base,
                tally,
                variant_called,
                variant_filter,
                gnomad_af,
            );

            if collect_reads {
                if let Some(details) = read_details.get(base) {
                    for detail in details {
                        read_records.push(locus.build_alt_read(&alt_allele, detail));
                    }
                }
            }
        }

        let (indels, indel_read_details) = tally_indels(
            &pileup,
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

            locus.push_indel_record(
                &mut records,
                indel,
                variant_called,
                variant_filter,
                gnomad_af,
            );
        }

        if collect_reads {
            for (alt_allele, details) in &indel_read_details {
                for detail in details {
                    read_records.push(locus.build_alt_read(alt_allele, detail));
                }
            }
        }
    }

    reporter.finish(start);
    Ok((records, read_records))
}

fn compute_input_sha256(path: &Path) -> Result<String> {
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
