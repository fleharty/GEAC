use anyhow::{Context, Result};
use rust_htslib::bam::Read;
use tracing::info;

use crate::bam::{gc_frac, open_bam, read_group_sample_id, RefCache};
use crate::cli::FragmentsArgs;
use crate::record::FragmentRecord;
use crate::region::{RegionInput, parse_region_input};
use crate::writer::parquet_fragments::FragmentsWriter;

pub fn collect_fragments(args: &FragmentsArgs, writer: &mut FragmentsWriter) -> Result<()> {
    let mut bam = open_bam(&args.input, &args.reference)?;
    let mut ref_cache = RefCache::new(&args.reference)?;

    let sample_id = match &args.sample_id {
        Some(id) => id.clone(),
        None => read_group_sample_id(bam.header())
            .context("--sample-id not provided and no SM tag in BAM/CRAM header")?,
    };

    let bam_contigs: Vec<(String, usize)> = {
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

    let region_input = parse_region_input(args.region.as_deref())?;
    let fetch_regions: Vec<Option<String>> = match &region_input {
        None => vec![None],
        Some(RegionInput::Single(s)) => vec![Some(s.clone())],
        Some(RegionInput::Intervals(ivs)) => ivs
            .named_intervals()
            .iter()
            .map(|iv| Some(format!("{}:{}-{}", iv.chrom, iv.start + 1, iv.end)))
            .collect(),
    };

    let mut fragments_written: u64 = 0;

    for fetch_region in &fetch_regions {
        if let Some(r) = fetch_region.as_deref() {
            bam.fetch(r)
                .with_context(|| format!("failed to fetch region '{r}'"))?;
        }

        for record in bam.records() {
            let record = record.context("error reading BAM record")?;

            // Only process R1 of proper pairs; skip duplicates, unmapped, secondary, supplementary
            if record.is_duplicate()
                || record.is_unmapped()
                || record.is_mate_unmapped()
                || !record.is_paired()
                || !record.is_proper_pair()
                || !record.is_first_in_template()
                || record.is_secondary()
                || record.is_supplementary()
                || (record.mapq() as u8) < args.min_map_qual
            {
                continue;
            }

            let insert_size = record.insert_size();
            if insert_size <= 0 {
                continue;
            }

            let frag_start = record.pos();
            let insert_size = insert_size as i32;
            let frag_end = frag_start + insert_size as i64;
            let midpoint = (frag_start + frag_end) / 2;
            let tid = record.tid() as usize;
            let map_qual = record.mapq() as i32;

            // Load reference sequence for this contig (cached per chromosome)
            let chrom = bam_contigs
                .get(tid)
                .map(|(name, _)| name.clone())
                .unwrap_or_default();

            // Trigger RefCache load for this tid
            ref_cache.get(&bam_contigs, tid, frag_start as usize)?;
            let seq = ref_cache.current_seq();

            let gc_content = gc_frac(seq, frag_start as usize, insert_size as usize);

            let end_motif_5p = extract_motif(seq, frag_start as usize, 4);
            let end_motif_3p = if frag_end as usize >= 4 {
                extract_motif(seq, frag_end as usize - 4, 4)
            } else {
                None
            };

            writer.push(FragmentRecord {
                sample_id: sample_id.clone(),
                chrom,
                frag_start,
                frag_end,
                midpoint,
                insert_size,
                gc_content,
                end_motif_5p,
                end_motif_3p,
                map_qual,
                read_type: args.read_type,
                pipeline: args.pipeline,
                batch: args.batch.clone(),
                label1: args.label1.clone(),
                label2: args.label2.clone(),
                label3: args.label3.clone(),
                timepoint: args.timepoint.clone(),
            })?;

            fragments_written += 1;
        }
    }

    info!(fragments_written, "fragment collection complete");
    Ok(())
}

fn extract_motif(seq: &[u8], start: usize, len: usize) -> Option<String> {
    let end = start + len;
    if end > seq.len() {
        return None;
    }
    let motif: String = seq[start..end].iter().map(|&b| b as char).collect();
    if motif.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T')) {
        Some(motif)
    } else {
        None
    }
}
