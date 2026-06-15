use anyhow::{Context, Result};
use rust_htslib::bam::Read;
use tracing::info;

use crate::bam::{gc_frac, open_bam, read_group_sample_id, RefCache};
use crate::cli::FragmentsArgs;
use crate::record::FragmentRecord;
use crate::region::{parse_region_input, RegionInput};
use crate::writer::parquet_fragments::FragmentsWriter;

pub fn collect_fragments(args: &FragmentsArgs, writer: &mut FragmentsWriter) -> Result<()> {
    let mut bam = open_bam(&args.input, &args.reference, args.index.as_deref())?;
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
        match fetch_region.as_deref() {
            Some(r) => bam
                .fetch(r)
                .with_context(|| format!("failed to fetch region '{r}'"))?,
            None => bam.fetch(".").context("failed to fetch all reads")?,
        }

        for record in bam.records() {
            let record = record.context("error reading BAM record")?;

            // Process the leftmost read of each proper pair (positive TLEN). By
            // SAM spec, exactly one read per pair has positive TLEN — the
            // leftmost one — and in FR libraries it is on the + strand
            // regardless of whether it is R1 or R2. Filtering on `is_first_in_template`
            // would drop pairs where R1 maps to the − strand (~half of all pairs).
            if record.is_duplicate()
                || record.is_unmapped()
                || record.is_mate_unmapped()
                || !record.is_paired()
                || !record.is_proper_pair()
                || record.is_secondary()
                || record.is_supplementary()
                || record.mapq() < args.min_map_qual
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

            // Centered 4-mer: [cut-2, cut+2) — 2 bp inside + 2 bp outside the fragment,
            // matching the standard cfDNA end-motif convention (Jiang et al. 2020).
            // The 5′ motif is read directly from the + strand.
            // The 3′ end is the 5′ end of the − strand, so its motif must be
            // reverse-complemented to be reported in the same 5′→3′ convention.
            let end_motif_5p = if frag_start as usize >= 2 {
                extract_motif(seq, frag_start as usize - 2, 4)
            } else {
                None
            };
            let end_motif_3p = if frag_end as usize >= 2 {
                extract_motif(seq, frag_end as usize - 2, 4).map(|m| reverse_complement(&m))
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
                pipeline: args.pipeline.clone(),
                subject_id: args.subject_id.clone(),
                sample_type: args.sample_type.clone(),
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

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_motif_basic() {
        let seq = b"ACGTACGTAC";
        assert_eq!(extract_motif(seq, 0, 4).as_deref(), Some("ACGT"));
        assert_eq!(extract_motif(seq, 2, 4).as_deref(), Some("GTAC"));
        assert_eq!(extract_motif(seq, 6, 4).as_deref(), Some("GTAC"));
    }

    #[test]
    fn extract_motif_out_of_bounds_returns_none() {
        let seq = b"ACGT";
        assert!(extract_motif(seq, 1, 4).is_none());
        assert!(extract_motif(seq, 4, 1).is_none());
    }

    #[test]
    fn extract_motif_rejects_non_acgt() {
        let seq = b"ACNT";
        assert!(extract_motif(seq, 0, 4).is_none());
        let lower = b"acgt";
        assert!(extract_motif(lower, 0, 4).is_none());
    }

    #[test]
    fn reverse_complement_basic() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("CCGG"), "CCGG");
    }

    #[test]
    fn reverse_complement_is_involution() {
        for motif in ["ACGT", "TTGA", "CCGA", "GATC", "AAAT"] {
            assert_eq!(reverse_complement(&reverse_complement(motif)), motif);
        }
    }

    #[test]
    fn end_motif_3p_orientation_matches_jiang_convention() {
        // Reference (+ strand) around a fragment 3' cut site: ...AATT|CGGC...
        // (cut between positions 3 and 4 on this snippet)
        // The + strand reads 5'-AATTCGGC-3'.
        // The − strand at the cut site reads 5'-GCCG AATT-3' (reverse complement).
        // The 3' end motif as the nuclease "sees" it on the − strand
        // (centered on cut, [cut-2, cut+2)) is the 4-mer at the 5' end of the
        // − strand fragment, i.e. revcomp of the + strand bases [cut-2, cut+2).
        let seq = b"AATTCGGC";
        // cut at index 4 → [2, 6) = b"TTCG" on + strand
        let plus = extract_motif(seq, 2, 4).unwrap();
        assert_eq!(plus, "TTCG");
        // 3' motif (reported in 5'→3' of − strand) = reverse complement
        assert_eq!(reverse_complement(&plus), "CGAA");
    }
}
