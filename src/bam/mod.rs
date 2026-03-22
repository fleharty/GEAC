use std::collections::HashMap;
use std::path::Path;
use std::sync::atomic::Ordering;
use std::time::Instant;

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::pileup::Indel;
use rust_htslib::faidx;

use crate::cli::CollectArgs;
use crate::gene_annotations::GeneAnnotations;
use crate::progress::ProgressReporter;
use crate::record::{AltBase, AltRead, VariantType};
use crate::repeat::compute_repeat_metrics;
use crate::targets::TargetIntervals;
use crate::vcf::VariantAnnotator;

/// Process a BAM/CRAM file and return all alt base records (and optionally per-read detail records).
///
/// When `args.reads_output` is true, the second element of the returned tuple contains one
/// `AltRead` record per read (fragment) that supports an alt base. When false, the second
/// element is always empty.
pub fn collect_alt_bases(args: &CollectArgs, annotator: Option<&dyn VariantAnnotator>, target_intervals: Option<&TargetIntervals>, gene_annots: Option<&GeneAnnotations>) -> Result<(Vec<AltBase>, Vec<AltRead>)> {
    let mut bam = open_bam(&args.input, &args.reference)?;
    let mut ref_cache = RefCache::new(&args.reference)?;

    // Resolve sample ID: CLI flag takes precedence, then SM tag from read group header.
    let sample_id = match &args.sample_id {
        Some(id) => id.clone(),
        None => read_group_sample_id(bam.header())
            .context("--sample-id was not provided and no SM tag found in BAM/CRAM header @RG line")?,
    };

    // Extract target info before the pileup loop to avoid conflicting borrows on `bam`.
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
        bam.fetch(region.as_str())
            .with_context(|| format!("failed to fetch region '{region}': check that the region is valid and the BAM is indexed"))?;
    }

    let start = Instant::now();
    let (reporter, progress) = ProgressReporter::start(args.progress_interval);

    /// Look up variant annotation for a given locus and allele.
    /// Returns (variant_called, variant_filter).
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
        } = tally_pileup(&pileup, args.min_base_qual, args.min_map_qual, ref_base, collect_reads);

        progress.positions_processed.fetch_add(1, Ordering::Relaxed);
        progress.reads_processed.fetch_add(total_depth as u64, Ordering::Relaxed);
        progress.update_locus(&chrom, pos);

        if total_depth == 0 {
            continue;
        }

        // On-target annotation — computed once per locus, shared by all alt records.
        let on_target: Option<bool> = target_intervals.map(|t| t.contains(&chrom, pos));

        // Gene annotation — computed once per locus, shared by all alt records.
        let gene: Option<String> = gene_annots.and_then(|g| g.get(&chrom, pos).map(str::to_owned));

        // Repetitiveness metrics — computed once per locus from the cached ref sequence.
        let repeat = compute_repeat_metrics(ref_cache.current_seq(), pos as usize, args.repeat_window);

        // Trinucleotide context — raw 3-mer centered on this position (null at chromosome edges).
        let trinuc_context: Option<String> = {
            let seq = ref_cache.current_seq();
            let p = pos as usize;
            if p > 0 && p + 1 < seq.len() {
                Some(format!("{}{}{}", seq[p - 1] as char, seq[p] as char, seq[p + 1] as char))
            } else {
                None
            }
        };

        // Extract ref base counts once for use in every alt record at this locus
        let ref_tally = bases.get(&ref_base);
        let ref_count = ref_tally.map_or(0, |t| t.total);
        let fwd_ref_count = ref_tally.map_or(0, |t| t.fwd);
        let rev_ref_count = ref_tally.map_or(0, |t| t.rev);
        let overlap_ref_agree = ref_tally.map_or(0, |t| t.overlap_alt_agree);

        // ── SNV records ───────────────────────────────────────────────────────
        for (base, tally) in &bases {
            if *base == ref_base || *base == 'N' {
                continue;
            }
            if tally.total == 0 {
                continue;
            }

            progress.alt_bases_found.fetch_add(1, Ordering::Relaxed);

            let alt_allele = base.to_string();
            let (variant_called, variant_filter) =
                vcf_annotation(annotator, &chrom, pos, &alt_allele);

            records.push(AltBase {
                sample_id: sample_id.clone(),
                chrom: chrom.clone(),
                pos,
                ref_allele: ref_base.to_string(),
                alt_allele: alt_allele.clone(),
                variant_type: VariantType::Snv,
                total_depth,
                alt_count: tally.total,
                ref_count,
                fwd_depth,
                rev_depth,
                fwd_alt_count: tally.fwd,
                rev_alt_count: tally.rev,
                fwd_ref_count,
                rev_ref_count,
                overlap_depth,
                overlap_alt_agree: tally.overlap_alt_agree,
                overlap_alt_disagree: tally.overlap_alt_disagree,
                overlap_ref_agree,
                read_type: args.read_type,
                pipeline: args.pipeline,
                batch: args.batch.clone(),
                variant_called,
                variant_filter,
                on_target,
                gene: gene.clone(),
                homopolymer_len: repeat.homopolymer_len,
                str_period:      repeat.str_period,
                str_len:         repeat.str_len,
                trinuc_context:  trinuc_context.clone(),
            });

            if collect_reads {
                if let Some(details) = read_details.get(base) {
                    for d in details {
                        read_records.push(AltRead {
                            sample_id: sample_id.clone(),
                            chrom: chrom.clone(),
                            pos,
                            alt_allele: alt_allele.clone(),
                            dist_from_read_start: d.qpos as i32,
                            dist_from_read_end: (d.read_len as i32) - (d.qpos as i32) - 1,
                            read_length: d.read_len as i32,
                            ab_count: d.ab_count,
                            ba_count: d.ba_count,
                            family_size: d.family_size,
                            base_qual: d.base_qual as i32,
                            map_qual: d.map_qual as i32,
                            insert_size: d.insert_size,
                        });
                    }
                }
            }
        }

        // ── Indel records ─────────────────────────────────────────────────────
        let indels = tally_indels(&pileup, pos, ref_cache.current_seq(), args.min_map_qual);

        for (_, indel) in &indels {
            if indel.total == 0 {
                continue;
            }

            progress.alt_bases_found.fetch_add(1, Ordering::Relaxed);

            let (variant_called, variant_filter) =
                vcf_annotation(annotator, &chrom, pos, &indel.alt_allele);

            records.push(AltBase {
                sample_id: sample_id.clone(),
                chrom: chrom.clone(),
                pos,
                ref_allele: indel.ref_allele.clone(),
                alt_allele: indel.alt_allele.clone(),
                variant_type: indel.variant_type,
                total_depth,
                alt_count: indel.total,
                ref_count,
                fwd_depth,
                rev_depth,
                fwd_alt_count: indel.fwd,
                rev_alt_count: indel.rev,
                fwd_ref_count,
                rev_ref_count,
                overlap_depth,
                overlap_alt_agree: indel.overlap_alt_agree,
                overlap_alt_disagree: indel.overlap_alt_disagree,
                overlap_ref_agree,
                read_type: args.read_type,
                pipeline: args.pipeline,
                batch: args.batch.clone(),
                variant_called,
                variant_filter,
                on_target,
                gene: gene.clone(),
                homopolymer_len: repeat.homopolymer_len,
                str_period:      repeat.str_period,
                str_len:         repeat.str_len,
                trinuc_context:  None,
            });
        }
    }

    reporter.finish(start);
    Ok((records, read_records))
}

// ── Pileup tallying ───────────────────────────────────────────────────────────

/// Per-read detail collected during tallying, used to build AltRead records.
struct ReadDetail {
    qpos: usize,
    read_len: usize,
    base_qual: u8,
    map_qual: u8,
    ab_count: Option<i32>,
    ba_count: Option<i32>,
    family_size: Option<i32>,
    insert_size: Option<i32>,
}

/// Position-level summary returned by `tally_pileup`.
struct PileupResult {
    /// Per-alt-base tallies (only populated for bases that passed filters)
    bases: HashMap<char, BaseTally>,
    total_depth: i32,
    fwd_depth: i32,
    rev_depth: i32,
    /// Number of overlapping fragment pairs at this position (pair count, not read count)
    overlap_depth: i32,
    /// Per-read details keyed by base (only populated when collect_reads is true;
    /// only non-ref, non-N bases are included)
    read_details: HashMap<char, Vec<ReadDetail>>,
}

/// Per-base tally at a pileup position.
#[derive(Default)]
struct BaseTally {
    total: i32,
    fwd: i32,
    rev: i32,
    /// Overlapping pairs where both reads agree on this base
    overlap_alt_agree: i32,
    /// Overlapping pairs where one read sees this base and the other sees something different
    overlap_alt_disagree: i32,
}

/// Full per-read data collected during the first pileup pass.
struct LocusRead {
    base: char,
    is_reverse: bool,
    is_first_in_pair: bool,
    qpos: usize,
    read_len: usize,
    base_qual: u8,
    map_qual: u8,
    /// fgbio aD tag: AB (top-strand) raw read count
    ab_count: Option<i32>,
    /// fgbio bD tag: BA (bottom-strand) raw read count
    ba_count: Option<i32>,
    /// fgbio cD tag: total family size (aD + bD for duplex; sole count for simplex)
    family_size: Option<i32>,
    /// SAM TLEN (insert size); None when 0 (unpaired / mate unmapped)
    insert_size: Option<i32>,
}

/// Tally each observed base at a pileup column with overlap detection.
///
/// Overlap is detected by grouping reads by query name. A query name appearing
/// twice at the same position means both reads of the fragment cover that position.
///
/// Depth is counted at the **fragment level** — each overlapping pair contributes 1 to
/// `total_depth`. Strand attribution for overlapping pairs uses the R1 read's orientation
/// (BAM flag 0x40). Non-overlapping singleton N bases are excluded from all tallies.
///
/// Overlapping pair classification rules:
///
/// | Pair (read 1 + read 2)          | total_depth | base tally      | agree / disagree  | overlap_depth |
/// |---------------------------------|-------------|-----------------|-------------------|---------------|
/// | same base + same base (non-N)   | +1          | that base +1    | agree +1          | +1            |
/// | alt + ref                       | +1          | alt +1          | disagree +1       | +1            |
/// | alt₁ + alt₂ (different alts)    | +1          | both +1         | disagree +1 each  | +1            |
/// | alt + N                         | +1          | alt +1          | —                 | +1            |
/// | ref + N                         | +1          | ref +1          | —                 | +1            |
/// | N + N                           | +1          | —               | —                 | +1            |
///
/// For `alt + ref` pairs, `overlap_alt_disagree` is still incremented even though the
/// fragment is classified as alt — the disagreement is a useful quality signal.
///
/// When `collect_reads` is true, `PileupResult.read_details` is populated with per-read
/// detail for every non-ref, non-N base. When false, `read_details` is always empty.
fn tally_pileup(
    pileup: &rust_htslib::bam::pileup::Pileup,
    min_base_qual: u8,
    min_map_qual: u8,
    ref_base: char,
    collect_reads: bool,
) -> PileupResult {
    // First pass: collect LocusRead per query name.
    let mut by_qname: HashMap<Vec<u8>, Vec<LocusRead>> = HashMap::new();

    for alignment in pileup.alignments() {
        if alignment.is_del() || alignment.is_refskip() {
            continue;
        }

        let record = alignment.record();

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
        let is_reverse = record.is_reverse();
        let is_first_in_pair = record.flags() & 0x40 != 0;
        let map_qual = record.mapq();
        let read_len = record.seq_len();

        let (ab_count, ba_count, family_size, insert_size) = if collect_reads {
            let ab = aux_i32(&record, b"aD");
            let ba = aux_i32(&record, b"bD");
            let fs = aux_i32(&record, b"cD");
            let tlen = record.insert_size();
            let ins = if tlen == 0 { None } else { Some(tlen.unsigned_abs() as i32) };
            (ab, ba, fs, ins)
        } else {
            (None, None, None, None)
        };

        by_qname
            .entry(record.qname().to_vec())
            .or_default()
            .push(LocusRead { base, is_reverse, is_first_in_pair, qpos, read_len, base_qual, map_qual, ab_count, ba_count, family_size, insert_size });
    }

    // Second pass: tally with overlap detection.
    let mut bases: HashMap<char, BaseTally> = HashMap::new();
    let mut read_details: HashMap<char, Vec<ReadDetail>> = HashMap::new();
    let mut total_depth: i32 = 0;
    let mut fwd_depth: i32 = 0;
    let mut rev_depth: i32 = 0;
    let mut overlap_depth: i32 = 0;

    // Helper macro: push a ReadDetail for base `b` from LocusRead `r` (only when collect_reads).
    macro_rules! push_detail {
        ($b:expr, $r:expr) => {
            if collect_reads && $b != ref_base && $b != 'N' {
                read_details.entry($b).or_default().push(ReadDetail {
                    qpos: $r.qpos,
                    read_len: $r.read_len,
                    base_qual: $r.base_qual,
                    map_qual: $r.map_qual,
                    ab_count: $r.ab_count,
                    ba_count: $r.ba_count,
                    family_size: $r.family_size,
                    insert_size: $r.insert_size,
                });
            }
        };
    }

    for reads in by_qname.values() {
        match reads.as_slice() {
            [r] => {
                // Non-overlapping read; skip uninformative N bases entirely.
                if r.base == 'N' {
                    continue;
                }
                total_depth += 1;
                if r.is_reverse { rev_depth += 1; } else { fwd_depth += 1; }

                let t = bases.entry(r.base).or_default();
                t.total += 1;
                if r.is_reverse { t.rev += 1; } else { t.fwd += 1; }

                push_detail!(r.base, r);
            }
            [r1, r2] => {
                // Overlapping fragment: fragment-level depth — always counts as 1.
                overlap_depth += 1;
                total_depth += 1;

                // Strand is attributed using R1's orientation (BAM flag 0x40).
                let r1_is_rev = if r1.is_first_in_pair { r1.is_reverse } else { r2.is_reverse };
                if r1_is_rev { rev_depth += 1; } else { fwd_depth += 1; }

                let b1 = r1.base;
                let b2 = r2.base;

                if b1 == 'N' && b2 == 'N' {
                    // N + N: fragment counted in depth but no base tally.
                } else if b1 == 'N' {
                    // N + informative: tally the informative base.
                    let t = bases.entry(b2).or_default();
                    t.total += 1;
                    if r1_is_rev { t.rev += 1; } else { t.fwd += 1; }
                    push_detail!(b2, r2);
                } else if b2 == 'N' {
                    // informative + N: tally the informative base.
                    let t = bases.entry(b1).or_default();
                    t.total += 1;
                    if r1_is_rev { t.rev += 1; } else { t.fwd += 1; }
                    push_detail!(b1, r1);
                } else if b1 == b2 {
                    // Both reads agree on the same base; use R1 for read details.
                    let t = bases.entry(b1).or_default();
                    t.total += 1;
                    if r1_is_rev { t.rev += 1; } else { t.fwd += 1; }
                    t.overlap_alt_agree += 1;
                    let det = if r1.is_first_in_pair { r1 } else { r2 };
                    push_detail!(b1, det);
                } else if b1 == ref_base {
                    // b2 is alt, b1 is ref — alt wins; classify as alt.
                    let t = bases.entry(b2).or_default();
                    t.total += 1;
                    if r1_is_rev { t.rev += 1; } else { t.fwd += 1; }
                    t.overlap_alt_disagree += 1;
                    push_detail!(b2, r2);
                } else if b2 == ref_base {
                    // b1 is alt, b2 is ref — alt wins; classify as alt.
                    let t = bases.entry(b1).or_default();
                    t.total += 1;
                    if r1_is_rev { t.rev += 1; } else { t.fwd += 1; }
                    t.overlap_alt_disagree += 1;
                    push_detail!(b1, r1);
                } else {
                    // alt₁ + alt₂ — two different non-ref bases; tally both.
                    // Note: sum of base counts will exceed total_depth for this fragment.
                    let t1 = bases.entry(b1).or_default();
                    t1.total += 1;
                    if r1_is_rev { t1.rev += 1; } else { t1.fwd += 1; }
                    t1.overlap_alt_disagree += 1;
                    push_detail!(b1, r1);

                    let t2 = bases.entry(b2).or_default();
                    t2.total += 1;
                    if r1_is_rev { t2.rev += 1; } else { t2.fwd += 1; }
                    t2.overlap_alt_disagree += 1;
                    push_detail!(b2, r2);
                }
            }
            _ => {
                // More than 2 reads with the same query name: shouldn't happen in
                // practice but handle gracefully by treating as non-overlapping.
                for r in reads {
                    if r.base == 'N' {
                        continue;
                    }
                    total_depth += 1;
                    if r.is_reverse { rev_depth += 1; } else { fwd_depth += 1; }
                    let t = bases.entry(r.base).or_default();
                    t.total += 1;
                    if r.is_reverse { t.rev += 1; } else { t.fwd += 1; }
                    push_detail!(r.base, r);
                }
            }
        }
    }

    PileupResult { bases, total_depth, fwd_depth, rev_depth, overlap_depth, read_details }
}

/// Read an integer auxiliary tag from a BAM record, returning None if absent or wrong type.
/// Read an integer auxiliary tag from a BAM record, returning None if absent or non-integer.
///
/// Float tags are intentionally excluded. For fgbio simplex consensus reads, `cE` is the
/// per-base error rate (float), not a family size count. For duplex consensus reads, `cE`
/// is the BA strand count (integer) and will be read correctly here.
fn aux_i32(record: &bam::Record, tag: &[u8; 2]) -> Option<i32> {
    match record.aux(tag) {
        Ok(bam::record::Aux::I8(v))  => Some(v as i32),
        Ok(bam::record::Aux::U8(v))  => Some(v as i32),
        Ok(bam::record::Aux::I16(v)) => Some(v as i32),
        Ok(bam::record::Aux::U16(v)) => Some(v as i32),
        Ok(bam::record::Aux::I32(v)) => Some(v),
        Ok(bam::record::Aux::U32(v)) => Some(v as i32),
        _ => None,
    }
}

// ── Indel tallying ────────────────────────────────────────────────────────────

/// Per-indel-allele tally at a pileup position.
struct IndelCount {
    ref_allele: String,
    alt_allele: String,
    variant_type: VariantType,
    total: i32,
    fwd: i32,
    rev: i32,
    /// Overlapping pairs where both reads agree on this indel allele
    overlap_alt_agree: i32,
    /// Overlapping pairs where reads disagree (one has this indel, the other differs or has none)
    overlap_alt_disagree: i32,
}

/// Decoded indel allele for one read at the anchor position, or None if no indel.
type IndelAllele = Option<(String, String, VariantType)>; // (alt_allele, ref_allele, variant_type)

/// Tally indel alleles at a pileup column with overlap detection.
///
/// Groups reads by query name. A name appearing twice means both reads of the
/// fragment are at the anchor position — their indel alleles are compared to
/// determine agreement or disagreement.
///
/// Every read that passes the mapping quality filter pushes an entry (Some or None)
/// so that overlapping pairs are correctly identified even when one read has no indel.
fn tally_indels(
    pileup: &rust_htslib::bam::pileup::Pileup,
    pos: i64,
    chrom_seq: &[u8],
    min_map_qual: u8,
) -> HashMap<String, IndelCount> {
    // First pass: collect (indel_allele_or_none, is_reverse, is_first_in_pair) per query name.
    let mut by_qname: HashMap<Vec<u8>, Vec<(IndelAllele, bool, bool)>> = HashMap::new();

    for alignment in pileup.alignments() {
        if alignment.is_refskip() {
            continue;
        }
        let record = alignment.record();
        if record.mapq() < min_map_qual {
            continue;
        }

        let is_reverse = record.is_reverse();
        let is_first_in_pair = record.flags() & 0x40 != 0;
        let qname = record.qname().to_vec();

        let allele: IndelAllele = match alignment.indel() {
            Indel::Ins(len) => {
                let Some(qpos) = alignment.qpos() else {
                    by_qname.entry(qname).or_default().push((None, is_reverse, is_first_in_pair));
                    continue;
                };
                let seq = record.seq();
                let len = len as usize;
                if qpos + len >= seq.len() {
                    by_qname.entry(qname).or_default().push((None, is_reverse, is_first_in_pair));
                    continue;
                }
                let inserted: String = (1..=len)
                    .map(|i| seq[qpos + i].to_ascii_uppercase() as char)
                    .collect();
                let ref_allele = chrom_seq
                    .get(pos as usize)
                    .map(|&b| (b as char).to_string())
                    .unwrap_or_default();
                Some((format!("+{inserted}"), ref_allele, VariantType::Insertion))
            }

            Indel::Del(len) => {
                let start = pos as usize + 1;
                let end = start + len as usize;
                let deleted: String = chrom_seq
                    .get(start..end)
                    .unwrap_or(&[])
                    .iter()
                    .map(|&b| b as char)
                    .collect();
                if deleted.len() != len as usize {
                    by_qname.entry(qname).or_default().push((None, is_reverse, is_first_in_pair));
                    continue;
                }
                Some((format!("-{deleted}"), deleted.clone(), VariantType::Deletion))
            }

            Indel::None => None,
        };

        by_qname.entry(qname).or_default().push((allele, is_reverse, is_first_in_pair));
    }

    // Second pass: tally with overlap detection.
    let mut indels: HashMap<String, IndelCount> = HashMap::new();

    for reads in by_qname.values() {
        match reads.as_slice() {
            [(allele, is_rev, _)] => {
                // Non-overlapping read
                if let Some((alt, ref_a, vt)) = allele {
                    let e = indels.entry(alt.clone()).or_insert_with(|| IndelCount {
                        ref_allele: ref_a.clone(), alt_allele: alt.clone(),
                        variant_type: *vt, total: 0, fwd: 0, rev: 0,
                        overlap_alt_agree: 0, overlap_alt_disagree: 0,
                    });
                    e.total += 1;
                    if *is_rev { e.rev += 1; } else { e.fwd += 1; }
                }
            }

            [(allele1, is_rev1, is_first1), (allele2, is_rev2, _is_first2)] => {
                // Overlapping fragment — fragment-level counting.
                // Strand is attributed using R1's orientation (BAM flag 0x40).
                let r1_is_rev = if *is_first1 { *is_rev1 } else { *is_rev2 };

                match (allele1, allele2) {
                    (None, None) => {
                        // Neither read has an indel — nothing to tally.
                    }
                    (Some((alt, ref_a, vt)), None) | (None, Some((alt, ref_a, vt))) => {
                        // Indel + no indel: alt wins; classify as indel with disagree.
                        let e = indels.entry(alt.clone()).or_insert_with(|| IndelCount {
                            ref_allele: ref_a.clone(), alt_allele: alt.clone(),
                            variant_type: *vt, total: 0, fwd: 0, rev: 0,
                            overlap_alt_agree: 0, overlap_alt_disagree: 0,
                        });
                        e.total += 1;
                        if r1_is_rev { e.rev += 1; } else { e.fwd += 1; }
                        e.overlap_alt_disagree += 1;
                    }
                    (Some((alt1, ref_a1, vt1)), Some((alt2, ref_a2, vt2))) => {
                        if alt1 == alt2 {
                            // Same indel: count once with agree.
                            let e = indels.entry(alt1.clone()).or_insert_with(|| IndelCount {
                                ref_allele: ref_a1.clone(), alt_allele: alt1.clone(),
                                variant_type: *vt1, total: 0, fwd: 0, rev: 0,
                                overlap_alt_agree: 0, overlap_alt_disagree: 0,
                            });
                            e.total += 1;
                            if r1_is_rev { e.rev += 1; } else { e.fwd += 1; }
                            e.overlap_alt_agree += 1;
                        } else {
                            // Different indels: tally both with disagree (edge case).
                            let e1 = indels.entry(alt1.clone()).or_insert_with(|| IndelCount {
                                ref_allele: ref_a1.clone(), alt_allele: alt1.clone(),
                                variant_type: *vt1, total: 0, fwd: 0, rev: 0,
                                overlap_alt_agree: 0, overlap_alt_disagree: 0,
                            });
                            e1.total += 1;
                            if r1_is_rev { e1.rev += 1; } else { e1.fwd += 1; }
                            e1.overlap_alt_disagree += 1;

                            let e2 = indels.entry(alt2.clone()).or_insert_with(|| IndelCount {
                                ref_allele: ref_a2.clone(), alt_allele: alt2.clone(),
                                variant_type: *vt2, total: 0, fwd: 0, rev: 0,
                                overlap_alt_agree: 0, overlap_alt_disagree: 0,
                            });
                            e2.total += 1;
                            if r1_is_rev { e2.rev += 1; } else { e2.fwd += 1; }
                            e2.overlap_alt_disagree += 1;
                        }
                    }
                }
            }

            _ => {
                // More than 2 reads with same name — treat as non-overlapping.
                for (allele, is_rev, _) in reads {
                    if let Some((alt, ref_a, vt)) = allele {
                        let e = indels.entry(alt.clone()).or_insert_with(|| IndelCount {
                            ref_allele: ref_a.clone(), alt_allele: alt.clone(),
                            variant_type: *vt, total: 0, fwd: 0, rev: 0,
                            overlap_alt_agree: 0, overlap_alt_disagree: 0,
                        });
                        e.total += 1;
                        if *is_rev { e.rev += 1; } else { e.fwd += 1; }
                    }
                }
            }
        }
    }

    indels
}

// ── Reference cache ───────────────────────────────────────────────────────────

/// Caches one chromosome sequence at a time to avoid per-base faidx seeks.
/// When the chromosome changes, the new sequence is fetched and the old one dropped.
struct RefCache {
    fai: faidx::Reader,
    current_tid: Option<usize>,
    current_chrom: String,
    /// Upper-cased sequence bytes for the cached chromosome (0-based)
    current_seq: Vec<u8>,
}

impl RefCache {
    fn new(reference: &Path) -> Result<Self> {
        let fai = faidx::Reader::from_path(reference)
            .with_context(|| format!("failed to open reference FASTA: {}", reference.display()))?;
        Ok(Self {
            fai,
            current_tid: None,
            current_chrom: String::new(),
            current_seq: Vec::new(),
        })
    }

    /// Returns (chrom_name, ref_base) for the given tid and 0-based position.
    /// Loads the chromosome sequence on first access or when the chromosome changes.
    /// `targets` is a pre-built slice of (chrom_name, chrom_len) indexed by tid.
    fn get(&mut self, targets: &[(String, usize)], tid: usize, pos: usize) -> Result<(String, char)> {
        if self.current_tid != Some(tid) {
            let (chrom, chrom_len) = targets
                .get(tid)
                .with_context(|| format!("tid {tid} not found in BAM header"))?;

            // fetch_seq end is 0-based inclusive
            self.current_seq = self
                .fai
                .fetch_seq(chrom, 0, chrom_len.saturating_sub(1))
                .with_context(|| format!("failed to fetch sequence for {chrom}"))?
                .to_vec();

            self.current_seq.make_ascii_uppercase();

            self.current_tid = Some(tid);
            self.current_chrom = chrom.clone();
        }

        let base = self
            .current_seq
            .get(pos)
            .map(|&b| b as char)
            .unwrap_or('N');

        Ok((self.current_chrom.clone(), base))
    }

    /// Returns the cached sequence for the current chromosome.
    /// Must be called after `get()` has loaded the chromosome.
    fn current_seq(&self) -> &[u8] {
        &self.current_seq
    }
}

// ── BAM/CRAM opener ───────────────────────────────────────────────────────────

/// Extract the SM (sample name) field from the first @RG line in the BAM header.
/// Returns an error if no @RG line exists or none has an SM tag.
fn read_group_sample_id(header: &bam::HeaderView) -> Result<String> {
    let header_text = std::str::from_utf8(header.as_bytes())
        .context("BAM header is not valid UTF-8")?;

    for line in header_text.lines() {
        if !line.starts_with("@RG") {
            continue;
        }
        for field in line.split('\t') {
            if let Some(sm) = field.strip_prefix("SM:") {
                return Ok(sm.to_string());
            }
        }
    }

    anyhow::bail!("no SM tag found in any @RG line of the BAM/CRAM header")
}

fn open_bam(input: &Path, reference: &Path) -> Result<bam::IndexedReader> {
    let mut reader = bam::IndexedReader::from_path(input)
        .with_context(|| format!("failed to open BAM/CRAM: {}", input.display()))?;
    reader
        .set_reference(reference)
        .with_context(|| format!("failed to set reference: {}", reference.display()))?;
    Ok(reader)
}
