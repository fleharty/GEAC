use std::collections::HashMap;
use std::path::Path;
use std::sync::atomic::Ordering;
use std::time::Instant;

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::pileup::Indel;
use rust_htslib::faidx;

use crate::cli::CollectArgs;
use crate::progress::ProgressReporter;
use crate::record::{AltBase, VariantType};
use crate::vcf::VcfIndex;

/// Process a BAM/CRAM file and return all alt base records.
pub fn collect_alt_bases(args: &CollectArgs, vcf_index: Option<&VcfIndex>) -> Result<Vec<AltBase>> {
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

    /// Look up VCF annotation for a given locus and allele.
    /// Returns (variant_called, variant_filter).
    fn vcf_annotation(
        vcf_index: Option<&VcfIndex>,
        chrom: &str,
        pos: i64,
        alt_allele: &str,
    ) -> (Option<bool>, Option<String>) {
        match vcf_index {
            None => (None, None),
            Some(idx) => match idx.get(chrom, pos, alt_allele) {
                Some(ann) => (Some(true), Some(ann.filter.clone())),
                None => (Some(false), None),
            },
        }
    }
    let mut records: Vec<AltBase> = Vec::new();

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
        } = tally_pileup(&pileup, args.min_base_qual, args.min_map_qual);

        progress.positions_processed.fetch_add(1, Ordering::Relaxed);
        progress.reads_processed.fetch_add(total_depth as u64, Ordering::Relaxed);
        progress.update_locus(&chrom, pos);

        if total_depth == 0 {
            continue;
        }

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
                vcf_annotation(vcf_index, &chrom, pos, &alt_allele);

            records.push(AltBase {
                sample_id: sample_id.clone(),
                chrom: chrom.clone(),
                pos,
                ref_allele: ref_base.to_string(),
                alt_allele,
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
                variant_called,
                variant_filter,
            });
        }

        // ── Indel records ─────────────────────────────────────────────────────
        let indels = tally_indels(&pileup, pos, ref_cache.current_seq(), args.min_map_qual);

        for (_, indel) in &indels {
            if indel.total == 0 {
                continue;
            }

            progress.alt_bases_found.fetch_add(1, Ordering::Relaxed);

            let (variant_called, variant_filter) =
                vcf_annotation(vcf_index, &chrom, pos, &indel.alt_allele);

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
                variant_called,
                variant_filter,
            });
        }
    }

    reporter.finish(start);
    Ok(records)
}

// ── Pileup tallying ───────────────────────────────────────────────────────────

/// Position-level summary returned by `tally_pileup`.
struct PileupResult {
    /// Per-alt-base tallies (only populated for bases that passed filters)
    bases: HashMap<char, BaseTally>,
    total_depth: i32,
    fwd_depth: i32,
    rev_depth: i32,
    /// Number of overlapping fragment pairs at this position (pair count, not read count)
    overlap_depth: i32,
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

/// Tally each observed base at a pileup column with overlap detection.
///
/// Overlap is detected by grouping reads by query name. A query name appearing
/// twice at the same position means both reads of the fragment cover that position.
fn tally_pileup(
    pileup: &rust_htslib::bam::pileup::Pileup,
    min_base_qual: u8,
    min_map_qual: u8,
) -> PileupResult {
    // First pass: collect (base, is_reverse) per query name.
    // Query names are stored as Vec<u8> to avoid UTF-8 allocation overhead.
    let mut by_qname: HashMap<Vec<u8>, Vec<(char, bool)>> = HashMap::new();

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

        by_qname
            .entry(record.qname().to_vec())
            .or_default()
            .push((base, is_reverse));
    }

    // Second pass: tally with overlap detection.
    let mut bases: HashMap<char, BaseTally> = HashMap::new();
    let mut total_depth: i32 = 0;
    let mut fwd_depth: i32 = 0;
    let mut rev_depth: i32 = 0;
    let mut overlap_depth: i32 = 0;

    for reads in by_qname.values() {
        match reads.as_slice() {
            [(base, is_rev)] => {
                // Non-overlapping read
                total_depth += 1;
                if *is_rev { rev_depth += 1; } else { fwd_depth += 1; }

                let t = bases.entry(*base).or_default();
                t.total += 1;
                if *is_rev { t.rev += 1; } else { t.fwd += 1; }
            }
            [(base1, is_rev1), (base2, is_rev2)] => {
                // Overlapping fragment: both reads cover this position
                overlap_depth += 1;
                total_depth += 2;
                if *is_rev1 { rev_depth += 1; } else { fwd_depth += 1; }
                if *is_rev2 { rev_depth += 1; } else { fwd_depth += 1; }

                let t1 = bases.entry(*base1).or_default();
                t1.total += 1;
                if *is_rev1 { t1.rev += 1; } else { t1.fwd += 1; }

                let t2 = bases.entry(*base2).or_default();
                t2.total += 1;
                if *is_rev2 { t2.rev += 1; } else { t2.fwd += 1; }

                if base1 == base2 {
                    bases.entry(*base1).or_default().overlap_alt_agree += 1;
                } else {
                    bases.entry(*base1).or_default().overlap_alt_disagree += 1;
                    bases.entry(*base2).or_default().overlap_alt_disagree += 1;
                }
            }
            _ => {
                // More than 2 reads with the same query name: shouldn't happen in
                // practice but handle gracefully by treating as non-overlapping.
                for &(base, is_rev) in reads {
                    total_depth += 1;
                    if is_rev { rev_depth += 1; } else { fwd_depth += 1; }
                    let t = bases.entry(base).or_default();
                    t.total += 1;
                    if is_rev { t.rev += 1; } else { t.fwd += 1; }
                }
            }
        }
    }

    PileupResult { bases, total_depth, fwd_depth, rev_depth, overlap_depth }
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
    // First pass: collect (indel_allele_or_none, is_reverse) per query name.
    let mut by_qname: HashMap<Vec<u8>, Vec<(IndelAllele, bool)>> = HashMap::new();

    for alignment in pileup.alignments() {
        if alignment.is_refskip() {
            continue;
        }
        let record = alignment.record();
        if record.mapq() < min_map_qual {
            continue;
        }

        let is_reverse = record.is_reverse();
        let qname = record.qname().to_vec();

        let allele: IndelAllele = match alignment.indel() {
            Indel::Ins(len) => {
                let Some(qpos) = alignment.qpos() else {
                    by_qname.entry(qname).or_default().push((None, is_reverse));
                    continue;
                };
                let seq = record.seq();
                let len = len as usize;
                if qpos + len >= seq.len() {
                    by_qname.entry(qname).or_default().push((None, is_reverse));
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
                    by_qname.entry(qname).or_default().push((None, is_reverse));
                    continue;
                }
                Some((format!("-{deleted}"), deleted.clone(), VariantType::Deletion))
            }

            Indel::None => None,
        };

        by_qname.entry(qname).or_default().push((allele, is_reverse));
    }

    // Second pass: tally with overlap detection.
    let mut indels: HashMap<String, IndelCount> = HashMap::new();

    for reads in by_qname.values() {
        match reads.as_slice() {
            [(allele, is_rev)] => {
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

            [(allele1, is_rev1), (allele2, is_rev2)] => {
                // Overlapping fragment — count both reads
                for (allele, is_rev) in [(allele1, is_rev1), (allele2, is_rev2)] {
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

                // Determine agreement
                let key1 = allele1.as_ref().map(|(a, _, _)| a.as_str());
                let key2 = allele2.as_ref().map(|(a, _, _)| a.as_str());
                match (key1, key2) {
                    (Some(a1), Some(a2)) if a1 == a2 => {
                        if let Some(e) = indels.get_mut(a1) {
                            e.overlap_alt_agree += 1;
                        }
                    }
                    (Some(a1), Some(a2)) => {
                        if let Some(e) = indels.get_mut(a1) { e.overlap_alt_disagree += 1; }
                        if let Some(e) = indels.get_mut(a2) { e.overlap_alt_disagree += 1; }
                    }
                    (Some(a1), None) => {
                        if let Some(e) = indels.get_mut(a1) { e.overlap_alt_disagree += 1; }
                    }
                    (None, Some(a2)) => {
                        if let Some(e) = indels.get_mut(a2) { e.overlap_alt_disagree += 1; }
                    }
                    (None, None) => {}
                }
            }

            _ => {
                // More than 2 reads with same name — treat as non-overlapping
                for (allele, is_rev) in reads {
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
