use std::collections::HashMap;
use std::path::Path;
use std::sync::atomic::Ordering;
use std::time::Instant;

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read};
use rust_htslib::faidx;

use crate::cli::CollectArgs;
use crate::progress::ProgressReporter;
use crate::record::{AltBase, VariantType};

/// Process a BAM/CRAM file and return all alt base records.
pub fn collect_alt_bases(args: &CollectArgs) -> Result<Vec<AltBase>> {
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
        // TODO: parse region string and set fetch interval
        let _ = region;
    }

    let start = Instant::now();
    let (reporter, progress) = ProgressReporter::start(args.progress_interval);
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

        for (base, tally) in &bases {
            if *base == ref_base || *base == 'N' {
                continue;
            }
            if tally.total == 0 {
                continue;
            }

            progress.alt_bases_found.fetch_add(1, Ordering::Relaxed);

            records.push(AltBase {
                sample_id: sample_id.clone(),
                chrom: chrom.clone(),
                pos,
                ref_allele: ref_base.to_string(),
                alt_allele: base.to_string(),
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
                variant_called: None,
                variant_filter: None,
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
                    // Both reads agree on this base
                    bases.entry(*base1).or_default().overlap_alt_agree += 1;
                } else {
                    // Reads disagree — flag both bases as discordant
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
