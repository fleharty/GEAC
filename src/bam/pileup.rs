use std::collections::HashMap;

use rust_htslib::bam;

/// Per-read detail collected during tallying, used to build AltRead records.
#[derive(Clone)]
pub(super) struct ReadDetail {
    pub(super) qpos: usize,
    pub(super) read_len: usize,
    pub(super) is_first_in_pair: bool,
    /// True when the read is on the reverse strand (BAM flag 0x10).
    pub(super) is_reverse: bool,
    /// Number of hard-clipped bases at the 5' end of the original synthesis.
    ///
    /// For forward-strand reads this equals the leading `H` operations in the CIGAR.
    /// For reverse-strand reads this equals the trailing `H` operations, because the
    /// stored sequence is the reverse complement and the CIGAR is written in reference
    /// order (trailing H = 5' end of the original molecule).
    pub(super) hard_clip_before: usize,
    pub(super) base_qual: u8,
    pub(super) map_qual: u8,
    pub(super) ab_count: Option<i32>,
    pub(super) ba_count: Option<i32>,
    pub(super) family_size: Option<i32>,
    pub(super) insert_size: Option<i32>,
    pub(super) n_before_alt: usize,
    pub(super) n_after_alt: usize,
    pub(super) n_n_before_alt: usize,
    pub(super) n_n_after_alt: usize,
    pub(super) leading_n_run_len: usize,
    pub(super) trailing_n_run_len: usize,
}

/// Position-level summary returned by `tally_pileup`.
pub(super) struct PileupResult {
    /// Per-alt-base tallies (only populated for bases that passed filters)
    pub(super) bases: HashMap<char, BaseTally>,
    pub(super) total_depth: i32,
    pub(super) fwd_depth: i32,
    pub(super) rev_depth: i32,
    /// Number of overlapping fragment pairs at this position (pair count, not read count)
    pub(super) overlap_depth: i32,
    /// Per-read details keyed by base (only populated when collect_reads is true;
    /// only non-ref, non-N bases are included)
    pub(super) read_details: HashMap<char, Vec<ReadDetail>>,
}

/// Per-base tally at a pileup position.
#[derive(Default)]
pub(super) struct BaseTally {
    pub(super) total: i32,
    pub(super) fwd: i32,
    pub(super) rev: i32,
    /// Overlapping pairs where both reads agree on this base
    pub(super) overlap_alt_agree: i32,
    /// Overlapping pairs where one read sees this base and the other sees something different
    pub(super) overlap_alt_disagree: i32,
}

/// Full per-read data collected during the first pileup pass.
struct LocusRead {
    base: char,
    is_reverse: bool,
    is_first_in_pair: bool,
    qpos: usize,
    read_len: usize,
    hard_clip_before: usize,
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
    n_before_alt: usize,
    n_after_alt: usize,
    n_n_before_alt: usize,
    n_n_after_alt: usize,
    leading_n_run_len: usize,
    trailing_n_run_len: usize,
}

#[derive(Clone, Copy)]
pub(super) struct ReadContextMetrics {
    pub(super) n_before_alt: usize,
    pub(super) n_after_alt: usize,
    pub(super) n_n_before_alt: usize,
    pub(super) n_n_after_alt: usize,
    pub(super) leading_n_run_len: usize,
    pub(super) trailing_n_run_len: usize,
}

pub(super) fn read_context_metrics(bases: &[u8], qpos: usize) -> ReadContextMetrics {
    let n_before_alt = qpos;
    let n_after_alt = bases.len().saturating_sub(qpos + 1);
    let n_n_before_alt = bases[..qpos].iter().filter(|&&b| b == b'N').count();
    let n_n_after_alt = bases[qpos.saturating_add(1)..]
        .iter()
        .filter(|&&b| b == b'N')
        .count();

    let mut leading_n_run_len = 0;
    for &b in bases[..qpos].iter().rev() {
        if b == b'N' {
            leading_n_run_len += 1;
        } else {
            break;
        }
    }

    let mut trailing_n_run_len = 0;
    for &b in &bases[qpos.saturating_add(1)..] {
        if b == b'N' {
            trailing_n_run_len += 1;
        } else {
            break;
        }
    }

    ReadContextMetrics {
        n_before_alt,
        n_after_alt,
        n_n_before_alt,
        n_n_after_alt,
        leading_n_run_len,
        trailing_n_run_len,
    }
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
pub(super) fn tally_pileup(
    pileup: &rust_htslib::bam::pileup::Pileup,
    min_base_qual: u8,
    min_map_qual: u8,
    include_duplicates: bool,
    include_secondary: bool,
    include_supplementary: bool,
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

        let base = record.seq()[qpos].to_ascii_uppercase() as char;
        let is_reverse = record.is_reverse();
        let is_first_in_pair = record.flags() & 0x40 != 0;
        let map_qual = record.mapq();
        let read_len = record.seq_len();
        let context = if collect_reads {
            let seq_bases = (0..read_len)
                .map(|i| record.seq()[i].to_ascii_uppercase())
                .collect::<Vec<_>>();
            read_context_metrics(&seq_bases, qpos)
        } else {
            ReadContextMetrics {
                n_before_alt: 0,
                n_after_alt: 0,
                n_n_before_alt: 0,
                n_n_after_alt: 0,
                leading_n_run_len: 0,
                trailing_n_run_len: 0,
            }
        };
        let (hc_leading, hc_trailing) = hard_clip_counts(&record);
        let hard_clip_before = if is_reverse { hc_trailing } else { hc_leading };

        let (ab_count, ba_count, family_size, insert_size) = if collect_reads {
            let ab = aux_i32(&record, b"aD");
            let ba = aux_i32(&record, b"bD");
            let fs = aux_i32(&record, b"cD");
            let tlen = record.insert_size();
            let ins = if tlen == 0 {
                None
            } else {
                Some(tlen.unsigned_abs() as i32)
            };
            (ab, ba, fs, ins)
        } else {
            (None, None, None, None)
        };

        by_qname
            .entry(record.qname().to_vec())
            .or_default()
            .push(LocusRead {
                base,
                is_reverse,
                is_first_in_pair,
                qpos,
                read_len,
                hard_clip_before,
                base_qual,
                map_qual,
                ab_count,
                ba_count,
                family_size,
                insert_size,
                n_before_alt: context.n_before_alt,
                n_after_alt: context.n_after_alt,
                n_n_before_alt: context.n_n_before_alt,
                n_n_after_alt: context.n_n_after_alt,
                leading_n_run_len: context.leading_n_run_len,
                trailing_n_run_len: context.trailing_n_run_len,
            });
    }

    // Second pass: tally with overlap detection.
    let mut bases: HashMap<char, BaseTally> = HashMap::new();
    let mut read_details: HashMap<char, Vec<ReadDetail>> = HashMap::new();
    let mut total_depth: i32 = 0;
    let mut fwd_depth: i32 = 0;
    let mut rev_depth: i32 = 0;
    let mut overlap_depth: i32 = 0;

    macro_rules! push_detail {
        ($b:expr, $r:expr) => {
            if collect_reads && $b != ref_base && $b != 'N' {
                read_details.entry($b).or_default().push(ReadDetail {
                    qpos: $r.qpos,
                    read_len: $r.read_len,
                    is_first_in_pair: $r.is_first_in_pair,
                    is_reverse: $r.is_reverse,
                    hard_clip_before: $r.hard_clip_before,
                    base_qual: $r.base_qual,
                    map_qual: $r.map_qual,
                    ab_count: $r.ab_count,
                    ba_count: $r.ba_count,
                    family_size: $r.family_size,
                    insert_size: $r.insert_size,
                    n_before_alt: $r.n_before_alt,
                    n_after_alt: $r.n_after_alt,
                    n_n_before_alt: $r.n_n_before_alt,
                    n_n_after_alt: $r.n_n_after_alt,
                    leading_n_run_len: $r.leading_n_run_len,
                    trailing_n_run_len: $r.trailing_n_run_len,
                });
            }
        };
    }

    for reads in by_qname.values() {
        match reads.as_slice() {
            [r] => {
                if r.base == 'N' {
                    continue;
                }
                total_depth += 1;
                if r.is_reverse {
                    rev_depth += 1;
                } else {
                    fwd_depth += 1;
                }

                let t = bases.entry(r.base).or_default();
                t.total += 1;
                if r.is_reverse {
                    t.rev += 1;
                } else {
                    t.fwd += 1;
                }

                push_detail!(r.base, r);
            }
            [r1, r2] => {
                overlap_depth += 1;
                total_depth += 1;

                let r1_is_rev = if r1.is_first_in_pair {
                    r1.is_reverse
                } else {
                    r2.is_reverse
                };
                if r1_is_rev {
                    rev_depth += 1;
                } else {
                    fwd_depth += 1;
                }

                let b1 = r1.base;
                let b2 = r2.base;

                if b1 == 'N' && b2 == 'N' {
                    // N + N: fragment counted in depth but no base tally.
                } else if b1 == 'N' {
                    let t = bases.entry(b2).or_default();
                    t.total += 1;
                    if r2.is_reverse {
                        t.rev += 1;
                    } else {
                        t.fwd += 1;
                    }
                    push_detail!(b2, r2);
                } else if b2 == 'N' {
                    let t = bases.entry(b1).or_default();
                    t.total += 1;
                    if r1.is_reverse {
                        t.rev += 1;
                    } else {
                        t.fwd += 1;
                    }
                    push_detail!(b1, r1);
                } else if b1 == b2 {
                    let t = bases.entry(b1).or_default();
                    t.total += 1;
                    if r1.is_reverse {
                        t.rev += 1;
                    } else {
                        t.fwd += 1;
                    }
                    if r2.is_reverse {
                        t.rev += 1;
                    } else {
                        t.fwd += 1;
                    }
                    t.overlap_alt_agree += 1;
                    let det = if r1.is_first_in_pair { r1 } else { r2 };
                    push_detail!(b1, det);
                } else if b1 == ref_base {
                    let t = bases.entry(b2).or_default();
                    t.total += 1;
                    if r2.is_reverse {
                        t.rev += 1;
                    } else {
                        t.fwd += 1;
                    }
                    t.overlap_alt_disagree += 1;
                    push_detail!(b2, r2);
                } else if b2 == ref_base {
                    let t = bases.entry(b1).or_default();
                    t.total += 1;
                    if r1.is_reverse {
                        t.rev += 1;
                    } else {
                        t.fwd += 1;
                    }
                    t.overlap_alt_disagree += 1;
                    push_detail!(b1, r1);
                } else {
                    let t1 = bases.entry(b1).or_default();
                    t1.total += 1;
                    if r1.is_reverse {
                        t1.rev += 1;
                    } else {
                        t1.fwd += 1;
                    }
                    t1.overlap_alt_disagree += 1;
                    push_detail!(b1, r1);

                    let t2 = bases.entry(b2).or_default();
                    t2.total += 1;
                    if r2.is_reverse {
                        t2.rev += 1;
                    } else {
                        t2.fwd += 1;
                    }
                    t2.overlap_alt_disagree += 1;
                    push_detail!(b2, r2);
                }
            }
            _ => {
                for r in reads {
                    if r.base == 'N' {
                        continue;
                    }
                    total_depth += 1;
                    if r.is_reverse {
                        rev_depth += 1;
                    } else {
                        fwd_depth += 1;
                    }
                    let t = bases.entry(r.base).or_default();
                    t.total += 1;
                    if r.is_reverse {
                        t.rev += 1;
                    } else {
                        t.fwd += 1;
                    }
                    push_detail!(r.base, r);
                }
            }
        }
    }

    PileupResult {
        bases,
        total_depth,
        fwd_depth,
        rev_depth,
        overlap_depth,
        read_details,
    }
}

/// Read an integer auxiliary tag from a BAM record, returning None if absent or non-integer.
///
/// Float tags are intentionally excluded. For fgbio simplex consensus reads, `cE` is the
/// per-base error rate (float), not a family size count. For duplex consensus reads, `cE`
/// is the BA strand count (integer) and will be read correctly here.
pub(super) fn aux_i32(record: &bam::Record, tag: &[u8; 2]) -> Option<i32> {
    match record.aux(tag) {
        Ok(bam::record::Aux::I8(v)) => Some(v as i32),
        Ok(bam::record::Aux::U8(v)) => Some(v as i32),
        Ok(bam::record::Aux::I16(v)) => Some(v as i32),
        Ok(bam::record::Aux::U16(v)) => Some(v as i32),
        Ok(bam::record::Aux::I32(v)) => Some(v),
        Ok(bam::record::Aux::U32(v)) => Some(v as i32),
        _ => None,
    }
}

/// Returns `(leading_hard_clips, trailing_hard_clips)` for the record's CIGAR string.
///
/// Hard-clipped bases are not stored in the BAM sequence, so `seq_len()` and `qpos` do
/// not account for them. Knowing how many bases were hard-clipped at each end is
/// required to compute the correct sequencing cycle number.
pub(super) fn hard_clip_counts(record: &bam::Record) -> (usize, usize) {
    use rust_htslib::bam::record::Cigar;

    let cigar = record.cigar();
    let cigar_slice: &[Cigar] = &cigar;
    let leading = match cigar_slice.first() {
        Some(Cigar::HardClip(n)) => *n as usize,
        _ => 0,
    };
    let trailing = if cigar_slice.len() > 1 {
        match cigar_slice.last() {
            Some(Cigar::HardClip(n)) => *n as usize,
            _ => 0,
        }
    } else {
        0
    };
    (leading, trailing)
}

/// Compute the 1-based sequencing cycle for an alt-supporting base.
///
/// The cycle number reflects the position in synthesis order (cycle 1 = first base
/// synthesized by the polymerase).
///
/// For **forward-strand** reads the stored sequence is in synthesis order, so:
///   `cycle = hard_clip_before + qpos + 1`
///
/// For **reverse-strand** reads BAM stores the reverse complement. `qpos = 0` is the
/// last synthesized base (3' end of synthesis), so synthesis order runs backwards:
///   `cycle = hard_clip_before + read_len − qpos`
///
/// `hard_clip_before` is the number of hard-clipped bases at the 5' end of the
/// original molecule (leading H for forward reads, trailing H for reverse reads).
#[inline]
pub(super) fn true_cycle(
    qpos: usize,
    read_len: usize,
    is_reverse: bool,
    hard_clip_before: usize,
) -> i32 {
    if is_reverse {
        (hard_clip_before + read_len - qpos) as i32
    } else {
        (hard_clip_before + qpos + 1) as i32
    }
}

#[cfg(test)]
mod tests {
    use super::read_context_metrics;

    #[test]
    fn read_context_metrics_handles_trailing_ns() {
        let m = read_context_metrics(b"ATGNNAA", 2);
        assert_eq!(m.n_before_alt, 2);
        assert_eq!(m.n_after_alt, 4);
        assert_eq!(m.n_n_before_alt, 0);
        assert_eq!(m.n_n_after_alt, 2);
        assert_eq!(m.leading_n_run_len, 0);
        assert_eq!(m.trailing_n_run_len, 2);
    }

    #[test]
    fn read_context_metrics_handles_leading_ns() {
        let m = read_context_metrics(b"ANNTAAA", 3);
        assert_eq!(m.n_before_alt, 3);
        assert_eq!(m.n_after_alt, 3);
        assert_eq!(m.n_n_before_alt, 2);
        assert_eq!(m.n_n_after_alt, 0);
        assert_eq!(m.leading_n_run_len, 2);
        assert_eq!(m.trailing_n_run_len, 0);
    }

    #[test]
    fn read_context_metrics_handles_non_adjacent_ns() {
        let m = read_context_metrics(b"NATANAA", 3);
        assert_eq!(m.n_n_before_alt, 1);
        assert_eq!(m.n_n_after_alt, 1);
        assert_eq!(m.leading_n_run_len, 0);
        assert_eq!(m.trailing_n_run_len, 1);
    }

    #[test]
    fn read_context_metrics_handles_alt_at_read_edges() {
        let start = read_context_metrics(b"TNNAA", 0);
        assert_eq!(start.n_before_alt, 0);
        assert_eq!(start.n_after_alt, 4);
        assert_eq!(start.leading_n_run_len, 0);
        assert_eq!(start.trailing_n_run_len, 2);

        let end = read_context_metrics(b"AANNT", 4);
        assert_eq!(end.n_before_alt, 4);
        assert_eq!(end.n_after_alt, 0);
        assert_eq!(end.leading_n_run_len, 2);
        assert_eq!(end.trailing_n_run_len, 0);
    }
}
