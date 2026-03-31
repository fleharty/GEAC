use std::collections::HashMap;

use rust_htslib::bam::pileup::Indel;

use crate::record::VariantType;

use super::pileup::{aux_i32, hard_clip_counts, ReadDetail};

/// Per-indel-allele tally at a pileup position.
pub(super) struct IndelCount {
    pub(super) ref_allele: String,
    pub(super) alt_allele: String,
    pub(super) variant_type: VariantType,
    pub(super) total: i32,
    pub(super) fwd: i32,
    pub(super) rev: i32,
    /// Overlapping pairs where both reads agree on this indel allele
    pub(super) overlap_alt_agree: i32,
    /// Overlapping pairs where reads disagree (one has this indel, the other differs or has none)
    pub(super) overlap_alt_disagree: i32,
}

/// Decoded indel allele for one read at the anchor position, or None if no indel.
type IndelAllele = Option<(String, String, VariantType)>; // (alt_allele, ref_allele, variant_type)

/// Compare two indel allele strings treating N (in either string) as a wildcard.
/// Returns true if the alleles are compatible (same length, every position matches or is N).
pub(super) fn indels_compatible(a: &str, b: &str) -> bool {
    if a == b {
        return true;
    }
    let ab = a.as_bytes();
    let bb = b.as_bytes();
    if ab.len() != bb.len() {
        return false;
    }
    ab.iter()
        .zip(bb.iter())
        .all(|(x, y)| *x == b'N' || *y == b'N' || x == y)
}

/// Tally indel alleles at a pileup column with overlap detection.
///
/// Groups reads by query name. A name appearing twice means both reads of the
/// fragment are at the anchor position — their indel alleles are compared to
/// determine agreement or disagreement.
///
/// Every read that passes the mapping quality filter pushes an entry (Some or None)
/// so that overlapping pairs are correctly identified even when one read has no indel.
pub(super) fn tally_indels(
    pileup: &rust_htslib::bam::pileup::Pileup,
    pos: i64,
    chrom_seq: &[u8],
    min_map_qual: u8,
    include_duplicates: bool,
    include_secondary: bool,
    include_supplementary: bool,
    collect_reads: bool,
) -> (
    HashMap<String, IndelCount>,
    HashMap<String, Vec<ReadDetail>>,
) {
    let mut by_qname: HashMap<Vec<u8>, Vec<(IndelAllele, bool, bool, Option<ReadDetail>)>> =
        HashMap::new();

    for alignment in pileup.alignments() {
        if alignment.is_refskip() {
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

        let is_reverse = record.is_reverse();
        let is_first_in_pair = record.flags() & 0x40 != 0;
        let qname = record.qname().to_vec();

        let allele: IndelAllele = match alignment.indel() {
            Indel::Ins(len) => {
                let Some(qpos) = alignment.qpos() else {
                    by_qname.entry(qname).or_default().push((
                        None,
                        is_reverse,
                        is_first_in_pair,
                        None,
                    ));
                    continue;
                };
                let seq = record.seq();
                let len = len as usize;
                if qpos + len >= seq.len() {
                    by_qname.entry(qname).or_default().push((
                        None,
                        is_reverse,
                        is_first_in_pair,
                        None,
                    ));
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
                    by_qname.entry(qname).or_default().push((
                        None,
                        is_reverse,
                        is_first_in_pair,
                        None,
                    ));
                    continue;
                }
                Some((
                    format!("-{deleted}"),
                    deleted.clone(),
                    VariantType::Deletion,
                ))
            }
            Indel::None => None,
        };

        let detail: Option<ReadDetail> = if collect_reads && allele.is_some() {
            let qpos = alignment.qpos().unwrap_or(0);
            let qual = record.qual();
            let tlen = record.insert_size();
            let (hc_leading, hc_trailing) = hard_clip_counts(&record);
            let hard_clip_before = if is_reverse { hc_trailing } else { hc_leading };
            Some(ReadDetail {
                qpos,
                read_len: record.seq_len(),
                is_first_in_pair,
                is_reverse,
                hard_clip_before,
                base_qual: qual.get(qpos).copied().unwrap_or(0),
                map_qual: record.mapq(),
                ab_count: aux_i32(&record, b"aD"),
                ba_count: aux_i32(&record, b"bD"),
                family_size: aux_i32(&record, b"cD"),
                insert_size: if tlen == 0 {
                    None
                } else {
                    Some(tlen.unsigned_abs() as i32)
                },
            })
        } else {
            None
        };

        by_qname
            .entry(qname)
            .or_default()
            .push((allele, is_reverse, is_first_in_pair, detail));
    }

    let mut indels: HashMap<String, IndelCount> = HashMap::new();
    let mut indel_details: HashMap<String, Vec<ReadDetail>> = HashMap::new();

    for reads in by_qname.values() {
        match reads.as_slice() {
            [(allele, is_rev, _, detail)] => {
                if let Some((alt, ref_a, vt)) = allele {
                    let e = indels.entry(alt.clone()).or_insert_with(|| IndelCount {
                        ref_allele: ref_a.clone(),
                        alt_allele: alt.clone(),
                        variant_type: *vt,
                        total: 0,
                        fwd: 0,
                        rev: 0,
                        overlap_alt_agree: 0,
                        overlap_alt_disagree: 0,
                    });
                    e.total += 1;
                    if *is_rev {
                        e.rev += 1;
                    } else {
                        e.fwd += 1;
                    }
                    if let Some(d) = detail {
                        indel_details
                            .entry(alt.clone())
                            .or_default()
                            .push(d.clone());
                    }
                }
            }
            [(allele1, is_rev1, is_first1, detail1), (allele2, is_rev2, _is_first2, detail2)] => {
                let r1_detail = if *is_first1 { detail1 } else { detail2 };

                match (allele1, allele2) {
                    (None, None) => {}
                    (Some((alt, ref_a, vt)), None) => {
                        let e = indels.entry(alt.clone()).or_insert_with(|| IndelCount {
                            ref_allele: ref_a.clone(),
                            alt_allele: alt.clone(),
                            variant_type: *vt,
                            total: 0,
                            fwd: 0,
                            rev: 0,
                            overlap_alt_agree: 0,
                            overlap_alt_disagree: 0,
                        });
                        e.total += 1;
                        if *is_rev1 {
                            e.rev += 1;
                        } else {
                            e.fwd += 1;
                        }
                        e.overlap_alt_disagree += 1;
                        if let Some(d) = detail1 {
                            indel_details
                                .entry(alt.clone())
                                .or_default()
                                .push(d.clone());
                        }
                    }
                    (None, Some((alt, ref_a, vt))) => {
                        let e = indels.entry(alt.clone()).or_insert_with(|| IndelCount {
                            ref_allele: ref_a.clone(),
                            alt_allele: alt.clone(),
                            variant_type: *vt,
                            total: 0,
                            fwd: 0,
                            rev: 0,
                            overlap_alt_agree: 0,
                            overlap_alt_disagree: 0,
                        });
                        e.total += 1;
                        if *is_rev2 {
                            e.rev += 1;
                        } else {
                            e.fwd += 1;
                        }
                        e.overlap_alt_disagree += 1;
                        if let Some(d) = detail2 {
                            indel_details
                                .entry(alt.clone())
                                .or_default()
                                .push(d.clone());
                        }
                    }
                    (Some((alt1, ref_a1, vt1)), Some((alt2, ref_a2, vt2))) => {
                        if indels_compatible(alt1, alt2) {
                            let (canon_alt, canon_ref, canon_vt) =
                                if alt1.contains('N') && !alt2.contains('N') {
                                    (&alt2, &ref_a2, vt2)
                                } else {
                                    (&alt1, &ref_a1, vt1)
                                };
                            let e =
                                indels
                                    .entry(canon_alt.to_string())
                                    .or_insert_with(|| IndelCount {
                                        ref_allele: canon_ref.to_string(),
                                        alt_allele: canon_alt.to_string(),
                                        variant_type: *canon_vt,
                                        total: 0,
                                        fwd: 0,
                                        rev: 0,
                                        overlap_alt_agree: 0,
                                        overlap_alt_disagree: 0,
                                    });
                            e.total += 1;
                            if *is_rev1 {
                                e.rev += 1;
                            } else {
                                e.fwd += 1;
                            }
                            if *is_rev2 {
                                e.rev += 1;
                            } else {
                                e.fwd += 1;
                            }
                            e.overlap_alt_agree += 1;
                            if let Some(d) = r1_detail {
                                indel_details
                                    .entry(canon_alt.to_string())
                                    .or_default()
                                    .push(d.clone());
                            }
                        } else {
                            let e1 = indels.entry(alt1.clone()).or_insert_with(|| IndelCount {
                                ref_allele: ref_a1.clone(),
                                alt_allele: alt1.clone(),
                                variant_type: *vt1,
                                total: 0,
                                fwd: 0,
                                rev: 0,
                                overlap_alt_agree: 0,
                                overlap_alt_disagree: 0,
                            });
                            e1.total += 1;
                            if *is_rev1 {
                                e1.rev += 1;
                            } else {
                                e1.fwd += 1;
                            }
                            e1.overlap_alt_disagree += 1;
                            if let Some(d) = detail1 {
                                indel_details
                                    .entry(alt1.clone())
                                    .or_default()
                                    .push(d.clone());
                            }

                            let e2 = indels.entry(alt2.clone()).or_insert_with(|| IndelCount {
                                ref_allele: ref_a2.clone(),
                                alt_allele: alt2.clone(),
                                variant_type: *vt2,
                                total: 0,
                                fwd: 0,
                                rev: 0,
                                overlap_alt_agree: 0,
                                overlap_alt_disagree: 0,
                            });
                            e2.total += 1;
                            if *is_rev2 {
                                e2.rev += 1;
                            } else {
                                e2.fwd += 1;
                            }
                            e2.overlap_alt_disagree += 1;
                            if let Some(d) = detail2 {
                                indel_details
                                    .entry(alt2.clone())
                                    .or_default()
                                    .push(d.clone());
                            }
                        }
                    }
                }
            }
            _ => {
                for (allele, is_rev, _, detail) in reads {
                    if let Some((alt, ref_a, vt)) = allele {
                        let e = indels.entry(alt.clone()).or_insert_with(|| IndelCount {
                            ref_allele: ref_a.clone(),
                            alt_allele: alt.clone(),
                            variant_type: *vt,
                            total: 0,
                            fwd: 0,
                            rev: 0,
                            overlap_alt_agree: 0,
                            overlap_alt_disagree: 0,
                        });
                        e.total += 1;
                        if *is_rev {
                            e.rev += 1;
                        } else {
                            e.fwd += 1;
                        }
                        if let Some(d) = detail {
                            indel_details
                                .entry(alt.clone())
                                .or_default()
                                .push(d.clone());
                        }
                    }
                }
            }
        }
    }

    (indels, indel_details)
}

#[cfg(test)]
mod tests {
    use super::indels_compatible;

    #[test]
    fn n_in_second_allele_is_compatible() {
        assert!(indels_compatible("GACTT", "GACTN"));
    }

    #[test]
    fn n_in_first_allele_is_compatible() {
        assert!(indels_compatible("GACTN", "GACTT"));
    }

    #[test]
    fn n_in_both_alleles_is_compatible() {
        assert!(indels_compatible("GNCTN", "GACTN"));
    }

    #[test]
    fn identical_alleles_are_compatible() {
        assert!(indels_compatible("GACTT", "GACTT"));
    }

    #[test]
    fn different_alleles_are_incompatible() {
        assert!(!indels_compatible("GACTT", "GACTA"));
    }

    #[test]
    fn different_lengths_are_incompatible() {
        assert!(!indels_compatible("GACTT", "GAC"));
    }

    #[test]
    fn all_n_second_is_compatible() {
        assert!(indels_compatible("GACTT", "NNNNN"));
    }
}
