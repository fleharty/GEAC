//! gnomAD VCF annotation — allele-frequency lookup via tabix-indexed VCF/BCF.
//!
//! The file must be bgzip-compressed and indexed (`.tbi` or `.csi`).
//! gnomAD distributes its VCFs in exactly this format, so no preprocessing is needed.
//!
//! Usage:
//! ```text
//! let mut index = GnomadIndex::open(path, None)?;
//! if let Some(af) = index.get("chr1", 12345, "A", "T")? {
//!     println!("gnomAD AF = {af}");
//! }
//! ```

use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bcf::{self, Read as BcfRead};

/// Tabix-indexed gnomAD VCF reader for per-locus allele-frequency lookup.
pub struct GnomadIndex {
    reader: bcf::IndexedReader,
    /// INFO field to extract as AF (e.g. `b"AF"`).
    af_field: Vec<u8>,
}

impl GnomadIndex {
    /// Open a bgzip+tabix-indexed gnomAD VCF or BCF file.
    ///
    /// `af_field` is the INFO key to use for allele frequency; defaults to `"AF"`.
    pub fn open(path: &Path, af_field: Option<&str>) -> Result<Self> {
        let reader = bcf::IndexedReader::from_path(path).with_context(|| {
            format!(
                "failed to open gnomAD VCF '{}' — file must be bgzip-compressed \
                 and tabix/CSI-indexed (.tbi or .csi index must exist alongside the file)",
                path.display()
            )
        })?;
        Ok(Self {
            reader,
            af_field: af_field.unwrap_or("AF").as_bytes().to_vec(),
        })
    }

    /// Look up the allele frequency for an alt allele at a locus.
    ///
    /// `pos` is **0-based** (matching GEAC's internal representation).
    /// Returns `Some(af)` when an exact-allele match is found and AF is available,
    /// `None` when the allele is absent from gnomAD or AF cannot be read.
    ///
    /// For indels (alt_allele starts with `+` or `-`), GEAC's representation is
    /// translated to VCF form (anchor + deleted/inserted bases) before matching.
    ///
    /// `ref_allele` must be the single anchor base at `pos` (not the deleted
    /// bases). Callers pass the reference base from the pileup position.
    pub fn get(
        &mut self,
        chrom: &str,
        pos: i64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<f32>> {
        // Clone header and af_field before the mutable borrow for records().
        let header = self.reader.header().clone();
        let af_field = self.af_field.clone();

        let rid = match Self::resolve_rid(&header, chrom) {
            Some(r) => r,
            None => return Ok(None),
        };

        // Fetch: 0-based half-open [pos, pos+1). Errors here mean the index can't
        // seek to this chromosome/position — treat as not found rather than fatal.
        if self
            .reader
            .fetch(rid, pos as u64, Some(pos as u64 + 1))
            .is_err()
        {
            return Ok(None);
        }

        let (expected_ref, expected_alt) = geac_to_vcf_alleles(ref_allele, alt_allele);

        for result in self.reader.records() {
            let record = result.context("error reading gnomAD VCF record")?;

            if record.pos() != pos {
                continue;
            }

            let alleles = record.alleles();
            if alleles.is_empty() {
                continue;
            }

            let vcf_ref = std::str::from_utf8(alleles[0]).unwrap_or(".");
            if vcf_ref != expected_ref {
                continue;
            }

            if let Some(alt_idx) = alleles[1..]
                .iter()
                .position(|a| std::str::from_utf8(a).unwrap_or(".") == expected_alt)
            {
                return Ok(Self::extract_af(&record, &af_field, alt_idx));
            }
        }

        Ok(None)
    }

    /// Try to resolve a chromosome name to a VCF RID, accounting for `chr` prefix mismatches.
    fn resolve_rid(header: &bcf::header::HeaderView, chrom: &str) -> Option<u32> {
        if let Ok(rid) = header.name2rid(chrom.as_bytes()) {
            return Some(rid);
        }
        // Toggle "chr" prefix to handle hg19/hg38 naming mismatches.
        let alt = if chrom.starts_with("chr") {
            chrom.trim_start_matches("chr").to_string()
        } else {
            format!("chr{}", chrom)
        };
        header.name2rid(alt.as_bytes()).ok()
    }

    /// Extract the AF value for `alt_idx` (0-based alt index) from a VCF record's INFO.
    fn extract_af(record: &bcf::Record, af_field: &[u8], alt_idx: usize) -> Option<f32> {
        match record.info(af_field).float() {
            Ok(Some(vals)) => vals.get(alt_idx).copied().filter(|v| !v.is_nan()),
            _ => None,
        }
    }
}

/// Translate GEAC's `(ref_allele, alt_allele)` pair into the `(VCF REF, VCF ALT)`
/// pair that should match the gnomAD record at the same anchor position.
///
/// `ref_allele` must be the single anchor base at the position (GEAC's
/// convention for every AltBase row). GEAC indel alt alleles are prefixed:
///
///   SNV        ref="A",  alt="T"    -> vcf REF="A",    ALT="T"
///   Insertion  ref="A",  alt="+AC"  -> vcf REF="A",    ALT="AAC"
///   Deletion   ref="A",  alt="-TCG" -> vcf REF="ATCG", ALT="A"
fn geac_to_vcf_alleles(ref_allele: &str, alt_allele: &str) -> (String, String) {
    if let Some(inserted) = alt_allele.strip_prefix('+') {
        (ref_allele.to_string(), format!("{ref_allele}{inserted}"))
    } else if let Some(deleted) = alt_allele.strip_prefix('-') {
        (format!("{ref_allele}{deleted}"), ref_allele.to_string())
    } else {
        (ref_allele.to_string(), alt_allele.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::geac_to_vcf_alleles;

    #[test]
    fn snv_passes_through_unchanged() {
        assert_eq!(geac_to_vcf_alleles("A", "T"), ("A".into(), "T".into()));
    }

    #[test]
    fn insertion_concatenates_anchor_with_inserted_bases() {
        // GEAC "+AC" at anchor "A" -> VCF REF="A", ALT="AAC"
        assert_eq!(geac_to_vcf_alleles("A", "+AC"), ("A".into(), "AAC".into()));
    }

    #[test]
    fn deletion_concatenates_anchor_with_deleted_bases() {
        // GEAC "-TCG" at anchor "A" -> VCF REF="ATCG", ALT="A"
        // This is the case previously broken: the caller used to pass the
        // deleted bases as ref_allele, producing "TCGTCG" for the VCF REF.
        assert_eq!(
            geac_to_vcf_alleles("A", "-TCG"),
            ("ATCG".into(), "A".into())
        );
    }

    #[test]
    fn single_base_deletion() {
        assert_eq!(geac_to_vcf_alleles("A", "-T"), ("AT".into(), "A".into()));
    }
}
