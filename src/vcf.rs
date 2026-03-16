use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bcf::{self, Read as BcfRead};
use rust_htslib::bcf::header::HeaderView;

/// Annotation attached to a called variant locus.
#[derive(Debug, Clone)]
pub struct VcfAnnotation {
    /// FILTER field: "PASS", semicolon-separated filter reasons, or "." if unknown.
    pub filter: String,
}

/// Common interface for variant annotation sources (VCF or TSV).
pub trait VariantAnnotator {
    /// Return annotation for the given locus and alt allele, or `None` if not found.
    fn get(&self, chrom: &str, pos: i64, alt_allele: &str) -> Option<&VcfAnnotation>;
}

/// Pre-loaded VCF lookup structures built before the pileup loop.
pub struct VcfIndex {
    /// Exact (chrom, 0-based pos, alt_allele) → annotation.
    /// Used for SNV matching.
    by_allele: HashMap<(String, i64, String), VcfAnnotation>,

    /// (chrom, 0-based pos) → annotation for the first record at that position.
    /// Used as a fallback for indels, whose allele representations differ from VCF.
    by_position: HashMap<(String, i64), VcfAnnotation>,
}

impl VcfIndex {
    /// Load a VCF or BCF file into memory.
    pub fn load(vcf_path: &Path) -> Result<Self> {
        let mut reader = bcf::Reader::from_path(vcf_path)
            .with_context(|| format!("failed to open VCF: {}", vcf_path.display()))?;

        let mut by_allele: HashMap<(String, i64, String), VcfAnnotation> = HashMap::new();
        let mut by_position: HashMap<(String, i64), VcfAnnotation> = HashMap::new();

        // Clone the HeaderView so we can use it after the reader is consumed by records()
        let header = reader.header().clone();

        for result in reader.records() {
            let record = result.context("error reading VCF record")?;

            let rid = record.rid().context("VCF record missing CHROM")?;
            let chrom = std::str::from_utf8(header.rid2name(rid)?)
                .context("VCF chromosome name is not valid UTF-8")?
                .to_string();

            // htslib stores VCF positions as 0-based internally
            let pos = record.pos();

            let filter = filter_string(&record, &header);
            let annotation = VcfAnnotation { filter };

            // Position-level fallback (first record at each position wins)
            by_position
                .entry((chrom.clone(), pos))
                .or_insert_with(|| annotation.clone());

            // Allele-level exact match for each alt
            for alt_bytes in record.alleles().iter().skip(1) {
                let alt = std::str::from_utf8(alt_bytes)
                    .unwrap_or(".")
                    .to_string();
                by_allele
                    .entry((chrom.clone(), pos, alt))
                    .or_insert_with(|| annotation.clone());
            }
        }

        Ok(Self { by_allele, by_position })
    }

}

impl VariantAnnotator for VcfIndex {
    /// For SNVs: tries exact (chrom, pos, alt_allele) match first, then position fallback.
    /// For indels (alt starts with '+' or '-'): position-only match because VCF and our
    /// indel representations differ.
    fn get(&self, chrom: &str, pos: i64, alt_allele: &str) -> Option<&VcfAnnotation> {
        let is_indel = alt_allele.starts_with('+') || alt_allele.starts_with('-');

        if is_indel {
            self.by_position.get(&(chrom.to_string(), pos))
        } else {
            self.by_allele
                .get(&(chrom.to_string(), pos, alt_allele.to_string()))
                .or_else(|| self.by_position.get(&(chrom.to_string(), pos)))
        }
    }
}

/// Build a human-readable filter string from a BCF record.
/// Returns "PASS" if the variant passed, or semicolon-joined filter names otherwise.
fn filter_string(record: &bcf::Record, header: &HeaderView) -> String {
    if record.has_filter("PASS".as_bytes()) {
        return "PASS".to_string();
    }

    let names: Vec<String> = record
        .filters()
        .filter_map(|id| {
            std::str::from_utf8(&header.id_to_name(id))
                .ok()
                .map(|s| s.to_string())
        })
        .filter(|s| s != "PASS")
        .collect();

    if names.is_empty() {
        ".".to_string()
    } else {
        names.join(";")
    }
}
