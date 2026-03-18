use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

use crate::vcf::{VcfAnnotation, VariantAnnotator};

/// Pre-loaded lookup built from a tab-separated variant list.
///
/// Expected columns (tab-separated, optional header line is auto-detected):
///   chrom   pos_start   pos_end   ref   var
///
/// Coordinates are 0-based half-open (BED convention): a SNV at chromosome
/// position N has pos_start = N, pos_end = N+1.  All matched loci are
/// annotated as variant_called = true, variant_filter = "PASS".
pub struct VariantsTsv {
    by_allele:   HashMap<(String, i64, String), VcfAnnotation>,
    by_position: HashMap<(String, i64), VcfAnnotation>,
}

impl VariantsTsv {
    pub fn load(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read variants TSV: {}", path.display()))?;

        let annotation = VcfAnnotation { filter: "PASS".to_string() };

        let mut by_allele:   HashMap<(String, i64, String), VcfAnnotation> = HashMap::new();
        let mut by_position: HashMap<(String, i64), VcfAnnotation>         = HashMap::new();

        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            let cols: Vec<&str> = line.splitn(6, '\t').collect();
            if cols.len() < 5 {
                continue;
            }

            // Auto-detect and skip a header line (pos_start column is not numeric)
            let pos_start: i64 = match cols[1].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };

            let chrom = cols[0].to_string();
            let var   = cols[4].to_string();

            by_position
                .entry((chrom.clone(), pos_start))
                .or_insert_with(|| annotation.clone());

            by_allele
                .entry((chrom, pos_start, var))
                .or_insert_with(|| annotation.clone());
        }

        Ok(Self { by_allele, by_position })
    }
}

impl VariantAnnotator for VariantsTsv {
    /// For SNVs: tries exact (chrom, pos_start, var) match, then position fallback.
    /// For indels: position-only match (indel representations differ from VCF/TSV notation).
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
