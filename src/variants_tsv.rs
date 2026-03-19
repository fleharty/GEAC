use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

use crate::vcf::{VcfAnnotation, VariantAnnotator};

/// Pre-loaded lookup built from a tab-separated variant list.
///
/// Expected columns (tab-separated):
///   chrom   pos_start   pos_end   ref   var   [...optional additional columns...]
///
/// Coordinates are 0-based half-open (BED convention): a SNV at chromosome
/// position N has pos_start = N, pos_end = N+1.  Matched loci are annotated
/// as variant_called = true.  If a header row is present and contains a
/// "post_filter_flag" column, that column's value is used as variant_filter;
/// "NoRules" is treated as "PASS".  If the column is absent or there is no
/// header, all variants are annotated as "PASS".
pub struct VariantsTsv {
    by_allele:   HashMap<(String, i64, String), VcfAnnotation>,
    by_position: HashMap<(String, i64), VcfAnnotation>,
}

impl VariantsTsv {
    pub fn load(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read variants TSV: {}", path.display()))?;

        let mut by_allele:   HashMap<(String, i64, String), VcfAnnotation> = HashMap::new();
        let mut by_position: HashMap<(String, i64), VcfAnnotation>         = HashMap::new();

        // Find the post_filter_flag column index from the header row, if present.
        // The header is identified by a non-numeric value in the pos_start column (col 1).
        let mut lines = content.lines().peekable();
        let filter_col: Option<usize> = match lines.peek() {
            Some(first) if first.split('\t').nth(1).map(|v| v.parse::<i64>().is_err()).unwrap_or(false) => {
                first.split('\t').position(|h| h.trim() == "post_filter_flag")
            }
            _ => None,
        };

        for line in lines {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 5 {
                continue;
            }

            // Skip the header row
            let pos_start: i64 = match cols[1].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };

            let chrom = cols[0].to_string();
            let var   = cols[4].to_string();

            // Use post_filter_flag column if we found it in the header.
            // "NoRules" means no filter rule fired — equivalent to PASS.
            let filter = match filter_col.and_then(|i| cols.get(i)).map(|s| s.trim()) {
                None | Some("") => "PASS".to_string(),
                Some("NoRules") => "PASS".to_string(),
                Some(f)         => f.to_string(),
            };

            let annotation = VcfAnnotation { filter };

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
