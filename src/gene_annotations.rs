use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

// ── Data structures ────────────────────────────────────────────────────────────

struct GeneRecord {
    start: u32,
    end:   u32,
    gene:  String,
}

/// Per-chromosome gene intervals, sorted by start with a prefix-max of end values
/// to support efficient overlap queries even when genes overlap.
struct ChromGenes {
    intervals:      Vec<GeneRecord>,
    /// prefix_max_end[i] = max(end[0], end[1], …, end[i])
    prefix_max_end: Vec<u32>,
}

impl ChromGenes {
    fn from_unsorted(mut records: Vec<GeneRecord>) -> Self {
        records.sort_unstable_by_key(|r| r.start);
        let prefix_max_end = {
            let mut running_max = 0u32;
            records
                .iter()
                .map(|r| {
                    running_max = running_max.max(r.end);
                    running_max
                })
                .collect()
        };
        Self { intervals: records, prefix_max_end }
    }

    /// Returns the gene name for the first interval covering `pos` (0-based), or None.
    fn get(&self, pos: u32) -> Option<&str> {
        // All intervals with start <= pos are at indices 0..idx.
        let idx = self.intervals.partition_point(|r| r.start <= pos);
        if idx == 0 {
            return None;
        }
        // Quick exit: if the maximum end of all candidate intervals is <= pos, none can cover it.
        if self.prefix_max_end[idx - 1] <= pos {
            return None;
        }
        // Scan backward; stop early when no earlier interval could cover pos.
        for i in (0..idx).rev() {
            if self.intervals[i].end > pos {
                return Some(&self.intervals[i].gene);
            }
            if self.prefix_max_end[i] <= pos {
                break;
            }
        }
        None
    }
}

// ── Public API ─────────────────────────────────────────────────────────────────

/// Pre-loaded gene interval lookup built from a GFF3 or GTF annotation file.
///
/// Only `gene`-level features are loaded (not transcript/exon), so the returned
/// name represents the gene body. For loci that fall in intergenic regions, or
/// when no annotation file is provided, the `gene` column is null.
///
/// **Format auto-detection:**
/// - A `##gff-version 3` header line → GFF3
/// - File extension `.gff3` or `.gff` → GFF3
/// - Otherwise → GTF
///
/// **Gene name extraction:**
/// - GTF:  `gene_name "…"` → `gene_id "…"` (fallback)
/// - GFF3: `Name=…`        → `gene_name=…` → `gene_id=…` (fallback)
pub struct GeneAnnotations {
    by_chrom: HashMap<String, ChromGenes>,
}

impl GeneAnnotations {
    pub fn load(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read gene annotations: {}", path.display()))?;

        let fmt = detect_format(&content, path);
        let mut raw: HashMap<String, Vec<GeneRecord>> = HashMap::new();

        for line in content.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.splitn(9, '\t').collect();
            if fields.len() < 9 {
                continue;
            }

            // Only load gene-level features.
            if fields[2] != "gene" {
                continue;
            }

            // Both GFF3 and GTF use 1-based inclusive coordinates.
            let start: u32 = fields[3]
                .parse::<u32>()
                .with_context(|| format!("invalid start coordinate: '{}'", fields[3]))?
                .saturating_sub(1); // → 0-based
            let end: u32 = fields[4]
                .parse()
                .with_context(|| format!("invalid end coordinate: '{}'", fields[4]))?; // → 0-based exclusive

            let gene = match fmt {
                Format::Gtf  => parse_gtf_gene_name(fields[8]),
                Format::Gff3 => parse_gff3_gene_name(fields[8]),
            };

            if let Some(gene) = gene {
                raw.entry(fields[0].to_string())
                    .or_default()
                    .push(GeneRecord { start, end, gene });
            }
        }

        let by_chrom = raw
            .into_iter()
            .map(|(chrom, records)| (chrom, ChromGenes::from_unsorted(records)))
            .collect();

        Ok(Self { by_chrom })
    }

    /// Returns the gene name for the given 0-based locus, or None.
    pub fn get(&self, chrom: &str, pos: i64) -> Option<&str> {
        self.by_chrom.get(chrom)?.get(pos as u32)
    }

    pub fn n_genes(&self) -> usize {
        self.by_chrom.values().map(|c| c.intervals.len()).sum()
    }
}

// ── Format detection ───────────────────────────────────────────────────────────

enum Format { Gtf, Gff3 }

fn detect_format(content: &str, path: &Path) -> Format {
    // Check for a GFF3 version header in the first few lines.
    for line in content.lines().take(10) {
        if line.starts_with("##gff-version") {
            return Format::Gff3;
        }
    }
    // Fall back to file extension.
    match path.extension().and_then(|e| e.to_str()).unwrap_or("") {
        "gff3" | "gff" => Format::Gff3,
        _ => Format::Gtf,
    }
}

// ── GTF attribute parsing ──────────────────────────────────────────────────────

fn parse_gtf_gene_name(attrs: &str) -> Option<String> {
    extract_gtf_attr(attrs, "gene_name")
        .or_else(|| extract_gtf_attr(attrs, "gene_id"))
}

/// Extract a quoted GTF attribute value: `key "value";`
fn extract_gtf_attr(attrs: &str, key: &str) -> Option<String> {
    let needle = format!("{} \"", key);
    let start = attrs.find(&needle)? + needle.len();
    let end = attrs[start..].find('"')? + start;
    Some(attrs[start..end].to_string())
}

// ── GFF3 attribute parsing ─────────────────────────────────────────────────────

fn parse_gff3_gene_name(attrs: &str) -> Option<String> {
    extract_gff3_attr(attrs, "Name")
        .or_else(|| extract_gff3_attr(attrs, "gene_name"))
        .or_else(|| extract_gff3_attr(attrs, "gene_id"))
}

/// Extract a GFF3 attribute value: `key=value;`
fn extract_gff3_attr(attrs: &str, key: &str) -> Option<String> {
    let needle = format!("{}=", key);
    for part in attrs.split(';') {
        if let Some(val) = part.trim().strip_prefix(&needle) {
            return Some(val.to_string());
        }
    }
    None
}
