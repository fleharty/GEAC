use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;

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
        let idx = self.intervals.partition_point(|r| r.start <= pos);
        if idx == 0 {
            return None;
        }
        if self.prefix_max_end[idx - 1] <= pos {
            return None;
        }
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

/// Pre-loaded gene interval lookup built from a gene annotation file.
///
/// **Format auto-detection** (by file extension, stripping `.gz` first):
/// - `.txt` / `.txt.gz`  → UCSC genePred (e.g. ncbiRefSeq.txt.gz from UCSC).
///                          One row per transcript; `name2` column used as gene symbol.
///                          Coordinates are already 0-based half-open.
/// - `.gff3` / `.gff`   → GFF3 (gene-level features only)
/// - `.gtf` / other     → GTF  (gene-level features only)
///
/// Both `.gz` and uncompressed files are supported for all formats.
pub struct GeneAnnotations {
    by_chrom: HashMap<String, ChromGenes>,
}

impl GeneAnnotations {
    pub fn load(path: &Path) -> Result<Self> {
        let fmt = detect_format(path);
        let reader = open_reader(path)
            .with_context(|| format!("failed to open gene annotations: {}", path.display()))?;
        Self::load_from_reader(reader, fmt)
    }

    fn load_from_reader(reader: impl BufRead, fmt: Format) -> Result<Self> {
        let mut raw: HashMap<String, Vec<GeneRecord>> = HashMap::new();

        for (lineno, line) in reader.lines().enumerate() {
            let line = line.with_context(|| {
                format!("failed to read line {}", lineno + 1)
            })?;

            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            match fmt {
                Format::GenePred => {
                    // UCSC genePred with leading bin column:
                    //  0:bin  1:name  2:chrom  3:strand  4:txStart  5:txEnd
                    //  6:cdsStart  7:cdsEnd  8:exonCount  9:exonStarts  10:exonEnds
                    //  11:score  12:name2  ...
                    //
                    // txStart/txEnd are 0-based half-open — no conversion needed.
                    let fields: Vec<&str> = line.splitn(14, '\t').collect();
                    if fields.len() < 13 {
                        continue;
                    }
                    let chrom = fields[2];
                    let start: u32 = match fields[4].parse() {
                        Ok(v) => v,
                        Err(_) => continue,
                    };
                    let end: u32 = match fields[5].parse() {
                        Ok(v) => v,
                        Err(_) => continue,
                    };
                    let gene = fields[12].trim().to_string();
                    if gene.is_empty() {
                        continue;
                    }
                    raw.entry(chrom.to_string())
                        .or_default()
                        .push(GeneRecord { start, end, gene });
                }

                Format::Gtf | Format::Gff3 => {
                    let fields: Vec<&str> = line.splitn(9, '\t').collect();
                    if fields.len() < 9 {
                        continue;
                    }
                    if fields[2] != "gene" {
                        continue;
                    }
                    // GTF/GFF3 are 1-based inclusive → convert to 0-based half-open.
                    let start: u32 = match fields[3].parse::<u32>() {
                        Ok(v) => v.saturating_sub(1),
                        Err(_) => continue,
                    };
                    let end: u32 = match fields[4].parse() {
                        Ok(v) => v,
                        Err(_) => continue,
                    };
                    let gene = match fmt {
                        Format::Gtf  => parse_gtf_gene_name(fields[8]),
                        Format::Gff3 => parse_gff3_gene_name(fields[8]),
                        _ => unreachable!(),
                    };
                    if let Some(gene) = gene {
                        raw.entry(fields[0].to_string())
                            .or_default()
                            .push(GeneRecord { start, end, gene });
                    }
                }
            }
        }

        let by_chrom = raw
            .into_iter()
            .map(|(chrom, records)| (chrom, ChromGenes::from_unsorted(records)))
            .collect();

        Ok(Self { by_chrom })
    }

    /// Returns the gene name for the given 0-based locus, or None.
    ///
    /// Tries the chromosome name as-is first, then bridges the common
    /// `chr1` ↔ `1` mismatch between UCSC and NCBI/Broad references.
    pub fn get(&self, chrom: &str, pos: i64) -> Option<&str> {
        let p = pos as u32;
        if let Some(cg) = self.by_chrom.get(chrom) {
            return cg.get(p);
        }
        // UCSC annotation (chr1) queried with NCBI-style chrom (1) → add prefix
        if let Some(cg) = self.by_chrom.get(&format!("chr{chrom}")) {
            return cg.get(p);
        }
        // NCBI annotation (1) queried with UCSC-style chrom (chr1) → strip prefix
        if let Some(stripped) = chrom.strip_prefix("chr") {
            if let Some(cg) = self.by_chrom.get(stripped) {
                return cg.get(p);
            }
        }
        None
    }

    pub fn n_genes(&self) -> usize {
        self.by_chrom.values().map(|c| c.intervals.len()).sum()
    }
}

// ── Format detection ───────────────────────────────────────────────────────────

enum Format { GenePred, Gtf, Gff3 }

fn detect_format(path: &Path) -> Format {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    // Strip a trailing .gz to reveal the real extension.
    let inner = name.strip_suffix(".gz").unwrap_or(name);
    match Path::new(inner).extension().and_then(|e| e.to_str()).unwrap_or("") {
        "txt"        => Format::GenePred,
        "gff3" | "gff" => Format::Gff3,
        _            => Format::Gtf,   // .gtf and anything else
    }
}

// ── File opening (transparent gzip) ───────────────────────────────────────────

fn open_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("cannot open {}", path.display()))?;
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

// ── GTF attribute parsing ──────────────────────────────────────────────────────

fn parse_gtf_gene_name(attrs: &str) -> Option<String> {
    extract_gtf_attr(attrs, "gene_name")
        .or_else(|| extract_gtf_attr(attrs, "gene_id"))
}

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

fn extract_gff3_attr(attrs: &str, key: &str) -> Option<String> {
    let needle = format!("{}=", key);
    for part in attrs.split(';') {
        if let Some(val) = part.trim().strip_prefix(&needle) {
            return Some(val.to_string());
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn load_gtf(s: &str) -> GeneAnnotations {
        GeneAnnotations::load_from_reader(Cursor::new(s), Format::Gtf).expect("gtf parse failed")
    }

    fn load_gff3(s: &str) -> GeneAnnotations {
        GeneAnnotations::load_from_reader(Cursor::new(s), Format::Gff3).expect("gff3 parse failed")
    }

    fn load_genepred(s: &str) -> GeneAnnotations {
        GeneAnnotations::load_from_reader(Cursor::new(s), Format::GenePred)
            .expect("genepred parse failed")
    }

    // ── GTF ───────────────────────────────────────────────────────────────────

    #[test]
    fn gtf_basic_lookup() {
        // GTF is 1-based inclusive: [100, 500] → 0-based half-open [99, 500)
        let ga = load_gtf(
            "chr1\t.\tgene\t100\t500\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
        );
        assert_eq!(ga.get("chr1", 99),  Some("GENE1")); // first base (0-based)
        assert_eq!(ga.get("chr1", 250), Some("GENE1")); // middle
        assert_eq!(ga.get("chr1", 499), Some("GENE1")); // last base
        assert_eq!(ga.get("chr1", 500), None);           // just outside
        assert_eq!(ga.get("chr1", 98),  None);           // just before
    }

    #[test]
    fn gtf_gene_id_fallback() {
        // No gene_name attribute → falls back to gene_id
        let ga = load_gtf(
            "chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_id \"FALLBACK\";\n",
        );
        assert_eq!(ga.get("chr1", 0), Some("FALLBACK"));
    }

    #[test]
    fn gtf_non_gene_features_skipped() {
        let ga = load_gtf(concat!(
            "chr1\t.\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"TX1\"; gene_name \"TX1\";\n",
            "chr1\t.\texon\t1\t200\t.\t+\t.\tgene_id \"EX1\"; gene_name \"EX1\";\n",
            "chr1\t.\tgene\t300\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"G1\";\n",
        ));
        // Transcript and exon ranges should not contribute
        assert_eq!(ga.get("chr1", 50),  None);
        assert_eq!(ga.get("chr1", 400), Some("G1"));
    }

    #[test]
    fn gtf_multiple_genes() {
        let ga = load_gtf(concat!(
            "chr1\t.\tgene\t100\t300\t.\t+\t.\tgene_name \"GENE_A\";\n",
            "chr1\t.\tgene\t500\t700\t.\t-\t.\tgene_name \"GENE_B\";\n",
            "chr2\t.\tgene\t1\t1000\t.\t+\t.\tgene_name \"GENE_C\";\n",
        ));
        assert_eq!(ga.get("chr1", 200), Some("GENE_A"));
        assert_eq!(ga.get("chr1", 400), None);
        assert_eq!(ga.get("chr1", 600), Some("GENE_B"));
        assert_eq!(ga.get("chr2", 500), Some("GENE_C"));
        assert_eq!(ga.get("chr3", 500), None);
    }

    #[test]
    fn gtf_comment_lines_skipped() {
        let ga = load_gtf(concat!(
            "# this is a comment\n",
            "##format: gtf\n",
            "chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_name \"G1\";\n",
        ));
        assert_eq!(ga.get("chr1", 0), Some("G1"));
    }

    // ── GFF3 ──────────────────────────────────────────────────────────────────

    #[test]
    fn gff3_basic_lookup() {
        let ga = load_gff3(
            "##gff-version 3\nchr1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=BRCA1\n",
        );
        assert_eq!(ga.get("chr1", 99),  Some("BRCA1"));
        assert_eq!(ga.get("chr1", 499), Some("BRCA1"));
        assert_eq!(ga.get("chr1", 500), None);
    }

    #[test]
    fn gff3_gene_name_fallbacks() {
        // Name= is missing → try gene_name= → try gene_id=
        let ga = load_gff3(
            "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;gene_id=FALLBACK_ID\n",
        );
        assert_eq!(ga.get("chr1", 0), Some("FALLBACK_ID"));
    }

    // ── genePred (UCSC ncbiRefSeq.txt) ────────────────────────────────────────

    #[test]
    fn genepred_basic_lookup() {
        // bin=0, name=NM_001, chrom=chr1, strand=+, txStart=99, txEnd=500, ..., name2=BRCA1
        // txStart/txEnd are 0-based half-open → [99, 500)
        let ga = load_genepred(
            "0\tNM_001\tchr1\t+\t99\t500\t99\t500\t1\t99,\t500,\t0\tBRCA1\tcmpl\tcmpl\t0\n",
        );
        assert_eq!(ga.get("chr1", 99),  Some("BRCA1"));
        assert_eq!(ga.get("chr1", 300), Some("BRCA1"));
        assert_eq!(ga.get("chr1", 499), Some("BRCA1"));
        assert_eq!(ga.get("chr1", 500), None);
        assert_eq!(ga.get("chr1", 98),  None);
    }

    #[test]
    fn genepred_multiple_transcripts_same_gene() {
        // Two isoforms of the same gene — both intervals should be searchable
        let ga = load_genepred(concat!(
            "0\tNM_001\tchr1\t+\t100\t300\t100\t300\t1\t100,\t300,\t0\tMYGENE\tcmpl\tcmpl\t0\n",
            "0\tNM_002\tchr1\t+\t200\t500\t200\t500\t1\t200,\t500,\t0\tMYGENE\tcmpl\tcmpl\t0\n",
        ));
        assert_eq!(ga.get("chr1", 150), Some("MYGENE"));
        assert_eq!(ga.get("chr1", 400), Some("MYGENE"));
    }

    #[test]
    fn genepred_skips_short_lines() {
        // Lines with fewer than 13 fields are silently ignored
        let ga = load_genepred("0\tNM_001\tchr1\t+\t99\t500\n");
        assert_eq!(ga.get("chr1", 100), None);
    }

    // ── chr prefix bridging ───────────────────────────────────────────────────

    #[test]
    fn chr_prefix_ucsc_annotation_ncbi_query() {
        // Annotation stored under "chr1", queried with "1" (NCBI/Broad style)
        let ga = load_gtf(
            "chr1\t.\tgene\t100\t500\t.\t+\t.\tgene_name \"TP53\";\n",
        );
        assert_eq!(ga.get("1",    200), Some("TP53")); // NCBI query → should match
        assert_eq!(ga.get("chr1", 200), Some("TP53")); // direct match still works
    }

    #[test]
    fn chr_prefix_ncbi_annotation_ucsc_query() {
        // Annotation stored under "1", queried with "chr1" (UCSC style)
        let ga = load_gtf(
            "1\t.\tgene\t100\t500\t.\t+\t.\tgene_name \"TP53\";\n",
        );
        assert_eq!(ga.get("chr1", 200), Some("TP53")); // UCSC query → should match
        assert_eq!(ga.get("1",    200), Some("TP53")); // direct match still works
    }

    #[test]
    fn chr_prefix_no_false_match_across_chroms() {
        // "chr10" should not accidentally match "chr1" via prefix stripping
        let ga = load_gtf(concat!(
            "chr1\t.\tgene\t100\t500\t.\t+\t.\tgene_name \"GENE1\";\n",
            "chr10\t.\tgene\t100\t500\t.\t+\t.\tgene_name \"GENE10\";\n",
        ));
        assert_eq!(ga.get("1",  200), Some("GENE1"));
        assert_eq!(ga.get("10", 200), Some("GENE10"));
    }

    // ── n_genes ───────────────────────────────────────────────────────────────

    #[test]
    fn n_genes_counts_all_intervals() {
        let ga = load_gtf(concat!(
            "chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_name \"A\";\n",
            "chr1\t.\tgene\t200\t300\t.\t+\t.\tgene_name \"B\";\n",
            "chr2\t.\tgene\t1\t100\t.\t+\t.\tgene_name \"C\";\n",
        ));
        assert_eq!(ga.n_genes(), 3);
    }
}
