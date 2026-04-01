use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;

// ── Public types ───────────────────────────────────────────────────────────────

/// The annotation returned for a genomic position.
pub struct GeneAnnotation {
    pub gene: String,
    /// Sub-gene feature type: "exon", "CDS", "5UTR", "3UTR", or "UTR".
    /// None when only a gene-level interval covers the position.
    pub feature_type: Option<String>,
    /// 1-based exon number within the transcript (from GTF/GFF3 `exon_number` attribute,
    /// or ordinal position within genePred exon arrays). None for CDS/UTR/gene features.
    pub exon_number: Option<i32>,
}

// ── Internal types ─────────────────────────────────────────────────────────────

/// Priority of overlapping features: higher wins when multiple intervals cover a position.
#[derive(Clone, Copy, PartialEq, Eq)]
enum FeatureKind {
    Gene,
    Utr,   // UTR where 5'/3' not determinable (e.g. GTF plain "UTR")
    Utr5,
    Utr3,
    Exon,
    Cds,
}

impl FeatureKind {
    fn priority(self) -> u8 {
        match self {
            Self::Cds             => 4,
            Self::Exon            => 3,
            Self::Utr5 | Self::Utr3 | Self::Utr => 2,
            Self::Gene            => 0,
        }
    }

    fn as_str(self) -> Option<&'static str> {
        match self {
            Self::Cds  => Some("CDS"),
            Self::Exon => Some("exon"),
            Self::Utr5 => Some("5UTR"),
            Self::Utr3 => Some("3UTR"),
            Self::Utr  => Some("UTR"),
            Self::Gene => None,
        }
    }
}

struct GeneRecord {
    start:       u32,
    end:         u32,
    gene:        String,
    kind:        FeatureKind,
    exon_number: Option<i32>,
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

    /// Returns the highest-priority interval covering `pos` (0-based), or None.
    fn get(&self, pos: u32) -> Option<&GeneRecord> {
        let idx = self.intervals.partition_point(|r| r.start <= pos);
        if idx == 0 {
            return None;
        }
        if self.prefix_max_end[idx - 1] <= pos {
            return None;
        }
        let mut best: Option<&GeneRecord> = None;
        for i in (0..idx).rev() {
            if self.intervals[i].end > pos {
                let candidate = &self.intervals[i];
                let better = match best {
                    None    => true,
                    Some(b) => candidate.kind.priority() > b.kind.priority(),
                };
                if better {
                    best = Some(candidate);
                }
            }
            if self.prefix_max_end[i] <= pos {
                break;
            }
        }
        best
    }
}

// ── Public API ─────────────────────────────────────────────────────────────────

/// Pre-loaded gene interval lookup built from a gene annotation file.
///
/// **Format auto-detection** (by file extension, stripping `.gz` first):
/// - `.txt` / `.txt.gz`  → UCSC genePred (e.g. ncbiRefSeq.txt.gz from UCSC).
///                          One row per transcript; `name2` column used as gene symbol.
///                          Exon intervals are derived from the exonStarts/exonEnds columns.
/// - `.gff3` / `.gff`   → GFF3 (gene, exon, CDS, UTR features)
/// - `.gtf` / other     → GTF  (gene, exon, CDS, UTR features)
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
            let line = line.with_context(|| format!("failed to read line {}", lineno + 1))?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            match fmt {
                Format::GenePred => parse_genepred_line(&line, &mut raw),
                Format::Gtf      => parse_gtf_line(&line, &mut raw),
                Format::Gff3     => parse_gff3_line(&line, &mut raw),
            }
        }

        let by_chrom = raw
            .into_iter()
            .map(|(chrom, records)| (chrom, ChromGenes::from_unsorted(records)))
            .collect();

        Ok(Self { by_chrom })
    }

    /// Returns the best annotation for the given 0-based locus, or None.
    ///
    /// "Best" means the highest-priority feature type (CDS > exon > UTR > gene).
    ///
    /// Tries the chromosome name as-is first, then bridges the common
    /// `chr1` ↔ `1` mismatch between UCSC and NCBI/Broad references.
    pub fn get(&self, chrom: &str, pos: i64) -> Option<GeneAnnotation> {
        let p = pos as u32;
        let record = self.lookup(chrom, p)?;
        Some(GeneAnnotation {
            gene:         record.gene.clone(),
            feature_type: record.kind.as_str().map(str::to_string),
            exon_number:  record.exon_number,
        })
    }

    fn lookup(&self, chrom: &str, pos: u32) -> Option<&GeneRecord> {
        if let Some(cg) = self.by_chrom.get(chrom) {
            if let Some(r) = cg.get(pos) {
                return Some(r);
            }
        }
        // UCSC annotation (chr1) queried with NCBI-style chrom (1) → add prefix
        if let Some(cg) = self.by_chrom.get(&format!("chr{chrom}")) {
            if let Some(r) = cg.get(pos) {
                return Some(r);
            }
        }
        // NCBI annotation (1) queried with UCSC-style chrom (chr1) → strip prefix
        if let Some(stripped) = chrom.strip_prefix("chr") {
            if let Some(cg) = self.by_chrom.get(stripped) {
                return cg.get(pos);
            }
        }
        None
    }

    pub fn n_intervals(&self) -> usize {
        self.by_chrom.values().map(|c| c.intervals.len()).sum()
    }
}

// ── Format detection ───────────────────────────────────────────────────────────

enum Format { GenePred, Gtf, Gff3 }

fn detect_format(path: &Path) -> Format {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    let inner = name.strip_suffix(".gz").unwrap_or(name);
    match Path::new(inner).extension().and_then(|e| e.to_str()).unwrap_or("") {
        "txt"          => Format::GenePred,
        "gff3" | "gff" => Format::Gff3,
        _              => Format::Gtf,
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

// ── GTF parsing ───────────────────────────────────────────────────────────────

fn parse_gtf_line(line: &str, raw: &mut HashMap<String, Vec<GeneRecord>>) {
    let fields: Vec<&str> = line.splitn(9, '\t').collect();
    if fields.len() < 9 {
        return;
    }
    let kind = match fields[2] {
        "gene"            => FeatureKind::Gene,
        "exon"            => FeatureKind::Exon,
        "CDS"             => FeatureKind::Cds,
        "five_prime_UTR"  | "five_prime_utr"  => FeatureKind::Utr5,
        "three_prime_UTR" | "three_prime_utr" => FeatureKind::Utr3,
        "UTR"             => FeatureKind::Utr,
        _                 => return,
    };

    // GTF is 1-based inclusive → convert to 0-based half-open.
    let start: u32 = match fields[3].parse::<u32>() {
        Ok(v) => v.saturating_sub(1),
        Err(_) => return,
    };
    let end: u32 = match fields[4].parse() {
        Ok(v) => v,
        Err(_) => return,
    };

    let attrs = fields[8];
    let gene = match parse_gtf_gene_name(attrs) {
        Some(g) => g,
        None    => return,
    };
    let exon_number = if matches!(kind, FeatureKind::Exon) {
        extract_gtf_attr(attrs, "exon_number")
            .and_then(|v| v.parse::<i32>().ok())
    } else {
        None
    };

    raw.entry(fields[0].to_string())
        .or_default()
        .push(GeneRecord { start, end, gene, kind, exon_number });
}

fn parse_gtf_gene_name(attrs: &str) -> Option<String> {
    extract_gtf_attr(attrs, "gene_name")
        .or_else(|| extract_gtf_attr(attrs, "gene_id"))
}

fn extract_gtf_attr(attrs: &str, key: &str) -> Option<String> {
    // Matches both quoted (`gene_name "X"`) and unquoted (`exon_number 1`) values.
    let needle_quoted = format!("{} \"", key);
    if let Some(start) = attrs.find(&needle_quoted) {
        let rest = &attrs[start + needle_quoted.len()..];
        if let Some(end) = rest.find('"') {
            return Some(rest[..end].to_string());
        }
    }
    // Unquoted: `key value;` — used for numeric attributes like exon_number in some GTFs.
    let needle_unquoted = format!("{} ", key);
    if let Some(start) = attrs.find(&needle_unquoted) {
        let rest = &attrs[start + needle_unquoted.len()..];
        let end = rest.find([';', ' ', '\t']).unwrap_or(rest.len());
        let val = rest[..end].trim();
        if !val.is_empty() {
            return Some(val.to_string());
        }
    }
    None
}

// ── GFF3 parsing ──────────────────────────────────────────────────────────────

fn parse_gff3_line(line: &str, raw: &mut HashMap<String, Vec<GeneRecord>>) {
    let fields: Vec<&str> = line.splitn(9, '\t').collect();
    if fields.len() < 9 {
        return;
    }
    let kind = match fields[2] {
        "gene"            => FeatureKind::Gene,
        "exon"            => FeatureKind::Exon,
        "CDS"             => FeatureKind::Cds,
        "five_prime_UTR"  => FeatureKind::Utr5,
        "three_prime_UTR" => FeatureKind::Utr3,
        _                 => return,
    };

    let start: u32 = match fields[3].parse::<u32>() {
        Ok(v) => v.saturating_sub(1),
        Err(_) => return,
    };
    let end: u32 = match fields[4].parse() {
        Ok(v) => v,
        Err(_) => return,
    };

    let attrs = fields[8];
    let gene = match parse_gff3_gene_name(attrs) {
        Some(g) => g,
        None    => return,
    };
    let exon_number = if matches!(kind, FeatureKind::Exon) {
        extract_gff3_attr(attrs, "exon_number")
            .and_then(|v| v.parse::<i32>().ok())
    } else {
        None
    };

    raw.entry(fields[0].to_string())
        .or_default()
        .push(GeneRecord { start, end, gene, kind, exon_number });
}

fn parse_gff3_gene_name(attrs: &str) -> Option<String> {
    extract_gff3_attr(attrs, "Name")
        .filter(|v| !v.starts_with("transcript:") && !v.starts_with("exon:"))
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

// ── genePred parsing ──────────────────────────────────────────────────────────

fn parse_genepred_line(line: &str, raw: &mut HashMap<String, Vec<GeneRecord>>) {
    // UCSC genePred with leading bin column:
    //  0:bin  1:name  2:chrom  3:strand  4:txStart  5:txEnd
    //  6:cdsStart  7:cdsEnd  8:exonCount  9:exonStarts  10:exonEnds
    //  11:score  12:name2  ...
    //
    // txStart/txEnd and exonStarts/exonEnds are 0-based half-open.
    let fields: Vec<&str> = line.splitn(14, '\t').collect();
    if fields.len() < 13 {
        return;
    }
    let chrom = fields[2];
    let tx_start: u32 = match fields[4].parse() {
        Ok(v) => v,
        Err(_) => return,
    };
    let tx_end: u32 = match fields[5].parse() {
        Ok(v) => v,
        Err(_) => return,
    };
    let gene = fields[12].trim().to_string();
    if gene.is_empty() {
        return;
    }

    // Emit the transcript-level gene interval.
    raw.entry(chrom.to_string())
        .or_default()
        .push(GeneRecord {
            start: tx_start, end: tx_end,
            gene: gene.clone(),
            kind: FeatureKind::Gene,
            exon_number: None,
        });

    // Emit per-exon intervals numbered 1..=exonCount.
    let exon_starts: Vec<u32> = fields[9]
        .split(',')
        .filter(|s| !s.is_empty())
        .filter_map(|s| s.parse().ok())
        .collect();
    let exon_ends: Vec<u32> = fields[10]
        .split(',')
        .filter(|s| !s.is_empty())
        .filter_map(|s| s.parse().ok())
        .collect();

    let chrom_vec = raw.entry(chrom.to_string()).or_default();
    for (i, (start, end)) in exon_starts.iter().zip(exon_ends.iter()).enumerate() {
        chrom_vec.push(GeneRecord {
            start: *start,
            end: *end,
            gene: gene.clone(),
            kind: FeatureKind::Exon,
            exon_number: Some((i + 1) as i32),
        });
    }
}

// ── Tests ──────────────────────────────────────────────────────────────────────

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

    fn gene(ga: &GeneAnnotations, chrom: &str, pos: i64) -> Option<String> {
        ga.get(chrom, pos).map(|a| a.gene)
    }

    fn feature(ga: &GeneAnnotations, chrom: &str, pos: i64) -> Option<String> {
        ga.get(chrom, pos).and_then(|a| a.feature_type)
    }

    fn exon_num(ga: &GeneAnnotations, chrom: &str, pos: i64) -> Option<i32> {
        ga.get(chrom, pos).and_then(|a| a.exon_number)
    }

    // ── GTF gene-only ─────────────────────────────────────────────────────────

    #[test]
    fn gtf_basic_lookup() {
        let ga = load_gtf(
            "chr1\t.\tgene\t100\t500\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
        );
        assert_eq!(gene(&ga, "chr1", 99),  Some("GENE1".into()));
        assert_eq!(gene(&ga, "chr1", 250), Some("GENE1".into()));
        assert_eq!(gene(&ga, "chr1", 499), Some("GENE1".into()));
        assert_eq!(gene(&ga, "chr1", 500), None);
        assert_eq!(gene(&ga, "chr1", 98),  None);
        // Gene-level feature_type is None
        assert_eq!(feature(&ga, "chr1", 250), None);
    }

    #[test]
    fn gtf_gene_id_fallback() {
        let ga = load_gtf(
            "chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_id \"FALLBACK\";\n",
        );
        assert_eq!(gene(&ga, "chr1", 0), Some("FALLBACK".into()));
    }

    #[test]
    fn gtf_exon_feature_type_and_number() {
        let ga = load_gtf(concat!(
            "chr1\t.\tgene\t1\t1000\t.\t+\t.\tgene_name \"BRCA1\";\n",
            "chr1\t.\texon\t100\t300\t.\t+\t.\tgene_name \"BRCA1\"; exon_number \"2\";\n",
        ));
        // Inside exon → feature_type = "exon", exon_number = 2
        assert_eq!(gene(&ga, "chr1", 150),    Some("BRCA1".into()));
        assert_eq!(feature(&ga, "chr1", 150), Some("exon".into()));
        assert_eq!(exon_num(&ga, "chr1", 150), Some(2));
        // Inside gene but outside exon → feature_type = None
        assert_eq!(feature(&ga, "chr1", 50),  None);
        assert_eq!(exon_num(&ga, "chr1", 50), None);
    }

    #[test]
    fn gtf_cds_beats_exon() {
        // CDS nested inside exon — CDS should win (higher priority).
        let ga = load_gtf(concat!(
            "chr1\t.\tgene\t1\t1000\t.\t+\t.\tgene_name \"TP53\";\n",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tgene_name \"TP53\"; exon_number \"1\";\n",
            "chr1\t.\tCDS\t150\t400\t.\t+\t0\tgene_name \"TP53\";\n",
        ));
        assert_eq!(feature(&ga, "chr1", 200), Some("CDS".into()));
        // In exon but outside CDS
        assert_eq!(feature(&ga, "chr1", 110), Some("exon".into()));
    }

    #[test]
    fn gtf_utr_variants() {
        let ga = load_gtf(concat!(
            "chr1\t.\tgene\t1\t1000\t.\t+\t.\tgene_name \"MYC\";\n",
            "chr1\t.\tfive_prime_UTR\t1\t100\t.\t+\t.\tgene_name \"MYC\";\n",
            "chr1\t.\tthree_prime_UTR\t900\t1000\t.\t+\t.\tgene_name \"MYC\";\n",
        ));
        assert_eq!(feature(&ga, "chr1", 50),  Some("5UTR".into()));
        assert_eq!(feature(&ga, "chr1", 950), Some("3UTR".into()));
        assert_eq!(feature(&ga, "chr1", 500), None); // gene only
    }

    #[test]
    fn gtf_non_gene_features_skipped() {
        let ga = load_gtf(concat!(
            "chr1\t.\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"TX1\"; gene_name \"TX1\";\n",
            "chr1\t.\tgene\t300\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"G1\";\n",
        ));
        assert_eq!(gene(&ga, "chr1", 50),  None);
        assert_eq!(gene(&ga, "chr1", 400), Some("G1".into()));
    }

    #[test]
    fn gtf_multiple_genes() {
        let ga = load_gtf(concat!(
            "chr1\t.\tgene\t100\t300\t.\t+\t.\tgene_name \"GENE_A\";\n",
            "chr1\t.\tgene\t500\t700\t.\t-\t.\tgene_name \"GENE_B\";\n",
            "chr2\t.\tgene\t1\t1000\t.\t+\t.\tgene_name \"GENE_C\";\n",
        ));
        assert_eq!(gene(&ga, "chr1", 200), Some("GENE_A".into()));
        assert_eq!(gene(&ga, "chr1", 400), None);
        assert_eq!(gene(&ga, "chr1", 600), Some("GENE_B".into()));
        assert_eq!(gene(&ga, "chr2", 500), Some("GENE_C".into()));
        assert_eq!(gene(&ga, "chr3", 500), None);
    }

    #[test]
    fn gtf_comment_lines_skipped() {
        let ga = load_gtf(concat!(
            "# this is a comment\n",
            "##format: gtf\n",
            "chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_name \"G1\";\n",
        ));
        assert_eq!(gene(&ga, "chr1", 0), Some("G1".into()));
    }

    // ── GFF3 ──────────────────────────────────────────────────────────────────

    #[test]
    fn gff3_basic_lookup() {
        let ga = load_gff3(
            "##gff-version 3\nchr1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=BRCA1\n",
        );
        assert_eq!(gene(&ga, "chr1", 99),  Some("BRCA1".into()));
        assert_eq!(gene(&ga, "chr1", 499), Some("BRCA1".into()));
        assert_eq!(gene(&ga, "chr1", 500), None);
    }

    #[test]
    fn gff3_exon_feature_type_and_number() {
        let ga = load_gff3(concat!(
            "chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1;Name=EGFR\n",
            "chr1\t.\texon\t100\t400\t.\t+\t.\tgene_name=EGFR;exon_number=3\n",
        ));
        assert_eq!(feature(&ga, "chr1", 200), Some("exon".into()));
        assert_eq!(exon_num(&ga, "chr1", 200), Some(3));
    }

    #[test]
    fn gff3_gene_name_fallbacks() {
        let ga = load_gff3(
            "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;gene_id=FALLBACK_ID\n",
        );
        assert_eq!(gene(&ga, "chr1", 0), Some("FALLBACK_ID".into()));
    }

    #[test]
    fn gff3_utr_features() {
        let ga = load_gff3(concat!(
            "chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1;Name=MYC\n",
            "chr1\t.\tfive_prime_UTR\t1\t100\t.\t+\t.\tgene_name=MYC\n",
            "chr1\t.\tthree_prime_UTR\t900\t1000\t.\t+\t.\tgene_name=MYC\n",
        ));
        assert_eq!(feature(&ga, "chr1", 50),  Some("5UTR".into()));
        assert_eq!(feature(&ga, "chr1", 950), Some("3UTR".into()));
    }

    // ── genePred ──────────────────────────────────────────────────────────────

    #[test]
    fn genepred_basic_lookup() {
        // bin=0, name=NM_001, chrom=chr1, strand=+, txStart=99, txEnd=500
        // exonCount=1, exonStarts=99, exonEnds=500, name2=BRCA1
        let ga = load_genepred(
            "0\tNM_001\tchr1\t+\t99\t500\t99\t500\t1\t99,\t500,\t0\tBRCA1\tcmpl\tcmpl\t0\n",
        );
        assert_eq!(gene(&ga, "chr1", 99),  Some("BRCA1".into()));
        assert_eq!(gene(&ga, "chr1", 300), Some("BRCA1".into()));
        assert_eq!(gene(&ga, "chr1", 499), Some("BRCA1".into()));
        assert_eq!(gene(&ga, "chr1", 500), None);
        assert_eq!(gene(&ga, "chr1", 98),  None);
        // Single exon → exon_number = 1
        assert_eq!(feature(&ga, "chr1", 200),   Some("exon".into()));
        assert_eq!(exon_num(&ga, "chr1", 200),  Some(1));
    }

    #[test]
    fn genepred_multi_exon_numbering() {
        // Two exons: [100,200) and [300,400)
        let ga = load_genepred(
            "0\tNM_001\tchr1\t+\t100\t400\t100\t400\t2\t100,300,\t200,400,\t0\tMYGENE\tcmpl\tcmpl\t0,0,\n",
        );
        assert_eq!(exon_num(&ga, "chr1", 150), Some(1));
        assert_eq!(exon_num(&ga, "chr1", 350), Some(2));
        // Intronic gap → gene-level only
        assert_eq!(feature(&ga, "chr1", 250), None);
        assert_eq!(gene(&ga, "chr1", 250),    Some("MYGENE".into()));
    }

    #[test]
    fn genepred_multiple_transcripts_same_gene() {
        let ga = load_genepred(concat!(
            "0\tNM_001\tchr1\t+\t100\t300\t100\t300\t1\t100,\t300,\t0\tMYGENE\tcmpl\tcmpl\t0\n",
            "0\tNM_002\tchr1\t+\t200\t500\t200\t500\t1\t200,\t500,\t0\tMYGENE\tcmpl\tcmpl\t0\n",
        ));
        assert_eq!(gene(&ga, "chr1", 150), Some("MYGENE".into()));
        assert_eq!(gene(&ga, "chr1", 400), Some("MYGENE".into()));
    }

    #[test]
    fn genepred_skips_short_lines() {
        let ga = load_genepred("0\tNM_001\tchr1\t+\t99\t500\n");
        assert_eq!(gene(&ga, "chr1", 100), None);
    }

    // ── chr prefix bridging ───────────────────────────────────────────────────

    #[test]
    fn chr_prefix_ucsc_annotation_ncbi_query() {
        let ga = load_gtf(
            "chr1\t.\tgene\t100\t500\t.\t+\t.\tgene_name \"TP53\";\n",
        );
        assert_eq!(gene(&ga, "1",    200), Some("TP53".into()));
        assert_eq!(gene(&ga, "chr1", 200), Some("TP53".into()));
    }

    #[test]
    fn chr_prefix_ncbi_annotation_ucsc_query() {
        let ga = load_gtf(
            "1\t.\tgene\t100\t500\t.\t+\t.\tgene_name \"TP53\";\n",
        );
        assert_eq!(gene(&ga, "chr1", 200), Some("TP53".into()));
        assert_eq!(gene(&ga, "1",    200), Some("TP53".into()));
    }

    #[test]
    fn chr_prefix_no_false_match_across_chroms() {
        let ga = load_gtf(concat!(
            "chr1\t.\tgene\t100\t500\t.\t+\t.\tgene_name \"GENE1\";\n",
            "chr10\t.\tgene\t100\t500\t.\t+\t.\tgene_name \"GENE10\";\n",
        ));
        assert_eq!(gene(&ga, "1",  200), Some("GENE1".into()));
        assert_eq!(gene(&ga, "10", 200), Some("GENE10".into()));
    }

    // ── n_intervals ───────────────────────────────────────────────────────────

    #[test]
    fn n_intervals_counts_all_records() {
        // 1 gene + 1 exon = 2 intervals
        let ga = load_gtf(concat!(
            "chr1\t.\tgene\t1\t1000\t.\t+\t.\tgene_name \"A\";\n",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tgene_name \"A\"; exon_number \"1\";\n",
            "chr2\t.\tgene\t1\t100\t.\t+\t.\tgene_name \"B\";\n",
        ));
        assert_eq!(ga.n_intervals(), 3);
    }
}
