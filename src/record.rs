use serde::{Deserialize, Serialize};

/// Type of variant at this locus
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum VariantType {
    Snv,
    Insertion,
    Deletion,
    Mnv,
}

impl std::fmt::Display for VariantType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VariantType::Snv => write!(f, "SNV"),
            VariantType::Insertion => write!(f, "insertion"),
            VariantType::Deletion => write!(f, "deletion"),
            VariantType::Mnv => write!(f, "MNV"),
        }
    }
}

/// Whether reads are raw, simplex consensus, or duplex consensus
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ReadType {
    Raw,
    Simplex,
    Duplex,
}

impl std::fmt::Display for ReadType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadType::Raw => write!(f, "raw"),
            ReadType::Simplex => write!(f, "simplex"),
            ReadType::Duplex => write!(f, "duplex"),
        }
    }
}

/// Which pipeline produced the BAM/CRAM
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Pipeline {
    Fgbio,
    Dragen,
    Raw,
}

impl std::fmt::Display for Pipeline {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Pipeline::Fgbio => write!(f, "fgbio"),
            Pipeline::Dragen => write!(f, "dragen"),
            Pipeline::Raw => write!(f, "raw"),
        }
    }
}

/// One record per alt allele observed at a locus in a sample.
/// Positions are 0-based.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AltBase {
    pub sample_id: String,
    pub chrom: String,
    /// 0-based position
    pub pos: i64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_type: VariantType,

    // Depth & alt support
    pub total_depth: i32,
    pub alt_count: i32,
    pub ref_count: i32,

    // Strand breakdown
    pub fwd_depth: i32,
    pub rev_depth: i32,
    pub fwd_alt_count: i32,
    pub rev_alt_count: i32,
    pub fwd_ref_count: i32,
    pub rev_ref_count: i32,

    // Fragment overlap
    /// Number of overlapping fragment pairs at this position
    pub overlap_depth: i32,
    /// Overlapping pairs where both reads support the alt
    pub overlap_alt_agree: i32,
    /// Overlapping pairs where reads disagree (one alt, one ref or other)
    pub overlap_alt_disagree: i32,
    /// Overlapping pairs where both reads support the ref
    pub overlap_ref_agree: i32,

    // Provenance
    pub read_type: ReadType,
    pub pipeline: Pipeline,

    // Variant calling annotation (None if no VCF was provided)
    pub variant_called: Option<bool>,
    /// None = not called; Some("PASS") = passing; Some("...") = filter reason(s)
    pub variant_filter: Option<String>,

    // Target annotation (None if no targets file was provided)
    /// true if the locus overlaps a target interval, false if not, null if no targets given
    pub on_target: Option<bool>,

    // Gene annotation (None if no annotation file was provided or locus is intergenic)
    pub gene: Option<String>,

    // Repetitiveness metrics (always populated; 0 = no repeat detected)
    /// Length of the longest homopolymer run overlapping this position (>= 1)
    pub homopolymer_len: i32,
    /// Period of the shortest tandem repeat unit whose tract includes this position
    /// (1 = homopolymer, 2 = dinucleotide, …). 0 if no repeat.
    pub str_period: i32,
    /// Total length (bp) of the STR tract. 0 if no repeat.
    pub str_len: i32,
}
