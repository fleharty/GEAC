use serde::{Deserialize, Serialize};

/// Type of variant at this locus
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum VariantType {
    Snv,
    Insertion,
    Deletion,
}

impl std::fmt::Display for VariantType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VariantType::Snv => write!(f, "SNV"),
            VariantType::Insertion => write!(f, "insertion"),
            VariantType::Deletion => write!(f, "deletion"),
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
    pub batch: Option<String>,
    /// Generic sample label 1 (user-defined; e.g. tissue type).
    pub label1: Option<String>,
    /// Generic sample label 2 (user-defined; e.g. library prep method).
    pub label2: Option<String>,
    /// Generic sample label 3 (user-defined; e.g. sequencer type).
    pub label3: Option<String>,
    /// SHA-256 of the input BAM/CRAM when collect was run with --input-checksum-sha256.
    pub input_checksum_sha256: Option<String>,

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

    // Trinucleotide context for SNVs (5'-base, ref-base, 3'-base from reference).
    // Null for indels or when the locus is at a chromosome boundary.
    pub trinuc_context: Option<String>,

    // gnomAD allele frequency annotation (None if --gnomad not provided or allele absent).
    /// Allele frequency from gnomAD INFO/AF field. Null when --gnomad is not supplied
    /// or when this exact allele is not present in the gnomAD VCF.
    pub gnomad_af: Option<f32>,
}

/// Cross-annotation of tumor alt-base loci against a paired normal sample.
///
/// Schema A1: one row per (tumor_sample_id, chrom, pos, tumor_alt_allele, normal_sample_id,
/// normal_alt_allele).  For every tumor locus an anchor row with `normal_alt_allele = None`
/// is always written so that `normal_depth` is available even when the normal is clean.
/// Additional rows are appended for each non-ref allele observed in the normal (SNV positions).
/// Positions are 0-based.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NormalEvidence {
    pub tumor_sample_id: String,
    pub chrom: String,
    /// 0-based position
    pub pos: i64,
    pub tumor_alt_allele: String,
    pub normal_sample_id: String,
    /// `None` = anchor row (no alt evidence in normal, or indel position).
    /// `Some(allele)` = a non-ref allele observed in the normal.
    pub normal_alt_allele: Option<String>,
    /// Total fragment depth in the normal at this position (reads passing filters).
    pub normal_depth: i32,
    /// Fragments in normal supporting `normal_alt_allele`; 0 for anchor rows.
    pub normal_alt_count: i32,
}

/// PoN cross-annotation for a tumor alt locus.
/// One row per (tumor_sample_id, chrom, pos, tumor_alt_allele).
/// Positions are 0-based.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PonEvidence {
    pub tumor_sample_id: String,
    pub chrom: String,
    /// 0-based position
    pub pos: i64,
    pub tumor_alt_allele: String,
    /// Number of PoN samples carrying this exact alt allele at this locus
    pub n_pon_samples: i64,
    /// Total samples in the PoN database (denominator for fraction)
    pub pon_total_samples: i64,
    /// Max(alt_count/total_depth) across PoN samples; None if locus absent from PoN
    pub max_pon_vaf: Option<f64>,
    /// Mean(alt_count/total_depth) across PoN samples; None if locus absent from PoN
    pub mean_pon_vaf: Option<f64>,
}

/// Per-position coverage record produced by `geac coverage`.
/// One row per position (or per bin when --bin-size > 1). Positions are 0-based.
#[derive(Debug, Clone)]
pub struct CoverageRecord {
    pub sample_id: String,
    pub chrom: String,
    /// 0-based start position
    pub pos: i64,
    /// Exclusive end: pos+1 normally
    pub end: i64,
    pub bin_n: i32,

    // ── Fragment depth ─────────────────────────────────────────────────────────
    // Unique fragments (duplicates excluded, overlapping pairs collapsed to 1)
    // passing --min-map-qual.
    pub total_depth: i32,
    pub min_depth: i32,
    pub max_depth: i32,
    pub fwd_depth: i32,
    pub rev_depth: i32,

    // ── Duplicate metrics ──────────────────────────────────────────────────────
    /// All reads at this position (excluding secondary/supplementary), including dups
    pub raw_read_depth: i32,
    /// Fraction of raw reads flagged BAM_FDUP (0x400)
    pub frac_dup: f32,

    // ── Overlap metrics ────────────────────────────────────────────────────────
    /// Fragment pairs where both reads cover this position
    pub overlap_depth: i32,
    /// overlap_depth / total_depth (0.0 when total_depth = 0)
    pub frac_overlap: f32,

    // ── Mappability signals (non-dup reads, before --min-map-qual filter) ──────
    pub mean_mapq: f32,
    pub frac_mapq0: f32,
    pub frac_low_mapq: f32,

    // ── Base quality signals (non-dup, passing --min-map-qual) ─────────────────
    pub mean_base_qual: f32,
    /// Lowest base quality observed (stored as i32 for Arrow Int32)
    pub min_base_qual_obs: i32,
    /// Highest base quality observed
    pub max_base_qual_obs: i32,
    /// Fraction of bases below --min-base-qual (default 20)
    pub frac_low_bq: f32,

    // ── Insert size (properly paired, non-dup, passing --min-map-qual, R1) ──────
    pub mean_insert_size: f32,
    pub min_insert_size: i32,
    pub max_insert_size: i32,
    pub n_insert_size_obs: i32,

    // ── GC content ─────────────────────────────────────────────────────────────
    /// Fraction of G+C bases in --gc-window bp window around this position
    pub gc_content: f32,

    // ── Annotations ────────────────────────────────────────────────────────────
    pub on_target: Option<bool>,
    pub gene: Option<String>,

    // ── Provenance ─────────────────────────────────────────────────────────────
    pub read_type: ReadType,
    pub pipeline: Pipeline,
    pub batch: Option<String>,
}

/// One record per alt-supporting read at a locus.
/// Linked to AltBase by (sample_id, chrom, pos, alt_allele).
/// Positions are 0-based.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AltRead {
    pub sample_id: String,
    pub chrom: String,
    /// 0-based position
    pub pos: i64,
    pub alt_allele: String,
    /// 1-based sequencing cycle (= query position + 1)
    pub cycle: i32,
    pub read_length: i32,
    /// true if this read is R1 (BAM flag 0x40), false if R2 or unpaired
    pub is_read1: bool,
    /// fgbio aD tag: AB (top-strand) raw read count; None if tag absent
    pub ab_count: Option<i32>,
    /// fgbio bD tag: BA (bottom-strand) raw read count; None if tag absent
    pub ba_count: Option<i32>,
    /// fgbio cD tag: total raw read count (aD + bD for duplex; sole count for simplex).
    /// This is the primary family size field — present for both simplex and duplex data.
    pub family_size: Option<i32>,
    pub base_qual: i32,
    pub map_qual: i32,
    /// SAM TLEN field (template length / insert size). None when TLEN is 0
    /// (unpaired reads or reads where the mate is unmapped).
    pub insert_size: Option<i32>,
    /// SHA-256 of the input BAM/CRAM when collect was run with --input-checksum-sha256.
    pub input_checksum_sha256: Option<String>,
}

/// Per-interval coverage summary produced by `geac coverage --intervals-output`.
/// One row per target interval. Coordinates are 0-based half-open [start, end).
#[derive(Debug, Clone)]
pub struct IntervalRecord {
    pub sample_id:    String,
    pub chrom:        String,
    /// 0-based interval start (from the targets file)
    pub start:        i64,
    /// 0-based interval end (exclusive)
    pub end:          i64,
    /// Name field from BED col 4 or Picard interval name col
    pub interval_name: Option<String>,
    pub gene:         Option<String>,

    // ── Depth summary ─────────────────────────────────────────────────────────
    /// Total number of positions in the interval (= end - start)
    pub n_bases:      i32,
    pub mean_depth:   f32,
    pub median_depth: f32,
    pub min_depth:    i32,
    pub max_depth:    i32,
    /// Fraction of positions with total_depth >= 1
    pub frac_at_1x:   f32,
    pub frac_at_10x:  f32,
    pub frac_at_20x:  f32,
    pub frac_at_30x:  f32,
    pub frac_at_50x:  f32,
    pub frac_at_100x: f32,

    // ── Aggregated QC signals (means across all positions in the interval) ────
    pub mean_gc_content:   f32,
    pub mean_mapq:         f32,
    pub mean_frac_mapq0:   f32,
    pub mean_frac_dup:     f32,
    pub mean_frac_overlap: f32,
    pub mean_base_qual:    f32,
    pub mean_insert_size:  f32,

    // ── Provenance ────────────────────────────────────────────────────────────
    pub read_type: ReadType,
    pub pipeline:  Pipeline,
    pub batch:     Option<String>,
}

