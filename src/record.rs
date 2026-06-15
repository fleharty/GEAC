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

/// How to extract family-size counts from a read's aux tags.
///
/// The `pipeline` label is now free-text (any string the user supplies), so the
/// behavioral part — which aux tags encode family size — lives here, resolved
/// from either a built-in preset (`fgbio`/`dragen`) or an explicit
/// `--family-size-tags` spec. A tag that is absent from a read yields `None`,
/// so pipelines that do not emit family size simply produce null counts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FamilySizeScheme {
    /// AB / top-strand raw read-count tag (fgbio `aD`, DRAGEN `XV`).
    pub ab_tag: Option<[u8; 2]>,
    /// BA / bottom-strand raw read-count tag (fgbio `bD`; DRAGEN has none).
    pub ba_tag: Option<[u8; 2]>,
    /// Total family-size tag (fgbio `cD`, DRAGEN `XW`).
    pub total_tag: Option<[u8; 2]>,
    /// What to do when `total_tag` is absent or non-informative.
    pub total_fallback: TotalFallback,
}

/// Fallback rule for the total family size when `total_tag` does not yield a value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TotalFallback {
    /// No fallback: total is exactly what `total_tag` reads (or `None`).
    None,
    /// fgbio: when `total_tag` (`cD`) is absent, use `ab + ba` (`aD + bD`).
    SumAbBa,
    /// DRAGEN: use `total_tag` (`XW`) when > 0, otherwise `ab_tag` (`XV`).
    PositiveOrAb,
}

impl FamilySizeScheme {
    /// fgbio consensus tags: `aD` / `bD` / `cD`, with `cD` falling back to `aD + bD`.
    pub fn fgbio() -> Self {
        Self {
            ab_tag: Some(*b"aD"),
            ba_tag: Some(*b"bD"),
            total_tag: Some(*b"cD"),
            total_fallback: TotalFallback::SumAbBa,
        }
    }

    /// DRAGEN family tags: `XV` (fragments) and `XW`, with `XW` preferred when > 0.
    pub fn dragen() -> Self {
        Self {
            ab_tag: Some(*b"XV"),
            ba_tag: None,
            total_tag: Some(*b"XW"),
            total_fallback: TotalFallback::PositiveOrAb,
        }
    }

    /// No family-size tags (raw reads, or a pipeline that does not emit them).
    pub fn none() -> Self {
        Self {
            ab_tag: None,
            ba_tag: None,
            total_tag: None,
            total_fallback: TotalFallback::None,
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
    pub pipeline: Option<String>,
    /// Biological subject identifier — groups samples (e.g. multiple timepoints
    /// or tissue types) that come from the same person/animal/study subject.
    /// Used by `geac compare` to scope longitudinal and replicate analyses.
    pub subject_id: Option<String>,
    /// Sample substrate type. Standard values: cfDNA, urine_cfDNA, tumor_tissue,
    /// normal_tissue, buffy_coat, lymphocytes. Free-text; unrecognized values
    /// are preserved but ignored by mode-aware downstream tooling.
    pub sample_type: Option<String>,
    pub batch: Option<String>,
    /// Generic sample label 1 (user-defined; e.g. tissue type).
    pub label1: Option<String>,
    /// Generic sample label 2 (user-defined; e.g. library prep method).
    pub label2: Option<String>,
    /// Generic sample label 3 (user-defined; e.g. sequencer type).
    pub label3: Option<String>,
    /// Timepoint label for longitudinal / serial-sample studies (user-defined; e.g. "T0", "week12").
    pub timepoint: Option<String>,
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

    // Per-locus N-context summary, aggregated from alt-supporting reads.
    // Only populated when collect was run with reads-output enabled; None otherwise.
    /// Number of alt-supporting reads contributing to the N-context summary.
    pub n_alt_reads_with_n_ctx: Option<i32>,
    /// Mean fraction of N bases in the window before the alt position, across alt reads.
    pub mean_frac_n_before: Option<f32>,
    /// Mean fraction of N bases in the window after the alt position, across alt reads.
    pub mean_frac_n_after: Option<f32>,
    /// mean_frac_n_after - mean_frac_n_before. Positive = trailing N bias.
    pub mean_delta_n_frac: Option<f32>,
    /// Fraction of alt reads where N context is asymmetric (high after, low/none before).
    pub frac_reads_asymmetric: Option<f32>,

    // Resource paths (stored for downstream use in the Explorer and IGV session generation).
    // When running on Terra, use --bam-uri / --variants-uri / --gnomad-uri to store the
    // original cloud URIs instead of the localized local paths.
    /// Path or URI of the input BAM/CRAM file used to generate this record.
    pub bam_path: Option<String>,
    /// Path or URI of the BAM/CRAM index file.
    pub bai_path: Option<String>,
    /// Path or URI of the variants file (VCF or TSV) used for annotation, if any.
    pub variants_path: Option<String>,
    /// Path or URI of the gnomAD VCF used for allele-frequency annotation, if any.
    pub gnomad_path: Option<String>,
    /// Path or URI of the gnomAD index (.tbi/.csi), if specified, for IGV sessions.
    pub gnomad_index_path: Option<String>,
    /// Path or URI of the target BED / interval list used for on-target annotation, if any.
    pub targets_path: Option<String>,
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

/// Sample-level target-depth and on-target burden metrics produced by `geac collect`
/// when `--targets` is provided.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleMetricsRecord {
    pub sample_id: String,
    pub subject_id: Option<String>,
    pub sample_type: Option<String>,
    pub batch: Option<String>,
    pub read_type: ReadType,
    pub pipeline: Option<String>,
    pub input_checksum_sha256: Option<String>,
    pub n_target_positions: i32,
    pub n_target_positions_covered: i32,
    pub mean_target_depth_covered: Option<f32>,
    pub mean_target_depth_all: Option<f32>,
    pub median_target_depth_covered: Option<f32>,
    pub median_target_depth_all: Option<f32>,
    pub pct_fragment_bases_on_target: Option<f32>,
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

    // ── Pre-computed annotation tracks ─────────────────────────────────────────
    /// One score per --track flag (order matches TrackSet::names()). Empty when no tracks given.
    pub track_values: Vec<Option<f32>>,

    // ── Annotations ────────────────────────────────────────────────────────────
    pub on_target: Option<bool>,
    pub gene: Option<String>,
    pub feature_type: Option<String>,
    pub exon_number: Option<i32>,

    // ── Provenance ─────────────────────────────────────────────────────────────
    pub read_type: ReadType,
    pub pipeline: Option<String>,
    pub subject_id: Option<String>,
    pub sample_type: Option<String>,
    pub batch: Option<String>,
    /// Generic sample label 1 (user-defined; e.g. tissue type).
    pub label1: Option<String>,
    /// Generic sample label 2 (user-defined; e.g. library prep method).
    pub label2: Option<String>,
    /// Generic sample label 3 (user-defined; e.g. sequencer type).
    pub label3: Option<String>,
    /// Timepoint label for longitudinal / serial-sample studies (user-defined; e.g. "T0", "week12").
    pub timepoint: Option<String>,
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
    pub read_type: ReadType,
    pub pipeline: Option<String>,
    pub subject_id: Option<String>,
    pub sample_type: Option<String>,
    pub batch: Option<String>,
    pub label1: Option<String>,
    pub label2: Option<String>,
    pub label3: Option<String>,
    pub timepoint: Option<String>,
    /// 1-based sequencing cycle (= query position + 1)
    pub cycle: i32,
    pub read_length: i32,
    /// true if this read is R1 (BAM flag 0x40), false if R2 or unpaired
    pub is_read1: bool,
    /// Pipeline-aware strand/family support count.
    /// fgbio: aD (AB/top-strand raw read count); DRAGEN: XV (family fragments).
    pub ab_count: Option<i32>,
    /// Pipeline-aware second-strand support count.
    /// fgbio: bD (BA/bottom-strand raw read count); DRAGEN: no equivalent tag exposed.
    pub ba_count: Option<i32>,
    /// Pipeline-aware total support count.
    /// fgbio: cD (with fallback to aD + bD when needed); DRAGEN: XW when present, otherwise XV.
    pub family_size: Option<i32>,
    pub base_qual: i32,
    pub map_qual: i32,
    /// SAM TLEN field (template length / insert size). None when TLEN is 0
    /// (unpaired reads or reads where the mate is unmapped).
    pub insert_size: Option<i32>,
    /// GC fraction of the fragment inferred from the reference sequence.
    /// Computed over the region [read_pos, read_pos + abs(TLEN)) on the reference.
    /// None when insert size is unavailable (unpaired / mate unmapped).
    pub frag_gc: Option<f32>,
    /// Number of stored read-sequence bases before the alt position.
    pub n_before_alt: i32,
    /// Number of stored read-sequence bases after the alt position.
    pub n_after_alt: i32,
    /// Count of `N` bases before the alt position.
    pub n_n_before_alt: i32,
    /// Count of `N` bases after the alt position.
    pub n_n_after_alt: i32,
    /// Contiguous run of `N` bases immediately before the alt.
    pub leading_n_run_len: i32,
    /// Contiguous run of `N` bases immediately after the alt.
    pub trailing_n_run_len: i32,
    /// SHA-256 of the input BAM/CRAM when collect was run with --input-checksum-sha256.
    pub input_checksum_sha256: Option<String>,
    /// FNV-1a 64-bit hash of the read's qname, cast to i64 for Parquet Int64 compatibility.
    /// Both reads of a pair share the same qname, so this is a stable per-fragment ID.
    /// Used for MNV detection: a cross-locus self-join on fragment_id finds reads that
    /// carry substitutions at two different positions.
    pub fragment_id: i64,
}

/// One record per fragment (read pair) produced by `geac fragments`.
/// Coordinates are 0-based. Only proper pairs with positive TLEN are emitted.
#[derive(Debug, Clone)]
pub struct FragmentRecord {
    pub sample_id: String,
    pub chrom: String,
    /// 0-based fragment start (leftmost position of the pair)
    pub frag_start: i64,
    /// 0-based fragment end (frag_start + abs(TLEN))
    pub frag_end: i64,
    /// Fragment midpoint: (frag_start + frag_end) / 2
    pub midpoint: i64,
    /// abs(TLEN) from the BAM record
    pub insert_size: i32,
    /// GC fraction of the fragment span computed from the reference sequence
    pub gc_content: Option<f32>,
    /// 4-mer reference sequence at the 5' fragment end
    pub end_motif_5p: Option<String>,
    /// 4-mer reference sequence at the 3' fragment end
    pub end_motif_3p: Option<String>,
    /// Mapping quality of R1
    pub map_qual: i32,
    pub read_type: ReadType,
    pub pipeline: Option<String>,
    pub subject_id: Option<String>,
    pub sample_type: Option<String>,
    pub batch: Option<String>,
    pub label1: Option<String>,
    pub label2: Option<String>,
    pub label3: Option<String>,
    pub timepoint: Option<String>,
}

/// Per-interval coverage summary produced by `geac coverage --intervals-output`.
/// One row per target interval. Coordinates are 0-based half-open [start, end).
#[derive(Debug, Clone)]
pub struct IntervalRecord {
    pub sample_id: String,
    pub chrom: String,
    /// 0-based interval start (from the targets file)
    pub start: i64,
    /// 0-based interval end (exclusive)
    pub end: i64,
    /// Name field from BED col 4 or Picard interval name col
    pub interval_name: Option<String>,
    pub gene: Option<String>,
    pub feature_type: Option<String>,
    pub exon_number: Option<i32>,

    // ── Depth summary ─────────────────────────────────────────────────────────
    /// Total number of positions in the interval (= end - start)
    pub n_bases: i32,
    pub mean_depth: f32,
    pub median_depth: f32,
    pub min_depth: i32,
    pub max_depth: i32,
    /// Fraction of positions with total_depth >= 1
    pub frac_at_1x: f32,
    pub frac_at_10x: f32,
    pub frac_at_20x: f32,
    pub frac_at_30x: f32,
    pub frac_at_50x: f32,
    pub frac_at_100x: f32,

    // ── Aggregated QC signals (means across all positions in the interval) ────
    pub mean_gc_content: f32,
    pub mean_mapq: f32,
    pub mean_frac_mapq0: f32,
    pub mean_frac_dup: f32,
    pub mean_frac_overlap: f32,
    pub mean_base_qual: f32,
    pub mean_insert_size: f32,

    // ── Provenance ────────────────────────────────────────────────────────────
    pub read_type: ReadType,
    pub pipeline: Option<String>,
    pub subject_id: Option<String>,
    pub sample_type: Option<String>,
    pub batch: Option<String>,
    pub label1: Option<String>,
    pub label2: Option<String>,
    pub label3: Option<String>,
    pub timepoint: Option<String>,
}
