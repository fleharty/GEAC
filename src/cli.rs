use std::path::PathBuf;

use clap::{Parser, Subcommand};

use crate::record::{FamilySizeScheme, ReadType, TotalFallback};

#[derive(Parser, Debug)]
#[command(
    name = "geac",
    about = "Genomic Evidence Atlas of Cohorts — collect alt base metrics from BAM/CRAM files",
    version
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Process a single BAM/CRAM file and write alt base records to Parquet
    Collect(CollectArgs),

    /// Merge per-sample Parquet files or existing DuckDB databases into a cohort DuckDB
    Merge(MergeArgs),

    /// Inspect a merged cohort DuckDB for schema, metadata, and resource-path issues
    Inspect(InspectArgs),

    /// Print a per-sample QC summary from one or more Parquet files
    Qc(QcArgs),

    /// Summarise recurrent loci across a cohort of Parquet files
    Cohort(CohortArgs),

    /// Find samples in a merged cohort DuckDB that appear to be from the same individual
    SampleIdentity(SampleIdentityArgs),

    /// Cross-annotate tumor alt-base loci against a paired normal BAM/CRAM
    AnnotateNormal(AnnotateNormalArgs),

    /// Cross-annotate tumor alt-base loci against a Panel of Normals DuckDB
    AnnotatePon(AnnotatePonArgs),

    /// Compute per-position coverage metrics from a BAM/CRAM file
    Coverage(CoverageArgs),

    /// Extract per-fragment metrics (insert size, GC content, end motifs) from a BAM/CRAM file
    Fragments(FragmentsArgs),

    /// Experimental tools (k-mer / fusion detection). APIs and outputs may
    /// change without notice; do not rely on these in production pipelines.
    Experimental(ExperimentalArgs),
}

#[derive(Parser, Debug)]
pub struct ExperimentalArgs {
    #[command(subcommand)]
    pub command: ExperimentalCommand,
}

#[derive(Subcommand, Debug)]
pub enum ExperimentalCommand {
    /// [EXPERIMENTAL] Build a gene-unique k-mer index from a reference FASTA + GTF
    BuildFusionIndex(BuildFusionIndexArgs),

    /// [EXPERIMENTAL] Detect gene fusions from a BAM/CRAM using a pre-built k-mer index
    Fusions(FusionsArgs),

    /// [EXPERIMENTAL] Extract read pairs containing k-mers unique to target genes
    ExtractGene(ExtractGeneArgs),

    /// [EXPERIMENTAL] Search a reference FASTA for all occurrences of a k-mer
    LocateKmer(LocateKmerArgs),

    /// [EXPERIMENTAL] Look up a k-mer sequence in a fusion index and report its gene + position
    LookupKmer(LookupKmerArgs),

    /// [EXPERIMENTAL] Scan a read sequence against a fusion index and report every k-mer hit
    ScanRead(ScanReadArgs),

    /// [EXPERIMENTAL] Aggregate per-sample kmer-hits TSVs from normal samples into a k-mer blacklist Parquet
    BuildFusionKmerBlacklist(BuildFusionKmerBlacklistArgs),

    /// [EXPERIMENTAL] For every genomic locus, compute the smallest k such that the k-mer at that position is unique genome-wide
    ComputeUniquenessMap(ComputeUniquenessMapArgs),

    /// [EXPERIMENTAL] List the k-mers shared between two genes' bodies (cross-gene k-mers that index building would discard)
    SharedKmers(SharedKmersArgs),

    /// [EXPERIMENTAL] Diagnose a suspected false-positive fusion call from a `fusions --reads-output` evidence BAM
    DiagnoseFusion(DiagnoseFusionArgs),
}

#[derive(Parser, Debug)]
pub struct BuildFusionKmerBlacklistArgs {
    /// One or more kmer-hits TSV files produced by `geac experimental fusions --kmer-hits-output`
    /// from normal (PoN) samples
    #[arg(long = "kmer-hits", required = true, num_args = 1..)]
    pub kmer_hits: Vec<PathBuf>,

    /// Output Parquet file (e.g. normals.fusion_kmer_blacklist.parquet)
    #[arg(short, long)]
    pub output: PathBuf,

    /// A k-mer must appear in at least this many distinct PoN samples to be included in the
    /// blacklist. Default 1 means any k-mer seen in any normal sample is blacklisted.
    #[arg(long, default_value_t = 1)]
    pub min_pon_samples: u32,
}

#[derive(Parser, Debug)]
pub struct CollectArgs {
    /// Input BAM or CRAM file
    #[arg(short, long)]
    pub input: PathBuf,

    /// Reference FASTA (required for CRAM and for ref allele lookup)
    #[arg(short = 'r', long)]
    pub reference: PathBuf,

    /// Sample identifier (used as sample_id in output records).
    /// If omitted, the SM tag from the BAM/CRAM read group header is used.
    /// Exits with an error if neither is provided.
    #[arg(short, long)]
    pub sample_id: Option<String>,

    /// Biological subject identifier — groups samples (e.g. multiple timepoints
    /// or tissue types) that come from the same person/animal/study subject.
    /// Used by `geac compare` to scope longitudinal and replicate analyses.
    #[arg(long)]
    pub subject_id: Option<String>,

    /// Sample substrate type (free-text). Standard values: cfDNA, urine_cfDNA,
    /// tumor_tissue, normal_tissue, buffy_coat, lymphocytes. Unrecognized
    /// values are preserved but skipped by mode-aware downstream tooling.
    #[arg(long)]
    pub sample_type: Option<String>,

    /// Optional batch label stored as a column in the output Parquet.
    /// Use this to tag samples with a processing group name so cohorts
    /// processed in separate runs can be filtered or compared in the Explorer.
    #[arg(long)]
    pub batch: Option<String>,

    /// Optional free-text label 1 for this sample (e.g. tissue type).
    /// Stored as a nullable string column `label1` in the output Parquet.
    #[arg(long)]
    pub label1: Option<String>,

    /// Optional free-text label 2 for this sample (e.g. library prep method).
    /// Stored as a nullable string column `label2` in the output Parquet.
    #[arg(long)]
    pub label2: Option<String>,

    /// Optional free-text label 3 for this sample (e.g. sequencer type).
    /// Stored as a nullable string column `label3` in the output Parquet.
    #[arg(long)]
    pub label3: Option<String>,

    /// Optional timepoint label for longitudinal / serial-sample studies.
    /// Stored as a nullable string column `timepoint` in the output Parquet.
    /// Values are arbitrary strings (e.g. "T0", "baseline", "week12") and are
    /// compared alphabetically in the Explorer; use zero-padded numbers or
    /// ISO-8601 dates for correct sort order.
    #[arg(long)]
    pub timepoint: Option<String>,

    /// Output Parquet file path.
    /// When --reads-output is also set, this path is used as a stem:
    /// e.g. "sample.parquet" → "sample.locus.parquet" + "sample.reads.parquet".
    #[arg(short, long)]
    pub output: PathBuf,

    /// Write per-read detail Parquet alongside the locus Parquet.
    /// The output path is derived from --output by replacing the extension:
    /// "sample.parquet" → "sample.locus.parquet" and "sample.reads.parquet".
    #[arg(long)]
    pub reads_output: bool,

    /// Compute SHA-256 for the input BAM/CRAM and store it in output Parquet provenance columns.
    /// Disabled by default because hashing large alignment files adds I/O and wall time.
    #[arg(long, default_value_t = false)]
    pub input_checksum_sha256: bool,

    /// Read type: raw, simplex, or duplex
    #[arg(long, default_value = "duplex")]
    pub read_type: ReadType,

    /// Pipeline that produced the BAM/CRAM. Free-text label stored verbatim in
    /// the output `pipeline` column. `fgbio` and `dragen` (case-insensitive) also
    /// select built-in family-size tag schemes (`aD`/`bD`/`cD` and `XV`/`XW`);
    /// any other value records the label and reads no family size unless
    /// `--family-size-tags` is given.
    #[arg(long)]
    pub pipeline: Option<String>,

    /// Override which aux tags encode family size, as comma-separated KEY=VALUE
    /// pairs. Keys: `ab`, `ba`, `total` (two-character BAM tags), and optional
    /// `fallback` (`sum` = use ab+ba when total is absent, like fgbio's cD; or
    /// `none`, the default). Example: `--family-size-tags ab=aD,ba=bD,total=cD,fallback=sum`.
    /// Takes precedence over the `--pipeline` preset. A tag absent from a read
    /// yields a null family size for that read.
    #[arg(long = "family-size-tags", value_name = "ab=XX,ba=YY,total=ZZ")]
    pub family_size_tags: Option<String>,

    /// Optional VCF/BCF file to annotate whether loci overlap called variants.
    /// Mutually exclusive with --variants-tsv.
    #[arg(long, conflicts_with = "variants_tsv")]
    pub vcf: Option<PathBuf>,

    /// Optional tab-separated variant list to use instead of a VCF.
    /// Expected columns: chrom  pos_start  pos_end  ref  var
    /// Coordinates are 0-based half-open (BED convention).
    /// Mutually exclusive with --vcf.
    #[arg(long, conflicts_with = "vcf")]
    pub variants_tsv: Option<PathBuf>,

    /// Optional GFF3 or GTF gene annotation file.
    /// When provided, each record is annotated with the gene name it overlaps.
    /// Only gene-level features are used (not exon/transcript).
    /// If omitted, or for intergenic loci, `gene` is null in the output Parquet.
    #[arg(long)]
    pub gene_annotations: Option<PathBuf>,

    /// Optional gnomAD VCF/BCF for allele-frequency annotation.
    /// The file must be bgzip-compressed and tabix/CSI-indexed.
    /// Each locus record gains a `gnomad_af` column with the INFO/AF value
    /// for the matching allele, or null if the allele is absent from gnomAD.
    #[arg(long)]
    pub gnomad: Option<PathBuf>,

    /// INFO field to use as the allele frequency from the gnomAD VCF.
    /// Defaults to "AF". Ignored when --gnomad is not set.
    #[arg(long, default_value = "AF")]
    pub gnomad_af_field: String,

    /// Optional BED file or Picard interval list of target regions.
    /// When provided, each record is annotated with `on_target = true/false`.
    /// If omitted, `on_target` is null in the output Parquet.
    /// Auto-detects format: files with `@` header lines are treated as Picard
    /// interval lists (1-based, end-inclusive); all others as BED (0-based, half-open).
    #[arg(long)]
    pub targets: Option<PathBuf>,

    /// Minimum base quality to consider a base
    #[arg(long, default_value_t = 1)]
    pub min_base_qual: u8,

    /// Minimum mapping quality to consider a read
    #[arg(long, default_value_t = 0)]
    pub min_map_qual: u8,

    /// Maximum reads per pileup column (0 = unlimited). htslib defaults to 8000,
    /// which silently downsamples at high-coverage loci and can drop rare-event
    /// reads (e.g. the single read supporting a long deletion).
    #[arg(long, default_value_t = 0)]
    pub max_pileup_depth: u32,

    /// Include PCR/optical duplicate reads (FLAG 0x400); excluded by default
    #[arg(long)]
    pub include_duplicates: bool,

    /// Include secondary alignments (FLAG 0x100); excluded by default
    #[arg(long)]
    pub include_secondary: bool,

    /// Include supplementary alignments (FLAG 0x800); excluded by default
    #[arg(long)]
    pub include_supplementary: bool,

    /// Exclude reads whose auxiliary tag equals a value, as `TAG:VALUE`
    /// (e.g. `--exclude-tag RX:bad`). The two-character TAG is matched against the
    /// read's aux value rendered as a string, so it works for string, char, and
    /// integer tags alike; the comparison is exact (`XY:3` drops `XY=3`, not
    /// `XY=30`). Reads lacking the tag are kept. Repeatable; a read is dropped if
    /// it matches any of the given filters.
    #[arg(long = "exclude-tag", value_name = "TAG:VALUE")]
    pub exclude_tag: Vec<String>,

    /// Restrict processing to a region. Accepts either a region string (e.g. "chr1:1000-2000")
    /// or a path to a BED file / Picard interval list. When a file path is given, the pileup
    /// is restricted to each interval in the file in turn.
    #[arg(long)]
    pub region: Option<String>,

    /// Window size (bases on each side of the locus) used to scan for homopolymers
    /// and short tandem repeats. Larger values detect longer repeat tracts but
    /// slightly increase per-locus cost.
    #[arg(long, default_value_t = 10)]
    pub repeat_window: usize,

    /// Canonical URI for the input BAM/CRAM (e.g. gs://bucket/sample.bam).
    /// When set, stored in output Parquet instead of the local --input path.
    /// Useful when running on Terra or other platforms that localize cloud files.
    #[arg(long)]
    pub bam_uri: Option<String>,

    /// Canonical URI for the BAM/CRAM index (e.g. gs://bucket/sample.bam.bai).
    /// When set, stored in output Parquet alongside bam_path for IGV session generation.
    #[arg(long)]
    pub bai_uri: Option<String>,

    /// Canonical URI for the variants file (VCF or TSV) passed via --vcf/--variants-tsv.
    /// When set, stored in output Parquet instead of the local file path.
    #[arg(long)]
    pub variants_uri: Option<String>,

    /// Canonical URI for the gnomAD VCF passed via --gnomad.
    /// When set, stored in output Parquet instead of the local file path.
    #[arg(long)]
    pub gnomad_uri: Option<String>,

    /// Canonical URI for the target BED / interval list passed via --targets.
    /// When set, stored in output Parquet instead of the local file path.
    #[arg(long)]
    pub targets_uri: Option<String>,

    /// Progress reporting interval in seconds (0 to disable)
    #[arg(long, default_value_t = 30)]
    pub progress_interval: u64,
}

#[derive(Parser, Debug)]
pub struct QcArgs {
    /// Input Parquet file(s)
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    /// Write a machine-readable TSV summary to this file (in addition to the stdout report)
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Restrict QC to on-target loci only (requires on_target column)
    #[arg(long)]
    pub on_target_only: bool,
}

#[derive(Parser, Debug)]
pub struct CohortArgs {
    /// Input Parquet files (one per sample)
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    /// Output file — written as Parquet if the extension is .parquet, TSV otherwise
    #[arg(short, long)]
    pub output: PathBuf,

    /// Minimum number of samples a locus must appear in to be reported (default: 2)
    #[arg(long, default_value_t = 2)]
    pub min_samples: u32,

    /// Minimum fraction of samples a locus must appear in, 0.0–1.0 (default: 0.0)
    #[arg(long, default_value_t = 0.0)]
    pub min_sample_fraction: f64,

    /// Restrict to on-target loci only (requires on_target column)
    #[arg(long)]
    pub on_target_only: bool,

    /// Number of top loci to print to stdout (by sample fraction, default: 20)
    #[arg(long, default_value_t = 20)]
    pub top_n: usize,
}

#[derive(Parser, Debug)]
pub struct SampleIdentityArgs {
    /// Merged cohort DuckDB produced by `geac merge`
    #[arg(short, long)]
    pub input: PathBuf,

    /// Output TSV of pairwise same-individual candidates
    #[arg(short, long)]
    pub output: PathBuf,

    /// Minimum total_depth for an SNV call to be usable as a fingerprint marker
    #[arg(long, default_value_t = 20)]
    pub min_depth: i32,

    /// Minimum VAF to treat an alt-base row as a germline non-reference call
    #[arg(long, default_value_t = 0.15)]
    pub min_vaf: f64,

    /// VAF at or above this threshold is treated as homozygous-alt; lower passing calls are heterozygous
    #[arg(long, default_value_t = 0.85)]
    pub hom_vaf: f64,

    /// Marker must be observed in at least this many samples
    #[arg(long, default_value_t = 2)]
    pub min_recurrence: i32,

    /// Maximum number of marker loci to use, ranked by recurrence
    #[arg(long, default_value_t = 5000)]
    pub max_markers: i32,

    /// Maximum number of samples to compare; set 0 to disable the guard
    #[arg(long, default_value_t = 5000)]
    pub max_samples: usize,

    /// Jaccard threshold for reporting likely same-individual pairs
    #[arg(long, default_value_t = 0.70)]
    pub min_jaccard: f64,

    /// Genotype concordance threshold at shared markers for reporting likely same-individual pairs
    #[arg(long, default_value_t = 0.95)]
    pub min_concordance: f64,

    /// Restrict marker panel to common gnomAD alleles when gnomad_af is present
    #[arg(long)]
    pub common_gnomad_only: bool,

    /// Lower gnomAD AF bound used with --common-gnomad-only
    #[arg(long, default_value_t = 0.05)]
    pub gnomad_min_af: f64,

    /// Upper gnomAD AF bound used with --common-gnomad-only
    #[arg(long, default_value_t = 0.95)]
    pub gnomad_max_af: f64,

    /// Write every pair sharing at least one marker instead of only threshold-passing pairs
    #[arg(long)]
    pub all_pairs: bool,
}

#[derive(Parser, Debug)]
pub struct MergeArgs {
    /// Input files: per-sample Parquet files and/or existing cohort.duckdb databases.
    /// Parquet files are routed by suffix (.reads.parquet, .coverage.parquet, etc.).
    /// DuckDB files (.duckdb) are attached and their data tables merged directly.
    /// Mix and match freely: e.g. new_sample.parquet existing_cohort.duckdb
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    /// Output DuckDB database file
    #[arg(short, long)]
    pub output: PathBuf,

    /// Write a TSV of on-target alt-base calls to this path.
    /// Columns mirror the main summary table in the GEAC explorer, ordered by
    /// chrom, pos, alt_allele, sample_id.  Silently skipped when alt_bases has
    /// no on_target column (i.e. no --targets BED was supplied to geac collect).
    #[arg(long)]
    pub on_target_tsv: Option<PathBuf>,
}

#[derive(Parser, Debug)]
pub struct InspectArgs {
    /// Merged cohort DuckDB produced by `geac merge`
    #[arg(short, long)]
    pub input: PathBuf,

    /// Exit non-zero when warnings are found, not only on structural errors
    #[arg(long)]
    pub strict: bool,
}

#[derive(Parser, Debug)]
pub struct AnnotateNormalArgs {
    /// Tumor locus Parquet produced by `geac collect`
    #[arg(long)]
    pub tumor_parquet: PathBuf,

    /// Normal BAM or CRAM file
    #[arg(long)]
    pub normal_bam: PathBuf,

    /// Reference FASTA (required for CRAM and for ref allele lookup)
    #[arg(short = 'r', long)]
    pub reference: PathBuf,

    /// Sample identifier for the normal sample.
    /// If omitted, the SM tag from the normal BAM/CRAM read group header is used.
    #[arg(long)]
    pub normal_sample_id: Option<String>,

    /// Output Parquet file path.  Should end in `.normal_evidence.parquet`
    /// so that `geac merge` routes it to the `normal_evidence` table.
    #[arg(short, long)]
    pub output: PathBuf,

    /// Minimum base quality to consider a base
    #[arg(long, default_value_t = 1)]
    pub min_base_qual: u8,

    /// Minimum mapping quality to consider a read
    #[arg(long, default_value_t = 0)]
    pub min_map_qual: u8,

    /// Maximum reads per pileup column (0 = unlimited). htslib defaults to 8000,
    /// which silently downsamples at high-coverage loci and can drop rare-event
    /// reads (e.g. the single read supporting a long deletion).
    #[arg(long, default_value_t = 0)]
    pub max_pileup_depth: u32,

    /// Include PCR/optical duplicate reads (FLAG 0x400); excluded by default
    #[arg(long)]
    pub include_duplicates: bool,

    /// Include secondary alignments (FLAG 0x100); excluded by default
    #[arg(long)]
    pub include_secondary: bool,

    /// Include supplementary alignments (FLAG 0x800); excluded by default
    #[arg(long)]
    pub include_supplementary: bool,
}

#[derive(Parser, Debug)]
pub struct AnnotatePonArgs {
    /// Tumor locus Parquet produced by `geac collect`
    #[arg(long)]
    pub tumor_parquet: PathBuf,

    /// PoN DuckDB database produced by `geac merge` from normal samples
    #[arg(long)]
    pub pon_db: PathBuf,

    /// Output Parquet file path.  Should end in `.pon_evidence.parquet`
    /// so that `geac merge` routes it to the `pon_evidence` table.
    #[arg(short, long)]
    pub output: PathBuf,
}

#[derive(Parser, Debug)]
pub struct CoverageArgs {
    /// Input BAM or CRAM file
    #[arg(short, long)]
    pub input: PathBuf,

    /// Reference FASTA (required for CRAM and GC content)
    #[arg(short = 'r', long)]
    pub reference: PathBuf,

    /// Sample identifier.
    /// If omitted, the SM tag from the BAM/CRAM read group header is used.
    #[arg(short, long)]
    pub sample_id: Option<String>,

    /// Biological subject identifier — groups samples (e.g. multiple timepoints
    /// or tissue types) that come from the same person/animal/study subject.
    #[arg(long)]
    pub subject_id: Option<String>,

    /// Sample substrate type (free-text). Standard values: cfDNA, urine_cfDNA,
    /// tumor_tissue, normal_tissue, buffy_coat, lymphocytes.
    #[arg(long)]
    pub sample_type: Option<String>,

    /// Optional batch label stored as a column in the output Parquet
    #[arg(long)]
    pub batch: Option<String>,

    /// Optional free-text label 1 for this sample (e.g. tissue type).
    /// Stored as a nullable string column `label1` in the output Parquet.
    #[arg(long)]
    pub label1: Option<String>,

    /// Optional free-text label 2 for this sample (e.g. library prep method).
    /// Stored as a nullable string column `label2` in the output Parquet.
    #[arg(long)]
    pub label2: Option<String>,

    /// Optional free-text label 3 for this sample (e.g. sequencer type).
    /// Stored as a nullable string column `label3` in the output Parquet.
    #[arg(long)]
    pub label3: Option<String>,

    /// Optional timepoint label for longitudinal / serial-sample studies.
    /// Stored as a nullable string column `timepoint` in the output Parquet.
    #[arg(long)]
    pub timepoint: Option<String>,

    /// Output Parquet file path (should end in .coverage.parquet)
    #[arg(short, long)]
    pub output: PathBuf,

    /// Read type: raw, simplex, or duplex
    #[arg(long, default_value = "duplex")]
    pub read_type: ReadType,

    /// Pipeline that produced the BAM/CRAM: fgbio, dragen, or raw
    #[arg(long)]
    pub pipeline: Option<String>,

    /// Optional BED file or Picard interval list of target regions.
    /// When provided, only positions within targets are emitted (including zero-depth positions).
    /// Without --targets, only positions with at least one read are emitted.
    #[arg(long)]
    pub targets: Option<PathBuf>,

    /// Restrict processing to a region. Accepts either a region string (e.g. "chr1:1000-2000")
    /// or a path to a BED file / Picard interval list. When a file path is given, the pileup
    /// is restricted to each interval in the file in turn.
    #[arg(long)]
    pub region: Option<String>,

    /// Optional GFF3 or GTF gene annotation file
    #[arg(long)]
    pub gene_annotations: Option<PathBuf>,

    /// Minimum mapping quality (used for total_depth; mapq stats computed from all non-dup reads)
    #[arg(long, default_value_t = 0)]
    pub min_map_qual: u8,

    /// Maximum reads per pileup column (0 = unlimited). htslib defaults to 8000,
    /// which silently downsamples at high-coverage loci.
    #[arg(long, default_value_t = 0)]
    pub max_pileup_depth: u32,

    /// Base quality threshold for frac_low_bq
    #[arg(long, default_value_t = 20)]
    pub min_base_qual: u8,

    /// Window size in bp (centred on position) for GC content computation
    #[arg(long, default_value_t = 100)]
    pub gc_window: usize,

    /// Suppress positions with total_depth strictly below this value (0 = keep all)
    #[arg(long, default_value_t = 0)]
    pub min_depth: i32,

    /// Aggregate consecutive positions into bins of this size (1 = per-position, no binning)
    #[arg(long, default_value_t = 1)]
    pub bin_size: i64,

    /// Progress reporting interval in seconds (0 to disable)
    #[arg(long, default_value_t = 20)]
    pub progress_interval: u64,

    /// Output path for per-interval summary Parquet (should end in .coverage.intervals.parquet).
    /// Requires --targets. When omitted, no interval file is written.
    #[arg(long)]
    pub intervals_output: Option<PathBuf>,

    /// Pre-computed annotation track in BEDGraph format, as NAME:path (repeatable).
    /// Each NAME becomes a nullable Float32 column in the output Parquet.
    /// Example: --track gem150:gem_150mer.bedgraph --track umap50:umap_k50.bedgraph
    #[arg(long = "track", value_name = "NAME:FILE")]
    pub tracks: Vec<String>,

    /// Fill in zero-depth positions across all reference contigs, even without --targets.
    /// Useful for WGS to detect dropout without a BED file.
    /// For whole-genome output, combine with --bin-size to keep output size manageable.
    /// Has no effect when --min-depth > 0 (zero-depth positions would always be filtered).
    #[arg(long)]
    pub fill_zeros: bool,
}

#[derive(Parser, Debug)]
pub struct FragmentsArgs {
    /// Input BAM or CRAM file
    #[arg(short, long)]
    pub input: PathBuf,

    /// Reference FASTA (required for GC content and end motifs)
    #[arg(short = 'r', long)]
    pub reference: PathBuf,

    /// Sample identifier.
    /// If omitted, the SM tag from the BAM/CRAM read group header is used.
    #[arg(short, long)]
    pub sample_id: Option<String>,

    /// Biological subject identifier — groups samples (e.g. multiple timepoints
    /// or tissue types) that come from the same person/animal/study subject.
    #[arg(long)]
    pub subject_id: Option<String>,

    /// Sample substrate type (free-text). Standard values: cfDNA, urine_cfDNA,
    /// tumor_tissue, normal_tissue, buffy_coat, lymphocytes.
    #[arg(long)]
    pub sample_type: Option<String>,

    /// Optional batch label stored as a column in the output Parquet
    #[arg(long)]
    pub batch: Option<String>,

    /// Optional free-text label 1 for this sample
    #[arg(long)]
    pub label1: Option<String>,

    /// Optional free-text label 2 for this sample
    #[arg(long)]
    pub label2: Option<String>,

    /// Optional free-text label 3 for this sample
    #[arg(long)]
    pub label3: Option<String>,

    /// Optional timepoint label for longitudinal studies
    #[arg(long)]
    pub timepoint: Option<String>,

    /// Output Parquet file path (should end in .fragments.parquet)
    #[arg(short, long)]
    pub output: PathBuf,

    /// Read type: raw, simplex, or duplex
    #[arg(long, default_value = "duplex")]
    pub read_type: ReadType,

    /// Pipeline that produced the BAM/CRAM: fgbio, dragen, or raw
    #[arg(long)]
    pub pipeline: Option<String>,

    /// Restrict processing to a region. Accepts either a region string (e.g. "chr1:1000-2000")
    /// or a path to a BED file / Picard interval list.
    #[arg(long)]
    pub region: Option<String>,

    /// Minimum mapping quality for R1 reads included in output
    #[arg(long, default_value_t = 0)]
    pub min_map_qual: u8,
}

#[derive(Parser, Debug)]
pub struct BuildFusionIndexArgs {
    /// Gene annotation file. Accepted formats:
    ///   GTF (.gtf, .gtf.gz)        — GENCODE, NCBI RefSeq, or UCSC GTF with 'gene' rows.
    ///   GenePred (.txt, .txt.gz)   — UCSC table-browser genePredExt+bin (e.g. ncbiRefSeq.txt.gz)
    ///                                or genePredExt without bin, or refFlat. Format is
    ///                                auto-detected from file extension and column layout.
    #[arg(long)]
    pub gene_annotation: PathBuf,

    /// Reference FASTA (must be indexed with samtools faidx — .fai required)
    #[arg(long)]
    pub fasta: PathBuf,

    /// K-mer length. Must match the value used later with `geac fusions`
    #[arg(long, default_value_t = 23)]
    pub kmer_size: u8,

    /// Drop genes with fewer than this many genome-unique k-mers; these genes
    /// have too little unique sequence to reliably detect fusions. Default 1
    /// keeps any gene with at least one usable k-mer; raise it to suppress
    /// repetitive genes with little unique sequence.
    #[arg(long, default_value_t = 1)]
    pub min_gene_kmers: u32,

    /// Output DuckDB index file (e.g. hg38_fusion_index.duckdb)
    #[arg(short, long)]
    pub output: PathBuf,

    /// Optional text file listing gene names to include (one per line).
    /// Lines starting with '#' and blank lines are ignored.
    /// When provided, only the listed genes are indexed; all others are skipped.
    /// Useful for building a small targeted index (e.g. a fusion gene panel).
    #[arg(long)]
    pub genes: Option<PathBuf>,

    /// Scan the full genome FASTA and remove k-mers that appear more than once
    /// genome-wide. Disabled by default because the candidate k-mer set (often
    /// >1 billion entries) requires tens of GB of RAM to hold in a HashMap.
    /// At k ≥ 19 cross-gene deduplication already eliminates nearly all
    /// non-unique k-mers; this flag adds a stricter guarantee at high memory cost.
    #[arg(long, default_value_t = false)]
    pub check_genome_uniqueness: bool,

    /// Retain k-mers occurring up to this many times genome-wide (instead of
    /// requiring exactly 1). A value > 1 keeps "near-unique" k-mers — useful for
    /// repetitive genes that would otherwise have little usable sequence. Only
    /// meaningful together with --check-genome-uniqueness; setting it > 1 without
    /// that flag is an error because no genome-wide counts are available.
    #[arg(long, default_value_t = 1)]
    pub max_genome_copies: u32,

    /// Optional TSV histogram of genome-wide k-mer copy number across the whole
    /// candidate set: two columns, `copies` and `n_kmers`. Requires
    /// --check-genome-uniqueness (the copy counts come from that pass).
    #[arg(long)]
    pub copy_histogram_output: Option<PathBuf>,

    /// Optional BED file of merged intervals covering all unique k-mer start positions.
    /// Adjacent k-mers are merged into a single interval, so the output is compact
    /// and suitable for use as a target BED with IGV, bedtools, or other tools.
    /// Chromosomes are emitted in FASTA index order.
    #[arg(long)]
    pub bed_output: Option<PathBuf>,

    /// Optional prefix for per-copy-tier BED files. Writes one BED per genome-wide
    /// copy number — `<prefix>.copies1.bed` (strictly unique), `<prefix>.copies2.bed`
    /// (occurs twice), … up to --max-genome-copies — so unique and near-unique regions
    /// can be loaded/colored separately in IGV or used as distinct target BEDs.
    /// Requires --check-genome-uniqueness (copy counts come from that pass).
    #[arg(long)]
    pub bed_output_by_copies: Option<PathBuf>,

    /// Optional per-gene uniqueness TSV. For each indexed gene reports the body
    /// length, count of unique k-mers, and two "% of gene that is k-mer unique"
    /// metrics: pct_unique_windows (unique k-mers / k-mer windows) and
    /// pct_unique_bases (gene-body bases covered by a unique k-mer). Computed
    /// before the --min-gene-kmers drop so even excluded genes report real values.
    #[arg(long)]
    pub gene_stats_output: Option<PathBuf>,

    /// Discard candidate k-mers that have any reference genome k-mer within
    /// Hamming distance N. Use --edit-distance-filter 1 to retain only k-mers
    /// where no single-base sequencing error in a reference read can produce the
    /// candidate k-mer. Triggers a full genome scan (similar overhead to
    /// --check-genome-uniqueness). Default 0 disables the filter.
    #[arg(long, default_value_t = 0)]
    pub edit_distance_filter: u32,

    /// Number of threads for the parallel genome scan used by
    /// --check-genome-uniqueness and --edit-distance-filter. 0 uses all logical
    /// CPUs on the machine. Has no effect when neither genome-scan flag is set.
    #[arg(long, default_value_t = 0)]
    pub threads: u32,
}

#[derive(Parser, Debug)]
pub struct FusionsArgs {
    /// Input BAM or CRAM file
    #[arg(long)]
    pub bam: PathBuf,

    /// Fusion k-mer index produced by `geac build-fusion-index`
    #[arg(long)]
    pub index: PathBuf,

    /// Reference FASTA (required for CRAM decoding; optional for BAM)
    #[arg(short = 'r', long)]
    pub reference: Option<PathBuf>,

    /// Sample identifier written to the output Parquet.
    /// If omitted, the SM tag from the BAM/CRAM read group header is used.
    #[arg(long)]
    pub sample_id: Option<String>,

    /// Output Parquet file (e.g. sample.fusions.parquet)
    #[arg(short, long)]
    pub output: PathBuf,

    /// K-mer length — must match the value used in `geac build-fusion-index`
    #[arg(long, default_value_t = 23)]
    pub kmer_size: u8,

    /// Minimum number of unique k-mer hits to a gene for a read to be assigned
    /// to that gene. Lower values increase sensitivity; higher values reduce noise
    #[arg(long, default_value_t = 3)]
    pub min_kmer_hits: u32,

    /// Minimum fragment pairs supporting a fusion to include it in the output
    #[arg(long, default_value_t = 2)]
    pub min_supporting_reads: u32,

    /// Minimum mapping quality for reads to be considered. Unmapped reads
    /// bypass this filter (k-mer matching is alignment-independent).
    #[arg(long, default_value_t = 0)]
    pub min_mapq: u8,

    /// Optional output BAM containing all reads from fusion-supporting fragments.
    /// Both reads of a pair are written even if only one read was gene-assigned.
    #[arg(long)]
    pub reads_output: Option<PathBuf>,

    /// Optional TSV file for human-readable fusion results.
    /// Columns: sample_id, gene_a, gene_b, chrom_a, chrom_b, supporting_reads, min_mapq,
    /// n_spanning_reads, n_coherent_fragments, n_pon_samples, pon_total_samples,
    /// max_pon_supporting_reads, filter
    #[arg(long)]
    pub tsv_output: Option<PathBuf>,

    /// Optional TSV with one row per k-mer hit from fusion-supporting reads.
    /// Columns: fusion, sample_id, read_name, read_end, chrom, pos, gene_matched,
    /// kmer_pos_in_read, kmer_hash, kmer_seq
    #[arg(long)]
    pub kmer_hits_output: Option<PathBuf>,

    /// Optional TSV with one breakpoint row per detected fusion.
    /// Columns: fusion, gene_a, chrom_a, breakpoint_a, bp_a_n, bp_a_std,
    ///          gene_b, chrom_b, breakpoint_b, bp_b_n, bp_b_std, n_spanning_reads
    #[arg(long)]
    pub breakpoints_output: Option<PathBuf>,

    /// Ignore k-mers occurring more than this many times genome-wide when
    /// assigning reads to genes. Requires an index built with
    /// --check-genome-uniqueness (which records per-k-mer copy number); using
    /// this flag against an older index without that data is an error.
    #[arg(long)]
    pub max_kmer_copies: Option<u32>,

    /// Optional fusion Panel-of-Normals DuckDB (a `geac merge` of normal-sample
    /// *.fusions.parquet files). Each output fusion is annotated with PoN columns
    /// (n_pon_samples, pon_total_samples, max_pon_supporting_reads).
    #[arg(long)]
    pub fusion_pon: Option<PathBuf>,

    /// Flag fusions seen in strictly more than this many PoN samples with
    /// filter="pon" (rows are kept, not dropped — downstream can exclude them).
    /// Requires --fusion-pon. Default: annotate only, every row stays filter="PASS".
    #[arg(long)]
    pub max_pon_samples: Option<u32>,

    /// Require at least this many spanning reads (single reads hitting both
    /// fusion genes) to show a coherent A-block→B-block k-mer partition.
    /// Spanning reads with interleaved gene k-mers indicate homology/paralog
    /// artifacts rather than a real junction. Default: 0 (disabled — report
    /// all candidates regardless of coherence). Set to 1 or higher to enable.
    #[arg(long, default_value_t = 0)]
    pub min_coherent_fragments: u32,

    /// Minimum k-mer hits from each gene for a spanning read to be considered
    /// anchored on both sides of the junction. Reads with fewer than this many
    /// k-mers from either gene are counted as spanning but not coherent.
    /// Only relevant when --min-coherent-fragments > 0.
    #[arg(long, default_value_t = 3)]
    pub min_anchor_kmers: u32,

    /// Optional k-mer blacklist Parquet produced by `geac experimental build-fusion-kmer-blacklist`.
    /// K-mers in this file are considered PoN-noisy: a read must have at least `--min-kmer-hits`
    /// *non-blacklisted* k-mer matches to contribute to fusion evidence. Blacklisted k-mers
    /// are not stripped — they still help once a read has already passed on clean evidence.
    #[arg(long)]
    pub fusion_kmer_blacklist: Option<PathBuf>,

    /// A k-mer must appear in at least this many PoN samples (n_pon_samples column in the
    /// blacklist Parquet) to be treated as blacklisted at call time.
    /// Allows reusing one blacklist file with different stringency thresholds.
    #[arg(long, default_value_t = 1)]
    pub min_kmer_blacklist_samples: u32,

    /// Warn instead of failing when the fusion index or PoN was built with a
    /// different GEAC version than the running binary. The warning names both
    /// versions. Omitting this flag is the safe default: a version mismatch is
    /// an error so schema or algorithm changes do not silently corrupt results.
    #[arg(long, default_value_t = false)]
    pub skip_version_check: bool,
}

#[derive(Parser, Debug)]
pub struct ExtractGeneArgs {
    /// Input BAM or CRAM file
    #[arg(long)]
    pub bam: PathBuf,

    /// Fusion k-mer index produced by `geac build-fusion-index`
    #[arg(long)]
    pub index: PathBuf,

    /// Reference FASTA (required for CRAM decoding; optional for BAM)
    #[arg(short = 'r', long)]
    pub reference: Option<PathBuf>,

    /// Sample identifier written to TSV output.
    /// If omitted, the SM tag from the BAM/CRAM read group header is used.
    #[arg(long)]
    pub sample_id: Option<String>,

    /// Gene name(s) to extract reads for. Repeatable: --gene PAX3 --gene FOXO1
    #[arg(long = "gene", required = true)]
    pub genes: Vec<String>,

    /// Output BAM containing all reads from fragments with a k-mer hit to any target gene.
    /// Both reads of a pair are written even if only one read carries a k-mer hit.
    #[arg(short, long)]
    pub output: PathBuf,

    /// K-mer length — must match the value used in `geac build-fusion-index`
    #[arg(long, default_value_t = 23)]
    pub kmer_size: u8,

    /// Minimum number of unique k-mer hits to a target gene for a read to count as matching.
    /// Defaults to 1: any read with a single unique k-mer from the target gene is captured.
    /// At k=19 a single k-mer match is already highly specific. Raise this only if you see
    /// false-positive reads from repetitive regions adjacent to the target gene.
    #[arg(long, default_value_t = 1)]
    pub min_kmer_hits: u32,

    /// Minimum mapping quality for reads to be considered during gene assignment.
    /// Defaults to 0 so that multimapping reads are not silently dropped. Reads with
    /// mapq=0 that contain unique k-mers are still highly likely to belong to the target
    /// gene; raise this only to restrict output to confidently-placed alignments.
    #[arg(long, default_value_t = 0)]
    pub min_mapq: u8,

    /// Optional TSV with one row per k-mer hit from matching reads.
    /// Columns: gene_matched, sample_id, read_name, read_end, chrom, pos, kmer_hash, kmer_seq
    #[arg(long)]
    pub kmer_hits_output: Option<PathBuf>,

    /// Optional BAM for reads from fragments with no k-mer match to any target gene.
    /// When provided, every record in the input BAM appears in exactly one of
    /// --output (fragment has a k-mer hit) or --complement-output (fragment has none).
    /// A BAI index is created automatically alongside this file.
    #[arg(long)]
    pub complement_output: Option<PathBuf>,

    /// Diagnostic: dump per-read k-mer-matching detail for reads overlapping
    /// this region. Format: "chr:start-end" (1-based, inclusive). When set,
    /// every primary read whose alignment overlaps the region is logged with
    /// its k-mer extraction and match counts. Use this to debug "expected
    /// reads missing" scenarios.
    #[arg(long)]
    pub debug_region: Option<String>,
}

#[derive(Parser, Debug)]
pub struct LocateKmerArgs {
    /// Reference FASTA (must be indexed with samtools faidx — .fai required)
    #[arg(long)]
    pub fasta: PathBuf,

    /// K-mer sequence to search for (e.g. ACGTACGTACGTACGTACG).
    /// Both the k-mer and its reverse complement are reported; palindromic
    /// k-mers (equal to their own RC) are reported once with strand "+".
    /// Length must be 1–31.
    #[arg(long)]
    pub kmer: String,

    /// Output TSV file. When omitted, results are written to stdout.
    /// Columns: chrom, start (0-based), end, strand (+/-)
    /// When --gene-annotations is provided, additional columns: gene, region
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Optional gene annotation file (GTF, GFF3, or genePred) used to label
    /// each k-mer hit with the gene name and genomic region (exon/CDS/UTR/intron/intergenic).
    #[arg(long)]
    pub gene_annotations: Option<PathBuf>,
}

#[derive(Parser, Debug)]
pub struct LookupKmerArgs {
    /// Fusion index built by `geac experimental build-fusion-index` (DuckDB)
    #[arg(long)]
    pub index: PathBuf,

    /// K-mer sequence to look up (e.g. ACGTACGTACGTACGTACG). The canonical form
    /// (min of forward and reverse complement) is computed and queried against
    /// the index. Length must be 1–31.
    #[arg(long)]
    pub kmer: String,
}

#[derive(Parser, Debug)]
pub struct ScanReadArgs {
    /// Fusion index built by `geac experimental build-fusion-index` (DuckDB)
    #[arg(long)]
    pub index: PathBuf,

    /// Read sequence to scan (e.g. the SEQ field from a BAM record). Must contain
    /// only ACGT/acgt bases; N bases reset the sliding window just as in fusions.
    #[arg(long)]
    pub read: String,

    /// K-mer size used when building the index.
    #[arg(long, default_value = "23")]
    pub kmer_size: u8,
}

#[derive(Parser, Debug)]
pub struct ComputeUniquenessMapArgs {
    /// Reference FASTA (must be indexed with samtools faidx — .fai required)
    #[arg(long)]
    pub fasta: PathBuf,

    /// Output bedgraph file. Each position records the smallest k for which its
    /// k-mer is genome-unique. Adjacent positions with the same value are merged
    /// into a single interval (run-length encoding). The bedgraph is full base-
    /// pair resolution; most intervals will be 1 bp wide in practice.
    #[arg(short, long)]
    pub output: PathBuf,

    /// Minimum k-mer size to test. K-mers shorter than this are never checked.
    #[arg(long, default_value_t = 15)]
    pub min_k: usize,

    /// Maximum k-mer size to test. Positions with no unique k-mer in [min_k,
    /// max_k] are assigned max_k + 1. Must be ≤ 31 (hardware limit of the 2-bit
    /// rolling k-mer encoder).
    #[arg(long, default_value_t = 31)]
    pub max_k: usize,

    /// Optional BED file of regions to include in the output. The genome-wide
    /// count pass always scans the full FASTA (uniqueness is a global property),
    /// but only positions overlapping these regions are written to the output.
    /// When omitted, every genomic position is written.
    #[arg(long)]
    pub regions: Option<PathBuf>,

    /// Write one output line per base rather than merging adjacent positions
    /// with the same value. Produces a larger file but guarantees a fixed-step
    /// layout. Equivalent in content to the default run-length-encoded output.
    #[arg(long, default_value_t = false)]
    pub no_merge: bool,
}

#[derive(Parser, Debug)]
pub struct SharedKmersArgs {
    /// Gene annotation file (GTF or GenePred), same formats as build-fusion-index.
    /// Gene bodies are taken exactly as build-fusion-index would see them, so the
    /// shared set here is what cross-gene dedup removes from the index.
    #[arg(long)]
    pub gene_annotation: PathBuf,

    /// Reference FASTA (must be indexed with samtools faidx — .fai required)
    #[arg(long)]
    pub fasta: PathBuf,

    /// First gene name (must match a gene_name in the annotation)
    #[arg(long)]
    pub gene_a: String,

    /// Second gene name (must match a gene_name in the annotation)
    #[arg(long)]
    pub gene_b: String,

    /// K-mer length. Use the same value as the index you are debugging.
    #[arg(long, default_value_t = 23)]
    pub kmer_size: u8,

    /// Output TSV path. Writes to stdout when omitted.
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Maximum substitution (Hamming) distance allowed between a gene-A k-mer and
    /// a gene-B k-mer for them to be reported as a match. 0 (default) reports only
    /// exact shared k-mers; N > 0 additionally pairs *near*-identical k-mers that
    /// differ by up to N bases — useful for diverged paralogs that share little
    /// exact sequence. Each output row carries both sequences and their ab_dist.
    #[arg(long, default_value_t = 0)]
    pub edit_distance: u32,

    /// Scan the full reference FASTA and report, for each matched k-mer, how many
    /// times its canonical form occurs genome-wide (`ref_copies_a` / `ref_copies_b`
    /// columns). A value of 1 means the k-mer is unique to its gene body; higher
    /// counts mean it also occurs elsewhere (other genes, paralogs, repeats).
    /// Triggers a full-genome scan.
    #[arg(long, default_value_t = false)]
    pub check_reference: bool,

    /// Threads for the --check-reference genome scan. 0 uses all logical CPUs.
    /// No effect unless --check-reference is set.
    #[arg(long, default_value_t = 0)]
    pub threads: u32,

    /// Optional built fusion index (DuckDB). When given, each matched k-mer is
    /// annotated with whether it is an actual *index* k-mer for its gene
    /// (`in_index_a` / `in_index_b` columns) — i.e. one the fusion caller can vote
    /// on. The false-positive-driving pairs are those where `ab_dist > 0` and the
    /// voted gene's side is in the index. Must have the same k-mer size.
    #[arg(long)]
    pub index: Option<PathBuf>,
}

#[derive(Parser, Debug)]
pub struct DiagnoseFusionArgs {
    /// Evidence BAM produced by `geac experimental fusions --reads-output`
    /// (reads tagged `FX:Z:GENEA::GENEB`). May also be a CRAM with --reference.
    #[arg(long)]
    pub reads: PathBuf,

    /// Fusion k-mer index (DuckDB) used to make the call. Must share the k-mer size.
    #[arg(long)]
    pub index: PathBuf,

    /// Reference FASTA (required only if --reads is a CRAM).
    #[arg(short = 'r', long)]
    pub reference: Option<PathBuf>,

    /// First gene of the suspected fusion pair.
    #[arg(long)]
    pub gene_a: String,

    /// Second gene of the suspected fusion pair.
    #[arg(long)]
    pub gene_b: String,

    /// K-mer length — must match the index and the original `fusions` run.
    #[arg(long, default_value_t = 23)]
    pub kmer_size: u8,

    /// Minimum k-mers from each gene for a spanning read to count as anchored on
    /// both sides (matches the `fusions --min-anchor-kmers` default).
    #[arg(long, default_value_t = 3)]
    pub min_anchor: u32,

    /// Human-readable diagnostic report. Writes to stdout when omitted.
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Optional per-read TSV with the gene-A/gene-B k-mer layout track for every
    /// evidence read. When omitted, a short preview is included in the report.
    #[arg(long)]
    pub per_read_output: Option<PathBuf>,
}

// Allow clap to parse ReadType and Pipeline from strings

impl std::str::FromStr for ReadType {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "raw" => Ok(ReadType::Raw),
            "simplex" => Ok(ReadType::Simplex),
            "duplex" => Ok(ReadType::Duplex),
            _ => anyhow::bail!(
                "invalid read-type '{}': expected raw, simplex, or duplex",
                s
            ),
        }
    }
}

impl CollectArgs {
    /// Resolve the family-size tag scheme for this run: an explicit
    /// `--family-size-tags` spec takes precedence; otherwise the `--pipeline`
    /// label selects a built-in preset (`fgbio`/`dragen`), falling back to no
    /// family-size extraction for any other (or absent) pipeline.
    pub fn family_size_scheme(&self) -> anyhow::Result<FamilySizeScheme> {
        resolve_family_size_scheme(self.pipeline.as_deref(), self.family_size_tags.as_deref())
    }
}

/// Pick a [`FamilySizeScheme`] from a pipeline label and/or an explicit tag spec.
fn resolve_family_size_scheme(
    pipeline: Option<&str>,
    spec: Option<&str>,
) -> anyhow::Result<FamilySizeScheme> {
    if let Some(spec) = spec {
        return parse_family_size_spec(spec);
    }
    Ok(match pipeline.map(str::to_lowercase).as_deref() {
        Some("fgbio") => FamilySizeScheme::fgbio(),
        Some("dragen") => FamilySizeScheme::dragen(),
        _ => FamilySizeScheme::none(),
    })
}

/// Parse `--family-size-tags ab=aD,ba=bD,total=cD,fallback=sum` into a scheme.
fn parse_family_size_spec(spec: &str) -> anyhow::Result<FamilySizeScheme> {
    let mut scheme = FamilySizeScheme::none();
    for entry in spec.split(',') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        let (key, value) = entry.split_once('=').ok_or_else(|| {
            anyhow::anyhow!("--family-size-tags entry '{entry}' must be KEY=VALUE, e.g. total=cD")
        })?;
        match key.trim() {
            "ab" => scheme.ab_tag = Some(parse_aux_tag(value.trim())?),
            "ba" => scheme.ba_tag = Some(parse_aux_tag(value.trim())?),
            "total" => scheme.total_tag = Some(parse_aux_tag(value.trim())?),
            "fallback" => {
                scheme.total_fallback = match value.trim().to_lowercase().as_str() {
                    "sum" => TotalFallback::SumAbBa,
                    "none" => TotalFallback::None,
                    other => anyhow::bail!(
                        "--family-size-tags fallback '{other}' must be 'sum' or 'none'"
                    ),
                }
            }
            other => anyhow::bail!(
                "--family-size-tags key '{other}' must be one of ab, ba, total, fallback"
            ),
        }
    }
    if scheme.ab_tag.is_none() && scheme.ba_tag.is_none() && scheme.total_tag.is_none() {
        anyhow::bail!("--family-size-tags must set at least one of ab, ba, total");
    }
    Ok(scheme)
}

/// Parse a two-character BAM auxiliary tag name into a byte pair.
fn parse_aux_tag(s: &str) -> anyhow::Result<[u8; 2]> {
    let bytes = s.as_bytes();
    if bytes.len() != 2 {
        anyhow::bail!("--family-size-tags tag '{s}' must be exactly two characters (e.g. cD)");
    }
    Ok([bytes[0], bytes[1]])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pipeline_presets_resolve_to_builtin_schemes() {
        assert_eq!(
            resolve_family_size_scheme(Some("fgbio"), None).unwrap(),
            FamilySizeScheme::fgbio()
        );
        // Case-insensitive.
        assert_eq!(
            resolve_family_size_scheme(Some("DRAGEN"), None).unwrap(),
            FamilySizeScheme::dragen()
        );
    }

    #[test]
    fn unknown_pipeline_resolves_to_no_family_size() {
        assert_eq!(
            resolve_family_size_scheme(Some("myprep_v3"), None).unwrap(),
            FamilySizeScheme::none()
        );
        assert_eq!(
            resolve_family_size_scheme(None, None).unwrap(),
            FamilySizeScheme::none()
        );
    }

    #[test]
    fn explicit_spec_overrides_pipeline_preset() {
        // Even with pipeline=fgbio, an explicit spec wins.
        let scheme = resolve_family_size_scheme(Some("fgbio"), Some("total=fz")).unwrap();
        assert_eq!(scheme.total_tag, Some(*b"fz"));
        assert_eq!(scheme.ab_tag, None);
        assert_eq!(scheme.total_fallback, TotalFallback::None);
    }

    #[test]
    fn spec_parses_all_keys_including_sum_fallback() {
        let scheme = parse_family_size_spec("ab=aD, ba=bD, total=cD, fallback=sum").unwrap();
        assert_eq!(
            scheme,
            FamilySizeScheme {
                ab_tag: Some(*b"aD"),
                ba_tag: Some(*b"bD"),
                total_tag: Some(*b"cD"),
                total_fallback: TotalFallback::SumAbBa,
            }
        );
    }

    #[test]
    fn spec_rejects_bad_tag_length_and_unknown_keys() {
        assert!(parse_family_size_spec("total=ABC").is_err());
        assert!(parse_family_size_spec("foo=aD").is_err());
        assert!(parse_family_size_spec("total").is_err());
        assert!(parse_family_size_spec("fallback=bogus").is_err());
        // No tag keys set at all.
        assert!(parse_family_size_spec("fallback=sum").is_err());
    }
}
