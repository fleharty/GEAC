use std::path::PathBuf;

use clap::{Parser, Subcommand};

use crate::record::{Pipeline, ReadType};

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

    /// Print a per-sample QC summary from one or more Parquet files
    Qc(QcArgs),

    /// Summarise recurrent loci across a cohort of Parquet files
    Cohort(CohortArgs),

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

    /// Pipeline that produced the BAM/CRAM: fgbio, dragen, or raw
    #[arg(long)]
    pub pipeline: Option<Pipeline>,

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

    /// Canonical URI for the variants file (VCF or TSV) passed via --vcf/--variants-tsv.
    /// When set, stored in output Parquet instead of the local file path.
    #[arg(long)]
    pub variants_uri: Option<String>,

    /// Canonical URI for the gnomAD VCF passed via --gnomad.
    /// When set, stored in output Parquet instead of the local file path.
    #[arg(long)]
    pub gnomad_uri: Option<String>,

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
    pub pipeline: Option<Pipeline>,

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
    pub pipeline: Option<Pipeline>,

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
    /// Gene annotation file in GTF format (.gtf or .gtf.gz)
    #[arg(long)]
    pub gtf: PathBuf,

    /// Reference FASTA (must be indexed with samtools faidx — .fai required)
    #[arg(long)]
    pub fasta: PathBuf,

    /// K-mer length. Must match the value used later with `geac fusions`
    #[arg(long, default_value_t = 23)]
    pub kmer_size: u8,

    /// Drop genes with fewer than this many genome-unique k-mers; these genes
    /// have too little unique sequence to reliably detect fusions (default: 100)
    #[arg(long, default_value_t = 100)]
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

    /// Optional BED file of merged intervals covering all unique k-mer start positions.
    /// Adjacent k-mers are merged into a single interval, so the output is compact
    /// and suitable for use as a target BED with IGV, bedtools, or other tools.
    /// Chromosomes are emitted in FASTA index order.
    #[arg(long)]
    pub bed_output: Option<PathBuf>,
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
    /// Columns: sample_id, gene_a, gene_b, chrom_a, chrom_b, supporting_reads, min_mapq
    #[arg(long)]
    pub tsv_output: Option<PathBuf>,

    /// Optional TSV with one row per k-mer hit from fusion-supporting reads.
    /// Columns: fusion, sample_id, read_name, gene_matched, kmer_hash, kmer_seq
    #[arg(long)]
    pub kmer_hits_output: Option<PathBuf>,
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

impl std::str::FromStr for Pipeline {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "fgbio" => Ok(Pipeline::Fgbio),
            "dragen" => Ok(Pipeline::Dragen),
            "raw" => Ok(Pipeline::Raw),
            _ => anyhow::bail!("invalid pipeline '{}': expected fgbio, dragen, or raw", s),
        }
    }
}
