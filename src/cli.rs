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

    /// Merge per-sample Parquet files into a cohort DuckDB database
    Merge(MergeArgs),

    /// Print a per-sample QC summary from one or more Parquet files
    Qc(QcArgs),

    /// Summarise recurrent loci across a cohort of Parquet files
    Cohort(CohortArgs),

    /// Cross-annotate tumor alt-base loci against a paired normal BAM/CRAM
    AnnotateNormal(AnnotateNormalArgs),

    /// Cross-annotate tumor alt-base loci against a Panel of Normals DuckDB
    AnnotatePon(AnnotatePonArgs),
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

    /// Optional batch label stored as a column in the output Parquet.
    /// Use this to tag samples with a processing group name so cohorts
    /// processed in separate runs can be filtered or compared in the Explorer.
    #[arg(long)]
    pub batch: Option<String>,

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

    /// Read type: raw, simplex, or duplex
    #[arg(long, default_value = "duplex")]
    pub read_type: ReadType,

    /// Pipeline that produced the BAM/CRAM: fgbio, dragen, or raw
    #[arg(long, default_value = "fgbio")]
    pub pipeline: Pipeline,

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

    /// Include PCR/optical duplicate reads (FLAG 0x400); excluded by default
    #[arg(long)]
    pub include_duplicates: bool,

    /// Include secondary alignments (FLAG 0x100); excluded by default
    #[arg(long)]
    pub include_secondary: bool,

    /// Include supplementary alignments (FLAG 0x800); excluded by default
    #[arg(long)]
    pub include_supplementary: bool,

    /// Restrict processing to this region (e.g. "chr1:1000-2000")
    #[arg(long)]
    pub region: Option<String>,

    /// Window size (bases on each side of the locus) used to scan for homopolymers
    /// and short tandem repeats. Larger values detect longer repeat tracts but
    /// slightly increase per-locus cost.
    #[arg(long, default_value_t = 10)]
    pub repeat_window: usize,

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
    /// Input Parquet files or a glob pattern (e.g. "samples/*.parquet")
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

// Allow clap to parse ReadType and Pipeline from strings

impl std::str::FromStr for ReadType {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "raw" => Ok(ReadType::Raw),
            "simplex" => Ok(ReadType::Simplex),
            "duplex" => Ok(ReadType::Duplex),
            _ => anyhow::bail!("invalid read-type '{}': expected raw, simplex, or duplex", s),
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
