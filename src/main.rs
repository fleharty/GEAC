mod bam;
mod cli;
mod cohort;
mod coverage;
mod gene_annotations;
mod gnomad;
mod merge;
mod normal;
mod pon;
mod progress;
mod qc;
mod record;
mod repeat;
mod targets;
mod variants_tsv;
mod vcf;
mod writer;

use anyhow::Result;
use clap::Parser;
use tracing::info;

use cli::{Cli, Command};

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::from_default_env()
                .add_directive(tracing::Level::INFO.into()),
        )
        .init();

    let cli = Cli::parse();

    match cli.command {
        Command::Collect(args) => {
            info!(
                sample_id = %args.sample_id.as_deref().unwrap_or("<from SM tag>"),
                input = %args.input.display(),
                output = %args.output.display(),
                "collecting alt bases"
            );

            // Load whichever annotation source was provided (at most one).
            let vcf_index = args
                .vcf
                .as_ref()
                .map(|p| {
                    info!(vcf = %p.display(), "loading VCF annotations");
                    vcf::VcfIndex::load(p)
                })
                .transpose()?;

            let tsv_index = args
                .variants_tsv
                .as_ref()
                .map(|p| {
                    info!(variants_tsv = %p.display(), "loading variants TSV annotations");
                    variants_tsv::VariantsTsv::load(p)
                })
                .transpose()?;

            let annotator: Option<&dyn vcf::VariantAnnotator> = match (&vcf_index, &tsv_index) {
                (Some(v), _) => Some(v),
                (_, Some(t)) => Some(t),
                _ => None,
            };

            let target_intervals = args
                .targets
                .as_ref()
                .map(|p| {
                    info!(targets = %p.display(), "loading target intervals");
                    targets::TargetIntervals::load(p)
                })
                .transpose()?;

            if let Some(ref ti) = target_intervals {
                info!(n_targets = ti.n_targets(), "target intervals loaded");
            }

            let gene_annots = args
                .gene_annotations
                .as_ref()
                .map(|p| {
                    info!(gene_annotations = %p.display(), "loading gene annotations");
                    gene_annotations::GeneAnnotations::load(p)
                })
                .transpose()?;

            if let Some(ref ga) = gene_annots {
                info!(n_genes = ga.n_genes(), "gene annotations loaded");
            }

            let mut gnomad_index = args
                .gnomad
                .as_ref()
                .map(|p| {
                    info!(gnomad = %p.display(), af_field = %args.gnomad_af_field, "loading gnomAD VCF");
                    gnomad::GnomadIndex::open(p, Some(&args.gnomad_af_field))
                })
                .transpose()?;

            let (records, read_records) = bam::collect_alt_bases(
                &args,
                annotator,
                target_intervals.as_ref().map(|t| t as &_),
                gene_annots.as_ref(),
                gnomad_index.as_mut(),
            )?;

            let (locus_output, reads_output_path) = if args.reads_output {
                let stem = args.output.file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("output");
                let parent = args.output.parent().unwrap_or(std::path::Path::new("."));
                (
                    parent.join(format!("{stem}.locus.parquet")),
                    Some(parent.join(format!("{stem}.reads.parquet"))),
                )
            } else {
                (args.output.clone(), None)
            };

            info!(n_records = records.len(), output = %locus_output.display(), "writing locus Parquet");
            writer::parquet::write_parquet(&records, &locus_output)?;

            if let Some(ref reads_path) = reads_output_path {
                info!(n_records = read_records.len(), output = %reads_path.display(), "writing reads Parquet");
                writer::parquet_reads::write_parquet(&read_records, reads_path)?;
            }

            info!("done");
        }

        Command::Merge(args) => {
            merge::merge(&args)?;
        }

        Command::Qc(args) => {
            qc::run_qc(&args)?;
        }

        Command::Cohort(args) => {
            cohort::run_cohort(&args)?;
        }

        Command::AnnotatePon(args) => {
            info!(
                tumor_parquet = %args.tumor_parquet.display(),
                pon_db        = %args.pon_db.display(),
                output        = %args.output.display(),
                "annotating PoN evidence"
            );

            let records = pon::annotate_pon(&args)?;

            info!(n_records = records.len(), output = %args.output.display(), "writing PoN evidence Parquet");
            writer::parquet_pon::write_parquet(&records, &args.output)?;

            info!("done");
        }

        Command::AnnotateNormal(args) => {
            info!(
                tumor_parquet = %args.tumor_parquet.display(),
                normal_bam    = %args.normal_bam.display(),
                output        = %args.output.display(),
                "annotating normal evidence"
            );

            let records = normal::annotate_normal(&args)?;

            info!(n_records = records.len(), output = %args.output.display(), "writing normal evidence Parquet");
            writer::parquet_normal::write_parquet(&records, &args.output)?;

            info!("done");
        }

        Command::Coverage(args) => {
            info!(
                input  = %args.input.display(),
                output = %args.output.display(),
                sample_id = %args.sample_id.as_deref().unwrap_or("<from SM tag>"),
                "collecting coverage"
            );

            let target_intervals = args
                .targets
                .as_ref()
                .map(|p| {
                    info!(targets = %p.display(), "loading target intervals");
                    targets::TargetIntervals::load(p)
                })
                .transpose()?;

            if let Some(ref ti) = target_intervals {
                info!(n_targets = ti.n_targets(), total_bases = ti.total_bases(), "target intervals loaded");
            }

            let gene_annots = args
                .gene_annotations
                .as_ref()
                .map(|p| {
                    info!(gene_annotations = %p.display(), "loading gene annotations");
                    gene_annotations::GeneAnnotations::load(p)
                })
                .transpose()?;

            if args.intervals_output.is_some() && target_intervals.is_none() {
                anyhow::bail!("--intervals-output requires --targets");
            }

            let mut cov_writer = writer::parquet_coverage::CoverageWriter::new(&args.output)?;

            let interval_records = coverage::collect_coverage(
                &args,
                target_intervals.as_ref().map(|t| t as &_),
                gene_annots.as_ref(),
                &mut cov_writer,
            )?;

            cov_writer.close()?;

            if let Some(ref intervals_path) = args.intervals_output {
                info!(
                    n_intervals = interval_records.len(),
                    output = %intervals_path.display(),
                    "writing per-interval summary"
                );
                let mut iv_writer =
                    writer::parquet_coverage_intervals::CoverageIntervalsWriter::new(intervals_path)?;
                for record in interval_records {
                    iv_writer.push(record)?;
                }
                iv_writer.close()?;
            }

            info!("done");
        }
    }

    Ok(())
}
