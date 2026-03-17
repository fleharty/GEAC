mod bam;
mod cli;
mod cohort;
mod gene_annotations;
mod merge;
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

            let records = bam::collect_alt_bases(
                &args,
                annotator,
                target_intervals.as_ref().map(|t| t as &_),
                gene_annots.as_ref(),
            )?;

            info!(n_records = records.len(), "writing Parquet output");
            writer::parquet::write_parquet(&records, &args.output)?;

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
    }

    Ok(())
}
