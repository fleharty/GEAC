mod bam;
mod cli;
mod merge;
mod progress;
mod record;
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

            let records = bam::collect_alt_bases(&args, annotator)?;

            info!(n_records = records.len(), "writing Parquet output");
            writer::parquet::write_parquet(&records, &args.output)?;

            info!("done");
        }

        Command::Merge(args) => {
            merge::merge(&args)?;
        }
    }

    Ok(())
}
