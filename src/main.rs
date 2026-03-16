mod bam;
mod cli;
mod progress;
mod record;
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

            let vcf_index = args
                .vcf
                .as_ref()
                .map(|p| {
                    info!(vcf = %p.display(), "loading VCF annotations");
                    vcf::VcfIndex::load(p)
                })
                .transpose()?;

            let records = bam::collect_alt_bases(&args, vcf_index.as_ref())?;

            info!(n_records = records.len(), "writing Parquet output");
            writer::parquet::write_parquet(&records, &args.output)?;

            info!("done");
        }

        Command::Merge(args) => {
            // TODO: implement DuckDB merge of per-sample Parquet files
            // For now, print the files that would be merged.
            info!(
                n_inputs = args.inputs.len(),
                output = %args.output.display(),
                "merge not yet implemented"
            );
            anyhow::bail!("merge subcommand not yet implemented");
        }
    }

    Ok(())
}
