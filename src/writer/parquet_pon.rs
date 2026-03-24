use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, Float64Array, Int64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::record::PonEvidence;

/// Write a slice of PonEvidence records to a Parquet file.
pub fn write_parquet(records: &[PonEvidence], output: &Path) -> Result<()> {
    let schema = pon_evidence_schema();
    let batch = records_to_batch(records, Arc::clone(&schema))?;

    let file = std::fs::File::create(output)
        .with_context(|| format!("failed to create output file: {}", output.display()))?;

    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();

    let mut writer = ArrowWriter::try_new(file, Arc::clone(&schema), Some(props))
        .context("failed to create Parquet writer")?;

    writer.write(&batch).context("failed to write record batch")?;
    writer.close().context("failed to finalize Parquet file")?;

    Ok(())
}

fn pon_evidence_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("tumor_sample_id",   DataType::Utf8,    false),
        Field::new("chrom",             DataType::Utf8,    false),
        Field::new("pos",               DataType::Int64,   false),
        Field::new("tumor_alt_allele",  DataType::Utf8,    false),
        Field::new("n_pon_samples",     DataType::Int64,   false),
        Field::new("pon_total_samples", DataType::Int64,   false),
        Field::new("max_pon_vaf",       DataType::Float64, true),  // nullable
        Field::new("mean_pon_vaf",      DataType::Float64, true),  // nullable
    ]))
}

fn records_to_batch(records: &[PonEvidence], schema: Arc<Schema>) -> Result<RecordBatch> {
    let tumor_sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.tumor_sample_id.as_str()),
    ));
    let chrom: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.chrom.as_str()),
    ));
    let pos: ArrayRef =
        Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.pos)));
    let tumor_alt_allele: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.tumor_alt_allele.as_str()),
    ));
    let n_pon_samples: ArrayRef =
        Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.n_pon_samples)));
    let pon_total_samples: ArrayRef =
        Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.pon_total_samples)));
    let max_pon_vaf: ArrayRef = Arc::new(Float64Array::from(
        records.iter().map(|r| r.max_pon_vaf).collect::<Vec<_>>(),
    ));
    let mean_pon_vaf: ArrayRef = Arc::new(Float64Array::from(
        records.iter().map(|r| r.mean_pon_vaf).collect::<Vec<_>>(),
    ));

    RecordBatch::try_new(
        schema,
        vec![
            tumor_sample_id,
            chrom,
            pos,
            tumor_alt_allele,
            n_pon_samples,
            pon_total_samples,
            max_pon_vaf,
            mean_pon_vaf,
        ],
    )
    .context("failed to create Arrow record batch")
}
