use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, Int32Array, Int64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::record::NormalEvidence;

/// Write a slice of NormalEvidence records to a Parquet file.
pub fn write_parquet(records: &[NormalEvidence], output: &Path) -> Result<()> {
    let schema = normal_evidence_schema();
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

fn normal_evidence_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("tumor_sample_id",   DataType::Utf8,  false),
        Field::new("chrom",             DataType::Utf8,  false),
        Field::new("pos",               DataType::Int64, false),
        Field::new("tumor_alt_allele",  DataType::Utf8,  false),
        Field::new("normal_sample_id",  DataType::Utf8,  false),
        Field::new("normal_alt_allele", DataType::Utf8,  true),  // nullable
        Field::new("normal_depth",      DataType::Int32, false),
        Field::new("normal_alt_count",  DataType::Int32, false),
    ]))
}

fn records_to_batch(records: &[NormalEvidence], schema: Arc<Schema>) -> Result<RecordBatch> {
    let tumor_sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.tumor_sample_id.as_str()),
    ));
    let chrom: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.chrom.as_str()),
    ));
    let pos: ArrayRef = Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.pos)));
    let tumor_alt_allele: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.tumor_alt_allele.as_str()),
    ));
    let normal_sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.normal_sample_id.as_str()),
    ));
    let normal_alt_allele: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.normal_alt_allele.as_deref())
            .collect::<Vec<_>>(),
    ));
    let normal_depth: ArrayRef =
        Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.normal_depth)));
    let normal_alt_count: ArrayRef =
        Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.normal_alt_count)));

    RecordBatch::try_new(
        schema,
        vec![
            tumor_sample_id,
            chrom,
            pos,
            tumor_alt_allele,
            normal_sample_id,
            normal_alt_allele,
            normal_depth,
            normal_alt_count,
        ],
    )
    .context("failed to create Arrow record batch")
}
