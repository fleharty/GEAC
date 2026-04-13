use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, Int32Array, Int64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::record::LocusDepthRecord;

/// Write a slice of LocusDepthRecord records to a Parquet file.
pub fn write_parquet(records: &[LocusDepthRecord], output: &Path) -> Result<()> {
    let schema = locus_depth_schema();
    let batch = records_to_batch(records, Arc::clone(&schema))?;

    let file = std::fs::File::create(output)
        .with_context(|| format!("failed to create output file: {}", output.display()))?;

    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();

    let mut writer = ArrowWriter::try_new(file, Arc::clone(&schema), Some(props))
        .context("failed to create Parquet writer")?;

    writer
        .write(&batch)
        .context("failed to write record batch")?;
    writer.close().context("failed to finalize Parquet file")?;

    Ok(())
}

fn locus_depth_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_id", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, false),
        Field::new("pos", DataType::Int64, false),
        Field::new("total_depth", DataType::Int32, false),
        Field::new("fwd_depth", DataType::Int32, false),
        Field::new("rev_depth", DataType::Int32, false),
    ]))
}

fn records_to_batch(records: &[LocusDepthRecord], schema: Arc<Schema>) -> Result<RecordBatch> {
    let sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.sample_id.as_str()),
    ));
    let chrom: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.chrom.as_str()),
    ));
    let pos: ArrayRef = Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.pos)));
    let total_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.total_depth),
    ));
    let fwd_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.fwd_depth),
    ));
    let rev_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.rev_depth),
    ));

    RecordBatch::try_new(
        schema,
        vec![sample_id, chrom, pos, total_depth, fwd_depth, rev_depth],
    )
    .context("failed to create Arrow record batch")
}
