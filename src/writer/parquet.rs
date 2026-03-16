use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{
    ArrayRef, BooleanArray, Int32Array, Int64Array, StringArray,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::record::AltBase;

/// Write a slice of AltBase records to a Parquet file.
pub fn write_parquet(records: &[AltBase], output: &Path) -> Result<()> {
    let schema = alt_base_schema();
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

fn alt_base_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_id", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, false),
        Field::new("pos", DataType::Int64, false),
        Field::new("ref_allele", DataType::Utf8, false),
        Field::new("alt_allele", DataType::Utf8, false),
        Field::new("variant_type", DataType::Utf8, false),
        Field::new("total_depth", DataType::Int32, false),
        Field::new("alt_count", DataType::Int32, false),
        Field::new("fwd_depth", DataType::Int32, false),
        Field::new("rev_depth", DataType::Int32, false),
        Field::new("fwd_alt_count", DataType::Int32, false),
        Field::new("rev_alt_count", DataType::Int32, false),
        Field::new("overlap_depth", DataType::Int32, false),
        Field::new("overlap_alt_agree", DataType::Int32, false),
        Field::new("overlap_alt_disagree", DataType::Int32, false),
        Field::new("read_type", DataType::Utf8, false),
        Field::new("pipeline", DataType::Utf8, false),
        Field::new("variant_called", DataType::Boolean, true),
        Field::new("variant_filter", DataType::Utf8, true),
    ]))
}

fn records_to_batch(records: &[AltBase], schema: Arc<Schema>) -> Result<RecordBatch> {
    let sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(records.iter().map(|r| r.sample_id.as_str())));
    let chrom: ArrayRef = Arc::new(StringArray::from_iter_values(records.iter().map(|r| r.chrom.as_str())));
    let pos: ArrayRef = Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.pos)));
    let ref_allele: ArrayRef = Arc::new(StringArray::from_iter_values(records.iter().map(|r| r.ref_allele.as_str())));
    let alt_allele: ArrayRef = Arc::new(StringArray::from_iter_values(records.iter().map(|r| r.alt_allele.as_str())));
    let variant_type: ArrayRef = Arc::new(StringArray::from_iter_values(records.iter().map(|r| r.variant_type.to_string())));
    let total_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.total_depth)));
    let alt_count: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.alt_count)));
    let fwd_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.fwd_depth)));
    let rev_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.rev_depth)));
    let fwd_alt_count: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.fwd_alt_count)));
    let rev_alt_count: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.rev_alt_count)));
    let overlap_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.overlap_depth)));
    let overlap_alt_agree: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.overlap_alt_agree)));
    let overlap_alt_disagree: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.overlap_alt_disagree)));
    let read_type: ArrayRef = Arc::new(StringArray::from_iter_values(records.iter().map(|r| r.read_type.to_string())));
    let pipeline: ArrayRef = Arc::new(StringArray::from_iter_values(records.iter().map(|r| r.pipeline.to_string())));
    let variant_called: ArrayRef = Arc::new(BooleanArray::from(
        records.iter().map(|r| r.variant_called).collect::<Vec<_>>(),
    ));
    let variant_filter: ArrayRef = Arc::new(StringArray::from(
        records.iter().map(|r| r.variant_filter.as_deref()).collect::<Vec<_>>(),
    ));

    RecordBatch::try_new(
        schema,
        vec![
            sample_id, chrom, pos, ref_allele, alt_allele, variant_type,
            total_depth, alt_count, fwd_depth, rev_depth, fwd_alt_count, rev_alt_count,
            overlap_depth, overlap_alt_agree, overlap_alt_disagree,
            read_type, pipeline, variant_called, variant_filter,
        ],
    )
    .context("failed to create Arrow record batch")
}
