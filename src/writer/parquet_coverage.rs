use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{
    ArrayRef, BooleanArray, Float32Array, Int32Array, Int64Array, StringArray,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::record::CoverageRecord;

/// Write a slice of CoverageRecord to a Parquet file.
pub fn write_parquet(records: &[CoverageRecord], output: &Path) -> Result<()> {
    let schema = coverage_schema();
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

fn coverage_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_id",          DataType::Utf8,    false),
        Field::new("chrom",              DataType::Utf8,    false),
        Field::new("pos",                DataType::Int64,   false),
        Field::new("end",                DataType::Int64,   false),
        // Fragment depth
        Field::new("total_depth",        DataType::Int32,   false),
        Field::new("fwd_depth",          DataType::Int32,   false),
        Field::new("rev_depth",          DataType::Int32,   false),
        // Duplicate metrics
        Field::new("raw_read_depth",     DataType::Int32,   false),
        Field::new("frac_dup",           DataType::Float32, false),
        // Overlap metrics
        Field::new("overlap_depth",      DataType::Int32,   false),
        Field::new("frac_overlap",       DataType::Float32, false),
        // Mappability signals
        Field::new("mean_mapq",          DataType::Float32, false),
        Field::new("frac_mapq0",         DataType::Float32, false),
        Field::new("frac_low_mapq",      DataType::Float32, false),
        // Base quality signals
        Field::new("mean_base_qual",     DataType::Float32, false),
        Field::new("min_base_qual_obs",  DataType::Int32,   false),
        Field::new("max_base_qual_obs",  DataType::Int32,   false),
        Field::new("frac_low_bq",        DataType::Float32, false),
        // Insert size
        Field::new("mean_insert_size",   DataType::Float32, false),
        Field::new("min_insert_size",    DataType::Int32,   false),
        Field::new("max_insert_size",    DataType::Int32,   false),
        Field::new("n_insert_size_obs",  DataType::Int32,   false),
        // GC content
        Field::new("gc_content",         DataType::Float32, false),
        // Annotations (nullable)
        Field::new("on_target",          DataType::Boolean, true),
        Field::new("gene",               DataType::Utf8,    true),
        // Provenance
        Field::new("read_type",          DataType::Utf8,    false),
        Field::new("pipeline",           DataType::Utf8,    false),
        Field::new("batch",              DataType::Utf8,    true),
    ]))
}

fn records_to_batch(records: &[CoverageRecord], schema: Arc<Schema>) -> Result<RecordBatch> {
    let sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.sample_id.as_str()),
    ));
    let chrom: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.chrom.as_str()),
    ));
    let pos: ArrayRef = Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.pos)));
    let end: ArrayRef = Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.end)));

    let total_depth:   ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.total_depth)));
    let fwd_depth:     ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.fwd_depth)));
    let rev_depth:     ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.rev_depth)));

    let raw_read_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.raw_read_depth)));
    let frac_dup:       ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.frac_dup)));

    let overlap_depth:  ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.overlap_depth)));
    let frac_overlap:   ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.frac_overlap)));

    let mean_mapq:     ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.mean_mapq)));
    let frac_mapq0:    ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.frac_mapq0)));
    let frac_low_mapq: ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.frac_low_mapq)));

    let mean_base_qual:    ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.mean_base_qual)));
    let min_base_qual_obs: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.min_base_qual_obs)));
    let max_base_qual_obs: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.max_base_qual_obs)));
    let frac_low_bq:       ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.frac_low_bq)));

    let mean_insert_size:  ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.mean_insert_size)));
    let min_insert_size:   ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.min_insert_size)));
    let max_insert_size:   ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.max_insert_size)));
    let n_insert_size_obs: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.n_insert_size_obs)));

    let gc_content: ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.gc_content)));

    let on_target: ArrayRef = Arc::new(BooleanArray::from(
        records.iter().map(|r| r.on_target).collect::<Vec<_>>(),
    ));
    let gene: ArrayRef = Arc::new(StringArray::from(
        records.iter().map(|r| r.gene.as_deref()).collect::<Vec<_>>(),
    ));

    let read_type: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.read_type.to_string()),
    ));
    let pipeline: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.pipeline.to_string()),
    ));
    let batch: ArrayRef = Arc::new(StringArray::from(
        records.iter().map(|r| r.batch.as_deref()).collect::<Vec<_>>(),
    ));

    RecordBatch::try_new(
        schema,
        vec![
            sample_id, chrom, pos, end,
            total_depth, fwd_depth, rev_depth,
            raw_read_depth, frac_dup,
            overlap_depth, frac_overlap,
            mean_mapq, frac_mapq0, frac_low_mapq,
            mean_base_qual, min_base_qual_obs, max_base_qual_obs, frac_low_bq,
            mean_insert_size, min_insert_size, max_insert_size, n_insert_size_obs,
            gc_content,
            on_target, gene,
            read_type, pipeline, batch,
        ],
    )
    .context("failed to create Arrow record batch")
}
