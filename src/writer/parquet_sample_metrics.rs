use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, Float32Array, Int32Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::record::SampleMetricsRecord;

pub fn write_parquet(records: &[SampleMetricsRecord], output: &Path) -> Result<()> {
    let schema = sample_metrics_schema();
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

fn sample_metrics_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_id", DataType::Utf8, false),
        Field::new("subject_id", DataType::Utf8, true),
        Field::new("sample_type", DataType::Utf8, true),
        Field::new("batch", DataType::Utf8, true),
        Field::new("read_type", DataType::Utf8, false),
        Field::new("pipeline", DataType::Utf8, true),
        Field::new("input_checksum_sha256", DataType::Utf8, true),
        Field::new("n_target_positions", DataType::Int32, false),
        Field::new("n_target_positions_covered", DataType::Int32, false),
        Field::new("mean_target_depth_covered", DataType::Float32, true),
        Field::new("mean_target_depth_all", DataType::Float32, true),
        Field::new("median_target_depth_covered", DataType::Float32, true),
        Field::new("median_target_depth_all", DataType::Float32, true),
        Field::new("pct_fragment_bases_on_target", DataType::Float32, true),
    ]))
}

fn records_to_batch(records: &[SampleMetricsRecord], schema: Arc<Schema>) -> Result<RecordBatch> {
    let sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.sample_id.as_str()),
    ));
    let subject_id: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.subject_id.as_deref())
            .collect::<Vec<_>>(),
    ));
    let sample_type: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.sample_type.as_deref())
            .collect::<Vec<_>>(),
    ));
    let batch: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.batch.as_deref())
            .collect::<Vec<_>>(),
    ));
    let read_type: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.read_type.to_string()),
    ));
    let pipeline: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.pipeline.map(|p| p.to_string()))
            .collect::<Vec<_>>(),
    ));
    let input_checksum_sha256: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.input_checksum_sha256.as_deref())
            .collect::<Vec<_>>(),
    ));
    let n_target_positions: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.n_target_positions),
    ));
    let n_target_positions_covered: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.n_target_positions_covered),
    ));
    let mean_target_depth_covered: ArrayRef = Arc::new(Float32Array::from(
        records
            .iter()
            .map(|r| r.mean_target_depth_covered)
            .collect::<Vec<_>>(),
    ));
    let mean_target_depth_all: ArrayRef = Arc::new(Float32Array::from(
        records
            .iter()
            .map(|r| r.mean_target_depth_all)
            .collect::<Vec<_>>(),
    ));
    let median_target_depth_covered: ArrayRef = Arc::new(Float32Array::from(
        records
            .iter()
            .map(|r| r.median_target_depth_covered)
            .collect::<Vec<_>>(),
    ));
    let median_target_depth_all: ArrayRef = Arc::new(Float32Array::from(
        records
            .iter()
            .map(|r| r.median_target_depth_all)
            .collect::<Vec<_>>(),
    ));
    let pct_fragment_bases_on_target: ArrayRef = Arc::new(Float32Array::from(
        records
            .iter()
            .map(|r| r.pct_fragment_bases_on_target)
            .collect::<Vec<_>>(),
    ));

    RecordBatch::try_new(
        schema,
        vec![
            sample_id,
            subject_id,
            sample_type,
            batch,
            read_type,
            pipeline,
            input_checksum_sha256,
            n_target_positions,
            n_target_positions_covered,
            mean_target_depth_covered,
            mean_target_depth_all,
            median_target_depth_covered,
            median_target_depth_all,
            pct_fragment_bases_on_target,
        ],
    )
    .context("failed to create Arrow record batch")
}
