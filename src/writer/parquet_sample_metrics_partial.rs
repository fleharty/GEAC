//! Writer for per-shard partial sample metrics (`*.sample_metrics_partial.parquet`).
//!
//! These hold combinable sufficient statistics — scalar sums plus a covered-depth
//! histogram stored as two `List` columns — that `geac aggregate-metrics` merges
//! into a final sample_metrics Parquet. This format is internal to the scatter
//! pipeline and must never be fed to `geac merge`.

use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, Int32Builder, Int64Array, Int64Builder, ListBuilder, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::record::SampleMetricsPartialRecord;

pub fn write_parquet(records: &[SampleMetricsPartialRecord], output: &Path) -> Result<()> {
    let schema = schema();
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

fn schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_id", DataType::Utf8, false),
        Field::new("subject_id", DataType::Utf8, true),
        Field::new("sample_type", DataType::Utf8, true),
        Field::new("batch", DataType::Utf8, true),
        Field::new("read_type", DataType::Utf8, false),
        Field::new("pipeline", DataType::Utf8, true),
        Field::new("input_checksum_sha256", DataType::Utf8, true),
        Field::new("n_target_positions", DataType::Int64, false),
        Field::new("total_fragment_bases", DataType::Int64, false),
        Field::new("on_target_fragment_bases", DataType::Int64, false),
        // List item nullability must match arrow's ListBuilder (nullable items);
        // the values themselves are never null.
        Field::new(
            "hist_depth",
            DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
            false,
        ),
        Field::new(
            "hist_count",
            DataType::List(Arc::new(Field::new("item", DataType::Int64, true))),
            false,
        ),
    ]))
}

fn records_to_batch(
    records: &[SampleMetricsPartialRecord],
    schema: Arc<Schema>,
) -> Result<RecordBatch> {
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
    let batch_col: ArrayRef = Arc::new(StringArray::from(
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
            .map(|r| r.pipeline.as_deref())
            .collect::<Vec<_>>(),
    ));
    let input_checksum_sha256: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.input_checksum_sha256.as_deref())
            .collect::<Vec<_>>(),
    ));
    let n_target_positions: ArrayRef = Arc::new(Int64Array::from_iter_values(
        records.iter().map(|r| r.n_target_positions),
    ));
    let total_fragment_bases: ArrayRef = Arc::new(Int64Array::from_iter_values(
        records.iter().map(|r| r.total_fragment_bases),
    ));
    let on_target_fragment_bases: ArrayRef = Arc::new(Int64Array::from_iter_values(
        records.iter().map(|r| r.on_target_fragment_bases),
    ));

    let mut depth_builder = ListBuilder::new(Int32Builder::new());
    for r in records {
        for &d in &r.hist_depth {
            depth_builder.values().append_value(d);
        }
        depth_builder.append(true);
    }
    let hist_depth: ArrayRef = Arc::new(depth_builder.finish());

    let mut count_builder = ListBuilder::new(Int64Builder::new());
    for r in records {
        for &c in &r.hist_count {
            count_builder.values().append_value(c);
        }
        count_builder.append(true);
    }
    let hist_count: ArrayRef = Arc::new(count_builder.finish());

    RecordBatch::try_new(
        schema,
        vec![
            sample_id,
            subject_id,
            sample_type,
            batch_col,
            read_type,
            pipeline,
            input_checksum_sha256,
            n_target_positions,
            total_fragment_bases,
            on_target_fragment_bases,
            hist_depth,
            hist_count,
        ],
    )
    .context("failed to create Arrow record batch")
}
