use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, BooleanArray, Float32Array, Int32Array, Int64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::record::AltRead;

/// Write a slice of AltRead records to a Parquet file.
pub fn write_parquet(records: &[AltRead], output: &Path) -> Result<()> {
    let schema = alt_read_schema();
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

fn alt_read_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_id", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, false),
        Field::new("pos", DataType::Int64, false),
        Field::new("alt_allele", DataType::Utf8, false),
        Field::new("read_type", DataType::Utf8, false),
        Field::new("pipeline", DataType::Utf8, true),
        Field::new("subject_id", DataType::Utf8, true),
        Field::new("sample_type", DataType::Utf8, true),
        Field::new("batch", DataType::Utf8, true),
        Field::new("label1", DataType::Utf8, true),
        Field::new("label2", DataType::Utf8, true),
        Field::new("label3", DataType::Utf8, true),
        Field::new("timepoint", DataType::Utf8, true),
        Field::new("cycle", DataType::Int32, false),
        Field::new("read_length", DataType::Int32, false),
        Field::new("is_read1", DataType::Boolean, false),
        Field::new("ab_count", DataType::Int32, true),
        Field::new("ba_count", DataType::Int32, true),
        Field::new("family_size", DataType::Int32, true),
        Field::new("base_qual", DataType::Int32, false),
        Field::new("map_qual", DataType::Int32, false),
        Field::new("insert_size", DataType::Int32, true),
        Field::new("frag_gc", DataType::Float32, true),
        Field::new("n_before_alt", DataType::Int32, false),
        Field::new("n_after_alt", DataType::Int32, false),
        Field::new("n_n_before_alt", DataType::Int32, false),
        Field::new("n_n_after_alt", DataType::Int32, false),
        Field::new("leading_n_run_len", DataType::Int32, false),
        Field::new("trailing_n_run_len", DataType::Int32, false),
        Field::new("input_checksum_sha256", DataType::Utf8, true),
        Field::new("fragment_id", DataType::Int64, false),
    ]))
}

fn records_to_batch(records: &[AltRead], schema: Arc<Schema>) -> Result<RecordBatch> {
    let sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.sample_id.as_str()),
    ));
    let chrom: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.chrom.as_str()),
    ));
    let pos: ArrayRef = Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.pos)));
    let alt_allele: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.alt_allele.as_str()),
    ));
    let read_type: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.read_type.to_string()),
    ));
    let pipeline: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.pipeline.clone())
            .collect::<Vec<_>>(),
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
    let label1: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.label1.as_deref())
            .collect::<Vec<_>>(),
    ));
    let label2: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.label2.as_deref())
            .collect::<Vec<_>>(),
    ));
    let label3: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.label3.as_deref())
            .collect::<Vec<_>>(),
    ));
    let timepoint: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.timepoint.as_deref())
            .collect::<Vec<_>>(),
    ));
    let cycle: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.cycle),
    ));
    let read_length: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.read_length),
    ));
    let is_read1: ArrayRef = Arc::new(BooleanArray::from(
        records.iter().map(|r| r.is_read1).collect::<Vec<_>>(),
    ));
    let ab_count: ArrayRef = Arc::new(Int32Array::from(
        records.iter().map(|r| r.ab_count).collect::<Vec<_>>(),
    ));
    let ba_count: ArrayRef = Arc::new(Int32Array::from(
        records.iter().map(|r| r.ba_count).collect::<Vec<_>>(),
    ));
    let family_size: ArrayRef = Arc::new(Int32Array::from(
        records.iter().map(|r| r.family_size).collect::<Vec<_>>(),
    ));
    let base_qual: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.base_qual),
    ));
    let map_qual: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.map_qual),
    ));
    let insert_size: ArrayRef = Arc::new(Int32Array::from(
        records.iter().map(|r| r.insert_size).collect::<Vec<_>>(),
    ));
    let frag_gc: ArrayRef = Arc::new(Float32Array::from(
        records.iter().map(|r| r.frag_gc).collect::<Vec<_>>(),
    ));
    let n_before_alt: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.n_before_alt),
    ));
    let n_after_alt: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.n_after_alt),
    ));
    let n_n_before_alt: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.n_n_before_alt),
    ));
    let n_n_after_alt: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.n_n_after_alt),
    ));
    let leading_n_run_len: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.leading_n_run_len),
    ));
    let trailing_n_run_len: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.trailing_n_run_len),
    ));
    let input_checksum_sha256: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.input_checksum_sha256.as_deref())
            .collect::<Vec<_>>(),
    ));
    let fragment_id: ArrayRef = Arc::new(Int64Array::from_iter_values(
        records.iter().map(|r| r.fragment_id),
    ));

    RecordBatch::try_new(
        schema,
        vec![
            sample_id,
            chrom,
            pos,
            alt_allele,
            read_type,
            pipeline,
            subject_id,
            sample_type,
            batch,
            label1,
            label2,
            label3,
            timepoint,
            cycle,
            read_length,
            is_read1,
            ab_count,
            ba_count,
            family_size,
            base_qual,
            map_qual,
            insert_size,
            frag_gc,
            n_before_alt,
            n_after_alt,
            n_n_before_alt,
            n_n_after_alt,
            leading_n_run_len,
            trailing_n_run_len,
            input_checksum_sha256,
            fragment_id,
        ],
    )
    .context("failed to create Arrow record batch")
}
