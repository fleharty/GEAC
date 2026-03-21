use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, Int32Array, Int64Array, StringArray};
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

    writer.write(&batch).context("failed to write record batch")?;
    writer.close().context("failed to finalize Parquet file")?;

    Ok(())
}

fn alt_read_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_id",            DataType::Utf8,  false),
        Field::new("chrom",                DataType::Utf8,  false),
        Field::new("pos",                  DataType::Int64, false),
        Field::new("alt_allele",           DataType::Utf8,  false),
        Field::new("dist_from_read_start", DataType::Int32, false),
        Field::new("dist_from_read_end",   DataType::Int32, false),
        Field::new("read_length",          DataType::Int32, false),
        Field::new("ab_count",             DataType::Int32, true),
        Field::new("ba_count",             DataType::Int32, true),
        Field::new("family_size",          DataType::Int32, true),
        Field::new("base_qual",            DataType::Int32, false),
        Field::new("map_qual",             DataType::Int32, false),
    ]))
}

fn records_to_batch(records: &[AltRead], schema: Arc<Schema>) -> Result<RecordBatch> {
    let sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(records.iter().map(|r| r.sample_id.as_str())));
    let chrom:     ArrayRef = Arc::new(StringArray::from_iter_values(records.iter().map(|r| r.chrom.as_str())));
    let pos:       ArrayRef = Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.pos)));
    let alt_allele: ArrayRef = Arc::new(StringArray::from_iter_values(records.iter().map(|r| r.alt_allele.as_str())));
    let dist_from_read_start: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.dist_from_read_start)));
    let dist_from_read_end:   ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.dist_from_read_end)));
    let read_length:          ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.read_length)));
    let ab_count:    ArrayRef = Arc::new(Int32Array::from(records.iter().map(|r| r.ab_count).collect::<Vec<_>>()));
    let ba_count:    ArrayRef = Arc::new(Int32Array::from(records.iter().map(|r| r.ba_count).collect::<Vec<_>>()));
    let family_size: ArrayRef = Arc::new(Int32Array::from(records.iter().map(|r| r.family_size).collect::<Vec<_>>()));
    let base_qual: ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.base_qual)));
    let map_qual:  ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.map_qual)));

    RecordBatch::try_new(
        schema,
        vec![
            sample_id, chrom, pos, alt_allele,
            dist_from_read_start, dist_from_read_end, read_length,
            ab_count, ba_count, family_size,
            base_qual, map_qual,
        ],
    )
    .context("failed to create Arrow record batch")
}
