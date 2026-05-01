use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, Float32Array, Int32Array, Int64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::record::FragmentRecord;

const BATCH_SIZE: usize = 100_000;

pub struct FragmentsWriter {
    writer: ArrowWriter<std::fs::File>,
    schema: Arc<Schema>,
    buf: Vec<FragmentRecord>,
}

impl FragmentsWriter {
    pub fn new(output: &Path) -> Result<Self> {
        let schema = fragments_schema();
        let file = std::fs::File::create(output)
            .with_context(|| format!("failed to create output file: {}", output.display()))?;
        let props = WriterProperties::builder()
            .set_compression(Compression::SNAPPY)
            .build();
        let writer = ArrowWriter::try_new(file, Arc::clone(&schema), Some(props))
            .context("failed to create Parquet writer")?;
        Ok(Self {
            writer,
            schema,
            buf: Vec::with_capacity(BATCH_SIZE),
        })
    }

    pub fn push(&mut self, record: FragmentRecord) -> Result<()> {
        self.buf.push(record);
        if self.buf.len() >= BATCH_SIZE {
            self.flush()?;
        }
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        let batch = records_to_batch(&self.buf, Arc::clone(&self.schema))?;
        self.writer
            .write(&batch)
            .context("failed to write record batch")?;
        self.buf.clear();
        Ok(())
    }

    pub fn close(mut self) -> Result<()> {
        self.flush()?;
        self.writer
            .close()
            .context("failed to finalize Parquet file")?;
        Ok(())
    }
}

fn fragments_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_id", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, false),
        Field::new("frag_start", DataType::Int64, false),
        Field::new("frag_end", DataType::Int64, false),
        Field::new("midpoint", DataType::Int64, false),
        Field::new("insert_size", DataType::Int32, false),
        Field::new("map_qual", DataType::Int32, false),
        Field::new("read_type", DataType::Utf8, false),
        Field::new("pipeline", DataType::Utf8, false),
        Field::new("gc_content", DataType::Float32, true),
        Field::new("end_motif_5p", DataType::Utf8, true),
        Field::new("end_motif_3p", DataType::Utf8, true),
        Field::new("subject_id", DataType::Utf8, true),
        Field::new("sample_type", DataType::Utf8, true),
        Field::new("batch", DataType::Utf8, true),
        Field::new("label1", DataType::Utf8, true),
        Field::new("label2", DataType::Utf8, true),
        Field::new("label3", DataType::Utf8, true),
        Field::new("timepoint", DataType::Utf8, true),
    ]))
}

fn records_to_batch(records: &[FragmentRecord], schema: Arc<Schema>) -> Result<RecordBatch> {
    let sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.sample_id.as_str()),
    ));
    let chrom: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.chrom.as_str()),
    ));
    let frag_start: ArrayRef = Arc::new(Int64Array::from_iter_values(
        records.iter().map(|r| r.frag_start),
    ));
    let frag_end: ArrayRef = Arc::new(Int64Array::from_iter_values(
        records.iter().map(|r| r.frag_end),
    ));
    let midpoint: ArrayRef = Arc::new(Int64Array::from_iter_values(
        records.iter().map(|r| r.midpoint),
    ));
    let insert_size: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.insert_size),
    ));
    let gc_content: ArrayRef = Arc::new(Float32Array::from(
        records.iter().map(|r| r.gc_content).collect::<Vec<_>>(),
    ));
    let end_motif_5p: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.end_motif_5p.as_deref())
            .collect::<Vec<_>>(),
    ));
    let end_motif_3p: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.end_motif_3p.as_deref())
            .collect::<Vec<_>>(),
    ));
    let map_qual: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.map_qual),
    ));
    let read_type: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.read_type.to_string()),
    ));
    let pipeline: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.pipeline.to_string()),
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
        records.iter().map(|r| r.batch.as_deref()).collect::<Vec<_>>(),
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

    RecordBatch::try_new(
        schema,
        vec![
            sample_id, chrom, frag_start, frag_end, midpoint, insert_size, map_qual, read_type,
            pipeline, gc_content, end_motif_5p, end_motif_3p, subject_id, sample_type,
            batch, label1, label2, label3, timepoint,
        ],
    )
    .context("failed to create Arrow record batch")
}
