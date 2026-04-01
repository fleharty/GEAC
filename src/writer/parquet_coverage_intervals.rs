use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, Float32Array, Int32Array, Int64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::record::IntervalRecord;

const BATCH_SIZE: usize = 10_000;

/// Streaming Parquet writer for IntervalRecord.
pub struct CoverageIntervalsWriter {
    writer: ArrowWriter<std::fs::File>,
    schema: Arc<Schema>,
    buf:    Vec<IntervalRecord>,
}

impl CoverageIntervalsWriter {
    pub fn new(output: &Path) -> Result<Self> {
        let schema = intervals_schema();
        let file = std::fs::File::create(output)
            .with_context(|| format!("failed to create output file: {}", output.display()))?;
        let props = WriterProperties::builder()
            .set_compression(Compression::SNAPPY)
            .build();
        let writer = ArrowWriter::try_new(file, Arc::clone(&schema), Some(props))
            .context("failed to create Parquet writer")?;
        Ok(Self { writer, schema, buf: Vec::with_capacity(BATCH_SIZE) })
    }

    pub fn push(&mut self, record: IntervalRecord) -> Result<()> {
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
        self.writer.write(&batch).context("failed to write record batch")?;
        self.buf.clear();
        Ok(())
    }

    pub fn close(mut self) -> Result<()> {
        self.flush()?;
        self.writer.close().context("failed to finalize Parquet file")?;
        Ok(())
    }
}

fn intervals_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_id",          DataType::Utf8,    false),
        Field::new("chrom",              DataType::Utf8,    false),
        Field::new("start",              DataType::Int64,   false),
        Field::new("end",                DataType::Int64,   false),
        Field::new("interval_name",      DataType::Utf8,    true),
        Field::new("gene",               DataType::Utf8,    true),
        // Depth summary
        Field::new("n_bases",            DataType::Int32,   false),
        Field::new("mean_depth",         DataType::Float32, false),
        Field::new("median_depth",       DataType::Float32, false),
        Field::new("min_depth",          DataType::Int32,   false),
        Field::new("max_depth",          DataType::Int32,   false),
        Field::new("frac_at_1x",         DataType::Float32, false),
        Field::new("frac_at_10x",        DataType::Float32, false),
        Field::new("frac_at_20x",        DataType::Float32, false),
        Field::new("frac_at_30x",        DataType::Float32, false),
        Field::new("frac_at_50x",        DataType::Float32, false),
        Field::new("frac_at_100x",       DataType::Float32, false),
        // Aggregated QC signals
        Field::new("mean_gc_content",    DataType::Float32, false),
        Field::new("mean_mapq",          DataType::Float32, false),
        Field::new("mean_frac_mapq0",    DataType::Float32, false),
        Field::new("mean_frac_dup",      DataType::Float32, false),
        Field::new("mean_frac_overlap",  DataType::Float32, false),
        Field::new("mean_base_qual",     DataType::Float32, false),
        Field::new("mean_insert_size",   DataType::Float32, false),
        // Provenance
        Field::new("read_type",          DataType::Utf8,    false),
        Field::new("pipeline",           DataType::Utf8,    false),
        Field::new("batch",              DataType::Utf8,    true),
    ]))
}

fn records_to_batch(records: &[IntervalRecord], schema: Arc<Schema>) -> Result<RecordBatch> {
    let sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.sample_id.as_str()),
    ));
    let chrom: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.chrom.as_str()),
    ));
    let start: ArrayRef = Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.start)));
    let end:   ArrayRef = Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.end)));

    let interval_name: ArrayRef = Arc::new(StringArray::from(
        records.iter().map(|r| r.interval_name.as_deref()).collect::<Vec<_>>(),
    ));
    let gene: ArrayRef = Arc::new(StringArray::from(
        records.iter().map(|r| r.gene.as_deref()).collect::<Vec<_>>(),
    ));

    let n_bases:      ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.n_bases)));
    let mean_depth:   ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.mean_depth)));
    let median_depth: ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.median_depth)));
    let min_depth:    ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.min_depth)));
    let max_depth:    ArrayRef = Arc::new(Int32Array::from_iter_values(records.iter().map(|r| r.max_depth)));
    let frac_at_1x:   ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.frac_at_1x)));
    let frac_at_10x:  ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.frac_at_10x)));
    let frac_at_20x:  ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.frac_at_20x)));
    let frac_at_30x:  ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.frac_at_30x)));
    let frac_at_50x:  ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.frac_at_50x)));
    let frac_at_100x: ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.frac_at_100x)));

    let mean_gc_content:   ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.mean_gc_content)));
    let mean_mapq:         ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.mean_mapq)));
    let mean_frac_mapq0:   ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.mean_frac_mapq0)));
    let mean_frac_dup:     ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.mean_frac_dup)));
    let mean_frac_overlap: ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.mean_frac_overlap)));
    let mean_base_qual:    ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.mean_base_qual)));
    let mean_insert_size:  ArrayRef = Arc::new(Float32Array::from_iter_values(records.iter().map(|r| r.mean_insert_size)));

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
            sample_id, chrom, start, end, interval_name, gene,
            n_bases, mean_depth, median_depth, min_depth, max_depth,
            frac_at_1x, frac_at_10x, frac_at_20x, frac_at_30x, frac_at_50x, frac_at_100x,
            mean_gc_content, mean_mapq, mean_frac_mapq0, mean_frac_dup,
            mean_frac_overlap, mean_base_qual, mean_insert_size,
            read_type, pipeline, batch,
        ],
    )
    .context("failed to create Arrow record batch")
}
