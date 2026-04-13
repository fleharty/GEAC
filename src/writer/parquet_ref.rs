use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, BooleanArray, Float32Array, Int32Array, Int64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::record::RefBase;

pub fn write_parquet(records: &[RefBase], output: &Path) -> Result<()> {
    let schema = ref_base_schema();
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

fn ref_base_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_id", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, false),
        Field::new("pos", DataType::Int64, false),
        Field::new("ref_allele", DataType::Utf8, false),
        Field::new("total_depth", DataType::Int32, false),
        Field::new("fwd_depth", DataType::Int32, false),
        Field::new("rev_depth", DataType::Int32, false),
        Field::new("ref_count", DataType::Int32, false),
        Field::new("fwd_ref_count", DataType::Int32, false),
        Field::new("rev_ref_count", DataType::Int32, false),
        Field::new("overlap_depth", DataType::Int32, false),
        Field::new("overlap_ref_agree", DataType::Int32, false),
        Field::new("read_type", DataType::Utf8, false),
        Field::new("pipeline", DataType::Utf8, false),
        Field::new("batch", DataType::Utf8, true),
        Field::new("label1", DataType::Utf8, true),
        Field::new("label2", DataType::Utf8, true),
        Field::new("label3", DataType::Utf8, true),
        Field::new("on_target", DataType::Boolean, true),
        Field::new("gene", DataType::Utf8, true),
        Field::new("homopolymer_len", DataType::Int32, true),
        Field::new("str_period", DataType::Int32, true),
        Field::new("str_len", DataType::Int32, true),
        Field::new("gnomad_af", DataType::Float32, true),
        Field::new("input_checksum_sha256", DataType::Utf8, true),
    ]))
}

fn records_to_batch(records: &[RefBase], schema: Arc<Schema>) -> Result<RecordBatch> {
    let sample_id: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.sample_id.as_str()),
    ));
    let chrom: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.chrom.as_str()),
    ));
    let pos: ArrayRef = Arc::new(Int64Array::from_iter_values(records.iter().map(|r| r.pos)));
    let ref_allele: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.ref_allele.as_str()),
    ));
    let total_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.total_depth),
    ));
    let fwd_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.fwd_depth),
    ));
    let rev_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.rev_depth),
    ));
    let ref_count: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.ref_count),
    ));
    let fwd_ref_count: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.fwd_ref_count),
    ));
    let rev_ref_count: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.rev_ref_count),
    ));
    let overlap_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.overlap_depth),
    ));
    let overlap_ref_agree: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.overlap_ref_agree),
    ));
    let read_type: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.read_type.to_string()),
    ));
    let pipeline: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.pipeline.to_string()),
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
    let on_target: ArrayRef = Arc::new(BooleanArray::from(
        records.iter().map(|r| r.on_target).collect::<Vec<_>>(),
    ));
    let gene: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.gene.as_deref())
            .collect::<Vec<_>>(),
    ));
    let homopolymer_len: ArrayRef = Arc::new(Int32Array::from(
        records
            .iter()
            .map(|r| r.homopolymer_len)
            .collect::<Vec<_>>(),
    ));
    let str_period: ArrayRef = Arc::new(Int32Array::from(
        records.iter().map(|r| r.str_period).collect::<Vec<_>>(),
    ));
    let str_len: ArrayRef = Arc::new(Int32Array::from(
        records.iter().map(|r| r.str_len).collect::<Vec<_>>(),
    ));
    let gnomad_af: ArrayRef = Arc::new(Float32Array::from(
        records.iter().map(|r| r.gnomad_af).collect::<Vec<_>>(),
    ));
    let input_checksum_sha256: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.input_checksum_sha256.as_deref())
            .collect::<Vec<_>>(),
    ));

    RecordBatch::try_new(
        schema,
        vec![
            sample_id,
            chrom,
            pos,
            ref_allele,
            total_depth,
            fwd_depth,
            rev_depth,
            ref_count,
            fwd_ref_count,
            rev_ref_count,
            overlap_depth,
            overlap_ref_agree,
            read_type,
            pipeline,
            batch,
            label1,
            label2,
            label3,
            on_target,
            gene,
            homopolymer_len,
            str_period,
            str_len,
            gnomad_af,
            input_checksum_sha256,
        ],
    )
    .context("failed to create Arrow record batch")
}
