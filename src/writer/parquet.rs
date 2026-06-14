use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, BooleanArray, Float32Array, Int32Array, Int64Array, StringArray};
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

    writer
        .write(&batch)
        .context("failed to write record batch")?;
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
        Field::new("ref_count", DataType::Int32, false),
        Field::new("fwd_depth", DataType::Int32, false),
        Field::new("rev_depth", DataType::Int32, false),
        Field::new("fwd_alt_count", DataType::Int32, false),
        Field::new("rev_alt_count", DataType::Int32, false),
        Field::new("fwd_ref_count", DataType::Int32, false),
        Field::new("rev_ref_count", DataType::Int32, false),
        Field::new("overlap_depth", DataType::Int32, false),
        Field::new("overlap_alt_agree", DataType::Int32, false),
        Field::new("overlap_alt_disagree", DataType::Int32, false),
        Field::new("overlap_ref_agree", DataType::Int32, false),
        Field::new("read_type", DataType::Utf8, false),
        Field::new("pipeline", DataType::Utf8, true),
        Field::new("subject_id", DataType::Utf8, true),
        Field::new("sample_type", DataType::Utf8, true),
        Field::new("batch", DataType::Utf8, true),
        Field::new("label1", DataType::Utf8, true),
        Field::new("label2", DataType::Utf8, true),
        Field::new("label3", DataType::Utf8, true),
        Field::new("timepoint", DataType::Utf8, true),
        Field::new("input_checksum_sha256", DataType::Utf8, true),
        Field::new("variant_called", DataType::Boolean, true),
        Field::new("variant_filter", DataType::Utf8, true),
        Field::new("on_target", DataType::Boolean, true),
        Field::new("gene", DataType::Utf8, true),
        Field::new("homopolymer_len", DataType::Int32, false),
        Field::new("str_period", DataType::Int32, false),
        Field::new("str_len", DataType::Int32, false),
        Field::new("trinuc_context", DataType::Utf8, true),
        Field::new("gnomad_af", DataType::Float32, true),
        Field::new("n_alt_reads_with_n_ctx", DataType::Int32, true),
        Field::new("mean_frac_n_before", DataType::Float32, true),
        Field::new("mean_frac_n_after", DataType::Float32, true),
        Field::new("mean_delta_n_frac", DataType::Float32, true),
        Field::new("frac_reads_asymmetric", DataType::Float32, true),
        Field::new("bam_path", DataType::Utf8, true),
        Field::new("bai_path", DataType::Utf8, true),
        Field::new("variants_path", DataType::Utf8, true),
        Field::new("gnomad_path", DataType::Utf8, true),
        Field::new("targets_path", DataType::Utf8, true),
    ]))
}

fn records_to_batch(records: &[AltBase], schema: Arc<Schema>) -> Result<RecordBatch> {
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
    let alt_allele: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.alt_allele.as_str()),
    ));
    let variant_type: ArrayRef = Arc::new(StringArray::from_iter_values(
        records.iter().map(|r| r.variant_type.to_string()),
    ));
    let total_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.total_depth),
    ));
    let alt_count: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.alt_count),
    ));
    let ref_count: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.ref_count),
    ));
    let fwd_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.fwd_depth),
    ));
    let rev_depth: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.rev_depth),
    ));
    let fwd_alt_count: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.fwd_alt_count),
    ));
    let rev_alt_count: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.rev_alt_count),
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
    let overlap_alt_agree: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.overlap_alt_agree),
    ));
    let overlap_alt_disagree: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.overlap_alt_disagree),
    ));
    let overlap_ref_agree: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.overlap_ref_agree),
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
    let input_checksum_sha256: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.input_checksum_sha256.as_deref())
            .collect::<Vec<_>>(),
    ));
    let variant_called: ArrayRef = Arc::new(BooleanArray::from(
        records.iter().map(|r| r.variant_called).collect::<Vec<_>>(),
    ));
    let variant_filter: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.variant_filter.as_deref())
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
    let homopolymer_len: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.homopolymer_len),
    ));
    let str_period: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.str_period),
    ));
    let str_len: ArrayRef = Arc::new(Int32Array::from_iter_values(
        records.iter().map(|r| r.str_len),
    ));
    let trinuc_context: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.trinuc_context.as_deref())
            .collect::<Vec<_>>(),
    ));
    let gnomad_af: ArrayRef = Arc::new(Float32Array::from(
        records.iter().map(|r| r.gnomad_af).collect::<Vec<_>>(),
    ));
    let n_alt_reads_with_n_ctx: ArrayRef = Arc::new(Int32Array::from(
        records
            .iter()
            .map(|r| r.n_alt_reads_with_n_ctx)
            .collect::<Vec<_>>(),
    ));
    let mean_frac_n_before: ArrayRef = Arc::new(Float32Array::from(
        records
            .iter()
            .map(|r| r.mean_frac_n_before)
            .collect::<Vec<_>>(),
    ));
    let mean_frac_n_after: ArrayRef = Arc::new(Float32Array::from(
        records
            .iter()
            .map(|r| r.mean_frac_n_after)
            .collect::<Vec<_>>(),
    ));
    let mean_delta_n_frac: ArrayRef = Arc::new(Float32Array::from(
        records
            .iter()
            .map(|r| r.mean_delta_n_frac)
            .collect::<Vec<_>>(),
    ));
    let frac_reads_asymmetric: ArrayRef = Arc::new(Float32Array::from(
        records
            .iter()
            .map(|r| r.frac_reads_asymmetric)
            .collect::<Vec<_>>(),
    ));
    let bam_path: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.bam_path.as_deref())
            .collect::<Vec<_>>(),
    ));
    let bai_path: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.bai_path.as_deref())
            .collect::<Vec<_>>(),
    ));
    let variants_path: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.variants_path.as_deref())
            .collect::<Vec<_>>(),
    ));
    let gnomad_path: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.gnomad_path.as_deref())
            .collect::<Vec<_>>(),
    ));
    let targets_path: ArrayRef = Arc::new(StringArray::from(
        records
            .iter()
            .map(|r| r.targets_path.as_deref())
            .collect::<Vec<_>>(),
    ));

    RecordBatch::try_new(
        schema,
        vec![
            sample_id,
            chrom,
            pos,
            ref_allele,
            alt_allele,
            variant_type,
            total_depth,
            alt_count,
            ref_count,
            fwd_depth,
            rev_depth,
            fwd_alt_count,
            rev_alt_count,
            fwd_ref_count,
            rev_ref_count,
            overlap_depth,
            overlap_alt_agree,
            overlap_alt_disagree,
            overlap_ref_agree,
            read_type,
            pipeline,
            subject_id,
            sample_type,
            batch,
            label1,
            label2,
            label3,
            timepoint,
            input_checksum_sha256,
            variant_called,
            variant_filter,
            on_target,
            gene,
            homopolymer_len,
            str_period,
            str_len,
            trinuc_context,
            gnomad_af,
            n_alt_reads_with_n_ctx,
            mean_frac_n_before,
            mean_frac_n_after,
            mean_delta_n_frac,
            frac_reads_asymmetric,
            bam_path,
            bai_path,
            variants_path,
            gnomad_path,
            targets_path,
        ],
    )
    .context("failed to create Arrow record batch")
}
