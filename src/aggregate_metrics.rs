//! `geac aggregate-metrics` — combine per-shard partial sample-metrics Parquets
//! into one final sample_metrics Parquet.
//!
//! Each input row carries combinable sufficient statistics for one interval shard
//! (see [`crate::record::SampleMetricsPartialRecord`]). Rows are grouped by sample
//! identity, their scalars summed and depth histograms merged, then the final
//! [`SampleMetricsRecord`] is reconstructed via [`crate::sample_metrics_math`].
//! Because the shards are disjoint and the histogram is lossless, the result is
//! identical to collecting the whole sample in one pass — medians included.

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, Result};
use arrow::array::{Array, Int32Array, Int64Array, ListArray, StringArray};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use crate::cli::AggregateMetricsArgs;
use crate::record::{ReadType, SampleMetricsPartialRecord, SampleMetricsRecord};
use crate::sample_metrics_math;
use crate::writer;

pub fn run(args: &AggregateMetricsArgs) -> Result<()> {
    let mut partials = Vec::new();
    for path in &args.inputs {
        partials.extend(read_partial_parquet(path)?);
    }
    let records = aggregate(&partials)?;
    writer::parquet_sample_metrics::write_parquet(&records, &args.output)?;
    Ok(())
}

/// Group partials by sample identity and reconstruct one final record per group.
/// Groups are emitted in deterministic (sorted) identity order.
pub fn aggregate(partials: &[SampleMetricsPartialRecord]) -> Result<Vec<SampleMetricsRecord>> {
    // Key by the identity fields; read_type is carried as a label and validated
    // when the final record is built.
    type Key = (
        String,
        Option<String>,
        Option<String>,
        Option<String>,
        String,
        Option<String>,
        Option<String>,
    );

    struct Acc {
        n_target_positions: i64,
        total_fragment_bases: i64,
        on_target_fragment_bases: i64,
        hist: BTreeMap<i32, i64>,
    }

    let mut groups: BTreeMap<Key, Acc> = BTreeMap::new();
    for p in partials {
        let key = (
            p.sample_id.clone(),
            p.subject_id.clone(),
            p.sample_type.clone(),
            p.batch.clone(),
            p.read_type.to_string(),
            p.pipeline.clone(),
            p.input_checksum_sha256.clone(),
        );
        let acc = groups.entry(key).or_insert_with(|| Acc {
            n_target_positions: 0,
            total_fragment_bases: 0,
            on_target_fragment_bases: 0,
            hist: BTreeMap::new(),
        });
        acc.n_target_positions += p.n_target_positions;
        acc.total_fragment_bases += p.total_fragment_bases;
        acc.on_target_fragment_bases += p.on_target_fragment_bases;
        anyhow::ensure!(
            p.hist_depth.len() == p.hist_count.len(),
            "partial metrics for {} has mismatched histogram arrays",
            p.sample_id
        );
        for (&d, &c) in p.hist_depth.iter().zip(&p.hist_count) {
            *acc.hist.entry(d).or_insert(0) += c;
        }
    }

    let mut out = Vec::with_capacity(groups.len());
    for (key, acc) in groups {
        let (sample_id, subject_id, sample_type, batch, read_type_label, pipeline, checksum) = key;
        let read_type = ReadType::from_label(&read_type_label).with_context(|| {
            format!("partial metrics for {sample_id} has unknown read_type '{read_type_label}'")
        })?;
        let stats = sample_metrics_math::depth_stats(
            &acc.hist,
            acc.n_target_positions,
            acc.total_fragment_bases,
            acc.on_target_fragment_bases,
        );
        out.push(SampleMetricsRecord {
            sample_id,
            subject_id,
            sample_type,
            batch,
            read_type,
            pipeline,
            input_checksum_sha256: checksum,
            n_target_positions: acc.n_target_positions,
            n_target_positions_covered: stats.n_target_positions_covered,
            mean_target_depth_covered: stats.mean_target_depth_covered,
            mean_target_depth_all: stats.mean_target_depth_all,
            median_target_depth_covered: stats.median_target_depth_covered,
            median_target_depth_all: stats.median_target_depth_all,
            pct_fragment_bases_on_target: stats.pct_fragment_bases_on_target,
        });
    }
    Ok(out)
}

fn read_partial_parquet(path: &Path) -> Result<Vec<SampleMetricsPartialRecord>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("failed to open partial metrics: {}", path.display()))?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .with_context(|| format!("failed to read parquet: {}", path.display()))?
        .build()
        .with_context(|| format!("failed to build parquet reader: {}", path.display()))?;

    let mut out = Vec::new();
    for batch in reader {
        let batch = batch.context("failed to read parquet batch")?;
        let sample_id = str_col(&batch, "sample_id")?;
        let subject_id = str_col(&batch, "subject_id")?;
        let sample_type = str_col(&batch, "sample_type")?;
        let batch_col = str_col(&batch, "batch")?;
        let read_type = str_col(&batch, "read_type")?;
        let pipeline = str_col(&batch, "pipeline")?;
        let checksum = str_col(&batch, "input_checksum_sha256")?;
        let n_target_positions = i64_col(&batch, "n_target_positions")?;
        let total_fragment_bases = i64_col(&batch, "total_fragment_bases")?;
        let on_target_fragment_bases = i64_col(&batch, "on_target_fragment_bases")?;
        let hist_depth = list_col(&batch, "hist_depth")?;
        let hist_count = list_col(&batch, "hist_count")?;

        for row in 0..batch.num_rows() {
            let depth_arr = hist_depth.value(row);
            let depth_arr = depth_arr
                .as_any()
                .downcast_ref::<Int32Array>()
                .context("hist_depth values are not Int32")?;
            let count_arr = hist_count.value(row);
            let count_arr = count_arr
                .as_any()
                .downcast_ref::<Int64Array>()
                .context("hist_count values are not Int64")?;

            out.push(SampleMetricsPartialRecord {
                sample_id: req_str(sample_id, row, "sample_id")?,
                subject_id: opt_str(subject_id, row),
                sample_type: opt_str(sample_type, row),
                batch: opt_str(batch_col, row),
                read_type: ReadType::from_label(&req_str(read_type, row, "read_type")?)
                    .with_context(|| format!("unknown read_type in {}", path.display()))?,
                pipeline: opt_str(pipeline, row),
                input_checksum_sha256: opt_str(checksum, row),
                n_target_positions: n_target_positions.value(row),
                total_fragment_bases: total_fragment_bases.value(row),
                on_target_fragment_bases: on_target_fragment_bases.value(row),
                hist_depth: depth_arr.values().to_vec(),
                hist_count: count_arr.values().to_vec(),
            });
        }
    }
    Ok(out)
}

fn str_col<'a>(batch: &'a arrow::record_batch::RecordBatch, name: &str) -> Result<&'a StringArray> {
    batch
        .column_by_name(name)
        .with_context(|| format!("missing column {name}"))?
        .as_any()
        .downcast_ref::<StringArray>()
        .with_context(|| format!("column {name} is not Utf8"))
}

fn i64_col<'a>(batch: &'a arrow::record_batch::RecordBatch, name: &str) -> Result<&'a Int64Array> {
    batch
        .column_by_name(name)
        .with_context(|| format!("missing column {name}"))?
        .as_any()
        .downcast_ref::<Int64Array>()
        .with_context(|| format!("column {name} is not Int64"))
}

fn list_col<'a>(batch: &'a arrow::record_batch::RecordBatch, name: &str) -> Result<&'a ListArray> {
    batch
        .column_by_name(name)
        .with_context(|| format!("missing column {name}"))?
        .as_any()
        .downcast_ref::<ListArray>()
        .with_context(|| format!("column {name} is not a List"))
}

fn opt_str(arr: &StringArray, row: usize) -> Option<String> {
    if arr.is_null(row) {
        None
    } else {
        Some(arr.value(row).to_string())
    }
}

fn req_str(arr: &StringArray, row: usize, name: &str) -> Result<String> {
    anyhow::ensure!(!arr.is_null(row), "required column {name} is null");
    Ok(arr.value(row).to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn partial(
        sample: &str,
        n_target: i64,
        total_frag: i64,
        on_target_frag: i64,
        hist: &[(i32, i64)],
    ) -> SampleMetricsPartialRecord {
        SampleMetricsPartialRecord {
            sample_id: sample.to_string(),
            subject_id: None,
            sample_type: None,
            batch: None,
            read_type: ReadType::Duplex,
            pipeline: None,
            input_checksum_sha256: None,
            n_target_positions: n_target,
            total_fragment_bases: total_frag,
            on_target_fragment_bases: on_target_frag,
            hist_depth: hist.iter().map(|&(d, _)| d).collect(),
            hist_count: hist.iter().map(|&(_, c)| c).collect(),
        }
    }

    #[test]
    fn aggregates_two_shards_exactly() {
        // Shard A: covered depths [10, 20], 3 target positions (1 zero-depth).
        // Shard B: covered depth [30], 2 target positions (1 zero-depth).
        // Combined: covered [10,20,30], 5 target positions, 2 zero-depth.
        let a = partial("s1", 3, 100, 80, &[(10, 1), (20, 1)]);
        let b = partial("s1", 2, 50, 50, &[(30, 1)]);
        let out = aggregate(&[a, b]).unwrap();
        assert_eq!(out.len(), 1);
        let r = &out[0];
        assert_eq!(r.n_target_positions, 5);
        assert_eq!(r.n_target_positions_covered, 3);
        // covered median of [10,20,30] = 20
        assert_eq!(r.median_target_depth_covered, Some(20.0));
        // all median of [0,0,10,20,30] = 10
        assert_eq!(r.median_target_depth_all, Some(10.0));
        // mean covered = 60/3 = 20; mean all = 60/5 = 12
        assert_eq!(r.mean_target_depth_covered, Some(20.0));
        assert_eq!(r.mean_target_depth_all, Some(12.0));
        // pct on target = 130/150
        assert!((r.pct_fragment_bases_on_target.unwrap() - 130.0 / 150.0).abs() < 1e-6);
    }

    #[test]
    fn separates_distinct_samples() {
        let a = partial("s1", 1, 10, 10, &[(5, 1)]);
        let b = partial("s2", 1, 10, 10, &[(7, 1)]);
        let out = aggregate(&[a, b]).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].sample_id, "s1");
        assert_eq!(out[1].sample_id, "s2");
    }
}
