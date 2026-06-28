# wdl/adhoc/ — disposable analysis WDLs

One-off WDLs tied to a specific investigation. These are **not supported pipelines** and may
be deleted once the investigation they serve is done. Anything here that proves durably useful
should be promoted to `wdl/` or `wdl/experimental/` with proper docs and tests.

## `fusion_unique_anchor_eval.wdl`

Re-runs the fusion caller on the **evidence BAMs** from a prior `geac experimental fusions`
cohort run (the small `*.fusions.bam` outputs) with `--emit-unique-anchor`, so each call gains
the `n_unique_anchored` column, then merges the tiny per-sample Parquets into a single cohort
DuckDB.

**Purpose:** evaluate the unique-anchor specificity filter (the `n_unique_anchored / supporting_reads`
N-sweep) on a cohort whose original run predates the column — **without** re-processing the
multi-GB input BAMs or localizing ~71 GB of evidence BAMs. The job runs entirely in the cloud and
returns one small DuckDB; the N-sweep + truth scoring then run offline on that single file.

Requires a geac image with the `n_unique_anchored` column (>= 0.4.60) and the **same copy-labeled
index** (`--check-genome-uniqueness`) used for the original run.

> **⚠ Index must match the evidence BAMs.** Re-running with a different index changes k-mer matching
> and the spanning-read set, so the breakpoint-geometry tags (`chimera`/`samelocus`) will diverge from
> the original full-BAM run. `n_unique_anchored` is robust to this; the geometry filters are not. Use
> the exact index that produced the evidence BAMs (for the panel edit-0 cohort, the copy-labeled
> edit-0 panel index).

## `fusion_truth_eval.wdl`

Scores the per-sample calls (with `n_unique_anchored`) against the Terra `truth_fusion` column and
sweeps the unique-anchor threshold N, emitting a cohort sensitivity/precision curve. Pure-stdlib
Python on the small fusions TSVs (no duckdb/pandas; runs on `python:3.11-slim`).

Scatter a `fusion_sample_set`; Terra fills the parallel arrays from
`this.fusion_samples.{fusions_tsv, fusion_sample_id, truth_fusion, truth_dilution}`. At each
threshold N a call is kept when `filter == "PASS"` and `n_unique_anchored >= N`; for a positive
sample, detection = a kept call whose gene pair equals the truth pair (order-insensitive), and every
other kept call is a false positive; for a negative, every kept call is a false positive.

Outputs `cohort_sweep.tsv` (per-N `sensitivity` / `precision` / `total_fp`) and `cohort_detail.tsv`
(per sample × threshold). Validated locally on the CIC/DUX4 fixtures: the true CIC::DUX4 sample had
`n_unique_anchored = 1350` vs a negative control's `4`, so any N in [5, 1350] separates them cleanly.

The truth column was loaded onto the Terra `fusion_sample` entities from the investigation's
`truth.tsv` (`truth_fusion`, `truth_expected`, `truth_dilution`); see the investigation docs.

### Inputs (example — replace placeholders with your own paths)

```json
{
  "FusionUniqueAnchorEval.evidence_bams": [
    "gs://BUCKET/RUN/SAMPLE_1.fusions.bam",
    "gs://BUCKET/RUN/SAMPLE_2.fusions.bam"
  ],
  "FusionUniqueAnchorEval.fusion_index": "gs://BUCKET/resources/INDEX_copylabeled.duckdb",
  "FusionUniqueAnchorEval.docker_image": "ghcr.io/OWNER/geac:0.4.60"
}
```

The caller-config inputs default to the original investigation's report config
(`min_supporting_reads=1`, `min_kmer_hits=1`, `min_mapq=0`, `max_breakpoint_std=100`,
`min_breakpoint_distance=10000`, `kmer_size=23`). Override them only to keep the reproduced
PASS / chimera / samelocus set faithful to the run being evaluated.

### Output

- `cohort_duckdb` — one DuckDB with a `fusions` table for all samples, carrying
  `n_unique_anchored`, `supporting_reads`, and `filter`. Pull this single file (KB–MB) and run the
  threshold sweep offline, e.g.:

  ```sql
  -- unique-anchored fraction per call; sweep the threshold N offline
  SELECT sample_id, gene_a, gene_b, supporting_reads, n_unique_anchored, filter
  FROM fusions
  WHERE filter = 'PASS';
  ```
