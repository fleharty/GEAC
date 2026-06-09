# WDL / Terra

WDL 1.0 workflows live in `wdl/`.

| Workflow | Status | Purpose |
|---|---|---|
| `geac_cohort.wdl` | Tested | Scatter `geac collect`, then merge into a cohort DuckDB |
| `geac_coverage.wdl` | Tested | Scatter `geac coverage`, then merge coverage outputs |
| `geac_fragments.wdl` | Untested | Scatter `geac fragments` across a cohort |
| `geac_cohort_loci.wdl` | Untested | Run `geac cohort` on per-sample Parquets |
| `geac_collect.wdl` | Untested | Single-sample `geac collect` wrapper |
| `geac_merge.wdl` | Untested | Standalone merge into DuckDB |
| `geac_qc.wdl` | Untested | Run `geac qc` on Parquet inputs |
| `geac_annotate_normal.wdl` | Untested | Cross-check tumor loci against a paired normal BAM |
| `geac_annotate_pon.wdl` | Untested | Cross-check tumor loci against a PoN DuckDB |
| `experimental/geac_build_fusion_index.wdl` | Experimental | Build a fusion k-mer index |
| `experimental/geac_fusions.wdl` | Experimental | Run `geac experimental fusions` |

## Docker Image

Use the versioned image that matches the GEAC release:

```text
ghcr.io/fleharty/geac:<version>
```

For example:

```text
ghcr.io/fleharty/geac:0.4.39
```

## Common Inputs

Most workflows include runtime controls:

| Input | Description |
|---|---|
| `docker_image` | GEAC Docker image |
| `memory_gb` | Runtime memory |
| `disk_gb` | Runtime disk |
| `preemptible` | Terra preemptible retry count |
| `threads` | Threads where supported |

## `geac_collect.wdl`

Required inputs:

| Input | Type | Description |
|---|---|---|
| `input_bam` | File | BAM or CRAM |
| `input_bam_index` | File | `.bai` or `.crai` |
| `reference_fasta` | File | Reference FASTA |
| `reference_fasta_index` | File | `.fai` |
| `read_type` | String | `duplex`, `simplex`, or `raw` |
| `pipeline` | String | `fgbio`, `dragen`, or `raw` |

Common optional inputs include `sample_id`, `subject_id`, `sample_type`, `batch`,
`label1`, `label2`, `label3`, `timepoint`, `vcf`, `vcf_index`, `variants_tsv`,
`targets`, `gene_annotations`, `region`, `region_bed`, `gnomad`, `gnomad_index`,
`gnomad_af_field`, `bam_uri`, `bai_uri`, `variants_uri`, `gnomad_uri`,
`targets_uri`,
`reads_output`, `input_checksum_sha256`, pileup quality thresholds, and inclusion
flags for duplicates/secondary/supplementary reads.

On Terra, bind `bam_uri`, `bai_uri`, `gnomad_uri`, and `targets_uri` to stable
string columns containing the original `gs://` URIs. The corresponding `File`
inputs are still localized for compute, but the string URI inputs are what GEAC
stores in Parquet metadata for IGV session generation.

Outputs:

| Output | Description |
|---|---|
| `locus_parquet` | Per-sample locus Parquet |
| `reads_parquets` | Empty or one `.reads.parquet` when `reads_output=true` |
| `sample_metrics_parquets` | Empty or one `.sample_metrics.parquet` when `targets` is set |

## `geac_cohort.wdl`

Accepts parallel arrays of BAMs, BAM indices, optional sample IDs, optional per-sample
VCFs/TSVs, and optional metadata arrays. Shared inputs include reference, targets,
gene annotations, region, pileup thresholds, gnomAD, runtime resources, and optional
fragment collection.

For Terra IGV support, provide optional `Array[String] bam_uris` and
`Array[String] bai_uris` in parallel with `input_bams` and `input_bam_indices`,
plus shared `String? gnomad_uri` and `String? targets_uri` when those resources
are used. These should come from stable data-table string columns containing the
original `gs://` paths. GEAC uses those string URIs for embedded resource metadata
and for the emitted `cohort_manifest`; the `File` inputs remain the localized
inputs used by Cromwell for compute.

Outputs:

| Output | Description |
|---|---|
| `cohort_db` | Merged DuckDB |
| `cohort_manifest` | Resource URI manifest plus sample metadata |
| `cohort_on_target_tsv` | On-target alt-base TSV when targets are provided |

## Other Workflows

| Workflow | Key output |
|---|---|
| `geac_merge.wdl` | `cohort_db` |
| `geac_coverage.wdl` | Coverage cohort DuckDB |
| `geac_fragments.wdl` | Fragment Parquets suitable for merge |
| `geac_qc.wdl` | QC summary TSV/report |
| `geac_annotate_normal.wdl` | `.normal_evidence.parquet` |
| `geac_annotate_pon.wdl` | `.pon_evidence.parquet` |

## Terra Notes

1. Import the desired WDL into the workspace.
2. Set `docker_image` to a versioned GEAC image, for example `ghcr.io/fleharty/geac:0.4.39`.
3. For single-sample workflows, bind BAM/CRAM, index, reference FASTA, and reference index columns.
4. For cohort workflows, provide parallel arrays or workspace table expressions.
5. Use `geac_merge.wdl` to combine existing Parquet artifacts.
