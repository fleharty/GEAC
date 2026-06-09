# GEAC Schema Reference

Coordinates are 0-based unless noted. `geac merge` stores routed Parquet files in
DuckDB tables with the names shown below.

## Table Routing

| File suffix | DuckDB table |
|---|---|
| `.locus.parquet` or ordinary `.parquet` | `alt_bases` |
| `.reads.parquet` | `alt_reads` |
| `.normal_evidence.parquet` | `normal_evidence` |
| `.pon_evidence.parquet` | `pon_evidence` |
| `.coverage.parquet` | `coverage` |
| `.coverage.intervals.parquet` | `coverage_intervals` |
| `.sample_metrics.parquet` | `sample_metrics` |
| `.fragments.parquet` | `fragments` |
| `.fusions.parquet` | `fusions` |

## `alt_bases`

One row per observed alt allele at a locus.

Core columns:

| Column | Type | Description |
|---|---|---|
| `sample_id` | string | Sample identifier |
| `chrom` | string | Chromosome |
| `pos` | int64 | 0-based position |
| `ref_allele` | string | Reference allele |
| `alt_allele` | string | Alt allele, e.g. `T`, `+ACG`, `-2` |
| `variant_type` | string | `SNV`, `insertion`, or `deletion` |
| `total_depth` | int32 | Fragment depth; overlapping mates count once |
| `alt_count` | int32 | Fragments supporting the alt allele |
| `ref_count` | int32 | Fragments supporting the reference allele |
| `fwd_depth`, `rev_depth` | int32 | Strand-specific fragment depth |
| `fwd_alt_count`, `rev_alt_count` | int32 | Strand-specific alt support |
| `fwd_ref_count`, `rev_ref_count` | int32 | Strand-specific reference support |
| `overlap_depth` | int32 | Overlapping fragment pairs at the locus |
| `overlap_alt_agree` | int32 | Overlaps where both mates support the alt |
| `overlap_alt_disagree` | int32 | Overlaps where mates disagree |
| `overlap_ref_agree` | int32 | Overlaps where both mates support reference |

Metadata and provenance:

| Column | Type | Description |
|---|---|---|
| `read_type` | string | `raw`, `simplex`, or `duplex` |
| `pipeline` | string? | `fgbio`, `dragen`, or `raw` |
| `subject_id` | string? | Biological subject identifier |
| `sample_type` | string? | Sample substrate type |
| `batch` | string? | Batch label |
| `label1`, `label2`, `label3` | string? | User-defined labels |
| `timepoint` | string? | Longitudinal timepoint label |
| `input_checksum_sha256` | string? | BAM/CRAM SHA-256 when requested |
| `bam_path`, `bai_path` | string? | Input alignment paths or canonical URIs |
| `variants_path` | string? | VCF/TSV annotation path or URI |
| `gnomad_path` | string? | gnomAD path or URI |

Annotations:

| Column | Type | Description |
|---|---|---|
| `variant_called` | bool? | Whether a VCF/TSV call overlaps the locus |
| `variant_filter` | string? | VCF FILTER value |
| `on_target` | bool? | Whether the locus overlaps targets |
| `gene` | string? | Gene annotation |
| `homopolymer_len` | int32 | Longest homopolymer run overlapping the locus |
| `str_period` | int32 | Shortest tandem-repeat period; `0` when absent |
| `str_len` | int32 | STR tract length; `0` when absent |
| `trinuc_context` | string? | SNV trinucleotide context |
| `gnomad_af` | float32? | gnomAD allele frequency |

Read-context summaries appear when `--reads-output` is used:

| Column | Type | Description |
|---|---|---|
| `n_alt_reads_with_n_ctx` | int32? | Alt-supporting reads with N-context data |
| `mean_frac_n_before`, `mean_frac_n_after` | float32? | Mean N fraction before/after the alt in read sequence order |
| `mean_delta_n_frac` | float32? | `after - before` N fraction |
| `frac_reads_asymmetric` | float32? | Fraction of strongly asymmetric N-context reads |

## `alt_reads`

Produced by `geac collect --reads-output`. Linked to `alt_bases` by
`sample_id`, `chrom`, `pos`, and `alt_allele`.

| Column | Type | Description |
|---|---|---|
| `sample_id`, `chrom`, `pos`, `alt_allele` | mixed | Locus key |
| `read_type`, `pipeline` | string | Read and pipeline metadata |
| `subject_id`, `sample_type`, `batch`, `label1`, `label2`, `label3`, `timepoint` | string? | Sample metadata |
| `cycle` | int32 | 1-based sequencing cycle at the alt |
| `read_length` | int32 | Stored read length |
| `is_read1` | bool | Whether the read is R1 |
| `ab_count`, `ba_count`, `family_size` | int32? | Pipeline-aware consensus/family support tags |
| `base_qual`, `map_qual` | int32 | Base and mapping quality |
| `insert_size` | int32? | SAM TLEN when available |
| `frag_gc` | float32? | Reference GC fraction across inferred fragment span |
| `n_before_alt`, `n_after_alt` | int32 | Read bases before/after the alt |
| `n_n_before_alt`, `n_n_after_alt` | int32 | N bases before/after the alt |
| `leading_n_run_len`, `trailing_n_run_len` | int32 | Adjacent N-run lengths |
| `input_checksum_sha256` | string? | Input checksum when requested |
| `fragment_id` | int64 | FNV-1a hash of read name for fragment-level joins |

## `sample_metrics`

Produced by `geac collect` when `--targets` is provided.

| Column | Type | Description |
|---|---|---|
| `sample_id` | string | Sample identifier |
| `subject_id`, `sample_type`, `batch` | string? | Optional sample metadata |
| `read_type`, `pipeline` | string? | Processing metadata |
| `input_checksum_sha256` | string? | Input checksum when requested |
| `n_target_positions` | int32 | Total target positions |
| `n_target_positions_covered` | int32 | Target positions with depth |
| `mean_target_depth_covered`, `mean_target_depth_all` | float32? | Mean target depth |
| `median_target_depth_covered`, `median_target_depth_all` | float32? | Median target depth |
| `pct_fragment_bases_on_target` | float32? | On-target fragment-base percentage |

## `normal_evidence`

Produced by `geac annotate-normal`.

| Column | Type | Description |
|---|---|---|
| `tumor_sample_id` | string | Tumor sample |
| `chrom`, `pos`, `tumor_alt_allele` | mixed | Tumor locus key |
| `normal_sample_id` | string | Normal sample |
| `normal_alt_allele` | string? | Normal alt allele; null anchor row records depth only |
| `normal_depth` | int32 | Normal depth at the locus |
| `normal_alt_count` | int32 | Normal support for `normal_alt_allele` |

## `pon_evidence`

Produced by `geac annotate-pon`.

| Column | Type | Description |
|---|---|---|
| `tumor_sample_id`, `chrom`, `pos`, `tumor_alt_allele` | mixed | Tumor locus key |
| `n_pon_samples` | int64 | PoN samples carrying the allele |
| `pon_total_samples` | int64 | Total PoN samples |
| `max_pon_vaf`, `mean_pon_vaf` | float64? | PoN VAF summaries |

## `coverage`

Produced by `geac coverage`.

| Column | Type | Description |
|---|---|---|
| `sample_id`, `chrom`, `pos`, `end`, `bin_n` | mixed | Position or bin key |
| `total_depth`, `min_depth`, `max_depth`, `fwd_depth`, `rev_depth` | int32 | Depth metrics |
| `raw_read_depth`, `frac_dup` | int32/float32 | Duplicate-aware raw read metrics |
| `overlap_depth`, `frac_overlap` | int32/float32 | Mate-overlap metrics |
| `mean_mapq`, `frac_mapq0`, `frac_low_mapq` | float32 | Mapping quality metrics |
| `mean_base_qual`, `min_base_qual_obs`, `max_base_qual_obs`, `frac_low_bq` | mixed | Base quality metrics |
| `mean_insert_size`, `min_insert_size`, `max_insert_size`, `n_insert_size_obs` | mixed | Insert-size metrics |
| `gc_content` | float32 | Reference GC fraction |
| `on_target`, `gene`, `feature_type`, `exon_number` | mixed | Annotation columns |
| `read_type`, `pipeline`, `subject_id`, `sample_type`, `batch`, `label1`, `label2`, `label3`, `timepoint` | mixed | Metadata |

Additional nullable float columns are added for each `--track NAME:FILE`.

## `coverage_intervals`

Produced by `geac coverage --intervals-output`.

| Column | Type | Description |
|---|---|---|
| `sample_id`, `chrom`, `start`, `end`, `interval_name` | mixed | Target interval key |
| `gene`, `feature_type`, `exon_number` | mixed | Annotation columns |
| `n_bases`, `mean_depth`, `median_depth`, `min_depth`, `max_depth` | mixed | Depth summaries |
| `frac_at_1x`, `frac_at_10x`, `frac_at_20x`, `frac_at_30x`, `frac_at_50x`, `frac_at_100x` | float32 | Breadth summaries |
| `mean_gc_content`, `mean_mapq`, `mean_frac_mapq0`, `mean_frac_dup`, `mean_frac_overlap`, `mean_base_qual`, `mean_insert_size` | float32 | Mean QC metrics |
| `read_type`, `pipeline`, `subject_id`, `sample_type`, `batch`, `label1`, `label2`, `label3`, `timepoint` | mixed | Metadata |

## `fragments`

Produced by `geac fragments`.

| Column | Type | Description |
|---|---|---|
| `sample_id`, `chrom`, `frag_start`, `frag_end`, `midpoint` | mixed | Fragment coordinates |
| `insert_size` | int32 | Absolute TLEN |
| `gc_content` | float32? | Reference GC fraction over fragment span |
| `end_motif_5p`, `end_motif_3p` | string? | 4-mer end motifs |
| `map_qual` | int32 | R1 mapping quality |
| `read_type`, `pipeline`, `subject_id`, `sample_type`, `batch`, `label1`, `label2`, `label3`, `timepoint` | mixed | Metadata |

## `fusions`

Produced by `geac experimental fusions`. See [EXPERIMENTAL.md](EXPERIMENTAL.md)
for the evolving schema.

## Provenance Tables

`geac merge` writes:

| Table | Purpose |
|---|---|
| `geac_metadata` | One-row database header with GEAC version, schema version, command line, platform, input counts, and row counts |
| `geac_inputs` | One row per input artifact with source metadata |

See [provenance.md](provenance.md) for details.
