# Provenance

GEAC records provenance at two levels:

- Per-sample collect outputs (`*.parquet`, `*.locus.parquet`, `*.reads.parquet`)
- Merged cohort DuckDB databases (`cohort.duckdb`)

## Collect-time provenance

`geac collect` writes provenance columns into the per-sample Parquet outputs:

- `read_type`
- `pipeline`
- `batch`
- `label1`
- `label2`
- `label3`
- `input_checksum_sha256`

`input_checksum_sha256` is enabled by default. It can be disabled with:

```bash
geac collect ... # --input-checksum-sha256 is on by default; pass false in the WDL to disable
```

This computes a SHA-256 of the input BAM/CRAM once during collection and stores the
same value on every output row for that sample. Hashing large alignment files adds
I/O and wall time, but the checksum is useful for provenance tracking and deduplication.

The WDLs expose the same behavior as:

- `input_checksum_sha256 = true` in [wdl/geac_collect.wdl](../wdl/geac_collect.wdl)
- `input_checksum_sha256 = true` in [wdl/geac_cohort.wdl](../wdl/geac_cohort.wdl)

## Merge-time provenance

`geac merge` writes two provenance tables to the output DuckDB.

### `geac_metadata`

`geac_metadata` is a one-row database header with:

- `schema_version`
- `geac_version`
- `created_at`
- `command_line`
- `output_path`
- `platform_os`
- `platform_arch`
- `platform_family`
- `n_alt_bases_inputs`
- `n_alt_reads_inputs`
- `n_normal_evidence_inputs`
- `n_pon_evidence_inputs`
- `n_coverage_inputs`
- `n_duckdb_inputs`
- `n_samples`
- `alt_bases_rows`
- `alt_reads_rows`
- `normal_evidence_rows`
- `pon_evidence_rows`
- `coverage_rows`
- `samples_rows`

### `geac_inputs`

`geac_inputs` contains one row per source artifact merged into the database:

- `input_path`
- `input_kind`
- `source_kind`
- `file_size_bytes`
- `modified_at`
- `checksum_sha256`
- `sample_count`
- `row_count`

`checksum_sha256` is currently reserved and written as null for merge inputs.

## Explorer visibility

The Streamlit explorers expose a compact view of both tables in the sidebar under:

- `Advanced` -> database metadata

This makes it possible to inspect provenance without opening DuckDB manually.
