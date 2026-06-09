# GEAC — Genomic Evidence Atlas of Cohorts

GEAC is a Rust command-line toolkit for building cohort-scale sequencing evidence
databases from BAM/CRAM files. It collects alt-base, coverage, read-level, and
fragment-level metrics into Parquet files, then merges them into DuckDB databases
for SQL queries and Streamlit exploration.

The stable CLI includes:

```text
collect          Process one BAM/CRAM into alt-base Parquet records
merge            Merge Parquet files or DuckDB databases into one cohort DuckDB
inspect          Check a merged cohort DuckDB for schema, metadata, and resource issues
qc               Print per-sample QC summaries
cohort           Summarize recurrent loci across samples
sample-identity  Find likely duplicate / same-individual samples
annotate-normal  Cross-check tumor loci against a paired normal BAM/CRAM
annotate-pon     Cross-check tumor loci against a Panel of Normals DuckDB
coverage         Emit per-position or binned coverage metrics
fragments        Emit per-fragment insert-size, GC, and end-motif metrics
experimental     K-mer / fusion tools; APIs may change
```

## Install

### Homebrew, macOS arm64

```bash
brew install fleharty/geac/geac
```

This installs `geac`, `geac-cohort`, and `geac-coverage-explorer`.

### Docker, linux/amd64

```bash
docker pull ghcr.io/fleharty/geac:latest
```

Use versioned images for reproducible workflows:

```bash
docker pull ghcr.io/fleharty/geac:0.4.39
```

The Docker image contains the `geac` binary for batch/cloud execution. The Streamlit
Explorer launchers are installed by Homebrew/source installs, not by the Docker image.

### From Source

Requires Rust, htslib, and pkg-config.

```bash
brew install htslib pkg-config
cargo build --release
```

## Quickstart

Collect one sample:

```bash
geac collect \
  --input sample.bam \
  --reference hg38.fa \
  --output SAMPLE_001.parquet \
  --read-type duplex \
  --pipeline fgbio
```

Merge samples into a cohort database:

```bash
geac merge --output cohort.duckdb samples/*.parquet
```

Find likely duplicate or same-individual samples:

```bash
geac sample-identity \
  --input cohort.duckdb \
  --output sample_identity_pairs.tsv
```

Run the main Explorer:

```bash
geac-cohort
```

Run the coverage Explorer:

```bash
geac-coverage-explorer
```

## Common Workflows

### Alt-Base Cohort

```text
geac collect  →  per-sample .parquet / .locus.parquet
geac merge    →  cohort.duckdb
geac-cohort   →  interactive review
```

Use `geac collect --reads-output` when you need read-level plots, per-read filters,
MNV candidate discovery, family-size filtering, or read-context diagnostics.

Use `geac collect --targets` to annotate on-target status and emit
`.sample_metrics.parquet` files with target-depth summaries.

### Tumor/Normal And PoN Annotation

```text
geac collect          →  TUMOR.locus.parquet
geac annotate-normal  →  TUMOR.normal_evidence.parquet
geac annotate-pon     →  TUMOR.pon_evidence.parquet
geac merge            →  cohort.duckdb with evidence tables
```

### Coverage

```bash
geac coverage \
  --input SAMPLE.bam \
  --reference hg38.fa \
  --output SAMPLE.coverage.parquet \
  --targets capture_targets.bed

geac merge --output coverage.duckdb *.coverage.parquet
geac-coverage-explorer
```

### Fragmentomics

```bash
geac fragments \
  --input SAMPLE.bam \
  --reference hg38.fa \
  --output SAMPLE.fragments.parquet \
  --read-type duplex \
  --pipeline fgbio

geac merge --output fragments.duckdb *.fragments.parquet
```

## Data Model

`geac merge` routes Parquet files into DuckDB tables by suffix:

| Suffix | DuckDB table |
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

Existing `.duckdb` inputs can also be merged directly.

## Query Example

```sql
SELECT chrom, pos, ref_allele, alt_allele,
       COUNT(DISTINCT sample_id) AS n_samples,
       SUM(alt_count) AS total_alt_reads
FROM alt_bases
GROUP BY chrom, pos, ref_allele, alt_allele
ORDER BY n_samples DESC;
```

## Documentation

| Topic | Link |
|---|---|
| CLI commands and options | [docs/cli.md](docs/cli.md) |
| DuckDB and Parquet schemas | [docs/schema.md](docs/schema.md) |
| Streamlit Explorer apps | [docs/explorer.md](docs/explorer.md) |
| WDL / Terra workflows | [docs/wdl.md](docs/wdl.md) |
| Release checklist | [docs/release.md](docs/release.md) |
| Provenance tables | [docs/provenance.md](docs/provenance.md) |
| Experimental fusion tools | [docs/EXPERIMENTAL.md](docs/EXPERIMENTAL.md) |
| Code audit log | [docs/CODE_AUDIT.md](docs/CODE_AUDIT.md) |

## Development

```bash
cargo check -q
cargo test -q
cargo build --release
```

For Explorer helper tests:

```bash
python -m pytest app/tests
```

Release notes and historical engineering notes live in
[docs/DEVELOPMENT_LOG.md](docs/DEVELOPMENT_LOG.md).
