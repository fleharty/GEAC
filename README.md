# GEAC — Genomic Evidence Atlas of Cohorts

Collect alt base metrics from duplex/simplex BAM/CRAM files across a sequencing cohort.
Each sample is processed independently into a Parquet file; samples can then be merged
into a DuckDB database for cohort-level queries.

## Overview

GEAC is a Rust command-line tool designed for large-scale sequencing cohorts (thousands
of samples, ~2 MB panels). It performs pileup-based analysis of consensus-called BAM/CRAM
files (fgbio, DRAGEN, or raw reads) and records every position where an alt allele is
observed, along with rich per-locus metrics:

- Forward/reverse strand counts for alt and reference alleles
- Overlapping fragment pair agreement (mate-overlap concordance)
- Variant type classification (SNV, insertion, deletion, MNV)
- Optional VCF annotation — whether a variant was called and its filter status

Each sample produces a Parquet file (~MB scale). Samples are then merged into a DuckDB
cohort database for efficient SQL queries across thousands of samples.

## Setup

Requires [Miniforge](https://github.com/conda-forge/miniforge) (or Anaconda) and internet access.

```bash
bash scripts/setup.sh
conda activate geac
cargo build --release
```

The setup script:
1. Creates the `geac` conda environment with `htslib`, Python dependencies, and Streamlit
2. Installs conda activate hooks so `LIBRARY_PATH`, `PKG_CONFIG_PATH`, and `RUSTFLAGS`
   (macOS rpath) are set automatically on `conda activate geac`
3. Installs or updates Rust via `rustup`

After building, optionally add a shell alias so you can run `geac` from anywhere:

```bash
echo 'alias geac="/path/to/GEAC/target/release/geac"' >> ~/.zshrc
source ~/.zshrc
```

## Usage

### Collect — process a single sample

Runs a pileup over the BAM/CRAM and writes one row per alt allele per locus to a Parquet file.

```bash
geac collect \
  --input sample.bam \
  --reference hg38.fa \
  --output SAMPLE_001.parquet \
  --read-type duplex \
  --pipeline fgbio
```

`--sample-id` is optional. If omitted, the SM tag is read from the BAM `@RG` header line.
If no SM tag is present, the command exits with an error.

Optional flags:

| Flag | Default | Description |
|---|---|---|
| `--sample-id` | from SM tag | Override the sample identifier |
| `--vcf` | — | Annotate loci with variant calling status from a VCF/BCF. Mutually exclusive with `--variants-tsv` |
| `--variants-tsv` | — | TSV variant list (columns: chrom, pos_start, pos_end, ref, var; 0-based). Alternative to `--vcf` |
| `--targets` | — | BED or Picard interval list of target regions; adds `on_target` bool column |
| `--gene-annotations` | — | GFF3, GTF, or UCSC genePred file; adds `gene` string column |
| `--repeat-window` | 10 | Bases on each side of locus to scan for homopolymers and STRs |
| `--min-base-qual` | 1 | Minimum base quality to count a read |
| `--min-map-qual` | 20 | Minimum mapping quality |
| `--region` | whole genome | Restrict to a genomic region (e.g. `chr1:1-1000000`) |
| `--threads` | 1 | Parallel processing threads |
| `--progress-interval` | 30 | Seconds between progress reports to stderr |

#### Read types and pipelines

`--read-type`: `duplex` | `simplex` | `raw`

`--pipeline`: `fgbio` | `dragen` | `raw`

These values are stored as metadata in the Parquet file and do not change processing behaviour —
they allow downstream filtering by sequencing strategy.

#### VCF annotation

When `--vcf` is provided, each alt allele record is annotated with:
- `variant_called` — `true` if a variant overlapping this locus was called, `false` if the
  locus was covered but no variant called, `null` if no VCF was provided
- `variant_filter` — the VCF FILTER value (`PASS`, a filter reason, or `null`)

SNVs are matched exactly by chrom/pos/alt allele. Indels are matched by position only,
since VCF left-aligned representation differs from GEAC's `+seq`/`-seq` notation.

#### Soft-clipped bases

Soft-clipped bases are **not** counted. `rust-htslib` pileup only yields bases that are
aligned to the reference at each position; soft clips are excluded by design.

#### Overlap detection

For paired-end reads where the two mates overlap the same locus, GEAC detects the overlap
by grouping pileup reads at a position by query name. If both mates are present at the
same locus, the pair is counted in `overlap_depth`. If both agree on the alt allele,
`overlap_alt_agree` is incremented; if they disagree, `overlap_alt_disagree` is incremented.
This is an integer count (not boolean) because multiple overlapping pairs can cover the same
locus in a deep pileup.

### Merge — combine samples into a cohort DuckDB

```bash
geac merge --output cohort.duckdb samples/*.parquet
```

Creates a DuckDB database with:
- `alt_bases` — all per-sample alt base records
- `samples` — one-row-per-sample summary (n_alt_loci, total_alt_reads, n_positions, etc.)
- Indices on `(chrom, pos)` and `sample_id` for fast queries

The output file must not already exist (use a new path or delete the old file first).

### Query the cohort (DuckDB)

You can query either a merged DuckDB or raw Parquet files directly.

```sql
-- Cohort frequency of each alt allele
SELECT chrom, pos, ref_allele, alt_allele,
       COUNT(DISTINCT sample_id) AS n_samples,
       SUM(alt_count)            AS total_alt_reads
FROM read_parquet('samples/*.parquet')
GROUP BY chrom, pos, ref_allele, alt_allele
ORDER BY n_samples DESC;

-- Positions with high strand bias in a specific sample
SELECT chrom, pos, ref_allele, alt_allele,
       fwd_alt_count, rev_alt_count,
       ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
FROM alt_bases
WHERE sample_id = 'SAMPLE_001'
  AND (fwd_alt_count = 0 OR rev_alt_count = 0)
  AND alt_count >= 5;

-- Overlap-discordant alt calls (potential errors)
SELECT chrom, pos, alt_allele, overlap_alt_disagree, overlap_alt_agree
FROM alt_bases
WHERE overlap_alt_disagree > overlap_alt_agree;
```

## Explorer UI

An interactive Streamlit web app for browsing and visualising alt base data:

```bash
conda activate geac
streamlit run app/geac_explorer.py
```

Then open `http://localhost:8501` in your browser. Enter a Parquet or DuckDB
file path in the text box to load data.

Features:
- **Summary statistics** — alt records, samples, total alt bases, mean VAF, mean depth,
  % variant called
- **Sidebar filters** — chromosome, samples, variant type, VAF range, min alt count,
  variant called status, min/max depth, on-target, gene name (partial match), homopolymer
  length, STR length; **Clear all filters** button resets all filters at once
- **Data table** — sortable, all schema columns, IGV session download button
- **Tabbed plots**
  - *VAF distribution* — separate histograms for SNV, insertion, deletion; click a bar
    to see matching records and download an IGV session
  - *SNV error spectrum* — click a substitution bar to see matching records and IGV session
  - *Strand bias* — forward vs. reverse alt reads scatter with 95% CI boundary lines;
    log1p or linear axis toggle; color by variant type, sample, on-target, or called status;
    click/shift-click to select points and view a drill-down table + IGV session
  - *SNV Trinucleotide Spectrum (SBS96)* — 3×2 grid of per-mutation-type panels with shared
    y-axis; shift-click to select multiple contexts; drill-down table and IGV session
  - *COSMIC Signature Decomposition* — NNLS fit of the SBS96 spectrum against COSMIC
    reference signatures; shows top N signatures with weights, etiology annotations,
    cosine similarity, and residual percentage
  - *Overlap agreement* — histogram of overlap concordance fractions
  - *Cohort* (DuckDB only) — per-sample summary table, VAF distribution overlay, strand
    balance scatter, SNV count bar chart stacked by SBS6 substitution type, and SBS96
    heatmap (samples × 96 contexts, normalised by sample); click a sample row to focus
    all other views
- **IGV integration** — provide a manifest TSV (`sample_id`, `bam_path`) in the sidebar
  to enable "Download IGV session" buttons throughout the app. Downloads a zip containing
  `session.xml` (BAM tracks + BED track) and `positions.bed` (one row per unique locus).
  Sessions are capped at 5 samples by default with an override option.

### Manifest format

```tsv
sample_id	bam_path	bai_path
SAMPLE_001	gs://my-bucket/bams/SAMPLE_001.bam	gs://my-bucket/bams/SAMPLE_001.bam.bai
SAMPLE_002	gs://my-bucket/bams/SAMPLE_002.bam	gs://my-bucket/bams/SAMPLE_002.bam.bai
SAMPLE_003	/local/path/to/SAMPLE_003.bam	/local/path/to/SAMPLE_003.bam.bai
```

`bai_path` is optional — if omitted or left blank, IGV will attempt to find the index automatically.

## Schema

Each Parquet file contains one row per alt allele observed at a locus.

| Column | Type | Description |
|---|---|---|
| `sample_id` | string | Sample identifier (from `--sample-id` or BAM SM tag) |
| `chrom` | string | Chromosome |
| `pos` | int64 | 0-based position |
| `ref_allele` | string | Reference allele |
| `alt_allele` | string | Alt allele (e.g. `T`, `+ACG`, `-2`) |
| `variant_type` | string | `SNV` / `insertion` / `deletion` / `MNV` |
| `total_depth` | int32 | Total read depth at position |
| `alt_count` | int32 | Reads supporting the alt allele |
| `ref_count` | int32 | Reads supporting the reference allele |
| `fwd_depth` | int32 | Forward strand depth |
| `rev_depth` | int32 | Reverse strand depth |
| `fwd_alt_count` | int32 | Forward strand alt reads |
| `rev_alt_count` | int32 | Reverse strand alt reads |
| `fwd_ref_count` | int32 | Forward strand reference reads |
| `rev_ref_count` | int32 | Reverse strand reference reads |
| `overlap_depth` | int32 | Number of reads at this locus that have an overlapping mate |
| `overlap_alt_agree` | int32 | Overlapping pairs where both mates support the alt |
| `overlap_alt_disagree` | int32 | Overlapping pairs where mates disagree on the alt |
| `overlap_ref_agree` | int32 | Overlapping pairs where both mates support the reference |
| `read_type` | string | `raw` / `simplex` / `duplex` |
| `pipeline` | string | `fgbio` / `dragen` / `raw` |
| `variant_called` | bool? | Whether a variant was called here (null if no VCF/TSV provided) |
| `variant_filter` | string? | VCF FILTER value (`PASS`, filter reason, or null) |
| `on_target` | bool? | Whether locus overlaps a target region (null if no `--targets` provided) |
| `gene` | string? | Gene name at locus (null if no `--gene-annotations` provided or intergenic) |
| `homopolymer_len` | int32? | Length of longest homopolymer overlapping locus within `--repeat-window` |
| `str_period` | int32? | Period of shortest tandem repeat unit at locus (null if no STR detected) |
| `str_len` | int32? | Total length of STR tract at locus (null if no STR detected) |
| `trinuc_context` | string? | Trinucleotide context for SNVs, e.g. `A[C>T]G` (null for indels/MNVs) |

## Docker

Two Dockerfiles are provided:

| File | When to use |
|---|---|
| `docker/Dockerfile` | Linux/CI — compiles geac inside Docker (two-stage build) |
| `docker/Dockerfile.prebuilt` | Apple Silicon Mac — copies a pre-compiled binary into a minimal image |

On Apple Silicon, `docker/Dockerfile` fails because compiling Rust + DuckDB C++ under
QEMU emulation (`--platform linux/amd64`) causes a SIGSEGV in `rustc`. The prebuilt
approach cross-compiles the binary natively on the Mac and copies it into the image,
avoiding QEMU entirely.

### Building on Apple Silicon (recommended)

**Prerequisites:** Podman, `gcloud` CLI, Zig, and `cargo-zigbuild` installed.

```bash
# Install build tools (once)
brew install zig
cargo install cargo-zigbuild
rustup target add x86_64-unknown-linux-gnu

# Authenticate Podman to GCR (once per session)
gcloud auth print-access-token | podman login -u oauth2accesstoken --password-stdin gcr.io
```

> **Note:** `gcloud auth configure-docker gcr.io` configures the Docker credential helper
> but Podman does not use it. Always use the `podman login` command above instead.

```bash
# 1. Cross-compile the geac binary for Linux/AMD64 on your Mac
cargo zigbuild --release --target x86_64-unknown-linux-gnu

# 2. Build the Docker image (no --platform flag needed — binary is already AMD64)
podman build -t gcr.io/my-gcp-project/geac:$(cat VERSION) -f docker/Dockerfile.prebuilt .

# 3. Push
podman push gcr.io/my-gcp-project/geac:$(cat VERSION)
```

### Building on Linux / CI

```bash
# Authenticate Podman to GCR (once per session)
gcloud auth print-access-token | podman login -u oauth2accesstoken --password-stdin gcr.io

# Build and push
bash scripts/build_docker.sh my-gcp-project --push
```

### Authenticating Podman to GCR

`gcloud auth configure-docker` only configures the Docker credential helper — Podman
does not read it. Authenticate Podman directly each session:

```bash
gcloud auth print-access-token | podman login -u oauth2accesstoken --password-stdin gcr.io
```

### Updating the version

```bash
echo "0.3.0" > VERSION
# also update version = "..." in Cargo.toml
# rebuild and push using the Apple Silicon steps above
git add VERSION Cargo.toml && git commit -m "Bump version to 0.3.0"
```

## WDL / Terra

Three WDL 1.0 workflows are provided in `wdl/`:

| Workflow | Purpose |
|---|---|
| `geac_collect.wdl` | Single-sample wrapper around `geac collect`; use this to scatter across a sample table |
| `geac_cohort.wdl` | Full cohort workflow: scatters `geac collect` then gathers with `geac merge` |
| `geac_merge.wdl` | Standalone merge — takes existing Parquets and builds a DuckDB |

### `geac_collect.wdl` inputs

| Input | Type | Description |
|---|---|---|
| `input_bam` | File | BAM or CRAM file |
| `input_bam_index` | File | `.bai` or `.crai` index |
| `reference_fasta` | File | Reference FASTA |
| `reference_fasta_index` | File | `.fai` index |
| `read_type` | String | `duplex` / `simplex` / `raw` |
| `pipeline` | String | `fgbio` / `dragen` / `raw` |
| `docker_image` | String | e.g. `gcr.io/my-project/geac:0.1.0` |
| `sample_id` | String? | Override sample ID (default: BAM SM tag) |
| `vcf` | File? | VCF/BCF for variant call annotation |
| `vcf_index` | File? | `.tbi` or `.csi` index for VCF |
| `variants_tsv` | File? | TSV variant list (alternative to `--vcf`) |
| `targets` | File? | BED or Picard interval list for on-target annotation |
| `gene_annotations` | File? | GFF3, GTF, or UCSC genePred for gene annotation |
| `region` | String? | Restrict to a region, e.g. `chr1:1-1000000` |
| `repeat_window` | Int | Bases each side of locus for homopolymer/STR scan (default: 10) |
| `min_base_qual` | Int | Default: 1 |
| `min_map_qual` | Int | Default: 20 |
| `threads` | Int | Default: 4 |
| `memory_gb` | Int | Default: 8 |
| `disk_gb` | Int | Default: 100 |
| `preemptible` | Int | Default: 2 |

Output: `parquet` (File) — per-sample alt base Parquet file.

### `geac_cohort.wdl` inputs

Per-sample parallel arrays: `input_bams`, `input_bam_indices`, optional `sample_ids`,
optional `variants_tsvs`, optional `vcfs` + `vcf_indices` (per-sample VCF annotation).
Shared inputs applied to all samples: `reference_fasta`, `targets`, `gene_annotations`,
`region`, `repeat_window`, `read_type`, `pipeline`, `min_base_qual`, `min_map_qual`, `threads`.

Outputs: `parquets` (Array[File]) and `cohort_db` (File, the merged DuckDB).

### `geac_merge.wdl` inputs

| Input | Type | Description |
|---|---|---|
| `parquets` | Array[File] | Per-sample Parquet files |
| `cohort_name` | String | Base name for the output DuckDB (default: `cohort`) |
| `docker_image` | String | geac Docker image |
| `memory_gb` | Int | Default: 16 |
| `disk_gb` | Int | Default: 50 |
| `preemptible` | Int | Default: 2 |

Output: `cohort_db` (File) — merged DuckDB database.

### Running on Terra

1. Import the desired WDL into your Terra workspace.
2. Set `docker_image` to your pushed GCR image (e.g. `gcr.io/my-project/geac:0.1.0`).
3. For `geac_collect.wdl`: link `input_bam`, `input_bam_index`, `reference_fasta`, and
   `reference_fasta_index` to your workspace data table columns; Terra will scatter automatically.
4. For `geac_cohort.wdl`: provide parallel arrays directly and let the workflow scatter and merge.
5. To merge existing Parquets, use `geac_merge.wdl` with the list of Parquet files.

## Architecture

```
geac collect  →  per-sample .parquet
geac merge    →  cohort .duckdb  (alt_bases + samples tables)
streamlit run app/geac_explorer.py  →  interactive browser
```

- **Rust + rust-htslib** for BAM/CRAM pileup processing
- **Apache Arrow + Parquet** for columnar per-sample storage
- **DuckDB (bundled)** for cohort-level SQL with no external database server
- **Streamlit + Altair** for the interactive explorer
- **conda + bioconda** for the htslib C library dependency
