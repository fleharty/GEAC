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
| `--vcf calls.vcf.gz` | — | Annotate loci with variant calling status from a VCF |
| `--min-base-qual` | 1 | Minimum base quality to count a read |
| `--min-map-qual` | 20 | Minimum mapping quality |
| `--region chr1:1-1000` | whole genome | Restrict to a genomic region |
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
- **Summary statistics** — record count, sample count, chromosome count, position count,
  total alt reads, mean VAF
- **Sidebar filters** — chromosome, samples (multi-select), variant type, VAF range,
  min alt count
- **Data table** — paginated, sortable, all schema columns
- **Tabbed plots**
  - VAF distribution (histogram by variant type)
  - SNV error spectrum (substitution counts)
  - Strand bias (forward vs reverse alt read scatter)
  - Overlap agreement fraction (histogram)

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
| `variant_called` | bool? | Whether a variant was called here (null if no VCF provided) |
| `variant_filter` | string? | VCF FILTER value (`PASS`, filter reason, or null) |

## Docker

The `docker/Dockerfile` performs a two-stage build: htslib and the `geac` binary are
compiled in a builder stage, then only the runtime libraries and binary are copied into
the final image, keeping it small. Images are built with [Podman](https://podman.io)
(no daemon required, drop-in compatible with Dockerfile syntax).

### Build and push to GCR

**Prerequisites:** Podman installed, and `gcloud` CLI authenticated.

```bash
# Authenticate Podman to GCR (once per session)
gcloud auth print-access-token | podman login -u oauth2accesstoken --password-stdin gcr.io

# Build locally (tagged with VERSION and latest)
bash scripts/build_docker.sh my-gcp-project

# Build and push in one step
bash scripts/build_docker.sh my-gcp-project --push
```

The version is read from the `VERSION` file in the repo root. Images are tagged as:
- `gcr.io/my-gcp-project/geac:0.1.0`
- `gcr.io/my-gcp-project/geac:latest`

When releasing a new version, update `VERSION`, rebuild, and push. Commit the `VERSION`
change so the git history and image tags stay in sync.

### Updating the version

```bash
echo "0.2.0" > VERSION
# also update version = "..." in Cargo.toml
bash scripts/build_docker.sh my-gcp-project --push
git add VERSION Cargo.toml && git commit -m "Bump version to 0.2.0"
```

## WDL / Terra

`wdl/geac_collect.wdl` is a single-sample WDL 1.0 workflow that wraps `geac collect`.
Run it on a sample set in Terra to scatter across all samples in parallel.

### Inputs

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
| `vcf` | File? | VCF for variant call annotation |
| `vcf_index` | File? | `.tbi` or `.csi` index |
| `region` | String? | Restrict to a region, e.g. `chr1:1-1000000` |
| `min_base_qual` | Int | Default: 1 |
| `min_map_qual` | Int | Default: 20 |
| `threads` | Int | Default: 4 |
| `memory_gb` | Int | Default: 8 |
| `disk_gb` | Int | Default: 100 |
| `preemptible` | Int | Default: 2 |

### Output

| Output | Type | Description |
|---|---|---|
| `parquet` | File | Per-sample alt base Parquet file |

### Running on Terra

1. Import `wdl/geac_collect.wdl` into your Terra workspace.
2. Set `docker_image` to your pushed GCR image (e.g. `gcr.io/my-project/geac:0.1.0`).
3. Link `input_bam`, `input_bam_index`, `reference_fasta`, and `reference_fasta_index`
   to your workspace data table columns.
4. Run on a sample set — Terra will scatter across samples automatically.
5. Collect the output `parquet` files and run `geac merge` to build the cohort DuckDB.

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
