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
| `--reads-output` | off | Also write per-read detail Parquet (see below) |

#### Per-read detail output (`--reads-output`)

When `--reads-output` is set, `geac collect` writes two files instead of one:

- `{stem}.locus.parquet` — the standard locus table (same schema as a regular run)
- `{stem}.reads.parquet` — one row per alt-supporting read (fragment) at each locus

For example, `--output SAMPLE_001.parquet --reads-output` produces:
- `SAMPLE_001.locus.parquet`
- `SAMPLE_001.reads.parquet`

The reads table is linked to the locus table by `(sample_id, chrom, pos, alt_allele)`.

**When to use:** filtering by family size (fgbio duplex reads), diagnosing end-of-read
artefacts via `dist_from_read_end`, or read-level phasing (e.g. MNV detection).

When `geac merge` is given a mix of `.locus.parquet` and `.reads.parquet` files, it routes
them automatically: locus files → `alt_bases` table; reads files → `alt_reads` table.

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

#### Multi-nucleotide variants (MNVs)

GEAC processes one reference position at a time (standard pileup model). MNVs — adjacent
substitutions on the same haplotype, e.g. `AG→TC` — are therefore split into individual
SNV records, one per position. There is no way to distinguish a true MNV from two
independent SNVs at neighbouring positions using only the locus table. Identifying MNVs
requires read-level phasing: checking whether both substitutions appear on the same read.
This is not currently implemented in the Explorer. The per-read detail table (produced by
`--reads-output`) provides the data needed: join the locus table to the reads table and
check whether the same read supports substitutions at adjacent positions.

#### Soft-clipped bases

Soft-clipped bases are **not** counted. `rust-htslib` pileup only yields bases that are
aligned to the reference at each position; soft clips are excluded by design.

#### Overlap detection

For paired-end reads where the two mates overlap the same locus, GEAC detects the overlap
by grouping pileup reads at a position by query name. Depth is counted at the **fragment
level** — each overlapping pair contributes 1 to `total_depth` regardless of how many
reads cover the position. Strand (`fwd_depth` / `rev_depth`) is attributed using the R1
read's orientation (BAM flag `0x40`).

The following rules govern how each overlapping pair is tallied:

| Pair (read 1 + read 2) | `total_depth` | base tally | `overlap_alt_agree` / `overlap_alt_disagree` | `overlap_depth` |
|---|---|---|---|---|
| same base + same base (non-N) | +1 | that base +1 | agree +1 | +1 |
| alt + ref | +1 | alt +1 | disagree +1 | +1 |
| alt₁ + alt₂ (two different alts) | +1 | both +1 | disagree +1 each | +1 |
| alt + N | +1 | alt +1 | — | +1 |
| ref + N | +1 | ref +1 | — | +1 |
| N + N | +1 | — | — | +1 |

For non-overlapping singleton reads, N bases are excluded from all tallies entirely.
`overlap_alt_disagree` for `alt + ref` pairs records that the two mates disagreed, even
though the fragment is classified as alt — this is intentional, as the disagreement is
itself a useful quality signal.

### Merge — combine samples into a cohort DuckDB

```bash
geac merge --output cohort.duckdb samples/*.parquet
```

Creates a DuckDB database with:
- `alt_bases` — all per-sample locus records
- `samples` — one-row-per-sample summary (n_alt_loci, total_alt_reads, n_positions, etc.)
- Indices on `(chrom, pos)` and `sample_id` for fast queries

If any `.reads.parquet` files are present in the input list (produced by `--reads-output`),
`geac merge` also creates:
- `alt_reads` — all per-read detail records, linked to `alt_bases` by `(sample_id, chrom, pos, alt_allele)`
- Index on `(sample_id, chrom, pos, alt_allele)` for efficient joins

Files are routed automatically by filename suffix — no extra flag is needed.

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
  - *Cohort* (DuckDB only) — per-sample summary table; VAF distribution overlay; strand
    balance scatter; alt loci count vs mean base quality scatter (outlier detection); SNV
    count bar chart stacked by SBS6 substitution type; SBS96 heatmap (samples × 96
    contexts, normalised by sample); click a sample row to focus all other views
  - *Reads* (DuckDB only, requires `--reads-output`) — family size histogram; read position
    bias line chart; mean base quality by distance from read end; insert size distribution;
    family size vs VAF scatter; mapping quality distribution; cohort artefact family size
    comparison (boxplot of family size by cohort frequency); all plots support aggregate /
    sample / batch color-by options
- **Per-read filters** (DuckDB only, requires `--reads-output`) — when an `alt_reads`
  table is present, a "Per-read filters" section appears in the sidebar with three
  range sliders, each with an include/exclude toggle:
  - *Family size* — filter by fgbio `cD` tag (total molecules per consensus read).
    Raising the minimum excludes singleton families that are likely PCR or sequencing
    errors. If a locus's alt count drops to zero after filtering, the locus is removed
    from the table entirely. This is the most useful filter for error-corrected data:
    a variant that disappears when singletons are excluded is almost certainly noise;
    one that holds up at family size ≥ 2 or 3 has stronger support.
  - *Dist from read end* — filter by distance of the alt base from the end of the read.
    Variants clustered near read ends are a common alignment artefact; raising the
    lower bound removes these reads.
  - *Mapping quality* — filter by per-read MAPQ. Excluding low-MAPQ reads removes
    potential multi-mapping artefacts at repetitive loci.

  When any per-read filter is active, `alt_count` and `vaf` are re-aggregated from
  the reads table using only reads that pass the filters. An `original_vaf` column
  is shown alongside `vaf` so you can see the pre-filter allele frequency for
  comparison. Note that `ref_count`, `total_depth`, and strand/overlap columns are
  not recomputed — they always reflect the full pileup. Per-read filters are best
  used as an exploratory tool: do variants hold up under quality thresholds?

- **IGV integration** — provide a manifest TSV (`sample_id`, `bam_path`) in the sidebar
  to enable "Download IGV session" buttons throughout the app. Downloads a zip containing
  `session.xml` (BAM tracks + BED track) and `positions.bed` (one row per unique locus).
  Sessions are capped at 5 samples by default with an override option.

### Project config (geac.toml)

Place a `geac.toml` file in the directory where you run Streamlit (or pass `--config /path/to/geac.toml` after `--` on the command line) to pre-populate sidebar fields:

```toml
data             = "/path/to/cohort.duckdb"         # pre-fill the data file path
manifest         = "/path/to/manifest.tsv"           # pre-fill the manifest path
cosmic           = "/path/to/COSMIC_v3.4_SBS_GRCh37.txt"
genome_build     = "hg19"                            # hg19 | hg38 | mm10 | mm39 | <any IGV ID>
auto_launch_igv  = false                             # auto-load sessions into running IGV
```

All keys are optional. `genome_build` accepts any IGV genome identifier — known values (`hg19`, `hg38`, `mm10`, `mm39`) are selected directly from the dropdown; anything else selects "other" and pre-fills the custom genome ID text box. `auto_launch_igv = true` checks the "Auto-launch IGV" checkbox by default, so every session is automatically sent to IGV via its REST API (port 60151) or launched as a subprocess if IGV is not already running.

### Manifest format

```tsv
sample_id	bam_path	bai_path
SAMPLE_001	gs://my-bucket/bams/SAMPLE_001.bam	gs://my-bucket/bams/SAMPLE_001.bam.bai
SAMPLE_002	gs://my-bucket/bams/SAMPLE_002.bam	gs://my-bucket/bams/SAMPLE_002.bam.bai
SAMPLE_003	/local/path/to/SAMPLE_003.bam	/local/path/to/SAMPLE_003.bam.bai
```

`bai_path` is optional — if omitted or left blank, IGV will attempt to find the index automatically.

## Schema

### Locus table (`*.locus.parquet` / `*.parquet`)

Each file contains one row per alt allele observed at a locus.

| Column | Type | Description |
|---|---|---|
| `sample_id` | string | Sample identifier (from `--sample-id` or BAM SM tag) |
| `chrom` | string | Chromosome |
| `pos` | int64 | 0-based position |
| `ref_allele` | string | Reference allele |
| `alt_allele` | string | Alt allele (e.g. `T`, `+ACG`, `-2`) |
| `variant_type` | string | `SNV` / `insertion` / `deletion` |
| `total_depth` | int32 | Fragment depth at position (each overlapping pair counts as 1) |
| `alt_count` | int32 | Fragments supporting the alt allele |
| `ref_count` | int32 | Fragments supporting the reference allele |
| `fwd_depth` | int32 | Forward strand fragment depth (R1 orientation for overlapping pairs) |
| `rev_depth` | int32 | Reverse strand fragment depth |
| `fwd_alt_count` | int32 | Forward strand alt fragments |
| `rev_alt_count` | int32 | Reverse strand alt fragments |
| `fwd_ref_count` | int32 | Forward strand reference fragments |
| `rev_ref_count` | int32 | Reverse strand reference fragments |
| `overlap_depth` | int32 | Number of overlapping fragment pairs at this locus |
| `overlap_alt_agree` | int32 | Overlapping pairs where both mates support the alt |
| `overlap_alt_disagree` | int32 | Overlapping pairs where mates disagree (one alt, one ref or different alt) |
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

### Reads table (`*.reads.parquet`)

Produced when `--reads-output` is passed to `geac collect`. One row per alt-supporting
fragment at a locus. Linked to the locus table by `(sample_id, chrom, pos, alt_allele)`.

| Column | Type | Description |
|---|---|---|
| `sample_id` | string | Sample identifier |
| `chrom` | string | Chromosome |
| `pos` | int64 | 0-based position |
| `alt_allele` | string | Alt allele (links to locus table) |
| `dist_from_read_start` | int32 | 0-based index of the alt base within the read |
| `dist_from_read_end` | int32 | Bases from the alt to the end of the read (`read_length - dist_from_read_start - 1`) |
| `read_length` | int32 | Total length of the read in bases |
| `ab_count` | int32? | fgbio `aD` tag: AB (top-strand) raw read count; null if tag absent |
| `ba_count` | int32? | fgbio `bD` tag: BA (bottom-strand) raw read count; null if tag absent |
| `family_size` | int32? | fgbio `cD` tag: total raw read count (`aD + bD` for duplex; sole count for simplex); null if tag absent |
| `base_qual` | int32 | Base quality at the alt position |
| `map_qual` | int32 | Mapping quality of the read |
| `insert_size` | int32? | SAM TLEN (template length / insert size); null when 0 (unpaired or mate unmapped) |

## Docker

Multi-platform images (linux/amd64 + linux/arm64) are built automatically by the GitHub
Actions release workflow (`.github/workflows/release.yml`) when a `v*.*.*` tag is pushed.
Images are published to the GitHub Container Registry:

```
ghcr.io/fleharty/geac:<version>
ghcr.io/fleharty/geac:latest
```

### Pulling the image

```bash
podman pull ghcr.io/fleharty/geac:latest
# or a specific version:
podman pull ghcr.io/fleharty/geac:0.3.2
```

### Running the Explorer

Use the provided `run_explorer.sh` script (requires Podman):

```bash
./run_explorer.sh /path/to/data/directory
```

Then open `http://localhost:8501` in your browser. The data directory is mounted at `/data`
inside the container — place your `cohort.duckdb` and optional `geac.toml` there.

### Running geac collect / merge directly

```bash
podman run --rm \
    -v /path/to/data:/data \
    ghcr.io/fleharty/geac:latest \
    collect --input /data/sample.bam --reference /data/ref.fa --output /data/sample.parquet \
    --read-type duplex --pipeline fgbio
```

### Cutting a release

Releases are gated on changes that require rebuilding `cohort.duckdb` (i.e. Rust / data
pipeline changes). To release:

```bash
echo "0.X.Y" > VERSION
# also update version = "..." in Cargo.toml
git add VERSION Cargo.toml && git commit -m "Bump version to 0.X.Y"
git tag vX.Y.Z && git push origin vX.Y.Z
```

The GitHub Actions workflow builds native amd64 and arm64 images and merges them into a
multi-platform manifest automatically.

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
| `reads_output` | Boolean | Also write per-read detail Parquet (default: false) |
| `threads` | Int | Default: 4 |
| `memory_gb` | Int | Default: 8 |
| `disk_gb` | Int | Default: 100 |
| `preemptible` | Int | Default: 2 |

Outputs: `locus_parquet` (File) — per-sample locus Parquet; `reads_parquets` (Array[File]) — per-read Parquet (one element when `reads_output=true`, empty otherwise).

### `geac_cohort.wdl` inputs

Per-sample parallel arrays: `input_bams`, `input_bam_indices`, optional `sample_ids`,
optional `variants_tsvs`, optional `vcfs` + `vcf_indices` (per-sample VCF annotation).
Shared inputs applied to all samples: `reference_fasta`, `targets`, `gene_annotations`,
`region`, `repeat_window`, `read_type`, `pipeline`, `min_base_qual`, `min_map_qual`, `threads`.

Outputs: `locus_parquets` (Array[File]), `reads_parquets` (Array[File], empty when `reads_output=false`), and `cohort_db` (File, the merged DuckDB).

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
