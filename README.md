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

### Homebrew (macOS arm64 — recommended)

```bash
brew install fleharty/geac/geac
```

This installs the `geac` binary plus `geac-cohort` and `geac-coverage-explorer` Streamlit launchers.

### Docker (linux/amd64 — Terra / cloud)

```bash
docker pull ghcr.io/fleharty/geac:latest
```

The Docker image contains only the `geac` binary (no Streamlit). It is intended for running `geac collect` on Terra or other cloud compute platforms.

### Local development (from source)

Requires Rust and htslib:

```bash
# macOS
brew install htslib pkg-config
cargo build --release

# Linux
# Build htslib from source (see .github/workflows/release.yml for the exact steps)
cargo build --release
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
| `--batch` | — | Batch/group label stored as a `batch` column (e.g. processing run name) |
| `--label1` | — | Free-text sample label 1 stored as `label1` column (e.g. tissue type) |
| `--label2` | — | Free-text sample label 2 stored as `label2` column (e.g. library prep method) |
| `--label3` | — | Free-text sample label 3 stored as `label3` column (e.g. sequencer type) |
| `--vcf` | — | Annotate loci with variant calling status from a VCF/BCF. Mutually exclusive with `--variants-tsv` |
| `--variants-tsv` | — | TSV variant list (columns: chrom, pos_start, pos_end, ref, var; 0-based). Alternative to `--vcf` |
| `--gnomad` | — | bgzip+tabix-indexed gnomAD VCF/BCF; adds `gnomad_af` float column (null = not in gnomAD) |
| `--gnomad-af-field` | `AF` | INFO field to use as allele frequency from the gnomAD VCF (e.g. `AF_joint`) |
| `--targets` | — | BED or Picard interval list of target regions; adds `on_target` bool column |
| `--gene-annotations` | — | GFF3, GTF, or UCSC genePred file; adds `gene` string column |
| `--repeat-window` | 10 | Bases on each side of locus to scan for homopolymers and STRs |
| `--min-base-qual` | 1 | Minimum base quality to count a read |
| `--min-map-qual` | 0 | Minimum mapping quality |
| `--include-duplicates` | off | Count PCR/optical duplicate reads (FLAG 0x400) |
| `--include-secondary` | off | Count secondary alignments (FLAG 0x100) |
| `--include-supplementary` | off | Count supplementary alignments (FLAG 0x800) |
| `--region` | whole genome | Restrict to a genomic region (e.g. `chr1:1-1000000`) |
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
artefacts via cycle number, or read-level phasing (e.g. MNV detection).

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

#### gnomAD allele frequency annotation

When `--gnomad` is provided, each alt allele record is annotated with:
- `gnomad_af` — the allele frequency from the gnomAD VCF's INFO/AF field (`null` if the
  exact allele is absent from gnomAD)

The file must be bgzip-compressed with a `.tbi` or `.csi` tabix index alongside it —
exactly the format gnomAD distributes. Chr-prefix mismatches between the BAM and the
gnomAD VCF (e.g. BAM uses `1`, gnomAD uses `chr1`) are handled automatically.

```bash
geac collect \
  --input sample.bam \
  --reference hg38.fa \
  --output SAMPLE_001.parquet \
  --gnomad gnomad.genomes.v4.vcf.gz
```

Use `--gnomad-af-field AF_joint` (or any other INFO key) to override the default `AF` field.

#### Sample labels

Three generic label columns (`label1`, `label2`, `label3`) let you attach free-text
metadata to every record at collection time:

```bash
geac collect \
  --input sample.bam \
  --reference hg38.fa \
  --output SAMPLE_001.parquet \
  --label1 "lung" \
  --label2 "KAPA HyperPrep" \
  --label3 "NovaSeq X"
```

Labels are stored as nullable strings alongside the existing `batch` column. They are
completely user-defined — use them for tissue type, library prep, sequencer, timepoint,
or any other per-sample dimension you want to filter or group by in the Explorer.

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

**Parquet files** are routed automatically by filename suffix — no extra flag is needed:

| Suffix | Table created |
|---|---|
| `.reads.parquet` | `alt_reads` — per-read detail records (from `--reads-output`) |
| `.normal_evidence.parquet` | `normal_evidence` — per-locus normal pileup evidence (from `geac annotate-normal`) |
| `.pon_evidence.parquet` | `pon_evidence` — per-locus PoN hit counts and VAFs (from `geac annotate-pon`) |
| `.coverage.parquet` | `coverage` — per-position coverage records (from `geac coverage`) |
| anything else | `alt_bases` — standard locus records |

Indices are created on each optional table for efficient joins back to `alt_bases`.

**DuckDB files** (`.duckdb`) can be passed directly alongside or instead of Parquet files.
Each known data table (`alt_bases`, `alt_reads`, `normal_evidence`, `pon_evidence`,
`coverage`) is copied from the source database into the output.  Inputs can be freely mixed:

```bash
# Combine two existing cohort databases
geac merge --output combined.duckdb cohort_a.duckdb cohort_b.duckdb

# Add new samples from Parquet into an existing database
geac merge --output updated.duckdb existing_cohort.duckdb new_sample.parquet

# Mix Parquet and DuckDB freely
geac merge --output cohort.duckdb batch1.duckdb batch2/*.parquet
```

The `samples` summary table is always rebuilt from the merged `alt_bases` at the end —
it is never copied from source DuckDB files so counts are always accurate.

A `geac_metadata` table is always written as a one-row database header, and a
`geac_inputs` table records one row per merged source artifact.  Together they
capture the merge tool version, schema version, command line, platform, input
counts, output row counts, and per-input file metadata:

```sql
SELECT * FROM geac_metadata;

SELECT * FROM geac_inputs;
```

The Explorer checks the database version at load time and warns if it differs from the
version it was built alongside. The sidebar `Advanced` expander also shows the current
`geac_metadata` header and the `geac_inputs` table for quick inspection without a manual
SQL query.

See [docs/provenance.md](docs/provenance.md) for the full provenance schema.

The output file must not already exist (use a new path or delete the old file first).

### Annotate Normal — cross-check tumor loci against a paired normal BAM

For each alt locus in the tumor Parquet, `geac annotate-normal` piles up the paired normal
BAM at that position and records how many fragments support each allele.  The result is a
`normal_evidence` Parquet that can be passed to `geac merge` so the Explorer can classify
loci as somatic candidates, germline-like, or artefacts.

```bash
geac annotate-normal \
  --tumor-parquet TUMOR.locus.parquet \
  --normal-bam    NORMAL.bam \
  --reference     hg38.fa \
  --output        TUMOR.normal_evidence.parquet
```

Optional flags:

| Flag | Default | Description |
|---|---|---|
| `--normal-sample-id` | from SM tag | Override the normal sample identifier |
| `--min-base-qual` | 1 | Minimum base quality to count a base in the normal pileup |
| `--min-map-qual` | 0 | Minimum mapping quality |
| `--include-duplicates` | off | Count duplicate reads in the normal |
| `--include-secondary` | off | Count secondary alignments |
| `--include-supplementary` | off | Count supplementary alignments |

Output naming convention: use `.normal_evidence.parquet` so `geac merge` routes the file
to the `normal_evidence` table automatically.

### Annotate PoN — cross-check tumor loci against a Panel of Normals DuckDB

`geac annotate-pon` queries a pre-built PoN DuckDB (produced by running `geac collect` on
each normal sample and then `geac merge`) to find how many PoN samples carry each tumor
alt allele and at what VAF — all via DuckDB analytics, with no BAM re-pileup.

```bash
geac annotate-pon \
  --tumor-parquet TUMOR.locus.parquet \
  --pon-db        pon.duckdb \
  --output        TUMOR.pon_evidence.parquet
```

The PoN DuckDB must contain an `alt_bases` table (produced by a standard `geac merge` run
on the normal cohort).  The output Parquet records, for each tumor alt locus, the number
of PoN samples that carry the same allele and the maximum and mean PoN VAF.

Output naming convention: use `.pon_evidence.parquet` so `geac merge` routes the file
to the `pon_evidence` table automatically.

### Coverage

`geac coverage` pileups a BAM/CRAM and emits per-position depth and GC content as a
Parquet file (`.coverage.parquet`).  When a targets BED or Picard interval list is
supplied, every target position is always emitted even if depth is zero; without
`--targets` only positions with at least one covering read are written.

```bash
geac coverage \
  --input     SAMPLE.bam \
  --reference hg38.fa \
  --output    SAMPLE.coverage.parquet \
  --targets   capture_targets.bed \
  --sample-id SAMPLE_001
```

Key options:

| Flag | Default | Description |
|------|---------|-------------|
| `--targets` | — | BED or Picard interval list; forces all target positions to be emitted |
| `--region` | — | Restrict to a genomic region (e.g. `chr1:1000-2000`) |
| `--gene-annotations` | — | GTF, GFF3, or UCSC genePred for gene/transcript annotation |
| `--sample-id` | SM tag | Override the sample ID stored in the output |
| `--batch` | — | Batch/group label stored as a column |
| `--read-type` | `duplex` | `duplex`, `simplex`, or `raw` |
| `--pipeline` | `fgbio` | `fgbio`, `dragen`, or `raw` |
| `--min-map-qual` | `0` | Minimum mapping quality |
| `--min-base-qual` | `20` | Minimum base quality |
| `--gc-window` | `100` | Window size (bp) for GC-content calculation |
| `--min-depth` | `0` | Only emit positions with depth ≥ this value |
| `--bin-size` | `1` | Merge consecutive positions into bins of this size |
| `--adaptive-depth-threshold` | — | Positions with depth below this value are emitted at single-base resolution (`bin_size=1`) and split any in-progress bin, preserving precision in low-coverage regions |

The output Parquet is routed to the `coverage` table by `geac merge` when its filename
ends in `.coverage.parquet`.

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

Two interactive Streamlit apps are included:

```bash
geac-cohort            # alt base / cohort explorer
geac-coverage-explorer   # per-position coverage explorer
```

Both are installed by `brew install fleharty/geac/geac`. Run either command from any directory — Streamlit opens a browser tab automatically. Enter a Parquet or DuckDB file path in the sidebar text box to load data.

For local development, run directly:

```bash
streamlit run /path/to/GEAC/app/geac_explorer.py
```

Then open `http://localhost:8501` in your browser.

Features:
- **Summary statistics** — alt records, samples, total alt bases, mean VAF, mean depth,
  % variant called
- **Sidebar filters** — chromosome, samples, variant type, VAF range, min alt count,
  variant called status, variant filter value (PASS / filter reason), min/max depth,
  on-target, gene name (partial match), homopolymer length, STR length;
  **Clear all filters** button resets all filters at once
- **Data table** — sortable, all schema columns, IGV session download button
- **Tabbed plots**
  - *VAF distribution* — separate histograms for SNV, insertion, deletion; click a bar
    to see matching records and download an IGV session
  - *Error spectrum* — SNV trinucleotide spectrum (SBS96) as a 3×2 grid of per-mutation-type
    panels with shared y-axis and fraction/count toggle; shift-click to select multiple
    contexts; drill-down table and IGV session. Optional COSMIC decomposition: provide a
    COSMIC SBS matrix path to overlay a reconstruction (black dots), show top-N signature
    exposures with etiology annotations, cosine similarity, and residual percentage. Also
    includes: per-sample signature exposure heatmap (DuckDB only); Called vs Uncalled
    comparison (butterfly chart + grouped signature bar, requires VCF annotation);
    VAF-stratified spectra (germline VAF > 30% vs somatic VAF ≤ 30%); family-size
    stratified spectra (singleton vs multi-member families, requires `--reads-output`);
    SBS96 heatmap across samples (DuckDB only)
  - *Strand bias* — forward vs. reverse alt reads scatter with 95% CI boundary lines;
    log1p or linear axis toggle; color by variant type, sample, on-target, or called status;
    click/shift-click to select points and view a drill-down table + IGV session
  - *Overlap agreement* — histogram of overlap concordance fractions
  - *Cohort* (DuckDB only) — per-sample summary table; VAF distribution overlay; strand
    balance scatter; alt loci count vs mean base quality scatter (outlier detection); SNV
    count bar chart stacked by SBS6 substitution type; click a sample row to focus all
    other views
  - *Reads* (DuckDB only, requires `--reads-output`) — family size histogram; read position
    bias (cycle number); mean base quality by cycle; insert size distribution with
    gap-correction toggle; insert size by allele frequency class; family size vs VAF
    scatter; mapping quality distribution; cohort artefact family size comparison (boxplot
    of family size by cohort frequency); all plots support aggregate / sample / batch
    color-by options
  - *Duplex/Simplex* (DuckDB only, requires `--reads-output`) — analyses focused on
    error-corrected sequencing: AB/BA strand balance distribution (`aD`/`bD` fgbio tags),
    read position bias by cycle, base quality distribution, family size vs VAF scatter,
    and insert size distribution; all gated on fgbio tag availability
  - *Tumor/Normal* (DuckDB only, requires `normal_evidence` table) — per-locus normal
    pileup summary joined to tumor alt loci; loci classified as Somatic candidate, Germline-like
    (normal VAF ≥ 20%), Artifact-like (normal VAF > 0%), No normal coverage, or No normal data;
    classification bar chart, tumor VAF vs normal VAF scatter, normal depth histogram, and
    data table
  - *Panel of Normals* (DuckDB only, requires `pon_evidence` table) — per-locus PoN hit
    summary; loci classified as PoN clean (not seen), Rare in PoN (< 10% of samples), or
    Common in PoN (≥ 10%); classification bar chart, tumor VAF vs PoN sample fraction
    scatter, max PoN VAF histogram, and data table sorted by PoN sample fraction
- **Per-read filters** (DuckDB only, requires `--reads-output`) — when an `alt_reads`
  table is present, a "Per-read filters" section appears in the sidebar with three
  range sliders, each with an include/exclude toggle:
  - *Family size* — filter by fgbio `cD` tag (total molecules per consensus read).
    Raising the minimum excludes singleton families that are likely PCR or sequencing
    errors. If a locus's alt count drops to zero after filtering, the locus is removed
    from the table entirely. This is the most useful filter for error-corrected data:
    a variant that disappears when singletons are excluded is almost certainly noise;
    one that holds up at family size ≥ 2 or 3 has stronger support.
  - *Cycle number* — filter by 1-based sequencing cycle (position within the read).
    Variants clustered at high cycle numbers (near the read end) are a common
    alignment artefact; lowering the upper bound removes these reads.
  - *Mapping quality* — filter by per-read MAPQ. Excluding low-MAPQ reads removes
    potential multi-mapping artefacts at repetitive loci.

  Two filter modes are available (controlled by the "Recompute alt count" checkbox):
  - **Locus-inclusion mode** (default) — loci where no reads pass the filter are
    hidden entirely. `alt_count` and `vaf` are unchanged for loci that remain.
  - **Re-aggregation mode** — `alt_count` is recomputed from reads passing the
    filter; an `original_vaf` column is shown alongside `vaf` for comparison.
    Loci where all reads fail show `alt_count = 0` but remain visible.

  In both modes, `ref_count`, `total_depth`, and strand/overlap columns always
  reflect the full pileup. Per-read filters are best used as an exploratory tool:
  do variants hold up under quality thresholds?

- **IGV integration** — provide a manifest TSV (`sample_id`, `bam_path`) in the sidebar
  to enable "Download IGV session" buttons throughout the app. Downloads a zip containing
  `session.xml` (BAM tracks + BED track) and `positions.bed` (one row per unique locus).
  Sessions are capped at 5 samples by default with an override option.

### Project config (geac.toml)

Place a `geac.toml` file in the directory where you run Streamlit (or pass `--config /path/to/geac.toml` after `--` on the command line) to pre-populate sidebar fields:

```toml
data             = "/path/to/cohort.duckdb"         # pre-fill the data file path
manifest         = "/path/to/manifest.tsv"           # pre-fill the manifest path
cosmic           = "/path/to/COSMIC_v3.4_SBS_GRCh38.txt"
genome_build     = "hg38"                            # hg19 | hg38 | mm10 | mm39 | <any IGV ID>
auto_launch_igv  = false                             # auto-load sessions into running IGV
target_regions   = "/path/to/targets.bed"            # optional track added to IGV sessions
gnomad_track     = "/path/to/gnomad.vcf.gz"          # optional gnomAD VCF/BCF track for IGV sessions
```

All keys are optional. `genome_build` accepts any IGV genome identifier — known values (`hg19`, `hg38`, `mm10`, `mm39`) are selected directly from the dropdown; anything else selects "other" and pre-fills the custom genome ID text box. `auto_launch_igv = true` checks the "Auto-launch IGV" checkbox by default, so every session is automatically sent to IGV via its REST API (port 60151) or launched as a subprocess if IGV is not already running.

If `gnomad_track` points to a `*.vcf.gz`, `*.vcf.bgz`, or `*.bcf`, the Explorer infers the
matching index path automatically (`.tbi` for VCF, `.csi` for BCF). The index file must still
exist and be readable by IGV.

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
| `batch` | string? | Batch/processing-group label (null if `--batch` not provided) |
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
| `cycle` | int32 | 1-based sequencing cycle at the alt position. Forward reads: `hard_clips_5prime + qpos + 1`; reverse reads: `hard_clips_5prime + read_length − qpos`. Hard-clipped bases at the 5′ end of synthesis are included so cycle reflects true polymerase position. |
| `read_length` | int32 | Stored sequence length in bases (hard-clipped bases excluded, soft-clipped bases included) |
| `is_read1` | bool | `true` if R1 (BAM flag `0x40`), `false` if R2 or unpaired |
| `ab_count` | int32? | fgbio `aD` tag: AB (top-strand) raw read count; null if tag absent |
| `ba_count` | int32? | fgbio `bD` tag: BA (bottom-strand) raw read count; null if tag absent |
| `family_size` | int32? | fgbio `cD` tag: total raw read count (`aD + bD` for duplex; sole count for simplex); null if tag absent |
| `base_qual` | int32 | Base quality at the alt position |
| `map_qual` | int32 | Mapping quality of the read |
| `insert_size` | int32? | SAM TLEN (template length / insert size); null when 0 (unpaired or mate unmapped) |

### Normal evidence table (`*.normal_evidence.parquet`)

Produced by `geac annotate-normal`. One row per (tumor locus × normal allele observed).
Always includes a NULL-allele anchor row to record normal depth even when no alt is seen.

| Column | Type | Description |
|---|---|---|
| `tumor_sample_id` | string | Tumor sample identifier |
| `chrom` | string | Chromosome |
| `pos` | int64 | 0-based position |
| `tumor_alt_allele` | string | Tumor alt allele being annotated |
| `normal_sample_id` | string | Normal sample identifier |
| `normal_alt_allele` | string? | Alt allele observed in the normal at this position (null = anchor/depth-only row) |
| `normal_depth` | int32 | Total fragment depth in the normal at this position |
| `normal_alt_count` | int32 | Fragments supporting `normal_alt_allele` (0 for anchor row) |

For SNV positions, one NULL anchor row is always written (capturing `normal_depth`), plus
one additional row for each non-reference base observed in the normal pileup.  For indel
positions, only the NULL anchor row is written.

### PoN evidence table (`*.pon_evidence.parquet`)

Produced by `geac annotate-pon`. One row per (tumor alt locus) with Panel of Normals hit
statistics derived from the PoN DuckDB.

| Column | Type | Description |
|---|---|---|
| `tumor_sample_id` | string | Tumor sample identifier |
| `chrom` | string | Chromosome |
| `pos` | int64 | 0-based position |
| `tumor_alt_allele` | string | Tumor alt allele being annotated |
| `n_pon_samples` | int64 | Number of PoN samples that carry this allele |
| `pon_total_samples` | int64 | Total number of samples in the PoN |
| `max_pon_vaf` | float64? | Highest VAF seen for this allele across PoN samples (null if `n_pon_samples = 0`) |
| `mean_pon_vaf` | float64? | Mean VAF across PoN samples that carry the allele (null if `n_pon_samples = 0`) |

## Docker

`linux/amd64` images are built automatically when a `v*.*.*` tag is pushed.
The image contains only the `geac` binary — it is intended for Terra and other
cloud compute platforms, not for running the Explorer.
Images are published to the GitHub Container Registry:

```
ghcr.io/fleharty/geac:<version>
ghcr.io/fleharty/geac:latest
```

### Pulling the image

```bash
docker pull ghcr.io/fleharty/geac:latest
# or a specific version:
docker pull ghcr.io/fleharty/geac:0.3.12
```

### Running geac on Terra

Set `docker_image` in your WDL inputs to `ghcr.io/fleharty/geac:<version>`.

### Running geac locally via Docker

```bash
docker run --rm \
    -v /path/to/data:/data \
    ghcr.io/fleharty/geac:latest \
    collect --input /data/sample.bam --reference /data/ref.fa --output /data/sample.parquet \
    --read-type duplex --pipeline fgbio
```

### Cutting a release

```bash
# 1. Bump version in Cargo.toml
# 2. Update GEAC_VERSION in app/geac_explorer.py and app/geac_coverage_explorer.py
# 3. Commit, push, then tag:
git tag v0.X.Y && git push origin v0.X.Y
```

The GitHub Actions workflow will:
1. Build and push the `linux/amd64` Docker image to ghcr.io
2. Build a native `macos-arm64` binary and attach it to the GitHub release
3. Update the Homebrew tap formula automatically

## WDL / Terra

Three WDL 1.0 workflows are provided in `wdl/`:

| Workflow | Purpose |
|---|---|
| `geac_collect.wdl` | Single-sample wrapper around `geac collect`; use this to scatter across a sample table |
| `geac_cohort.wdl` | Full cohort workflow: scatters `geac collect` then gathers with `geac merge` |
| `geac_merge.wdl` | Standalone merge — takes existing Parquets and builds a DuckDB |
| `geac_annotate_normal.wdl` | Single-sample wrapper around `geac annotate-normal`; cross-checks tumor loci against a paired normal BAM |
| `geac_annotate_pon.wdl` | Single-sample wrapper around `geac annotate-pon`; cross-checks tumor loci against a pre-built PoN DuckDB |

### `geac_collect.wdl` inputs

| Input | Type | Description |
|---|---|---|
| `input_bam` | File | BAM or CRAM file |
| `input_bam_index` | File | `.bai` or `.crai` index |
| `reference_fasta` | File | Reference FASTA |
| `reference_fasta_index` | File | `.fai` index |
| `read_type` | String | `duplex` / `simplex` / `raw` |
| `pipeline` | String | `fgbio` / `dragen` / `raw` |
| `docker_image` | String | e.g. `ghcr.io/fleharty/geac:0.3.12` |
| `sample_id` | String? | Override sample ID (default: BAM SM tag) |
| `vcf` | File? | VCF/BCF for variant call annotation |
| `vcf_index` | File? | `.tbi` or `.csi` index for VCF |
| `variants_tsv` | File? | TSV variant list (alternative to `--vcf`) |
| `targets` | File? | BED or Picard interval list for on-target annotation |
| `gene_annotations` | File? | GFF3, GTF, or UCSC genePred for gene annotation |
| `region` | String? | Restrict to a region, e.g. `chr1:1-1000000` |
| `repeat_window` | Int | Bases each side of locus for homopolymer/STR scan (default: 10) |
| `min_base_qual` | Int | Default: 1 |
| `min_map_qual` | Int | Default: 0 |
| `include_duplicates` | Boolean | Count duplicate reads (default: false) |
| `include_secondary` | Boolean | Count secondary alignments (default: false) |
| `include_supplementary` | Boolean | Count supplementary alignments (default: false) |
| `batch` | String? | Optional batch label stored in the output Parquet |
| `reads_output` | Boolean | Also write per-read detail Parquet (default: false) |
| `threads` | Int | Default: 1 |
| `memory_gb` | Int | Default: 8 |
| `disk_gb` | Int | Default: 100 |
| `preemptible` | Int | Default: 2 |

Outputs: `locus_parquet` (File) — per-sample locus Parquet; `reads_parquets` (Array[File]) — per-read Parquet (one element when `reads_output=true`, empty otherwise).

### `geac_cohort.wdl` inputs

Per-sample parallel arrays: `input_bams`, `input_bam_indices`, optional `sample_ids`,
optional `variants_tsvs`, optional `vcfs` + `vcf_indices` (per-sample VCF annotation).
Shared inputs applied to all samples: `reference_fasta`, `targets`, `gene_annotations`,
`region`, `repeat_window`, `read_type`, `pipeline`, `min_base_qual`, `min_map_qual`,
`include_duplicates`, `include_secondary`, `include_supplementary`, `batches`
(optional parallel array of per-sample batch labels), `threads`.

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

### `geac_annotate_normal.wdl` inputs

| Input | Type | Description |
|---|---|---|
| `tumor_parquet` | File | Locus Parquet from `geac collect` for the tumor sample |
| `normal_bam` | File | Normal BAM or CRAM |
| `normal_bam_index` | File | `.bai` or `.crai` index |
| `reference_fasta` | File | Reference FASTA |
| `reference_fasta_index` | File | `.fai` index |
| `docker_image` | String | geac Docker image |
| `normal_sample_id` | String? | Override normal sample ID (default: BAM SM tag) |
| `min_base_qual` | Int | Default: 1 |
| `min_map_qual` | Int | Default: 0 |
| `include_duplicates` | Boolean | Count duplicate reads (default: false) |
| `include_secondary` | Boolean | Count secondary alignments (default: false) |
| `include_supplementary` | Boolean | Count supplementary alignments (default: false) |
| `memory_gb` | Int | Default: 8 |
| `disk_gb` | Int | Default: 100 |
| `preemptible` | Int | Default: 2 |

Output: `normal_evidence_parquet` (File) — `{tumor_stem}.normal_evidence.parquet`.

### `geac_annotate_pon.wdl` inputs

| Input | Type | Description |
|---|---|---|
| `tumor_parquet` | File | Locus Parquet from `geac collect` for the tumor sample |
| `pon_db` | File | PoN DuckDB from `geac merge` on normal samples |
| `docker_image` | String | geac Docker image |
| `memory_gb` | Int | Default: 4 |
| `disk_gb` | Int | Default: 50 |
| `preemptible` | Int | Default: 2 |

Output: `pon_evidence_parquet` (File) — `{tumor_stem}.pon_evidence.parquet`.

### Running on Terra

1. Import the desired WDL into your Terra workspace.
2. Set `docker_image` to `ghcr.io/fleharty/geac:<version>` (e.g. `ghcr.io/fleharty/geac:0.3.12`).
3. For `geac_collect.wdl`: link `input_bam`, `input_bam_index`, `reference_fasta`, and
   `reference_fasta_index` to your workspace data table columns; Terra will scatter automatically.
4. For `geac_cohort.wdl`: provide parallel arrays directly and let the workflow scatter and merge.
5. To merge existing Parquets, use `geac_merge.wdl` with the list of Parquet files.

## Architecture

```
geac collect  →  per-sample .locus.parquet  [+ .reads.parquet with --reads-output]
                       │
                       ├──► geac annotate-normal  (paired normal BAM)
                       │         →  .normal_evidence.parquet
                       │
                       └──► geac annotate-pon  (PoN DuckDB)
                                 →  .pon_evidence.parquet

geac merge  →  cohort .duckdb
    alt_bases         (locus Parquets or existing .duckdb files)
    samples           (one-row-per-sample summary, always rebuilt)
    alt_reads         (.reads.parquet files, optional)
    normal_evidence   (.normal_evidence.parquet files, optional)
    pon_evidence      (.pon_evidence.parquet files, optional)
    coverage          (.coverage.parquet files, optional)

    # inputs can be mixed: Parquet files, .duckdb files, or both

geac-cohort  →  interactive alt base / cohort browser
geac-coverage-explorer  →  interactive coverage browser
```

- **Rust + rust-htslib** for BAM/CRAM pileup processing
- **Apache Arrow + Parquet** for columnar per-sample storage
- **DuckDB (bundled)** for cohort-level SQL with no external database server
- **Streamlit + Altair** for the interactive explorers
- **Homebrew** for macOS installation; **ghcr.io** Docker image for Terra/cloud
