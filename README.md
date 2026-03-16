# GEAC — Genomic Evidence Atlas of Cohorts

Collect alt base metrics from duplex/simplex BAM/CRAM files across a sequencing cohort.
Each sample is processed independently into a Parquet file; samples can then be merged
into a DuckDB database for cohort-level queries.

## Explorer UI

An interactive web app for browsing and visualising alt base data:

```bash
conda activate geac
marimo run app/geac_explorer.py
```

Then open `http://localhost:2718` in your browser. Enter a Parquet or DuckDB
file path in the text box to load data. Provides:

- Summary statistics (records, samples, positions, mean VAF)
- Filters: chromosome, sample, variant type, VAF range, min alt count
- Paginated data table
- Tabbed plots: VAF distribution, error spectrum, strand bias, overlap agreement

For interactive/notebook mode (editable cells):

```bash
marimo edit app/geac_explorer.py
```

## Setup

Requires [Miniforge](https://github.com/conda-forge/miniforge) (or Anaconda) and internet access.

```bash
bash scripts/setup.sh
conda activate geac
cargo build --release
```

The setup script:
1. Creates the `geac` conda environment with `htslib`
2. Installs conda activate hooks so `LIBRARY_PATH`/`PKG_CONFIG_PATH` are set automatically
3. Installs or updates Rust via `rustup`

## Usage

### Collect — process a single sample

```bash
geac collect \
  --input sample.bam \
  --reference hg38.fa \
  --sample-id SAMPLE_001 \
  --output SAMPLE_001.parquet \
  --read-type duplex \
  --pipeline fgbio
```

Optional flags:
- `--vcf calls.vcf.gz` — annotate loci with variant calling status
- `--min-base-qual 1` — minimum base quality (default: 1)
- `--min-map-qual 20` — minimum mapping quality (default: 20)
- `--region chr1:1000-2000` — restrict to a genomic region
- `--threads 4` — parallel processing

### Merge — combine samples into a cohort DuckDB

```bash
geac merge --output cohort.duckdb samples/*.parquet
```

### Query the cohort (DuckDB)

```sql
-- Cohort frequency of each alt allele
SELECT chrom, pos, ref_allele, alt_allele,
       COUNT(DISTINCT sample_id) AS n_samples,
       SUM(alt_count) AS total_alt_reads
FROM read_parquet('samples/*.parquet')
GROUP BY chrom, pos, ref_allele, alt_allele
ORDER BY n_samples DESC;
```

## Schema

Each Parquet file contains one row per alt allele observed at a locus.

| Column | Type | Description |
|---|---|---|
| `sample_id` | string | Sample identifier |
| `chrom` | string | Chromosome |
| `pos` | int64 | 0-based position |
| `ref_allele` | string | Reference allele |
| `alt_allele` | string | Alt allele |
| `variant_type` | string | SNV / insertion / deletion / MNV |
| `total_depth` | int32 | Total read depth at position |
| `alt_count` | int32 | Reads supporting the alt |
| `fwd_depth` | int32 | Forward strand depth |
| `rev_depth` | int32 | Reverse strand depth |
| `fwd_alt_count` | int32 | Forward strand alt reads |
| `rev_alt_count` | int32 | Reverse strand alt reads |
| `overlap_depth` | int32 | Reads in overlapping fragment pairs |
| `overlap_alt_agree` | int32 | Overlapping pairs where both reads see alt |
| `overlap_alt_disagree` | int32 | Overlapping pairs where reads disagree |
| `read_type` | string | raw / simplex / duplex |
| `pipeline` | string | fgbio / dragen / raw |
| `variant_called` | bool? | Whether a variant was called here (null if no VCF) |
| `variant_filter` | string? | PASS, filter reason, or null if not called |
