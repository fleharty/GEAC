# GEAC Explorer Apps

GEAC ships two Streamlit launchers when installed through Homebrew:

```bash
geac-cohort
geac-coverage-explorer
```

For local development:

```bash
streamlit run app/geac_explorer.py
streamlit run app/geac_coverage_explorer.py
```

## Cohort Explorer

`geac-cohort` loads a merged DuckDB or compatible Parquet input and provides
interactive alt-base review.

Common sidebar filters include chromosome, samples, variant type, VAF, alt count,
depth, called status, filter value, target status, gene, repeat metrics, labels,
timepoint, and per-read filters when `alt_reads` is present.

Main tabs:

| Tab | Purpose |
|---|---|
| Summary | Cohort counts, sortable locus table, IGV session download, per-locus drill-down |
| VAF distribution | VAF histograms, depth views, and common-het depth retention |
| Error spectrum | SBS96/SBS6 spectra, COSMIC comparisons, NMF exploration |
| Strand bias | Forward/reverse alt support plots and drill-downs |
| Cohort | Per-sample summaries and cohort-level outlier views |
| Sample Identity | Germline-SNV fingerprinting for duplicates, swaps, and replicates |
| Reads | Read-level family size, cycle, mapping quality, insert-size, and N-context views |
| Duplex/Simplex | Error-corrected sequencing views gated on fgbio/DRAGEN tags |
| Tumor/Normal | Tumor loci joined to `normal_evidence` |
| Panel of Normals | Tumor loci joined to `pon_evidence` |
| Pipeline comparison | A/B comparison of the same sample under two `pipeline` values |
| Read-type comparison | Side-by-side comparison of read types in the same cohort |
| Overlap agreement | Mate-overlap concordance and discordance |
| MNV candidates | Cross-locus `fragment_id` joins for adjacent co-occurring SNVs |
| Fragmentomics | Insert-size, GC, end-motif, and nucleosome-style fragment views |

The Sample Identity tab uses the same core ideas as `geac sample-identity`: a common
SNV marker panel, pairwise fingerprint overlap, and genotype concordance. When
`subject_id` is populated it can separate expected same-subject matches from
unexpected same-individual candidates.

## Coverage Explorer

`geac-coverage-explorer` loads a merged DuckDB with a `coverage` table or a single
`.coverage.parquet` file.

Tabs:

| Tab | Purpose |
|---|---|
| Summary | Per-sample depth table, mean-depth bars, and QC fractions |
| Depth Distribution | Depth histograms and fraction-at-threshold summaries |
| GC Bias | Mean or normalized depth by GC content |
| Low Coverage | Loci below a configurable depth threshold across samples |
| Depth Profile | Linked depth, MAPQ, MAPQ0, and GC views across a gene or region |
| IGV | Embedded IGV.js viewer with GCS support through ADC token |
| Intervals | Per-interval and per-exon coverage summaries when `coverage_intervals` is present |

## `geac.toml`

Place `geac.toml` in the directory where Streamlit is launched, or pass a config path
after `--`:

```bash
streamlit run app/geac_explorer.py -- --config /path/to/geac.toml
```

Example:

```toml
data = "/path/to/cohort.duckdb"
manifest = "/path/to/manifest.tsv"
cosmic = "/path/to/COSMIC_v3.4_SBS_GRCh38.txt"
genome_build = "hg38"
auto_launch_igv = false
target_regions = "/path/to/targets.bed"
gnomad_track = "/path/to/gnomad.vcf.gz"
gnomad_track_index = "/path/to/gnomad.vcf.gz.tbi"
```

Relative paths are resolved against the config file directory. URI-style values such
as `gs://`, `http://`, and `https://` are preserved.

## IGV Manifest

```tsv
sample_id	bam_path	bai_path
SAMPLE_001	gs://bucket/SAMPLE_001.bam	gs://bucket/SAMPLE_001.bam.bai
SAMPLE_002	/local/SAMPLE_002.bam	/local/SAMPLE_002.bam.bai
```

`bai_path` is optional; if it is blank, IGV attempts to infer the index path.
