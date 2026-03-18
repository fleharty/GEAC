# GEAC Manual Test Checklist

Work through each item top to bottom. Check off items as verified, note failures inline.

---

## Testing Sessions

| Date | Version | Tester | Notes |
|---|---|---|---|
| 2026-03-18 | v0.1.1 | MF | Initial manual test session |

---

## CLI — `geac collect`

### Basic operation
- [x] Runs on a BAM with no optional flags; produces valid Parquet
- [x] DuckDB query against output Parquet returns valid rows
- [x] Reads sample ID from BAM SM tag when `--sample-id` is omitted
- [x] `--sample-id` overrides SM tag correctly
- [x] `--read-type` and `--pipeline` values stored correctly in output
- [ ] `--region` restricts output to the specified region only
- [ ] `--threads > 1` runs without error
- [ ] `--progress-interval 0` suppresses progress output

### Annotations
- [ ] `--vcf` populates `variant_called` and `variant_filter` columns
- [ ] `--variants-tsv` populates same columns (mutually exclusive with `--vcf`)
- [ ] `--targets` (BED) populates `on_target` correctly
- [ ] `--targets` (Picard interval list with `@` header) populates `on_target` correctly
- [ ] `--gene-annotations` (GTF) populates `gene` column; intergenic loci are null
- [ ] `--gene-annotations` (GFF3) works
- [ ] `--repeat-window` changes homopolymer/STR detection window
- [ ] `homopolymer_len`, `str_period`, `str_len` columns are populated
- [ ] `trinuc_context` populated for SNVs; null for indels/MNVs

### CRAM support
- [ ] Runs on a CRAM file with `--reference`

---

## CLI — `geac merge`

- [ ] Merges multiple Parquets into a `.duckdb` file
- [ ] Output DuckDB contains `alt_bases` and `samples` tables
- [ ] Fails with a clear error if output file already exists
- [ ] Handles Parquets with different optional annotation columns (`union_by_name`)

---

## CLI — `geac qc`

- [ ] Prints per-sample QC report to stdout
- [ ] `--output` writes TSV alongside stdout report
- [ ] `--on-target-only` restricts metrics to on-target loci

---

## CLI — `geac cohort`

- [ ] Produces TSV output listing recurrent loci
- [ ] Produces Parquet output when `--output` ends in `.parquet`
- [ ] `--min-samples` threshold filters correctly
- [ ] `--min-sample-fraction` threshold filters correctly
- [ ] `--top-n` controls number of loci printed to stdout
- [ ] `--on-target-only` restricts to on-target loci

---

## Explorer — File loading

- [ ] Loads a single `.parquet` file
- [ ] Loads a `.duckdb` file (cohort mode)
- [ ] Shows error on invalid path

---

## Explorer — Summary statistics

- [ ] Overall and filtered metric rows both display correctly
- [ ] Metrics update when filters change

---

## Explorer — Sidebar filters

- [ ] Chromosome filter restricts data
- [ ] Sample multiselect restricts data; blank = all samples
- [ ] Gene text filter (partial match) works; unavailable when column absent
- [ ] Variant type multiselect works
- [ ] VAF range slider works
- [ ] Min alt count filter works
- [ ] Min fwd alt count filter works
- [ ] Min rev alt count filter works
- [ ] Min overlap alt agree filter works
- [ ] Min overlap alt disagree filter works
- [ ] Variant called filter (Yes / No / Unknown) works
- [ ] Target bases filter (On target / Off target) works; unavailable when column absent
- [ ] Homopolymer length range filter works; unavailable when column absent
- [ ] STR length range filter works; unavailable when column absent
- [ ] Min depth / Max depth filters work
- [ ] Display limit changes row count in data table
- [ ] **Clear all** resets every filter to default

---

## Explorer — Data table

- [ ] Table displays filtered records with correct columns
- [ ] Clicking a row triggers position drill-down
- [ ] Drill-down shows all samples at that locus

---

## Explorer — IGV integration

- [ ] Manifest loads correctly from default path
- [ ] Manifest loads from custom path
- [ ] Warning appears when >5 samples in selection
- [ ] Sample picker appears and restricts BED to selected samples only
- [ ] "Load all N samples" override works
- [ ] Warning appears when sample not in manifest
- [ ] "Prepare IGV session" generates a downloadable zip
- [ ] Zip contains `session.xml` and `positions.bed`
- [ ] `positions.bed` only contains positions from selected samples
- [ ] Genome selector (hg19, hg38, mm10, mm39, other) sets correct genome in session XML

---

## Explorer — Tab 1: VAF Distribution

- [ ] Histograms render for each variant type present in data
- [ ] Clicking a bar shows matching records table
- [ ] IGV session download works from drill-down

---

## Explorer — Tab 2: Error Spectrum

- [ ] SNV trinucleotide spectrum (SBS96) renders as 3×2 grid
- [ ] All 96 bars shown (16 contexts × 6 mutation types)
- [ ] Shared y-axis across all 6 panels
- [ ] Labels show full context e.g. `A[C>A]A`
- [ ] Shift-click selects multiple contexts
- [ ] Drill-down table reflects all selected contexts
- [ ] IGV session download works from multi-select
- [ ] COSMIC signature decomposition loads when matrix path provided
- [ ] Cosine similarity and residual % displayed
- [ ] Top-N signatures slider changes chart (colors stable)
- [ ] Etiology annotations shown in table and tooltip

---

## Explorer — Tab 3: Strand Bias

- [ ] Scatter plot renders with diagonal and 95% CI boundary lines
- [ ] Linear axis mode works
- [ ] log1p axis mode works; tick labels show linear values
- [ ] Color by: Variant type works
- [ ] Color by: Sample works
- [ ] Color by: On target shows green/red (when column present)
- [ ] Color by: Called variant shows green/red (when column present)
- [ ] Click selects a point; drill-down table appears
- [ ] Shift-click selects multiple points
- [ ] IGV session download works from selection

---

## Explorer — Tab 4: Overlap Agreement

- [ ] Histogram renders when overlap data present
- [ ] "No overlapping fragments" message shown when absent

---

## Explorer — Tab 5: Cohort (DuckDB only)

- [ ] Tab only visible when loading `.duckdb` file
- [ ] Per-sample summary table shows all samples regardless of sample filter
- [ ] Clicking a sample row enables "Filter all tabs to {sample}" button
- [ ] Focus button sets sample filter for all other tabs
- [ ] VAF distribution overlay shows one curve per sample
- [ ] Strand balance scatter shows one dot per sample with dashed reference at 0.5

---

## WDL workflows *(Terra — future)*

- [ ] `geac_collect.wdl` runs on Terra with a single sample
- [ ] `geac_cohort.wdl` scatters across multiple samples and produces DuckDB
- [ ] `geac_merge.wdl` merges existing Parquets into DuckDB
- [ ] Optional inputs (`vcf`, `targets`, `gene_annotations`, `region`) work when provided
