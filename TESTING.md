# GEAC Manual Test Checklist

Work through each item top to bottom. Check off items as verified, note failures inline.

---

## High Priority — To Explore

- [ ] **IGV.js integration** — evaluate embedding [IGV.js](https://github.com/igvteam/igv.js) (the
  JavaScript port of IGV) directly inside the Explorer rather than generating session zip files for
  the desktop app. Could enable tighter integration: clicking a locus in any plot immediately renders
  the BAM reads inline without leaving the browser. Key questions to investigate: Streamlit
  compatibility (likely via `components.html` or a custom component), BAM/CRAM streaming
  requirements (igv.js can load remote files via URL or local files via a dev server), authentication
  for Terra-hosted BAMs, and whether the existing manifest/sample-picker workflow can feed into it.

---

## Known Issues / To Investigate

- [ ] **Insert size kink at ~250 bp** — both insert size plots show a kink (inflection/discontinuity)
  around 250 bp. Hypothesis: the kink occurs at exactly 2 × read_length = 2 × 131 = 262 bp (test
  data is 2×131 bp paired-end), which is the threshold where paired reads transition from overlapping
  (both reads cover the alt position → `[r1, r2]` branch in `tally_pileup`) to non-overlapping
  (only one read covers the position → `[r]` branch). This is a behavioral transition in the pileup
  logic, not necessarily a double-count, but worth confirming. To investigate: verify the kink is at
  exactly 262 bp, and inspect whether the per-read detail emission differs in any way between the two
  branches that could bias the insert size histogram near that threshold.

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

## Explorer — Sidebar: Per-read filters

- [ ] Per-read filters section only appears when `alt_reads` table is present
- [ ] Family size filter slider works; section hidden with caption when `cD` tag absent
- [ ] Family size include mode restricts reads to range (NULL rows excluded)
- [ ] Family size exclude mode removes reads within range (NULL rows kept)
- [ ] Dist from read end slider works in include and exclude modes
- [ ] Mapping quality (per-read) slider works in include and exclude modes
- [ ] Insert size slider appears only when `insert_size` data present
- [ ] At full range (20–500), caption reads "no filter active" and all reads pass through
- [ ] Narrowing the slider updates caption with active range description
- [ ] When filter is active, reads with `NULL` insert size (unpaired/mate-unmapped) are excluded
- [ ] Per-read filters re-compute `alt_count` and VAF; downstream plots and table reflect new values
- [ ] **Clear all** resets per-read filters to defaults

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

## Explorer — Reads plots

- [ ] Family size distribution histogram renders
- [ ] Read position bias plot renders
- [ ] Dist from read end histogram renders
- [ ] Base quality distribution renders
- [ ] Insert size distribution renders; color-by Sample/Batch works
- [ ] Insert size by AF class plot renders two lines (germline VAF > 30%, somatic VAF ≤ 30%)
- [ ] Insert size by AF class x-axis reflects current insert size slider bounds
- [ ] Family size vs VAF scatter renders
- [ ] Mapping quality distribution stacked bar renders (repetitive vs non-repetitive)
- [ ] Cohort artefact family size comparison renders (DuckDB only, ≥2 samples)

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
- [ ] SNV count bar chart renders with one bar per sample sorted by total count
- [ ] SBS6 stacking shows correct substitution type breakdown per sample
- [ ] SBS96 heatmap renders with samples as rows and 96 contexts as columns
- [ ] Heatmap color reflects normalized fraction (not raw count)
- [ ] Heatmap and SBS6 chart show "No SNVs" message when no data present
- [ ] Heatmap shows "Trinucleotide context unavailable" when column absent

---

## WDL workflows *(Terra — future)*

- [ ] `geac_collect.wdl` runs on Terra with a single sample
- [ ] `geac_cohort.wdl` scatters across multiple samples and produces DuckDB
- [ ] `geac_merge.wdl` merges existing Parquets into DuckDB
- [ ] Optional inputs (`vcf`, `targets`, `gene_annotations`, `region`) work when provided
