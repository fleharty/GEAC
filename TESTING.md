# GEAC Manual Test Checklist

Work through each item top to bottom. Check off items as verified, note failures inline.

---

## High Priority — To Explore

- [x] **VCF/TSV filter value sidebar filter** — add an Explorer sidebar control to filter loci by
  the `variant_filter` column (populated from VCF FILTER field or equivalent TSV column). The
  current `variant_called` selector only distinguishes called vs not-called; users also need to
  filter by specific filter values (e.g. show only `PASS`, or exclude `PASS` to focus on variants
  that failed a filter). Implementation: query the distinct `variant_filter` values present in the
  data and offer a multiselect (with an "All" default). Should be unavailable when `variant_filter`
  is entirely NULL (no VCF/TSV was provided at collect time), consistent with how other optional
  columns are handled.

- [ ] **Replicate and longitudinal comparison mode** — add a framework for comparing multiple
  Parquets derived from the same biological specimen, covering two related use cases:

  - *Technical replicates*: given N Parquets from the same sample (e.g. repeated library prep or
    sequencing runs), identify variants that appear consistently across replicates (high
    reproducibility, likely real signal) vs variably (likely artefact or stochastic noise). Key
    metrics: per-locus concordance rate, CV of VAF across replicates, fraction of replicates
    supporting the call.

  - *Longitudinal time points*: given Parquets from the same patient/specimen at multiple time
    points, track VAF trends over time. Key use cases: monitoring clonal evolution, treatment
    response, minimal residual disease (MRD), or clonal haematopoiesis dynamics. Would need a
    time-point ordering concept (e.g. `--timepoint` flag at collect time, stored as a column) and
    Explorer visualisations for VAF-over-time trajectories per locus.

  Possible implementation: a `geac compare` subcommand that takes a manifest of
  `(sample_id, replicate_id, timepoint)` groupings and produces a comparison Parquet/DuckDB with
  per-locus concordance and trend statistics. Explorer would add a "Comparison" tab.

- [ ] **Sensitivity / specificity evaluation mode** — add a benchmarking framework for evaluating
  pipeline performance when a truth set is available. Covers several truth scenarios:

  - *Well-characterised reference samples*: e.g. Genome in a Bottle (GIAB) or cell line standards
    where the true variant set is known. Compute TP/FP/FN/TN at each locus, stratified by VAF,
    depth, variant type, repeat context, etc.

  - *Spike-in experiments*: healthy donor samples with known germline variants spiked in at defined
    dilutions to simulate somatic variants at controlled VAFs. Truth = the spike-in variant list at
    the expected VAF tier.

  - *In silico dilutions*: mix two Parquets at defined ratios to simulate a somatic sample with
    ground truth.

  Key outputs: precision-recall curve, ROC curve, sensitivity/specificity at user-defined VAF
  thresholds, stratified by variant type / repeat context / family size / insert size. Possible
  implementation: a `geac benchmark` subcommand that takes a Parquet and a truth VCF/TSV and
  produces a benchmarking report. Explorer would add precision-recall and ROC plots to the
  comparison tab.

- [x] **Tumor-normal mode** — implemented as `geac annotate-normal`. Produces a
  `normal_evidence` Parquet that is routed to the `normal_evidence` table in `geac merge`.
  Explorer has a dedicated Tumor/Normal tab.

- [ ] **Replicate concordance mode** — given N Parquets from the same specimen (technical replicates
  or repeated library preps), measure per-locus reproducibility across replicates. Key metrics:
  fraction of replicates supporting each locus, CV of VAF across replicates, concordance rate.
  Variants that hold up across all replicates have much stronger support than those appearing in
  only one. Particularly powerful for error-corrected duplex data. Possible design: a `geac compare`
  subcommand that takes a manifest of `(sample_id, replicate_id)` groupings and produces a
  comparison Parquet with per-locus concordance statistics. Explorer would add a Replicates tab
  showing per-locus reproducibility plots and a concordance filter dimension.

- [ ] **Truth-set evaluation mode** — extends replicate mode when ground truth is known (spike-ins,
  cell line standards such as GIAB, or in silico dilutions). Given a truth VCF/TSV alongside one
  or more Parquets, compute TP/FP/FN/TN at each locus and generate precision-recall and ROC curves.
  Stratify sensitivity/specificity by VAF, family size, repeat context, trinucleotide context, and
  consequence type. Possible design: a `geac benchmark` subcommand that takes a Parquet and a truth
  VCF/TSV and produces a benchmarking report. Explorer would add precision-recall and ROC plots.
  Shares infrastructure with replicate mode — both involve multi-sample per-locus aggregation with
  a manifest-driven join. Truth evaluation is the natural extension once replicates are working.

- [x] **Duplex/Simplex consensus analysis tab** — implemented as the "Duplex/Simplex" Explorer
  tab. Shown only when `alt_reads` table is present. Includes AB/BA strand balance, read
  position bias by cycle, base quality distribution, family size vs VAF scatter, and insert
  size distribution.

- [ ] **IGV.js integration** — evaluate embedding [IGV.js](https://github.com/igvteam/igv.js) (the
  JavaScript port of IGV) directly inside the Explorer rather than generating session zip files for
  the desktop app. Could enable tighter integration: clicking a locus in any plot immediately renders
  the BAM reads inline without leaving the browser. Key questions to investigate: Streamlit
  compatibility (likely via `components.html` or a custom component), BAM/CRAM streaming
  requirements (igv.js can load remote files via URL or local files via a dev server), authentication
  for Terra-hosted BAMs, and whether the existing manifest/sample-picker workflow can feed into it.

---

## Known Issues / To Investigate

- [ ] **Insertions involving N bases in blacklist-filtered events** — when inspecting blacklist-filtered
  loci, insertions containing N bases appear in the output. Investigate whether this is expected
  behaviour (N-containing insertions are real soft-clipped or low-quality bases that happen to pass
  the current filters) or a bug (e.g. N bases should be excluded from alt alleles the same way N
  ref bases are skipped). Check the indel extraction logic in `tally_pileup` and confirm whether
  N-containing alt alleles should be suppressed at collect time.

- [ ] **SBS96 reconstruction overlay shifts between Fraction and Count** — the unified
  trinucleotide spectrum chart ("SNV Trinucleotide Spectrum — bars = observed, dots = reconstruction")
  changes subtly when toggling the Y-axis between Fraction and Count. The observed bars should be
  identical in shape (just rescaled), and the dots should track them exactly. Likely cause: the
  observed spectrum is normalised to `total_snvs` (sum of all observed counts) while the NNLS
  reconstruction is normalised to `recon_total` (sum of reconstructed counts), and these two
  denominators differ slightly — meaning the fraction-mode dot positions differ from the
  count-mode dot positions in relative terms. To investigate: confirm the denominators, then decide
  whether to normalise both to the same total (e.g. always use `total_snvs`) or accept the small
  discrepancy.

- [x] **Insert size kink at ~250 bp** — resolved by gap correction. The kink at 2 × read_length
  reflects the capture probability transition where fragments longer than 2R are only captured
  with probability 2R/L. The gap-corrected y-axis (default on) divides counts by this capture
  probability, removing the visual discontinuity.

---

## Rust CLI — To Explore / Future Features

- [ ] **Variant consequence annotation** — when `--gene-annotations` (ncbiRefSeq genePred) and
  `--reference` are both provided, classify each SNV by its coding consequence: synonymous,
  missense, nonsense (stop-gain), or start-loss. Implementation: extend `gene_annotations.rs`
  to store CDS exon intervals per transcript; for each SNV in a CDS, compute the cumulative
  CDS offset to find the codon position (offset % 3), fetch the codon from the reference FASTA,
  substitute the alt allele (reverse-complementing for minus-strand genes), translate both
  codons and compare. For multi-isoform genes, report the most severe consequence across
  transcripts. Store as a new `consequence` column in the Parquet (null for indels, UTR, and
  intergenic loci). The Explorer sidebar could then filter by consequence type, enabling direct
  comparison of missense vs synonymous error spectra.

- [ ] **Exon / intron feature type annotation** — when `--gene-annotations` is provided, also
  record the genomic feature type at each locus (exon, intron, UTR, CDS) as a new `feature_type`
  column in the Parquet. Would require storing per-feature-type intervals from the GTF/GFF3 at
  collect time rather than just gene name. A position exonic in any transcript of a gene should
  be classified as exonic (any-transcript priority rule). The Explorer sidebar could then offer a
  feature type filter alongside the existing gene filter, enabling direct comparison of exonic vs
  intronic error spectra — useful for separating true sequencing errors from mapping artefacts at
  intron/capture boundaries.

- [x] **Duplicate / secondary / supplementary read filtering** — implemented via
  `--include-duplicates`, `--include-secondary`, `--include-supplementary` flags on both
  `geac collect` and `geac annotate-normal`. All three classes are excluded by default.

- [ ] **`--min-alt-count` filter at collect time** — currently every position with even one alt read
  is written to the Parquet. For WGS data this produces very large outputs. A minimum alt count
  threshold at collect time would dramatically reduce output size with minimal loss of signal.
  Consider defaulting to 2 or 3 to suppress single-read noise.

- [ ] **Depth cap / downsampling** — at extremely deep positions (e.g. hotspots in amplicon data)
  pileup processing is slow and depth statistics become noisy. Consider a `--max-depth` flag that
  randomly downsamples reads at a position to the cap before tallying, consistent with how samtools
  and other tools handle this.

- [ ] **`--pipeline dragen` is currently cosmetic** — the pipeline flag is stored as metadata but
  processing always looks for fgbio tags (`aD`, `bD`, `cD`). DRAGEN consensus reads use different
  tags. If DRAGEN consensus support is needed, the tag lookup logic should branch on `pipeline`.

- [ ] **Read strand in alt_reads** — `is_reverse` and `is_first_in_pair` are not stored in the
  alt_reads Parquet. Locus-level strand counts exist (fwd/rev alt count), but per-read strand
  information would enable finer-grained strand bias analysis in the Explorer's reads tab.

- [ ] **NM tag (edit distance) in alt_reads** — storing the SAM `NM` tag (number of mismatches)
  per alt-supporting read would help distinguish reads with many mismatches (likely mapping
  artefacts) from clean alt-supporting reads. Useful as an additional per-read filter in the
  Explorer.

- [ ] **`qc` command enhancements** — current QC output covers SBS6 substitution breakdown but
  lacks: indel composition, on-target fraction, duplicate/secondary rate, and per-chromosome stats.
  These would make `geac qc` a more complete QC report for pipeline validation.

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
- [ ] `--batch` label stored in output Parquet `batch` column
- [ ] `--region` restricts output to the specified region only
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

### Read filtering flags
- [ ] Default behaviour excludes duplicate, secondary, and supplementary reads
- [ ] `--include-duplicates` causes duplicate reads (FLAG 0x400) to be included
- [ ] `--include-secondary` causes secondary alignments (FLAG 0x100) to be included
- [ ] `--include-supplementary` causes supplementary alignments (FLAG 0x800) to be included
- [ ] Combining `--include-duplicates --include-secondary --include-supplementary` includes all three classes
- [ ] Depth counts differ when flags are included vs excluded for a BAM known to contain such reads

### CRAM support
- [ ] Runs on a CRAM file with `--reference`

---

## CLI — `geac merge`

- [ ] Merges multiple Parquets into a `.duckdb` file
- [ ] Output DuckDB contains `alt_bases` and `samples` tables
- [ ] Fails with a clear error if output file already exists
- [ ] Handles Parquets with different optional annotation columns (`union_by_name`)
- [ ] `.reads.parquet` files routed to `alt_reads` table; locus Parquets still go to `alt_bases`
- [ ] `.normal_evidence.parquet` files routed to `normal_evidence` table
- [ ] `.pon_evidence.parquet` files routed to `pon_evidence` table
- [ ] Mix of all four suffix types in one `geac merge` call produces all four tables
- [ ] Fails with clear error when only special-suffix files provided and no locus Parquets

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

## CLI — `geac annotate-normal`

- [ ] Runs against a tumor Parquet and normal BAM; produces a `.normal_evidence.parquet` file
- [ ] Normal sample ID read from BAM SM tag when `--normal-sample-id` is omitted
- [ ] `--normal-sample-id` overrides SM tag correctly
- [ ] NULL anchor row always present for each tumor locus (captures normal depth)
- [ ] Non-reference bases observed in normal generate additional rows with `normal_alt_allele` set
- [ ] `normal_depth = 0` rows emitted for tumor loci not covered by normal BAM
- [ ] SNV positions produce NULL anchor + per-allele rows; indel positions produce anchor row only
- [ ] `--min-base-qual` and `--min-map-qual` affect which normal reads are counted
- [ ] `--include-duplicates` / `--include-secondary` / `--include-supplementary` change depth counts
- [ ] Output passed to `geac merge` lands in `normal_evidence` table

---

## CLI — `geac annotate-pon`

- [ ] Runs against a tumor Parquet and PoN DuckDB; produces a `.pon_evidence.parquet` file
- [ ] `n_pon_samples = 0` for tumor loci absent from the PoN
- [ ] `max_pon_vaf` and `mean_pon_vaf` are null when `n_pon_samples = 0`
- [ ] `pon_total_samples` matches the distinct sample count in the PoN's `alt_bases` table
- [ ] Output passed to `geac merge` lands in `pon_evidence` table
- [ ] Fails with clear error if PoN DuckDB does not contain an `alt_bases` table

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
- [ ] Cycle number slider works in include and exclude modes
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
- [ ] Read position bias (cycle number) plot renders
- [ ] Cycle number distribution histogram renders
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

## Explorer — Signatures: Future Ideas

- [x] **Count / fraction toggle on SBS96** — add a count/fraction y-axis toggle to the SBS96 spectrum,
  consistent with the insert size and other plots. Essential for comparing samples with different
  read depths without being misled by raw count differences.

- [x] **Reconstructed spectrum overlay** — after COSMIC NNLS fitting, overlay the reconstructed
  spectrum as dots on top of the observed SBS96 bars. Refit uses only the top-N signatures selected
  by the slider; cosine similarity and residual % reflect the refit quality.

- [ ] **Consolidate Error Spectrum and Sig. Comparison tabs** — the two tabs currently contain
  overlapping and duplicated plots (e.g. SBS96 spectrum appears in both). Merge them into a single
  tab with a clear internal layout: observed spectrum first, then COSMIC decomposition (signature
  exposures + reconstructed overlay), then the per-sample comparison plots currently living in Sig.
  Comparison. Eliminates duplication and makes the signature analysis workflow linear and
  self-contained.

- [ ] **Family-size stratified spectrum** — split the SBS96 into singleton reads (family_size = 1)
  vs multi-member families. Singletons are enriched for sequencing errors; the difference between
  the two spectra isolates the true variant signal from the error process. Particularly useful for
  error characterisation in error-corrected sequencing workflows.

- [ ] **VAF-stratified spectrum** — show germline (VAF > 30%) vs somatic (VAF ≤ 30%) SBS96 side by
  side, analogous to the insert size by AF class plot. If the error spectrum differs between the
  two, that reveals what is driving low-VAF calls.

- [ ] **Per-sample COSMIC decomposition (cohort stacked bar)** — run NNLS on each sample
  independently and display a stacked bar chart of signature exposures across samples (samples on
  x-axis, signatures stacked). This is the standard SigProfiler-style cohort plot and is especially
  valuable in DuckDB/cohort mode.

- [ ] **Reference trinucleotide frequency normalisation** — divide each context count by its
  expected frequency in the genome or target region to remove background trinucleotide composition
  bias. Without this, a C>T spike could simply reflect CpG abundance in the target rather than a
  real mutagenic process. Important for small-panel error spectrum work.

- [ ] **Indel spectrum (COSMIC ID signatures)** — extend the tab to show an indel context spectrum
  (repeat unit type, repeat length, indel length) and fit it against COSMIC ID signatures via NNLS.
  Would make the tab a full mutational signature suite covering SBS, and ID.

- [ ] **Strand asymmetry overlay** — add a per-context indicator on the SBS96 bars showing
  forward/reverse strand imbalance. Transcription-coupled repair and replication strand effects
  leave characteristic asymmetries that are invisible in the standard pyrimidine-collapsed SBS96.

- [ ] **Per-variant signature attribution** — after NNLS, compute each variant's most likely
  signature based on its SBS96 context weight vector and the fitted exposures. Add a
  `top_signature` column to the drill-down table so individual variants can be traced back to a
  mutagenic process.

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

## Explorer — Tab: Tumor/Normal

- [ ] Tab shows info message when `normal_evidence` table is absent
- [ ] Tab shows data when `normal_evidence` table is present
- [ ] Each locus correctly classified: Somatic candidate / Germline-like / Artifact-like / No normal coverage / No normal data
- [ ] Classification bar chart renders with correct counts per category
- [ ] Tumor VAF vs normal VAF scatter renders; color by classification works
- [ ] Normal depth histogram renders
- [ ] Data table expander shows joined data sorted correctly

---

## Explorer — Tab: Panel of Normals

- [ ] Tab shows info message when `pon_evidence` table is absent
- [ ] Tab shows data when `pon_evidence` table is present
- [ ] Summary metrics row: PoN total samples, clean count, rare count, common count, no-data count
- [ ] Each locus correctly classified: PoN clean / Rare in PoN / Common in PoN
- [ ] Classification bar chart renders
- [ ] Tumor VAF vs PoN sample fraction scatter renders; color by classification works
- [ ] Max PoN VAF histogram renders (only loci with PoN evidence)
- [ ] Data table expander sorted by `pon_sample_fraction` descending

---

## WDL workflows *(Terra — future)*

- [ ] `geac_collect.wdl` runs on Terra with a single sample
- [ ] `geac_cohort.wdl` scatters across multiple samples and produces DuckDB
- [ ] `geac_merge.wdl` merges existing Parquets into DuckDB
- [ ] Optional inputs (`vcf`, `targets`, `gene_annotations`, `region`) work when provided
- [ ] `geac_annotate_normal.wdl` runs on Terra with tumor Parquet + normal BAM
- [ ] `geac_annotate_pon.wdl` runs on Terra with tumor Parquet + PoN DuckDB
