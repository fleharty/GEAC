# GEAC — Development Log

Historical record of completed work and design notes that shaped the project.
The active backlog lives in [TODO.md](../TODO.md); the milestone roadmap lives in
[ROADMAP.md](../ROADMAP.md). Non-obvious bugs and multi-attempt fixes are tracked
separately in [CHALLENGES.md](../CHALLENGES.md).

This file is append-mostly. New entries go at the bottom of their section. The goal
is to preserve the *why* behind decisions that aren't obvious from the code.

---

## Original "Coverage complete" plan (v0.4.0)

Captured 2026-04-01. All three gating items shipped in v0.4.0.

1. **Gene annotation extensions** — `feature_type` and `exon_number` added to
   GTF/GFF3/genePred lookup; propagated through coverage records and Parquet schema.
2. **BEDGraph annotation tracks** — `src/track.rs` with `AnnotationTrack`/`TrackSet`;
   binary-search lookup; chr-prefix bridging; `--track NAME:FILE` repeatable flag;
   dynamic Arrow schema columns in coverage Parquet.
3. **Read-type comparison view** — `tab_read_type` in Explorer: locus concordance
   tiles + stacked bar, VAF density overlay, VAF correlation scatter, strand balance
   density, SBS96 side-by-side, unique-loci table.

---

## Rust / CLI — completed

- `geac qc` subcommand — per-sample summary of error rates by substitution type,
  strand bias metrics, and overlap concordance.
- `geac cohort` subcommand — per-locus artifact frequencies across samples (flag
  positions seen in N% of samples); outputs TSV or Parquet.
- On-target annotation — `--targets` BED/Picard interval list flag; records
  `on_target bool?` column in Parquet.
- Gene annotation — `--gene-annotations` flag; accepts GFF3, GTF, or UCSC genePred
  (.txt/.txt.gz); records `gene string?` column.
- Locus repetitiveness metrics — `homopolymer_len`, `str_period`, `str_len` columns;
  `--repeat-window` flag (default 10 bp).
- Trinucleotide context — `trinuc_context` column computed from reference at each
  SNV locus.
- Variant annotation — `--vcf` and `--variants-tsv` flags; annotates `variant_called`
  / `variant_filter` columns.
- Fragment overlap metrics — `overlap_alt_agree`, `overlap_alt_disagree`,
  `overlap_ref_agree` columns for read-pair concordance.
- Integration tests — synthetic-BAM end-to-end suite under `tests/integration.rs`.
- Audit `alt_count` double-counting — resolved: `total_depth`, `alt_count`, and
  `ref_count` are now fragment-level counts. Each overlapping pair contributes 1 to
  `total_depth` regardless of how many reads cover the position. See `tally_pileup`
  doc comment for full classification rules.
- Fix N-base handling in overlap tally — in `tally_pileup` (`src/bam/mod.rs`), if
  one read of an overlapping pair has an `N` at the position, `overlap_alt_disagree`
  is incorrectly incremented for the other read's alt base. An `N` is uninformative
  and should not count as a disagreement. Fix: skip the overlap agreement/disagreement
  logic when either base is `N`, and exclude `N` bases from `total_depth` and alt
  tallies entirely.
- Re-examine N-base handling before v0.3.0 — resolved as part of the fragment-level
  depth overhaul. See `tally_pileup` doc comment for the full classification table
  including N cases.

---

## Per-read detail table (two-table design) — shipped v0.3.2

**Motivation:** Read-end proximity and family size are inherently per-read properties.
Aggregating them into the locus table (e.g. `mean_dist_from_end`) loses distributional
information needed for principled artifact filtering. The most correct design stores
per-read detail in a second table linked to the existing locus table.

### Design

Two output files per sample from `geac collect`:
- `{sample}.parquet` — existing locus-level table (one row per alt locus per sample); unchanged
- `{sample}.reads.parquet` — new per-read table (one row per alt-supporting read); columns:
  - `sample_id`, `chrom`, `pos`, `alt_allele` — foreign key back to locus table
  - `dist_from_read_start` — 0-based position of the alt base within the read
  - `dist_from_read_end` — distance from the alt base to the 3' end of the read
  - `read_length` — total read length after soft-clipping
  - `ab_count` — top-strand family size (fgbio `cD` tag or equivalent)
  - `ba_count` — bottom-strand family size
  - `family_size` — ab_count + ba_count
  - `base_qual` — base quality at the alt position
  - `map_qual` — mapping quality of the read
  - `insert_size` — SAM TLEN (insert size); null when 0 (unpaired / mate unmapped)

### Pros
- Enables principled per-read filtering (e.g. "alt reads where family_size >= 3 AND dist_from_read_end > 10")
- Joint filters are natural SQL: no need to pre-compute every possible aggregate
- Supports future analyses not yet anticipated
- Locus table remains unchanged — no migration of existing Parquet files
- DuckDB handles multi-table joins efficiently with Parquet pushdown
- Per-read drill-down in Explorer at a specific locus becomes very rich
- Only alt-supporting reads are stored, so table is much smaller than total read count

### Cons
- `geac collect` must write two files instead of one — more complex output handling
- WDL workflows need to propagate and merge both file types
- `geac merge` needs a second merge step for the reads table
- Explorer must be aware of the optional reads table (graceful fallback if absent)
- Larger total storage footprint per sample
- Adds implementation complexity to the Rust BAM processing loop

### Implementation steps (all completed)
- Step 1: Define `AltRead` struct in `src/record.rs` with the columns listed above.
- Step 2: Populate `AltRead` records during BAM pileup in `src/bam/mod.rs`; parse
  fgbio/DRAGEN family size tags (`cD`/`cE` or `RX`/`MI`) per pipeline.
- Step 3: Add `src/writer/parquet_reads.rs` — write `AltRead` records to
  `{stem}.reads.parquet`.
- Step 4: Update `geac collect` CLI to emit both files; add `--reads-output` flag
  (optional; if omitted, reads table is not written).
- Step 5: Update `geac merge` to accept and merge reads Parquets into a second
  DuckDB table (`alt_reads`) alongside `alt_bases`.
- Step 6: Update WDL workflows to handle the optional reads Parquet output from
  `Collect` and pass it to `Merge`.
- Step 7: Explorer — in the position drill-down, JOIN `alt_reads` on
  `(sample_id, chrom, pos, alt_allele)` to show per-read detail (dist from end,
  family size, base qual) when reads table is present.
- Step 8: Explorer — sidebar filters for `family_size`, `dist_from_read_end`, and
  `map_qual` with include/exclude toggles; `alt_count` and `vaf` re-aggregated from
  reads when filters are active.
- `alt_reads` schema v2: `is_read1` + `cycle` rename — bundled two breaking schema
  changes into a single re-collect:
  - **Add `is_read1`** — `is_first_in_pair` (BAM flag 0x40) is tracked internally
    during pileup but never written to `AltRead`. Adding it as a boolean column
    enables R1/R2-stratified artefact analysis (e.g. R2-biased substitution patterns).
  - **Rename `dist_from_read_start` → `cycle`** (1-based: `qpos + 1`) and **drop
    `dist_from_read_end`** (derivable as `read_length - cycle`). Unifies the column
    name with the "Cycle number" label used in the Explorer sidebar and Reads tab.
  - **Invert sidebar filter direction**: change from "min dist_from_read_end" to
    "max cycle" slider (same artefact-rejection intent; cycle > threshold excludes
    end-of-read reads).
  - Changes spanned: `AltRead` struct (`src/record.rs`), Parquet schema
    (`src/writer/parquet_reads.rs`), both collection sites in `src/bam/mod.rs`
    (SNV and indel paths), and Explorer (`app/geac_explorer.py`) filter +
    visualization.
  - Required re-running `geac collect --reads-output` — shipped in v0.3.9.

---

## Intra-sample comparison (read-type) — shipped v0.4.0

Read-type comparison view in Explorer — locus concordance, VAF density overlay,
VAF correlation scatter, strand balance density, SBS96 side-by-side, unique-loci
table. DuckDB only; gated on ≥2 distinct `read_type` values.

---

## Per-read filter validation

- **Terra cohort test with reads output** — validated on Terra using `geac:0.3.2`
  with `reads_output = true`; both `.locus.parquet` and `.reads.parquet` outputs
  confirmed per sample; cohort DuckDB contains `alt_bases` + `alt_reads` tables;
  per-read filters work in Explorer.

---

## Reads tab (Explorer) — completed plots

Requires `alt_reads` table (data collected with `--reads-output`). Tab is hidden
when the table is absent. All plots are gated on the current locus-level filters
so they reflect only the records visible in the main table.

- **Family size histogram** — distribution of `family_size` across all alt-supporting
  reads for the current filtered locus set; overlay per-sample curves or show cohort
  aggregate. A true variant should have a family size distribution similar to
  background depth; artefacts are enriched in singletons (family_size = 1).
- **Read position bias histogram** — distribution of `dist_from_read_end` for
  alt-supporting reads; a spike near 0 is a red flag for alignment artefacts or
  damaged bases at read ends.
- **Base qual vs dist from read end scatter** — one point per alt read; low base
  qual + near read end = likely artefact. Color by family_size or sample.
- **Family size vs VAF scatter** — one point per locus; x = mean family_size of alt
  reads, y = VAF. True low-VAF variants should have reasonable family sizes;
  artefacts at low VAF tend to cluster at low family size.
- **Mapping quality distribution** — histogram of `map_qual` for alt reads, split
  by repetitiveness (homopolymer ≥ 5 or STR length ≥ 6 vs non-repetitive);
  characterises which artefact classes are driven by mapping uncertainty.
- **Cohort artefact family size comparison** — for loci seen in many samples vs
  few samples, compare mean family size of alt reads via boxplot; cohort artefacts
  (sequencing noise) should show lower mean family size than recurrent true
  variants.
- **N-rich supporting loci** — added "N-rich supporting loci" section to the Reads
  tab (row 3c, after N-asymmetry). Enabled by checkbox; computes per-locus metrics
  via GROUP BY on `alt_reads`: `frac_reads_with_any_n`, `frac_reads_with_trailing_n`,
  `mean_trailing_n_run_len`, `mean_total_n_frac`. Threshold slider filters to loci
  above a minimum N-read fraction. Table sorted descending by
  `frac_reads_with_any_n`; IGV buttons for the full set. Cached on filter +
  threshold strings.

---

## Per-read filter fixes (from audit) — completed

Audit document: [`docs/per-read-filter-audit.md`](per-read-filter-audit.md).

**Bugs:**
- Re-aggregation COALESCE gives original count instead of 0 (bug #1) — in
  re-aggregation mode (`recompute_vaf=True`), a SNV locus where ALL reads fail the
  filter has no row in the `ar_agg` subquery; `COALESCE(NULL, ab.alt_count)` falls
  back to the original count instead of 0. Fix: use
  `COUNT(*) FILTER (WHERE ...)` and `COUNT(*) AS has_reads` in a single pass so the
  code can distinguish "no reads in alt_reads" (indels → preserve original count)
  from "reads exist but none pass" (SNVs → show 0).
- Warning banner text is wrong in locus-inclusion mode (bug #2) — the warning
  always said "alt_count and VAF are re-aggregated from reads passing the filter"
  regardless of whether `recompute_vaf` is True or False. In the default
  locus-inclusion mode they are not re-aggregated.
- Insert size filter missing from warning banner (bug #3) — `_active_parts` in
  the warning construction omitted insert size; activating only the insert size
  filter produced "Per-read filters active ()".
- Family-size stratified spectrum bypasses per-read filters (bug #4) — the
  `locus_fs` CTE in the family-size stratified SBS96 spectrum queried `alt_reads`
  without `_reads_where`, so the singleton/multi classification ignored the active
  per-read filters.

**Polish / labelling:**
- "Cycle number" label mismatch (semantic #6) — bundled into the `alt_reads`
  schema v2 item above (`is_read1` + `cycle` rename). Renaming
  `dist_from_read_start` → `cycle` and dropping `dist_from_read_end` resolved this
  at the schema level.
- Insert size filter: add exclude mode and document NULL behaviour (pitfall #11) —
  `insert_size BETWEEN x AND y` silently dropped all unpaired reads
  (`insert_size IS NULL`). Added an exclude-mode toggle (consistent with family
  size / MAPQ) and a sidebar caption noting that activating the filter excludes
  unpaired reads.

**Efficiency:**
- Cache slider bound MAX queries (efficiency #8) — `_reads_maxes` was computed on
  every Streamlit rerun from a full `alt_reads` scan. Gated behind a session_state
  check so it only runs once per session (the database is read-only).

---

## Explorer (Streamlit) — completed features

- IGV session download — manifest-driven BAM tracks + BED positions zip, capped at
  5 samples.
- IGV sample picker — verified: BED file correctly contains only positions from
  selected samples. The SQL query always includes a `sample_id IN (...)` clause
  before fetching `igv_df`; `make_bed` and `make_igv_session` both operate on that
  already-filtered dataframe.
- BED file deletion-locus fix — deletion loci now span the full deleted region
  `[pos, pos+del_len)`.
- Position-level drill-down — click a locus and see all samples/alleles at that
  position.
- Export filtered data to CSV — handled by Streamlit's built-in dataframe toolbar
  download button.
- On-target filter — sidebar selectbox "Target bases": All / On target / Off target.
- Gene filter — sidebar text input with partial match (ILIKE); depends on `gene`
  column being populated.
- Repeat filter — sidebar range sliders for `homopolymer_len` and `str_len`.
- Strand bias plot — dashed y=x diagonal + 95% binomial CI band; gene name in
  hover tooltip.
- Strand bias click drill-down — click/shift-click to select points; shows table
  of selected loci and IGV session with correct BAMs and BED.
- Strand bias selection: confirmed `toggle="event.shiftKey"` works in the current
  Altair/Vega-Lite version; the only blocker was a `pos_display` KeyError in the
  drill-down query (fixed).
- Strand bias selection: cast `pos` to int before building SQL WHERE clauses
  (Altair returned float, causing silent failures).
- SNV trinucleotide spectrum (SBS96) — 3×2 grid of per-mutation-type panels with
  shared y-axis; click drill-down.
- Cohort comparison view (5 steps complete):
  - Per-sample summary table — one row per sample_id with n_snv, n_insertion,
    n_deletion, mean_depth, mean_vaf, strand_balance, overlap_concordance; clicking
    a row filters all other tabs to that sample.
  - VAF distribution overlay — all samples on one plot as density curves, colored
    by sample; highlights shifted VAF distributions.
  - Strand balance scatter — one dot per sample (x = mean strand balance, y = mean
    VAF); outliers immediately visible.
  - SNV count bar chart — n_snv per sample, stacked/colored by SBS6 substitution
    type breakdown.
  - SBS96 heatmap — samples as rows, 96 trinucleotide contexts as columns,
    color = normalized count; reveals samples with unusual mutational profiles.
- NMF decomposition — fit the per-sample SBS96 spectrum against COSMIC reference
  signatures using NNLS; show the largest contributing signatures and their
  weights.
- Save/load filter state — JSON export/import of the full sidebar filter state,
  covering all locus-level and per-read filters.
- Pipeline comparison tab (DuckDB only) (5 steps complete) — side-by-side analysis
  of the same sample processed through two different pipelines (e.g. fgbio vs
  dragen). Workflow: run `geac collect` twice with different `--pipeline` and the
  same `--sample-id`, then `geac merge` both Parquets into one DuckDB. Steps:
  - Locus concordance summary (counts unique-A / unique-B / shared, by variant type).
  - VAF correlation scatter with Pearson r.
  - Unique-to-pipeline loci table.
  - SBS96 spectrum side-by-side.
  - Depth comparison scatter (`total_depth` per locus, A vs B).
- Bait-bias analysis — replaced second-BAM-pass `--emit-ref-sites` / `ref_bases`
  approach with gnomAD ~50% AF het-site depth-retention analysis using
  `sample_metrics.median_target_depth_all` as per-sample baseline; IQR boxplot
  comparing SNVs, insertions, and deletions; no second pass required.

---

## Coverage Analysis — design and shipped work (v0.4.0)

**Motivation:** `geac collect` only records positions where an alt base was observed.
`geac coverage` fills the denominator — depth at every covered position — enabling
true per-base error rates and identification of systematically undercovered sites
across a cohort.

Coverage has three confounders that must be measured, not just controlled for:

1. **Mappability** — a region may appear undercovered because reads cannot be
   placed uniquely. Without a mappability signal, low coverage and multi-mapping
   are indistinguishable from genuine dropout (GC bias, FFPE degradation, probe
   failure).
2. **Duplicates** — PCR duplicates inflate raw read counts but represent the same
   original molecule. Per-region duplication rates reveal library complexity
   problems and can differ substantially between GC-rich and GC-poor targets.
3. **Fragment overlap** — when paired reads are longer than the insert, both reads
   cover the same bases but provide only one independent observation. High overlap
   inflates apparent depth while providing no additional evidence. The fraction of
   overlapping fragments is itself a useful QC signal (short inserts relative to
   read length).

### CLI

```
geac coverage \
  --input          sample.bam \
  --reference      ref.fa \
  --output         sample.coverage.parquet \
  [--targets       targets.bed]                  # BED or Picard interval list
  [--region        chr1:1-50000]                 # alternative to --targets for a single region
  [--track         NAME:file.bedgraph]           # pre-computed annotation track (repeatable)
  [--gene-annotations genes.gtf]                 # GTF or GFF3; annotates gene, feature_type, exon_number
  [--sample-id     override]
  [--read-type     raw|simplex|duplex]
  [--pipeline      fgbio|dragen|raw]
  [--min-map-qual  20]
  [--min-base-qual 20]                           # threshold for frac_low_bq; default 20
  [--gc-window     100]                          # bp window for GC content (centred on position)
  [--min-depth     0]                            # suppress positions with total_depth below this
  [--bin-size      1]                            # aggregate N bp into one row (1 = per-position)
  [--summarize-intervals]                        # emit one row per target interval instead of per position
  [--threads       1]
```

`--targets` is strongly recommended — it bounds output size and ensures zero-depth
positions are still recorded (complete dropout is important to capture). Without
`--targets` or `--region`, the whole BAM is scanned; fine for targeted panels,
impractical for WGS without `--bin-size`.

`--track` can be repeated for multiple annotation tracks (e.g. mappability at two
k-mer lengths, a CpG density track, a GC content track). Each `NAME` becomes a
column in the output Parquet. Common sources: ENCODE GEM tracks (150-mer), genmap,
Umap, custom BEDGraph from any tool.

### Output schema (`CoverageRecord`)

One row per position (or per bin). Positions are 0-based.

```
sample_id:        String
chrom:            String
pos:              i64          # 0-based start
end:              i64          # pos+1 normally; pos+bin_size when --bin-size > 1

# ── Fragment depth ────────────────────────────────────────────────────────────
# "Fragment depth" counts unique fragments, not raw reads:
#   - Duplicate reads (BAM flag 0x400) are excluded
#   - Overlapping read pairs (same qname, both covering this position) count as 1

total_depth:      i32          # unique fragments passing --min-map-qual
fwd_depth:        i32          # forward-strand fragments
rev_depth:        i32          # reverse-strand fragments

# ── Duplicate metrics ─────────────────────────────────────────────────────────
# Computed over all reads at this position before any quality filter.
# High frac_dup indicates PCR over-amplification or poor library complexity at this locus.

raw_read_depth:   i32          # all reads including duplicates and low-MAPQ
frac_dup:         f32          # fraction of raw reads marked BAM_FDUP (0x400)

# ── Overlap metrics ───────────────────────────────────────────────────────────
# Computed over non-duplicate reads passing --min-map-qual.
# High frac_overlap means inserts are shorter than 2× read length; depth is inflated.
# overlap_depth counts fragment pairs (not reads), so it is always <= total_depth / 2.

overlap_depth:    i32          # number of fragment pairs where both reads cover this position
frac_overlap:     f32          # overlap_depth / fragment_count at this position

# ── BAM-derived mappability signals ──────────────────────────────────────────
# Computed over all non-duplicate reads before the --min-map-qual filter.
# Costs nothing since we are already iterating reads for depth counting.

mean_mapq:        f32          # mean MAPQ of all (non-dup) reads at this position
frac_mapq0:       f32          # fraction with MAPQ = 0 (definitive multi-mappers)
frac_low_mapq:    f32          # fraction with MAPQ < --min-map-qual

# ── Base quality signals ──────────────────────────────────────────────────────
# Computed over bases at this position that pass the --min-map-qual filter.
# Systematically low base quality at a site reduces effective depth just as low
# read depth does — a site with total_depth=50 but frac_low_bq=0.8 has only ~10
# usable bases. min/max capture the spread; a wide range is a different problem
# from uniformly low quality.
# --min-base-qual defaults to 20 (same default as geac collect).

mean_base_qual:   f32          # mean base quality across all bases at this position
min_base_qual:    u8           # lowest base quality observed (Phred 0–93)
max_base_qual:    u8           # highest base quality observed
frac_low_bq:      f32          # fraction of bases below --min-base-qual (default 20)

# ── Soft-clipping signal ──────────────────────────────────────────────────────
# Computed over non-duplicate reads passing --min-map-qual.
# Heavy soft-clipping at a position indicates reads that partially align — a sign
# of structural variation, probe edge effects, or adapter contamination.
# frac_soft_clipped is the fraction of reads where the query position falls within
# a soft-clipped region of the CIGAR (i.e. the base is present in the read but
# not contributing to the alignment at this reference position).

frac_soft_clipped: f32        # fraction of reads soft-clipped at this position

# ── Insert size distribution ──────────────────────────────────────────────────
# Computed from properly paired, non-duplicate reads passing --min-map-qual.
# Insert size = TLEN (template length) from the BAM record; only meaningful for
# paired-end reads where both mates are mapped (FLAG: properly paired, 0x2).
# Short inserts relative to read length explain high frac_overlap and reduced
# effective depth. High variance indicates a heterogeneous library.
# Unpaired or single-end reads contribute 0 usable insert size observations;
# n_insert_size_obs records how many paired reads contributed to these stats.

mean_insert_size:     f32     # mean insert size across paired reads at this position
median_insert_size:   f32     # median insert size (requires buffering; see note below)
min_insert_size:      i32     # smallest insert size observed
max_insert_size:      i32     # largest insert size observed
n_insert_size_obs:    i32     # number of properly paired reads contributing

# Note on median_insert_size: computing a true median requires storing all insert
# sizes seen at each position, which is memory-intensive for deep coverage. Use
# reservoir sampling (e.g. keep up to 1000 values) or an approximate algorithm
# (e.g. t-digest) to keep memory bounded. The median is more robust to outliers
# (e.g. chimeric read pairs) than the mean.

# ── GC content ────────────────────────────────────────────────────────────────
# Computed directly from the reference FASTA (already required as --reference).
# No external track needed. Window size is configurable via --gc-window (default:
# 100 bp centred on the position). GC content is the primary explainer of
# amplification/capture dropout that is independent of mappability.

gc_content:       f32         # fraction of G+C bases in --gc-window around this position

# ── Pre-computed annotation tracks ───────────────────────────────────────────
# One column per --track NAME:file entry. Column name = NAME, type = Float32, nullable.
# Example: --track gem150:gem_150mer.bedgraph  →  column "gem150"
#          --track umap50:umap_k50.bedgraph    →  column "umap50"
# Requires dynamic Arrow schema construction at runtime (not a fixed Rust struct).

<track_name>:     Option<f32>  # 0.0–1.0 score from the named BEDGraph track

# ── Gene / feature annotation ─────────────────────────────────────────────────
# Populated when --gene-annotations is provided.
# feature_type and exon_number extend the existing GeneAnnotations infrastructure.

gene:             Option<String>   # gene name
feature_type:     Option<String>   # "exon", "intron", "5UTR", "3UTR", "CDS"
exon_number:      Option<i32>      # exon number within the transcript (from GTF/GFF3 attribute)

# ── Optional target annotation ────────────────────────────────────────────────
on_target:        Option<bool>     # populated when --targets is given

# ── Provenance ────────────────────────────────────────────────────────────────
read_type:        ReadType
pipeline:         Pipeline
```

### Read counting semantics

The three-layer decomposition at each position:

```
All reads
  └─ subtract BAM_FDUP reads       → raw_read_depth, frac_dup
       └─ subtract low-MAPQ reads  → mean_mapq, frac_mapq0, frac_low_mapq
            └─ collapse same-qname pairs as 1 fragment  → total_depth, overlap_depth, frac_overlap
```

This means `total_depth` is directly comparable to the `total_depth` in `alt_bases`
from `geac collect`, which uses the same duplicate-exclusion and overlap-collapsing
logic.

### Mappability diagnostic table

| `frac_mapq0` | track score | `total_depth` | Interpretation |
|---|---|---|---|
| High | Low | Low | Classic multi-mapping — expected, filter confidently |
| High | High | Low | Unexpected low MAPQ — SV, misassembly, or aligner artifact |
| Low | Low | Low | Genuine dropout (GC bias, FFPE, probe failure) |
| Low | Low | Normal | Mappability track k-mer length may not match read length |

### Per-interval summary mode (`--summarize-intervals`)

When `--summarize-intervals` is given alongside `--targets`, `geac coverage` emits
one row per target interval instead of one row per position. This is the natural
output format for the customer-facing Coverage Explorer — customers want to know
"exon 3 of BRCA1: mean depth 45x, 94% at ≥30x", not a table of 200 individual
positions.

The per-interval schema adds aggregated columns and drops the position-level ones:

```
sample_id:           String
chrom:               String
start:               i64          # 0-based interval start (from targets file)
end:                 i64          # 0-based interval end
interval_name:       Option<String>  # name field from BED col 4 / Picard interval name
gene:                Option<String>
feature_type:        Option<String>
exon_number:         Option<i32>

# Depth summary across all positions in the interval
mean_depth:          f32
median_depth:        f32
min_depth:           i32
max_depth:           i32
frac_at_1x:          f32          # fraction of bases with total_depth >= 1
frac_at_10x:         f32
frac_at_20x:         f32
frac_at_30x:         f32
frac_at_50x:         f32
frac_at_100x:        f32
n_bases:             i32          # total number of positions in the interval

# Aggregated signals (means across positions in the interval)
mean_gc_content:     f32
mean_mapq:           f32
mean_frac_mapq0:     f32
mean_frac_dup:       f32
mean_frac_overlap:   f32
mean_frac_soft_clipped: f32
mean_base_qual:      f32
mean_insert_size:    f32

read_type:           ReadType
pipeline:            Pipeline
```

The depth threshold columns (`frac_at_Nx`) use fixed thresholds rather than a
configurable value so that interval summaries from different runs are directly
comparable. The customer explorer can then filter by whichever threshold is
meaningful for that panel.

Per-interval and per-position outputs are written to separate Parquet files:
`{sample}.coverage.parquet` (per-position) and `{sample}.coverage.intervals.parquet`
(per-interval). `geac merge` inserts both into the DuckDB as `coverage` and
`coverage_intervals` tables respectively.

### Pre-computed annotation tracks (`--track`)

**Format**: BEDGraph (chrom, start, end, score). Covers ENCODE GEM, genmap, Umap,
and any custom track. bigWig support can be added later if needed.

**Multiple tracks** are useful in practice: e.g. mappability at the experiment's
read length alongside a GC-content track lets you disentangle GC bias from
repeat-element dropout.

**Implementation**: each track is loaded into a sorted `Vec<(i64, i64, f32)>` per
chromosome. Lookup at each pileup position uses binary search (O(log n)). For
targeted panels this fits comfortably in memory. For WGS without `--targets`, a
streaming approach (advance through the sorted track in lock-step with the sorted
pileup) avoids loading the full ~2 GB track.

Because the number of tracks is not known until runtime, `CoverageRecord` cannot
be a plain Rust struct with fixed fields. Instead, the Arrow `Schema` and
`RecordBatch` are constructed dynamically in `src/writer/parquet_coverage.rs`
based on the track names provided.

### Gene annotation extensions

`--gene-annotations` reuses `src/gene_annotations.rs` but extends it to also store
`feature_type` (exon / intron / UTR / CDS) and `exon_number` from the GTF/GFF3
attribute field. This enables Explorer queries like:

```sql
-- Coverage of BRCA1 exon 1 across all samples
SELECT sample_id, pos, total_depth, frac_mapq0, gem150
FROM coverage
WHERE gene = 'BRCA1' AND feature_type = 'exon' AND exon_number = 1
ORDER BY pos;
```

### DuckDB integration

`geac coverage` outputs `{sample}.coverage.parquet`. `geac merge` detects these
by schema (presence of `frac_dup`; absence of `alt_allele`) and inserts into a
`coverage` table in the cohort DuckDB alongside `alt_bases`.

```sql
-- Systematically undercovered positions across the cohort
SELECT chrom, pos, gene, exon_number,
       COUNT(DISTINCT sample_id)  AS n_samples,
       AVG(total_depth)           AS mean_depth,
       AVG(frac_dup)              AS mean_frac_dup,
       AVG(frac_mapq0)            AS mean_frac_mapq0,
       AVG(gem150)                AS mean_mappability   -- if track was provided
FROM coverage
GROUP BY chrom, pos, gene, exon_number
HAVING AVG(total_depth) < 20
ORDER BY mean_depth;

-- Alt bases in low-mappability, low-coverage context
SELECT a.*, c.mean_depth, c.frac_mapq0, c.gem150
FROM alt_bases a
JOIN (
    SELECT chrom, pos, AVG(total_depth) AS mean_depth,
           AVG(frac_mapq0) AS frac_mapq0, AVG(gem150) AS gem150
    FROM coverage GROUP BY chrom, pos
) c ON a.chrom = c.chrom AND a.pos = c.pos
WHERE c.frac_mapq0 > 0.3;
```

### Implementation steps (Steps 1–10 completed; Step 11 in TODO.md)

- Step 1: Extend `src/gene_annotations.rs` — add `feature_type` and `exon_number`
  to the annotation lookup result; update GTF/GFF3 parser to extract the
  `exon_number` attribute.
- Step 2: Add `src/track.rs` — `AnnotationTrack` struct; BEDGraph loader;
  binary-search lookup; `TrackSet` holding multiple named tracks; chr-prefix
  bridging; `--track NAME:FILE` repeatable flag in `CoverageArgs`; dynamic Arrow
  schema columns in `CoverageWriter`.
- Step 3: Add `CoverageArgs` to `src/cli.rs`; add `Command::Coverage` variant.
- Step 4: Add `src/coverage/mod.rs` — pileup loop with three-layer read counting
  (raw → de-dup → de-overlap → total_depth); all BAM-derived signals (mapq, base
  qual, insert size, GC content, overlap, dup fraction); zero-depth fill-in for
  target positions; `compute_gc_content` from reference cache.
- Step 5: Add per-interval aggregation pass in `src/coverage/mod.rs` — after the
  per-position pass, group positions by target interval and compute the interval
  summary schema; emit as a separate `Vec<IntervalRecord>`.
- Step 6: Add `src/writer/parquet_coverage.rs` — fixed Arrow schema matching
  `CoverageRecord`; Float32 columns for fractional signals.
- Step 7: Update `src/main.rs` to handle `Command::Coverage`.
- Step 8: Update `src/merge.rs` — detect `.coverage.parquet` by suffix; insert
  into `coverage` DuckDB table; index on `(sample_id, chrom, pos)`.
- Step 9: Add `wdl/geac_coverage.wdl` — scatter `geac coverage` over a cohort,
  merge all `.coverage.parquet` files into a `coverage` table in the cohort
  DuckDB.
- Step 10: Integration tests — 9 new tests covering all core coverage signals:
  `coverage_basic_depth`, `coverage_frac_dup_excludes_duplicates`,
  `coverage_mapq0_tracked_and_excluded`, `coverage_gc_content_computed_from_reference`,
  `coverage_gc_content_zero_for_all_a_reference`,
  `coverage_targets_emits_zero_depth_positions`,
  `coverage_no_targets_skips_zero_depth`,
  `coverage_insert_size_from_paired_reads`,
  `merge_routes_coverage_parquet_to_coverage_table` (all passing).

---

## Customer-facing Coverage Explorer (`app/geac_coverage_explorer.py`) — done

- App scaffold: file-path text input, `DataSource.open_coverage()`, version-mismatch warning.
- Sidebar filters: samples, chromosome, gene partial-match, on-target.
- Sidebar IGV integration: manifest path input, `load_manifest`, GCS OAuth token.
- Sidebar Advanced expander: `geac_metadata` and `geac_inputs` display.
- Tab 1 — Summary: per-sample depth table, mean-depth bar chart, QC-fractions grouped bar.
- Tab 2 — Depth distribution: per-sample depth histogram, fraction-at-threshold table.
- Tab 3 — GC bias: mean depth vs GC bin (5% bins) per sample; frac_mapq0 overlay expander.
- Tab 4 — Low coverage: undercovered positions table (depth + fraction-of-samples sliders),
  gene bar chart with click-to-drill-down, IGV locus link on row select.
- Tab 5 — IGV: embedded IGV.js viewer; BAM/CRAM tracks from manifest; per-track height slider; GCS OAuth.
- Tab 6 — Intervals (DuckDB + coverage_intervals only):
  - Undercovered intervals table (depth threshold + fraction-of-samples sliders).
  - GC bias scatter: mean_gc_content vs mean_depth per interval, colored by mean_frac_mapq0.
  - Per-exon heatmap: genes × intervals, color = frac_at_30x; sort and top-N controls.
- `geac_config.py`: TOML config loader for pre-populating paths (data, manifest, genome, etc.).
- `geac_metadata` and `geac_inputs` DuckDB tables (provenance) — written by `geac merge`.

---

## Fragmentomics — context

**Motivation:** Cell-free DNA (cfDNA) derived from tumour cells has distinct
fragmentation patterns compared to normal cfDNA. Capturing fragment-level features
enables nucleosome positioning analysis and cancer detection via fragmentation
signatures — without requiring additional sequencing.

All required infrastructure is already in place: `RefCache` for reference access,
`insert_size` in the reads schema, and read start/end positions. `frag_gc` is
already implemented in `alt_reads`.

### ~~Medium-term: `frag_midpoint` on `alt_reads`~~ — won't do

`dist_from_read_end` (cycle) already serves artefact detection on `alt_reads`.
Nucleosome positioning — the other use case for midpoint — requires the full
fragment population, not just alt-supporting reads. The `fragments` table already
has `midpoint` for every fragment, making this the correct foundation for
nucleosome positioning work. Adding `frag_midpoint` to `alt_reads` would be a
biased, low-value duplicate.

### Long-term: all-reads end-motif table (forward work, see TODO.md)

End motifs (the k-mer at the fragment cut site) are only meaningful as a
comparison between alt-supporting and reference-supporting reads. Since
`alt_reads` only covers alt-supporting reads, adding `frag_end_motif` there has
no reference baseline to compare against.

The correct approach is a new table (separate from `alt_reads`) that captures,
for every read at a pileup position: the 4-mer cut-site motif and whether the
read supports alt or ref. The motif should be **reference-based**
(`seq[frag_start-2..frag_start+2]`, centered on the cleavage point — 2 bases
outside the fragment, 2 inside) to avoid hard-clip ambiguity. The Explorer would
then compare motif distributions between alt-supporting and reference-supporting
reads as a bait-bias signal.

Forward work items live in `TODO.md` under **Fragmentomics**.

---

## CI / Release — completed

- GitHub Actions release workflow — on push of a `v*.*.*` tag, builds native
  amd64 (ubuntu-latest) and arm64 (ubuntu-22.04-arm) images, pushes by digest,
  then merges into a multi-platform manifest at `ghcr.io/fleharty/geac:<version>`
  and `:latest`. Uses `GITHUB_TOKEN` — no external credentials needed.

---

## WDL / Terra — completed

- WDL task wrapping `geac collect` — single-sample workflow in
  `wdl/geac_collect.wdl`.
- Terra-compatible Docker image — multi-stage `docker/Dockerfile` with htslib +
  geac binary.
- WDL workflow — scatter `geac collect` across a sample list, then gather with
  `geac merge` (`wdl/geac_cohort.wdl`).
- WDL task wrapping `geac merge` — standalone workflow in `wdl/geac_merge.wdl`.
- Test on Terra with a small cohort — v0.3.0 successfully run on a cohort of
  samples.

---

## MNV detection (Steps 1–4) — 2026-05-04

**Goal:** Identify multi-nucleotide variants (MNVs) — pairs of nearby substitutions
that co-occur on the same DNA fragment, distinguishing true MNVs from coincidental
independent SNVs at neighbouring positions.

**Approach:** Hash each read's `qname` (FNV-1a 64-bit, inline — no new dependency)
to produce a stable `fragment_id`. Both reads of a pair share the same qname, so
the hash is a stable per-fragment identifier. Store it in `alt_reads` as `Int64`.
The MNV candidate query is a self-join on `fragment_id` within `alt_reads`, restricted
to pairs within a configurable distance.

**Why FNV-1a inline rather than a crate:** The hash is 3 lines and has no dependencies.
The collision probability over typical cohort sizes is negligible. Adding a hash crate
for this would be premature.

**Key files changed:**
- `src/bam/pileup.rs` — `fnv1a_64`, `fragment_id` on `LocusRead` and `ReadDetail`
- `src/bam/indel.rs` — same for the indel `ReadDetail` construction path
- `src/record.rs` — `fragment_id: i64` on `AltRead`
- `src/bam/builders.rs` — propagated into `build_alt_read`
- `src/writer/parquet_reads.rs` — `fragment_id` column in Parquet schema and batch
- `schema/geac_schema.json` — added to `alt_reads.required_columns`
- `app/explorer/tabs/mnv_candidates.py` — new tab with distance/min-co-count sliders,
  candidates table, IGV drill-down for selected pairs, dinucleotide context bar chart

**MNV candidate query design:**
```sql
WITH filtered_reads AS (
    SELECT fragment_id, sample_id, chrom, pos, alt_allele
    FROM alt_reads [joined to loci passing sidebar filters]
),
co_occurring AS (
    SELECT a.sample_id, a.chrom, a.pos AS pos1, a.alt_allele AS alt_allele1,
           b.pos AS pos2, b.alt_allele AS alt_allele2, COUNT(*) AS co_count
    FROM filtered_reads a
    JOIN filtered_reads b
        ON a.fragment_id = b.fragment_id AND a.sample_id = b.sample_id
        AND a.chrom = b.chrom AND a.pos < b.pos AND (b.pos - a.pos) <= :max_dist
    GROUP BY ... HAVING COUNT(*) >= :min_co
)
-- join back to locus table for ref_allele, alt_count, frac_cooccurring
```

**Step 5 — Integration test (2026-05-04):** `reads_fragment_id_enables_mnv_detection`
in `tests/integration.rs`. Uses `write_mnv_bam` (added to `tests/common/mod.rs`): 3
fragments carry both T@pos50 and G@pos51, 2 carry only T@pos50, 2 carry only G@pos51,
5 carry ref. Asserts `co_count=3` and `frac_cooccurring=0.6` via a direct DuckDB
self-join on `fragment_id`; also verifies that single-substitution fragments do not
appear in the co-occurring set.

---

## `--max-pileup-depth` and Explorer indel-length filter (v0.4.27, 2026-05-11)

**Bug discovered:** A 42 bp deletion at chr1:1204427 in an amplicon BAM (14 126 raw
reads at the anchor; FLAG=163, MAPQ=60, primary, not duplicate, supported by exactly
one read) did **not** appear in the per-sample `alt_bases` Parquet. A nearby 4 bp
deletion at the same locus *did* appear, ruling out a general indel-tally bug. The
parquet's `total_depth` at neighbouring SNV positions clipped near ~6 984.

**Root cause:** `rust_htslib::bam::IndexedReader::pileup()` defaults `bam_plp_set_maxcnt`
to **8000**. When raw depth exceeds the cap, htslib silently downsamples — the lone
deletion-bearing read is discarded in the lottery. The default is a WGS-era heuristic
that's wrong for amplicon panels and any low-VAF / rare-event work.

**Fix:** Added `--max-pileup-depth` to `CollectArgs`, `CoverageArgs`, and
`AnnotateNormalArgs`. Sentinel `0` (the documented default) means *unlimited*; it
maps internally to `i32::MAX as u32` because htslib has no true "unlimited" and
literal `0` passed to `set_max_depth` would cap at zero. Helper
`resolve_max_pileup_depth` in `src/bam/ref_utils.rs`. All three `bam.pileup()` call
sites (`src/bam/mod.rs`, `src/coverage/mod.rs`, `src/normal.rs`) call
`set_max_depth` before iterating. Plumbed through `wdl/geac_collect.wdl`,
`wdl/geac_coverage.wdl`, `wdl/geac_annotate_normal.wdl`, `wdl/geac_cohort.wdl`. See
CHALLENGES.md entry for the full diagnostic timeline.

**Why default unlimited rather than a bounded large number:** picked 0=unlimited to
match the bug's framing (every dropped read at a panel locus is a possible missed
variant). On WGS this is a memory risk; if it ever bites, the default can move to
something like 1 000 000 without changing the CLI surface.

**Explorer indel-length filter (same release):** added a sidebar slider for
`indel_len_range` (default `(1, 500)`) that filters insertions/deletions by length
derived on-the-fly from `LENGTH(alt_allele) - 1`. SNVs are explicitly preserved
(`variant_type = 'SNV' OR ...`) so the length filter doesn't double-up with the
Variant-type multiselect. Considered adding a stored `indel_len` column to the
schema; rejected because `LENGTH(varchar)` is O(1) in DuckDB and a schema bump
would force every existing Parquet to be re-collected for no measurable speedup.
A new `mnv_min_dist` slider was added in the MNV candidates tab alongside the
existing max-distance slider (changed the join condition to `BETWEEN min AND max`).

## Fusion index — k-mer copy-number quantification & two-tier uniqueness (experimental)

The fusion k-mer index originally treated uniqueness as binary. Two refinements
were added to `geac experimental build-fusion-index` and `... fusions`.

**Two tiers already existed, only one was reported.** Step 1's cross-gene dedup
already yields *panel-wide* uniqueness (k-mers unique among the indexed gene
bodies); `--check-genome-uniqueness` is the strict genome-wide add-on. The work
was to quantify and expose both, not to invent panel-wide.

- **Copy counting.** The genome pass already counted each candidate k-mer's
  genome-wide occurrences (`genome_counts: HashMap<u64,u8>`, saturating at 255)
  but discarded everything except `== 1`. We now keep that map and (a) write a
  global `--copy-histogram-output` (`copies`,`n_kmers`), (b) add per-gene copy
  buckets + strict/relaxed base-coverage columns to `--gene-stats-output`, and
  (c) store per-k-mer `genome_copies` in the DuckDB `kmers` table (nullable;
  NULL when the genome pass didn't run — backward compatible with the other
  index readers, which all select explicit columns).
- **Relaxed build tier.** `--max-genome-copies N` (default 1) retains k-mers
  occurring `1..=N` times genome-wide. The genome-copy retain was moved to run
  *after* the stats/histogram writers so they see the full panel-unique survivor
  set together with every count.
- **Query-time filter.** `geac experimental fusions --max-kmer-copies N` ignores
  k-mers exceeding N copies (or with unknown/NULL copies). Errors clearly if the
  index lacks `genome_copies`. Default path is unchanged (no extra lookups).
- **Per-copy-tier BEDs.** `--bed-output-by-copies <prefix>` writes one merged BED
  per genome-wide copy number (`<prefix>.copies1.bed`, `.copies2.bed`, … up to
  `--max-genome-copies`), each merged only among k-mers of that exact copy count,
  so strictly-unique vs near-unique regions load as separate IGV tracks / target
  BEDs. The existing `--bed-output` still emits the combined set. Merge logic is
  shared via `write_merged_bed_intervals`. Adjacent tiers may overlap by up to
  k−1 bases at boundaries (independent k-mer sets with overlapping ±k windows).
  For tiers ≥2 a 4th BED column names every **full-GTF** gene each interval's
  k-mers occur in (comma-separated, sorted), or `intergenic` — so a near-unique
  region is labelled with the paralog(s) driving its non-uniqueness even when
  those paralogs are outside the `--genes` panel (e.g. `GENEB,GENEBP1`,
  `GENEC,intergenic`). Implemented by a `GeneIntervals` index built from the
  unfiltered GTF and a per-k-mer `GeneAnnot` (gene-index set + intergenic flag)
  accumulated during the genome scan; the 1× tier stays plain 3-column. The
  annotation map is only allocated when `--bed-output-by-copies` is set; for a
  whole-genome (non-panel) build it adds memory proportional to the candidate set.

Demonstrated on a synthetic gene with a duplicated 100 bp block: strict
genome-unique coverage = 81%, relaxed (≤2) = 100% — the same pattern as a real
tandem-repeat/paralog gene like DUX4 (~16% genome-wide).

**Known limitation (not fixed):** `--max-genome-copies` only relaxes *genome-wide*
paralog copies; it does not rescue *intra*-gene tandem repeats (those survive
step 1 already, assigned to their single gene). Also, the base-coverage metric
undercounts repeat genes because `kmer_to_pos` stores only each k-mer's first
position; switching to all-positions would fix it but was left out of scope.

**Post-release audit fixes (v0.4.25–30 review).** A code audit of the fusion work
surfaced and fixed: (1) cross-file label inconsistency — the Parquet `gene_a`/`gene_b`
ordered by gene index while the `FX` tag / kmer-hits TSV ordered alphabetically; both
now use one canonical gene-index order (see `CHALLENGES.md`). (2) The top-2 gene
selection was duplicated and depended on `HashMap` iteration order; it is now a single
`fragment_top_pair()` with a deterministic tie-break, and `assign_gene`'s per-read
winner breaks ties on lowest gene index. (3) Added the first integration test for the
experimental pipeline (`fusions_label_ordering_consistent_across_outputs`) — the
fusion code previously had no automated coverage beyond `kmer.rs` unit tests. Remaining
known gaps recorded in the audit: directionality (5'/3' partner) is intentionally
dropped by pair normalization; `chrom_a`/`chrom_b` are annotated gene loci, not
observed breakpoints; pre-0.4.30 indexes (no `genome_copies` column) are believed
compatible on the default path but unverified against a real old index.

## Fusion breakpoint localization & k-mer-hits memory fix (experimental)

Two related changes to `geac experimental fusions`, plus Terra WDLs.

**`--kmer-hits-output` memory fix.** The original implementation set
`collect_detail = kmer_hits_output.is_some()` and accumulated every matching
k-mer for *every* gene-assigned read in `ReadHit.kmer_hits` during the single
BAM scan — gigabytes of peak RSS on deep WGS BAMs, most of it for reads that
never belonged to a passing fusion. Fixed by mirroring the existing
`--reads-output` pattern: the first pass keeps only the minimal `ReadHit`
(gene_idx + mapq), and k-mer detail is re-derived in a *second BAM pass*
restricted to fusion-supporting qnames. `ReadHit` lost its `chrom`/`pos`/
`is_read1` fields (now read straight from the BAM record in the second pass).

**`kmer_pos_in_read` + `--breakpoints-output`.** `kmer_iter` already yields
`(kmer, start_pos_in_read)`, so the 0-based offset was threaded through
`assign_gene` into the kmer-hits TSV. On a junction-spanning read the offset
where Gene A k-mers give way to Gene B k-mers, projected onto the alignment
start, locates the breakpoint to within a few bp. `--breakpoints-output` does
this in-binary: it shares the *same* second BAM pass as `--kmer-hits-output`
(only one extra scan when both are set), votes per-gene chromosomes from
single-gene reads, collects spanning-read estimates, and writes one row per
fusion with median coordinate, contributing-read count, and std (spread).

**Why in-binary rather than a script.** A prototype `scripts/call_breakpoints.py`
post-processed the kmer-hits TSV, but folding it into the command keeps the
workflow single-pass and makes breakpoints a first-class Terra output. The
script was removed.

**Known limitation.** Reverse-strand reads report `kmer_pos_in_read` relative to
the stored SEQ (the reverse complement of the reference strand), which widens
`bp_*_std`. `scripts/reconstruct_fusions.sh` (de novo assembly + minimap2) is the
sequence-level cross-check that handles strand correctly; the two approaches are
complementary. WDLs: `wdl/experimental/geac_build_fusion_index.wdl` and
`geac_fusions.wdl`. Cromwell rejects the `None` literal, so optional task outputs
use `if cond then [path] else []` (0/1-element arrays) rather than `File?`.

---

## Fusion — junction-coherence / co-linearity filter (2026-06-01)

**Motivation.** The most common recurrent false-positive class is paralog/homology: a
single non-chimeric read from one locus carries k-mers that the index treats as unique
to *two* genes, so it looks like a junction. These reads are distinguishable from real
fusion-spanning reads by the *spatial pattern* of k-mer hits inside the read.

On a real junction-spanning read, Gene-A k-mers occupy a contiguous left block and
Gene-B k-mers occupy a disjoint right block (or vice versa). On a homology/paralog
read, k-mers from both genes are *interleaved* — the same bases match both genes.
`kmer_pos_in_read` already captures this signal.

**Implementation.** `assign_gene` was extended to track per-gene `(count, min_pos,
max_pos)` for all genes hit by a read (previously only the winner's count was
tracked). `ReadHit` now stores the top-2 gene indices and their k-mer position ranges.

During fragment aggregation, `read_coherence(rh, ga, gb, k, min_anchor)` classifies
each read as:
- **not spanning** — the read only hits one of the two fusion genes
- **spanning but not anchored** — hits both genes but fewer than `min_anchor` k-mers
  on one side
- **spanning and coherent** — hits both genes with enough anchor on both sides, and
  `max(A_positions) + k ≤ min(B_positions)` or the mirror — disjoint blocks

`FusionRecord` now carries `n_spanning_reads` and `n_coherent_reads` in both Parquet
and TSV output so users can inspect coherence without re-running.

**New flags.** `--min-coherent-fragments N` (default 0 = disabled): require at least
N coherent spanning reads. `--min-anchor-kmers N` (default 3): minimum k-mer hits per
gene to count as anchored. Discordant pairs (R1→GeneA, R2→GeneB, no read spanning
both) contribute `n_spanning_reads += 0`, so the filter never penalizes clean
discordant-pair evidence.

**Why default-off.** The filter requires at least one spanning read with a clean
A/B partition to pass. For low-AF fusions or short insert sizes, all fragments may
be discordant pairs, so hard-requiring ≥1 coherent read would silently drop real
calls. The user must opt in with `--min-coherent-fragments 1`. The `n_coherent_reads`
column is always written, so users can post-filter or inspect without re-running.
