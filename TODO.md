# GEAC — TODO

Active backlog. Only unchecked, actionable items live here. For shipped work, design
notes, and historical context, see [docs/DEVELOPMENT_LOG.md](docs/DEVELOPMENT_LOG.md).
For high-level milestones and release themes, see [ROADMAP.md](ROADMAP.md).
For targeted review findings that have not yet been scheduled, see
[docs/CODE_AUDIT.md](docs/CODE_AUDIT.md).

Conventions:
- Items are grouped by area, not chronology.
- Each top-level bullet should be self-contained enough to start work on. If it
  needs deep context, link to the relevant section of `docs/DEVELOPMENT_LOG.md`.
- When an item is completed, move it (with any useful context) into
  `docs/DEVELOPMENT_LOG.md` rather than leaving a `[x]` here.

---

## Per-read filter validation

- [ ] **Single-site read-level filter test** — run `geac collect --reads-output
  --region chr:start-end` on one sample restricted to the region of a known het
  variant. Single locus record with a known number of alt-supporting reads. Use
  the Explorer to manipulate family size and dist-from-read-end sliders; verify
  `alt_count` changes as expected. Confirm include vs exclude toggle behaviour
  matches intuition.

---

## Sample Identity tab (Explorer)

Audit findings captured 2026-06-08. The tab shipped as a DuckDB-only cohort
analysis; these are follow-up hardening items before treating it as production
QC.

- [ ] **Add a large-cohort execution guard** — `compute_pairwise_identity()` caps
  marker count (`max_markers`) but not sample count. The core self-join is still
  effectively `sum(marker_sample_count^2)`, so a cohort with hundreds/thousands
  of samples and common markers can make the first tab load very expensive before
  the heatmap's 60-sample display guard helps. Add a hard sample/pair budget, a
  required explicit "Run" action for large cohorts, or a more scalable pairwise
  overlap implementation.
- [ ] **Fix Sample Identity cache invalidation for `si_t_low`** —
  `app/explorer/tabs/sample_identity.py` builds `_si_all_loci` from candidate
  pairs using `p.t_low`, but `_si_cache_key` omits `p.t_low`. Changing the swap
  Jaccard threshold can reuse stale all-loci scores, so newly relevant flagged
  pairs may show missing `all_loci_jaccard` until another cache-keyed control
  changes.
- [ ] **Decide whether Sample Identity honors sidebar filters** — helper SQL
  queries `alt_bases` directly instead of `ctx.table_expr` / `ctx.where`, so
  batch, pipeline, subject, sample, timepoint, and other sidebar filters do not
  affect identity analysis. Either thread the filtered table expression into the
  helpers, or label the tab clearly as a global whole-cohort analysis.
- [ ] **Detect conflicting `subject_id` values per sample** — `subject_map()` uses
  `ANY_VALUE(subject_id)` grouped by `sample_id`. If merged inputs contain the
  same sample with conflicting subjects, classification becomes arbitrary and can
  hide the labeling problem the tab is meant to surface. Report such conflicts
  explicitly and avoid flagging dependent pairs until they are resolved.
- [ ] **Reconcile the CLI vs Explorer flag taxonomies — the headless CLI can't
  report a sample swap** (found 2026-06-13 auditing how the flag is consumed). The
  Rust CLI (`src/sample_identity.rs`) and the Explorer (`classify_flags` in
  `app/explorer/sample_identity_helpers.py`) are two *independent* implementations:
  the Explorer recomputes pairwise identity + flags straight from `alt_bases` and
  never reads the CLI's `--output` TSV, so the CLI `flag` column is consumed only by
  whoever reads that TSV directly (humans / batch / WDL — no WDL parses it today).
  They have drifted in both vocabulary and capability:

  | case | CLI (`sample_identity.rs`) | Explorer (`classify_flags`) |
  |------|----------------------------|-----------------------------|
  | same subject + high similarity | `EXPECTED_MATCH` | `EXPECTED_MATCH` |
  | diff/absent subject + identical | `SAME_INDIVIDUAL` | `UNKNOWN_DUPLICATE` |
  | **same subject + low similarity (a swap)** | **`""` or dropped** | **`POSSIBLE_SWAP`** |
  | thresholds | one band (`min_jaccard`/`min_concordance`) | two bands (`t_high`/`conc_high`, `t_low`/`conc_low`) |

  The material gap: the CLI **cannot emit a swap signal at all** — `compare_pair`
  drops non-passing pairs unless `--all-pairs` (`if !args.all_pairs && !passes`), and
  even then the flag match-arm only labels when `passes`, so a same-subject /
  genetically-divergent pair gets `""`. A Terra/WDL pipeline parsing the TSV thus
  gets a strictly weaker result than the Explorer. Also note: the
  `SAME_INDIVIDUAL`-when-`subject_id`-is-`None` overloading I originally flagged is
  *moot for the Explorer* (it ignores the column) but still confuses direct TSV
  readers. Two fix levels:
    1. *Minimum (doc):* state in `docs/cli.md` that swap detection requires the
       Explorer and that the CLI flag set is `{EXPECTED_MATCH, SAME_INDIVIDUAL}` on
       passing pairs only.
    2. *Proper (align CLI to the Explorer model):* give the CLI the two-band
       taxonomy — emit `POSSIBLE_SWAP` for same-subject pairs below
       `t_low`/`conc_low` (requires *not* filtering same-subject pairs before
       flagging), rename `SAME_INDIVIDUAL`→`UNKNOWN_DUPLICATE` for parity, add the
       threshold flags + expose in the WDL. Behavior change to a released command
       (new/renamed flag values) — needs a decision before implementing.

---

## WDL / Terra pipeline

- [ ] **Verify Terra call-caching after explicit resource URI inputs** — WDL now
  accepts `bam_uris` / `bai_uris`, `gnomad_uri`, and `targets_uri` string inputs
  so embedded IGV paths do not depend on Cromwell `File→String` coercion. Submit
  the same `geac_cohort.wdl` inputs twice with those values bound to stable
  `gs://` string columns and confirm `GeacCohort.Collect` shards are cache hits.
  If not, diff the two runs'
  `callCaching.hashes` via Cromwell's `/api/workflows/v1/callcaching/diff`
  endpoint to identify the remaining differing hash key.

---

## Reads tab (Explorer)

Requires `alt_reads` table. Plots already shipped are listed in
`docs/DEVELOPMENT_LOG.md`.

- [ ] **Reads tab review** — work through all plots with real cohort data; assess
  usefulness, make changes, remove plots that don't add value.
- [ ] **Investigate trailing-N pattern after alt-supporting bases** — in some
  reads, an alt base is followed by a run of `N` bases. Determine whether this
  is primarily a property of the input BAM/consensus pipeline, an alignment/context
  effect (e.g. indels, read-end failure, masking), or something GEAC should
  explicitly annotate or visualize. Start with real examples in IGV and the Reads
  tab; collect compact read-context metrics in `alt_reads` (`n_before_alt`,
  `n_after_alt`, `n_n_before_alt`, `n_n_after_alt`, `leading_n_run_len`,
  `trailing_n_run_len`) and compare the before/after N burden distributions.
  Compare by pipeline/read_type where possible; assess whether these loci are
  enriched for low-confidence or pipeline-unique calls.
- [ ] **Family size vs VAF click-through** — add click/shift-click selection with
  drill-down table and IGV buttons, same as the strand bias plot. Currently
  blocked: `selection_point` with `on_select="rerun"` returns `{"fsvaf_select":
  {}}` regardless of what is clicked; strand bias works identically so the root
  cause is unknown. Investigate when Altair/Streamlit version context is clearer.

---

## Strand bias tab (Explorer)

Audit findings captured 2026-06-08. The tab should remain, but it should behave
more like a QC workflow first and a locus drill-down second.

- [ ] **Add QC eligibility controls/defaults** — the scatter currently includes
  all alt loci, including single-read and very low-depth observations. Add
  minimum alt-count/depth controls or a default "QC-eligible loci" mode so the
  primary plot emphasizes actionable strand-balance signals.
- [ ] **Reshape the tab around sample-level QC summaries** — start with per-sample
  metrics such as median strand balance, fraction of loci outside the expected
  band, total evaluable loci, and optional batch/read-type grouping. Keep the
  existing locus scatter as a drill-down from suspicious samples.
- [ ] **Make large-cohort sampling stratified or explain its limits** — the fixed
  5,000-locus reservoir sample can underrepresent small problematic samples or
  batches. Prefer stratified sampling by sample/variant type/batch when columns
  are available, or surface a warning that the sampled scatter is exploratory.

---

## Pipeline comparison tab (Explorer)

The current tab is useful for comparing the same sample processed through two
`pipeline` values. Follow-up work should preserve the post-join threshold handling
that separates truly absent calls from sub-threshold evidence, but generalize the
workflow beyond pipeline-only comparisons.

- [ ] **Generalize Pipeline comparison into an A/B strata comparison tab** — allow
  users to choose the comparison field (`pipeline`, `read_type`, `batch`,
  `sample_type`, `label1`–`label3`, or another categorical column) and then select
  two values to compare. Keep `pipeline` as the default preset because it is the
  common current workflow.
- [ ] **Add A/B comparison preflight checks** — before running the expensive join,
  show record counts by selected stratum, overlapping `sample_id` values, samples
  present only on one side, and warnings when `sample_id` is not unique across
  subject/timepoint/batch dimensions that are not part of the join.
- [ ] **Make join keys explicit and safer** — default to
  `(sample_id, chrom, pos, alt_allele)`, but warn when selected data has repeated
  sample IDs across `subject_id`, `sample_type`, `timepoint`, `batch`, or
  `read_type`. Consider allowing advanced users to include extra metadata columns
  in the join key.
- [ ] **Refactor comparison SQL construction into helpers** — centralize WHERE
  and literal escaping for the A/B join, selected-point drilldowns, and IGV
  conditions. Preserve the current post-join threshold reclassification semantics.
- [ ] **Prioritize the tab flow around "what changed?" before "why?"** — first
  surface concordance, true-unique, and sub-threshold counts; then expose VAF/depth
  correlations, unique-loci characterization, SBS96 spectra, and IGV review as
  follow-up sections.

---

## Coverage Explorer (Customer-facing)

- [ ] **Make interval/exon review the default coverage workflow** — the current
  Depth Profile tab is useful for drill-down, but customer coverage review should
  start from `coverage_intervals`: low-coverage interval table, exon/interval
  heatmap, and interval-level QC summaries. Keep base/bin depth profiles as the
  detail view opened from a selected gene/exon/interval.
- [ ] **Clarify and harden coverage profile bin semantics** —
  `load_expanded_depth_profile()` groups by `pos` after width-weighting bins. That
  is reasonable when bins align across samples, but can be misleading if bin
  boundaries differ. Either document/enforce aligned-bin input or move toward
  interval-overlap aggregation for mixed-bin datasets.
- [ ] **Broaden coverage profile tests** — add cases for empty filters, sparse
  discontinuous targets, mixed bin widths, invalid/zero-width bins, and non-aligned
  bins across samples.
- [ ] Surface `feature_type` and `exon_number` in the intervals heatmap — available
  from v0.4.0; use `exon_number` as the x-axis instead of `interval_name` when
  present, so exons line up consistently across genes regardless of naming
  conventions.
- [ ] Surface BEDGraph track columns in the intervals GC scatter — detect any
  `Float32` columns beyond the fixed schema in `coverage_intervals` and offer them
  as a color-by selector.
- [ ] `config/panel_example.toml` — example config file for operators deploying
  the Coverage Explorer for a specific panel; documents all supported `geac.toml`
  keys with comments.
- [ ] **Longitudinal tracking** — add `run_id` and `run_date` provenance fields
  to `CoverageArgs` and `CoverageRecord`; surface a run-over-run depth trend line
  in the Summary tab.

---

## Coverage Analysis — Explorer integration

Steps 1–10 of `geac coverage` shipped in v0.4.0 (see `docs/DEVELOPMENT_LOG.md`).
Remaining:

- [ ] **Step 11: Explorer — "Coverage" tab (DuckDB mode only)**
  - Systematically undercovered intervals table (from `coverage_intervals`):
    configurable depth threshold and fraction-of-samples; columns include gene,
    exon_number, mean_gc, mean_mappability; explains whether dropout is GC,
    mappability, or other.
  - Scatter plot: `mean_gc_content` vs `mean_depth` per interval — GC bias
    visible as a U-shaped curve; color by `mean_frac_mapq0` to overlay
    mappability.
  - Per-sample depth distribution histogram; `frac_dup`, `frac_overlap`,
    `frac_soft_clipped` summary bars for QC overview.
  - Per-exon coverage heatmap (genes as rows, exons as columns, color = mean
    depth or `frac_at_30x`) — the primary customer-facing view.

---

## Fragmentomics

Long-term: end-motif table comparing alt-supporting vs reference-supporting reads
as a bait-bias signal. Design context in `docs/DEVELOPMENT_LOG.md`.

- [ ] **Rework the Fragmentomics Explorer tab around QC-first workflows** — the
  current tab combines insert-size motif trends, FFT periodicity, motif-by-GC,
  fragment GC distribution, sample selection, grouping, smoothing, and layout
  controls in one dense view. Split or reorganize it into clearer subviews:
  `Fragment QC` (fragment count, median/IQR insert size, mono/di-nucleosome
  fractions, GC median/IQR), `End motifs` (motif frequency by insert size/GC),
  and `Periodicity` (FFT/nucleosome signal).
- [ ] **Add guardrails to the FFT / nucleosome-periodicity view** — keep it as an
  exploratory analysis, but add minimum-count checks, a visible "exploratory"
  label, and summary metrics for peak strength/background instead of only drawing
  a 10.5 bp reference line. Avoid making the plot look like a validated clinical
  signal without validation.
- [ ] **Harden Fragmentomics SQL construction** — `app/explorer/tabs/fragmentomics.py`
  interpolates selected sample IDs and motif lists into SQL strings. Route these
  through shared SQL literal escaping or parameterized DuckDB queries so unusual
  sample IDs/motifs cannot break queries.
- [ ] **Reduce control density in the Fragmentomics tab** — replace the seven-widget
  row with compact controls grouped by task, and hide expert controls (smoothing,
  FFT insert-size range, facet columns) in an expander or per-view settings area.
- [ ] Design the new table schema (per-read, not per-alt-read).
- [ ] Implement collection in `geac collect` (opt-in flag, like `--reads-output`).
- [ ] Add Explorer plot comparing end-motif frequency for alt vs ref reads.

---

## Fusion caller (experimental)

Full design, rationale, and the niche this caller aims to own are in
`docs/FUSION_DEVELOPMENT.md`. The highest-value next items (Tier 1, specificity):

### ▶ Resume here — FP diagnostics (updated 2026-06-12)

**Where we are.** v0.4.45 shipped tooling to attribute false `GENEA::GENEB` calls to
their cause: `shared-kmers` (index/reference level — fuzzy `--edit-distance` matching,
`--check-reference` copy counts, `--index` to flag real index k-mers) and
`diagnose-fusion` (read level — from a `fusions --reads-output` evidence BAM, gives a
homology-vs-junction coherence verdict, original alignment map, suspicious-k-mer table
flagging 1-substitution error paths, and a per-read A/B layout track). Design notes in
`docs/DEVELOPMENT_LOG.md` (2026-06-12); usage/interpretation in `docs/EXPERIMENTAL.md`
("Diagnosing false-positive fusion calls").

**What I'm doing now.** Running `diagnose-fusion` on a labeled set of fusions I believe
are *real* vs *false*, to learn which signals separate them and calibrate thresholds for
this panel. Discriminators: coherent fraction (real high / artifact ~0), minority-side
median anchor (real ≥3 / artifact 1–2), alignment map (real = two loci / artifact = one),
and the suspicious-k-mer table (`1edit_from_other=yes` → error path; high `ref_copies` →
repeat).

**Next actions when I restart:**
- [ ] **Add `diagnose-fusion --summary-tsv`** — emit one machine-readable row per run
  (`gene_a, gene_b, n_fragments, n_spanning, coherent_frac, median_anchor_a/b,
  pct_minority_kmers_1edit, median_ref_copies, verdict`) so the labeled real/false set
  can be swept and the separation plotted instead of eyeballing each report. (I offered
  this; not yet built.)
- [ ] Use the swept table to pick per-panel thresholds, then feed them back into the
  **"Reduce FP fusion calls"** item below (and ultimately into `fusions` filter defaults).

- [ ] **Overlap / adjacency filter — gene-body form** (breakpoint form shipped). The
  breakpoint form is done: `--min-breakpoint-distance` tags `filter=samelocus` when both
  breakpoints are same-chromosome within X bp (catches single-locus paralog leakage, e.g.
  GNA12→GNA13::GNA11; see `docs/DEVELOPMENT_LOG.md` 2026-06-17). Still open: the
  complementary *annotated-gene-body-overlap* check — reject pairs whose gene bodies
  overlap or sit within X bp, using the gene coordinates already in the index, as an
  upstream filter that does not depend on having spanning reads / a second BAM pass.
- [ ] **Fusion filter `chimera`/`samelocus`: drop vs tag for production.** Both
  breakpoint filters currently *tag* `filter` (rows kept), matching the PoN convention.
  Decide whether production should drop these rows instead (or have `geac merge` / the
  Explorer exclude non-`PASS` by default), so downstream consumers don't have to filter
  manually. Low-risk either way; needs a product decision, not new mechanism.
- [ ] **Split-read vs discordant-pair separation** — report `split_reads` and
  `discordant_mates` as distinct columns.
- [ ] **VAF-like quantification** — supporting reads as a fraction of local depth;
  gates the longitudinal low-AF monitoring niche.
- [ ] **Reduce FP fusion calls — broaden validation & set defaults.** Two mechanism-level
  specificity filters now ship (`--max-breakpoint-std`/`--min-breakpoint-reads` →
  `chimera`; `--min-breakpoint-distance` → `samelocus`) and gave clean separation on the
  validation sample (1 real PASS / 533–1252 artifacts tagged; 0/1177 in a normal). No PoN.
  Root causes confirmed: chimeric molecules (scattered breakpoints) + unindexed-paralog
  leakage (GNA12; single-locus co-located breakpoints), *not* index k-mer leakage. Still
  open: (a) validate on more cell-line BAMs (raw/simplex/duplex) and read-types to confirm
  the thresholds (`std 100`, `reads 5`, `distance 10000`) generalize — a second tumor (a
  canonical EWSR1::FLI1) showed the original `std 10` was too tight: a real high-depth
  junction spreads over tens of bp, so the recommendation was raised to `100` (artifacts
  still sit at 10⁴–10⁷ bp; see `docs/DEVELOPMENT_LOG.md` / `CHALLENGES.md` 2026-06-17);
  (b) decide whether to make the filters **on-by-default** with those values vs. opt-in as
  now; (c) feed any `diagnose-fusion --summary-tsv` sweep results back in. This gates
  benchmarking.
- [ ] **Robust breakpoint-tightness metric (replace raw std).** Raw `bp_*_std` is
  outlier-sensitive: a handful of stray reads among thousands inflate it even when the bulk
  converges, and the real-call operating point moved an order of magnitude between samples
  (≈2 bp → tens of bp). *Largely addressed:* `cluster_around_median` (±1 kb) trims far-off
  reads before std/median (v0.4.52), and the `strong_support` tier (v0.4.53) lets a
  well-supported call exceed the std ceiling so a single threshold no longer splits real
  junctions. Still open: a fully principled per-side metric — e.g. *fraction of reads within
  ±W bp of the modal breakpoint*, or MAD instead of std — and justifying the 1 kb cluster
  window and the 25-read / 250-bp tier constants across read-types (promote to flags?).
- [ ] **Repeat-buried junctions at low input.** The concordant-pair fallback + strong-support
  tier recover repeat-buried fusions (e.g. EWSR1::ERG) only at high input; at low input the
  concordant support falls below the tier and they go undetected. The root-cause lever is the
  index, not the breakpoint filter — a higher index edit distance recovers some single-locus
  leakage at source (see the GNA12 case), so evaluate whether a higher-edit-distance index
  also restores low-input sensitivity for repeat-buried junctions.
- [ ] **(Stretch) Promote rescued `a`/`b` windows to genuine supporting evidence.** Today the
  rescue is rendering-only. Counting edit-1 windows toward anchoring / spanning evidence could
  improve sensitivity on SNP-dense reads, but it changes quantitative output and filters — a
  separate, opt-in design with its own validation, NOT folded into the diagnostic track.
- [ ] **(Design — needs more thought) Layout track: "matches-both" marker + explicit breakpoint
  boundaries.** Two requests that turn out to be one idea — characterizing the *ambiguous zone*
  at the junction, which the current `A`/`B`/`a`/`b`/`N`/`.` track can't describe. Capturing the
  design before building; representation question below must be settled first.
  - **(1) "matches both A and B" marker** (e.g. `*`). A window whose k-mer is present in *both*
    partners' reference (cross-gene homology / microhomology). The index can't answer this today:
    shared k-mers are marked `MULTI_GENE` and **dropped** at build time, so such a window renders
    `.` and the edit-1 rescue sees it as a *split*. To mark it honestly we need the shared k-mers
    back. Two paths:
    - **(a) Keep them at build time** — write dropped `MULTI_GENE` k-mers to a small side table
      (`kmer_hash → set of gene indices`); shared k-mers are a small fraction on a good index, so
      cheap, and FL/diagnose just look them up. Needed for the cohort `fusions` FL pass.
    - **(b) Recompute at render time** — the index now stores each gene's `chrom`/`tx_start`/
      `tx_end`, so for a *single* pair you can fetch both genes' reference windows from the FASTA,
      build their full (non-deduped) k-mer sets, and intersect. Clean for `diagnose-fusion`; too
      heavy per-pair across a 255-gene cohort scan.
    This directly surfaces the homology that drove the ADVL1823 chr1q false pairs (NTRK1/BCAN).
  - **(2) Explicit breakpoint boundaries** (`AAAA||BBBB`, always two — A's confident end and B's
    confident start). The two ideas merge here: the bases *between* the pipes are either
    microhomology (match both → the `*` marker) or inserted/novel (match neither → `.`), which
    distinguishes blunt vs microhomology-mediated vs templated-insertion junctions. The current
    `.` gap of width ~`k-1` between the A and B blocks already *is* the (k-mer-fuzzy) breakpoint
    signature; pipes make it base-precise.
  - **Key unresolved question — window-space vs base-space.** The FL track is window-indexed (one
    char per k-mer start, so `k` shorter than the read and offset); a `|` is a *base-level* concept.
    `||` likely pushes toward a **base-aligned** track (one char per read base), which would also be
    easier to read in IGV than today's k-shifted track. Decide this deliberately before building —
    it's the foundation for both markers. Don't wedge base-level pipes into a window-level string.
  - **Other constraints.** Only junction-*spanning* reads have a breakpoint to mark (one-sided reads
    just get their block). Diagnostic/rendering only — must not feed assignment/filters. Precise
    per-read base-level breakpoint overlaps with what `breakpoints.tsv` already computes at the
    fusion level, so much of the value is a better *per-read view*, not a new measurement — be
    clear-eyed about that.
  - **Recommended phasing.** Prototype in `diagnose-fusion` first (single pair; can afford path (b);
    already human-facing; natural place to iterate on the visual and settle window-vs-base). Then
    decide whether to port to the cohort `FL` tag — which would likely want path (a), the shared-
    k-mers side table.
- [ ] **Benchmarking harness — phases 2 & 3.** Phase 1 (local scoring) shipped:
  `experimental/fusion_benchmark/score_fusions.py` + `experimental/fusion_benchmark/README.md` score per-sample `*.fusions.tsv`
  against a truth manifest → TP/FN/FP per sample, per-tier, recall-by-dilution, and 5'->3'
  orientation accuracy (CI-tested via `app/tests/test_score_fusions.py`). Remaining:
  - **Phase 2 — Terra WDL**: a workflow that scatters `GeacFusions` over a cell-line set and
    runs the scorer as a gather step, emitting the rollup + per-call TSV as outputs (so a
    cohort benchmark is one submission, not a manual localize + score).
  - **Phase 3 — GitHub Action**: trigger the Terra benchmark on demand / on release tag and
    post the rollup as a workflow summary, so "did this change move recall/precision?" is
    visible in CI. (Truth manifests carry sample identifiers — keep them in a private location,
    not this repo.)
- [ ] **(low priority) Strand-backfill helper for legacy indexes.** A
  `build-fusion-index`-independent helper that adds `strand`/`tx_start`/`tx_end` +
  `gene_strand=1` to an existing index's `genes` table from a gene annotation, keyed by
  gene_name+chrom, without re-extracting k-mers — so a pre-0.4.55 index gains orientation
  without a full rebuild. Low priority: new indexes are strand-aware, and the cohort indexes
  are being rebuilt fresh, so legacy-index backfill mostly won't be needed. (The load-time
  WARNING for strandless indexes already shipped — `fusions.rs` warns once at load that
  orientation is disabled and `partner_order` will be `index` cohort-wide.)
- [ ] **Build a fresh edit-distance-1 fusion index (with padding) and upload it to Terra.**
  The index currently live on Terra was built with geac 0.4.44 — no strand → orientation
  abstains cohort-wide (see the strand-warning item above). Rebuild the edit-distance-1 k-mer
  list on a current geac (so the `genes` table carries `strand`/`tx_start`/`tx_end` +
  `gene_strand=1`), with `--gene-padding` applied so junction-adjacent k-mers are captured.
  **Encode the building geac version in the index filename** (e.g.
  `brightseq_fusion_index_edit_1_v0.4.58.duckdb`) so a localized index is never ambiguous about
  which version produced it — the same provenance lesson from the `cohort_v0.4.52`-vs-`0.4.58`
  confusion. Upload to the Terra workspace and repoint the fusion configs at the versioned file.
  Re-running the cohort against it should turn `partner_order` from `index` into real 5'->3'
  labels with no change to detection/PASS counts (strand is label-only).
- [ ] **Repeat-aware fusion detection for repeat-resident partners (e.g. DUX4).** Cell line E
  (CIC::DUX4) is detected 0/18 even at 100% input: DUX4 lives in the D4Z4 macrosatellite, so
  the genome-unique k-mer index has only sparse, mostly-non-unique k-mers there and can't
  anchor the partner. The genome-uniqueness/`--edit-distance-filter` passes that make the index
  trustworthy elsewhere are exactly what starve a repeat-resident gene. Investigate a
  repeat-aware path: e.g. allow a higher `--max-genome-copies` *per flagged gene*, a
  repeat-tolerant anchor mode that accepts multi-copy k-mers when the partner side is unique,
  or split-read/soft-clip evidence anchored only on the unique (CIC) side. Validate it recovers
  CIC::DUX4 without inflating FPs from the repeat.
  - *Fast-iteration fixture:* before touching the detector, carve a tiny experimentation BAM
    from a 100%-input CIC::DUX4 sample (cell line E) containing every read pair where **either**
    mate matches **any** CIC or DUX4 k-mer — regardless of that k-mer's genome-uniqueness (i.e.
    do *not* apply the genome-copy/edit-distance filters that the main index uses; we want the
    repeat-buried reads that those filters would discard). Pad the gene intervals by **300 bp**
    on each side so junction-spanning and soft-clipped reads aren't truncated. The result is a
    small BAM holding all the reads any repeat-aware approach could possibly use, enabling the
    ~5s evidence-BAM iteration loop instead of re-scanning the full BAM each experiment.
- [ ] **Copy-aware "unique-anchor" specificity filter — v1 shipped + validated (opt-in); depth-relative threshold + soft-clip pairing open.**
  The complement to the repeat-aware item above: it buys *specificity* where repeat-tolerance buys
  *sensitivity*. For a repeat-resident partner (DUX4, and the other low-uniqueness panel genes —
  pseudogene/segdup families like PMS2, SDHA, NOTCH2, BRCA1, and the handful with ~zero unique
  k-mers) the multi-copy k-mers match reads genome-wide and manufacture a high-support call in
  essentially *every* sample, regardless of whether the fusion is present.
  - **Shipped (v1):** `geac experimental fusions --min-unique-anchor-reads N` (default 0 = off; see
    `src/fusions.rs`, `docs/EXPERIMENTAL.md`). The runtime index now retains a genome-unique k-mer
    set + per-gene unique totals (only when the flag is on); each supporting fragment is "unique-
    anchored" if any read carries a `genome_copies==1` k-mer on the **higher-uniqueness partner**
    (asymmetric — chosen per call from per-gene unique totals). A call PASSes only if ≥ N fragments
    are unique-anchored, else `filter=no_unique_anchor` (tag, not drop; only touches calls still
    PASS). Requires `--check-genome-uniqueness` index. Unit-tested (`read_unique_anchored`).
  - **Validated (evidence-BAM cohort, geac 0.4.59):** sweeping `N` on a repeat-partner cohort, a
    sufficient threshold removes **all** repeat-partner false positives (100% precision) while the
    surviving true calls are the **highest-input** replicates in proper **dose-order** — turning a
    dose-blind, high-FP artifact into clean dose-ordered detection. N=1/2 did nothing (every sample
    has a few unique-anchored fragments from normal expression of the unique partner paired with
    promiscuous repeat k-mers); the discriminating signal is the *count* (true high-input has tens,
    an FP a handful), so the threshold must sit above that noise floor.
  - **Still open:**
    - **Depth-relative threshold.** A fixed absolute `N` is depth-dependent (which is *why* it gives
      dose-response, but isn't portable across coverage). Add a relative variant — fraction of
      supporting reads, or scale by local coverage — before considering on-by-default.
    - **Pair with the unique-side soft-clip** (repeat-aware item) for the *sensitivity* half — v1
      gives specificity only; low-input repeat-partner fusions still drop out (they need soft-clip
      to be detected at all).
    - **Full-panel safety check.** Validation so far is a 2-gene (CIC/DUX4) index; confirm the filter
      doesn't harm fusions between two well-anchored (high-uniqueness) genes (expected safe — both
      sides unique-anchor trivially — but unverified on the full panel).
    - **Alternative scoring** to evaluate: soft copy-weighting (weight a supporting read by
      `1/genome_copies`) vs the current hard gate; how many unique k-mers should constitute an anchor
      (currently ≥1); a per-gene uniqueness floor below which *neither* side can anchor (soft-clip-only).
    - **DONE — `n_unique_anchored` is now an output column** (`--emit-unique-anchor`, or implied by
      `--min-unique-anchor-reads > 0`). A single run records the raw count so the whole
      precision/recall curve is recoverable offline as `n_unique_anchored / supporting_reads`.
      See `docs/DEVELOPMENT_LOG.md`.
- [ ] **`genome_copies` saturates at 255 (`u8`).** `build-fusion-index` counts genome-wide k-mer
  occurrences in a `HashMap<u64, u8>`, so any k-mer occurring ≥255× is right-censored at 255 (the
  `kmers.genome_copies` column and `gene_stats.tsv` copy tiers all cap there). The `==1` unique
  count is unaffected, but true copy numbers for high-repeat k-mers (DUX4/D4Z4, satellite, etc.) are
  lost — which matters for copy-weighted anchoring and for honest per-gene copy reporting. Widen the
  counter to `u16`/`u32` (saturating), and the column accordingly.

---

## Code health / tech debt

Captured 2026-05-01 from a code review. Items here are not blocking any release;
they are organizational/structural improvements to take on between feature
milestones. Each "what to do" item below maps to a "why" weakness for context.

### Weaknesses observed

- **A few Python files are still oversized** — the original ~7,853-LOC
  `app/geac_explorer.py` has largely been decomposed into real tab modules under
  `app/explorer/tabs/` (the entrypoint is now ~2,000 LOC). The remaining offenders
  are `app/explorer/tabs/reads.py` (~1,345 LOC, a single `render()` of nested
  closures — see `docs/CODE_AUDIT.md` 2026-06-10) and `geac_coverage_explorer.py`
  / `error_spectrum.py` (~1,400–1,550 LOC each). Target: no single file >1,500 LOC.
- **No app-level Python tests** — `app/tests/` covers ~400 LOC of helper
  functions, but neither `geac_explorer.py` nor `geac_coverage_explorer.py` has
  integration tests.
- **No reproducibility pinning** — there is no `rust-toolchain.toml`, no
  `=`-pinned versions in `Cargo.toml` (e.g. `rust-htslib`, `parquet`, `duckdb`),
  and no `pyproject.toml`/`ruff`/`mypy` config.
- **CLI metadata-flag boilerplate** — `src/cli.rs` has ~200 lines of copy-paste
  for the optional `--label1` / `--label2` / `--label3` style metadata flags.
- **Sparse module-level docs in Rust** — most `.rs` files have function-level
  doc comments but no `//!` describing module intent, invariants, or how the
  module fits into the larger pipeline.
- **No lint enforcement in CI** — no `clippy`, `rustfmt`, or `ruff` step gates
  merges. Style and lint regressions can land silently.
- **Schema contract is one-way at runtime** — `app/explorer/schema.py` loads
  `schema/geac_schema.json` at startup, but nothing fails the build if a column
  added to the schema is never referenced by the app, or vice versa.

### Proposed work items

- [ ] **Decompose `app/geac_explorer.py`** — introduce a `BaseTab` abstract class
  with a `render()` method; move each tab's logic into its own file under
  `app/explorer/tabs/` (the directory already exists). Goal: no single Python
  file >1,500 LOC. Unlocks unit tests of tab logic without a Streamlit runtime.
- [ ] **Pin Rust toolchain and critical dependencies** — add `rust-toolchain.toml`
  with the channel currently used in CI; switch `rust-htslib`, `parquet`, `arrow`,
  and `duckdb` in `Cargo.toml` to `=`-pinned versions. Document the bump cadence
  in `CLAUDE.md`.
- [ ] **Generate schema bindings at build time** — add a `build.rs` that reads
  `schema/geac_schema.json` and emits a `const SCHEMA_COLUMNS: &[&str]` per
  table; add a Python test that fails if any schema column is unreferenced in
  the explorer code.
- [ ] **Add module-level `//!` docs to all public Rust modules** — at minimum:
  `src/coverage/mod.rs`, `src/bam/pileup.rs`, `src/merge.rs`, `src/track.rs`,
  `src/writer/*`. Each should describe module purpose, key invariants, and
  dependencies.
- [ ] **Add lint/format CI** — new GitHub Actions job that runs `cargo fmt --check`,
  `cargo clippy -- -D warnings`, and `ruff check app/` on every PR. Fix existing
  warnings before turning the gate on.
- [ ] **Collapse CLI metadata flag boilerplate** — replace `--label1`/`--label2`/
  `--label3` (and any parallel groups) in `src/cli.rs` with a single repeatable
  `--label KEY=VALUE` flag backed by `Vec<(String, String)>`. Update WDLs and
  README.
- [ ] **Add Streamlit integration smoke tests** — at minimum, a test that imports
  each tab module and instantiates its render function against a tiny synthetic
  Parquet fixture. Catches import-time and obvious-render regressions without
  requiring a full browser.
- [ ] **Add an end-to-end Rust→DuckDB→Python test** — run `collect` then `merge`, open the
  resulting DuckDB from Python, and assert a known locus's counts. The cross-layer contract
  is currently verified only piecewise (Rust output vs schema, and Python reads, but never
  the full chain). See `docs/CODE_AUDIT.md` 2026-06-10.
- [ ] **Extract a shared Parquet-writer helper** — the 8 modules under `src/writer/` repeat
  the same `File::create` → `WriterProperties(SNAPPY)` → `ArrowWriter::try_new` → `write` →
  `close` shell, varying only the schema and record-to-batch functions. Collapse into one
  `write_parquet_generic<T>(records, path, schema_fn, batch_fn)` (~150-200 LOC, single
  error/compression policy).
- [ ] **Migrate this backlog into GitHub Issues** — once the team grows beyond a
  single contributor. Keep `ROADMAP.md` for milestones; let issues hold
  individual work items with priority/milestone labels.
