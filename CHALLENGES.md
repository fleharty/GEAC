# Development Challenges

A running log of non-obvious problems encountered during GEAC development, how they
were diagnosed, and what the fix was. Updated as new issues arise.

---

## Rust / Core Tool

### Segfault on BAM contigs absent from FASTA index (e.g. hs37d5 decoy)
**Symptom:** `geac collect` (and `geac coverage`) crashed with a segfault after printing:
```
[E::faidx_adjust_position] The sequence "hs37d5" was not found
```
when run on a DRAGEN BAM aligned to GRCh37 with the `hs37d5` decoy contig, using a
FASTA that did not include that decoy sequence.
**Root cause:** `RefCache::get` called `faidx::Reader::fetch_seq` for a contig present
in the BAM header but absent from the FASTA index. htslib's internal `faidx_adjust_position`
printed the error message and then crashed instead of returning a recoverable error, so
the Rust `?` error-handling path was never reached.
**Fix:** `RefCache::new` now calls `fai.seq_names()` at startup to build a `HashSet` of
all sequences in the FASTA index. In `get()`, any contig not in that set returns `'N'`
and logs a one-time `WARN` — the existing `if ref_base == 'N' { continue; }` in the
pileup loop then skips those positions cleanly.
**Lesson:** Never call htslib `fetch_seq` for a sequence that might not be in the index.
Check existence first; htslib does not return a safe Rust error for missing sequences in
all code paths.

### Coverage memory growth on large genomes
**Symptom:** `geac coverage` on a whole-genome BAM would accumulate all
`CoverageRecord` objects in a `Vec` before writing, consuming multiple GB of RAM.
**Root cause:** The original design collected all records, sorted them, then wrote
the Parquet file in one shot.
**Fix:** Replaced with a streaming `CoverageWriter` that buffers 100k records at a
time and flushes to disk as it goes (`src/writer/parquet_coverage.rs`). The sort was
also dropped — pileup order is already genomic order.
**Lesson:** For genome-scale outputs, never buffer the full result set in memory.

### `fwd_alt_count`/`rev_alt_count` always 0 for reverse-strand reads in standard paired-end data
**Symptom:** `rev_alt_count` was always 0 for any paired-end BAM processed with the raw
pipeline. `fwd_alt_count` matched `alt_count`, making strand-bias detection useless for
overlapping-pair libraries.
**Root cause:** In `tally_pileup` and `tally_indels`, the overlapping-pair branches derived
a single `r1_is_rev` flag (from the biological R1's `is_reverse` bit) and used it to
assign the *entire fragment* to either `fwd` or `rev`. For standard paired-end data (R1
forward, R2 reverse), `r1_is_rev` is always `false`, so every overlapping pair went to
`fwd_*` regardless of which read actually carried the alt allele.
The design intent was: `alt_count`/`ref_count` are fragment-level (deduplicated molecules);
`fwd_alt_count`/`rev_alt_count` are read-level (each read counted on its own strand). The
implementation never matched this intent.
**Fix:** In each overlapping-pair branch of `tally_pileup` (and symmetrically in
`tally_indels`), replaced `r1_is_rev` with the per-read `is_reverse` flag:
- N + informative: use `r2.is_reverse`
- informative + N: use `r1.is_reverse`
- both agree (alt or ref): **two** increments — one for `r1.is_reverse`, one for `r2.is_reverse`
- b1=ref, b2=alt: use `r2.is_reverse` (R2 carries the alt)
- b1=alt, b2=ref: use `r1.is_reverse` (R1 carries the alt)
- two different alts: `r1.is_reverse` for t1, `r2.is_reverse` for t2
After the fix, `fwd + rev ≠ total` for overlapping pairs is expected and correct — for
"both agree alt", `fwd + rev = 2` while `total = 1`. The `r1_is_rev` variable is still
retained for `fwd_depth`/`rev_depth`, which remain fragment-level.
**Validation:** The R2-only artefact (`r2_artefact.bam`) is the most compelling check:
before the fix it showed `fwd_alt_count=2, rev_alt_count=0`; after the fix it correctly
shows `fwd_alt_count=0, rev_alt_count=2`, directly reflecting that only R2 reads carry
the alt allele.
**Lesson:** Fragment-level and read-level tallies require separate strand attribution logic.
When both concepts live in the same tally function, verify that each increment site uses
the per-read flag, not a single fragment-derived flag.

### `alt_reads` table missing all insertion and deletion records
**Symptom:** Family size (and other per-read) filters in the Explorer had no visible
effect on insertion and deletion VAF distributions. Debug output confirmed the
`alt_reads` table contained zero insertion or deletion rows, only SNVs.
**Root cause:** `tally_indels()` in `src/bam/mod.rs` returned only a
`HashMap<String, IndelCount>` (aggregate counts). No `AltRead` records were ever
pushed for indel-supporting reads. The code that populated `AltRead` records for
SNVs (using `read_details` from `tally_pileup`) had no equivalent for indels.
**Fix:** Extended `tally_indels()` to also accept a `collect_reads: bool` parameter
and return a second `HashMap<String, Vec<ReadDetail>>` keyed by alt allele. In the
first pass each alignment now also captures a `ReadDetail` (qpos, read length, base
qual, map qual, fgbio tags, insert size) when `collect_reads` is true and the allele
is non-None. In the second pass the detail is stored alongside the count increment,
split by case (single read / agree overlap / disagree overlap / multi-read). The
caller now iterates `indel_read_details` and pushes one `AltRead` per entry.
**Lesson:** When adding per-read output for SNVs, explicitly audit whether the same
path is needed for the indel tally — they are separate code paths.

### Rust import path: `ReadType` not in scope in submodule
**Symptom:** Compiler error when `BinAccumulator` in `src/coverage/mod.rs` tried to
reference `crate::cli::ReadType`.
**Root cause:** `ReadType` is defined in `crate::record`, not `crate::cli`. The CLI
module re-exports it for clap parsing but the canonical definition lives in `record`.
**Fix:** Changed import to `use crate::record::ReadType`.

### Edit tool conflict: ambiguous surrounding context
**Symptom:** An edit to add `adaptive_depth_threshold` to `CoverageArgs` in
`src/cli.rs` failed because two structs had identical surrounding lines, making
the match non-unique.
**Fix:** Extended the `old_string` context to include a nearby unique field
(`min_depth`) so the replacement target was unambiguous.

---

## DuckDB Query Engine

### Internal assertion: `inequal types (BIGINT != VARCHAR)` on duplicate complex subquery
**Symptom:** The "Cohort artefact vs rare variant: family size comparison" section threw
a DuckDB `InternalException: INTERNAL Error: Failed to bind column reference "pos" …
inequal types (BIGINT != VARCHAR)` at runtime, even though both `alt_bases.pos` and
`alt_reads.pos` are declared `Int64` in the Parquet schema.
**Root cause:** The same complex `table_expr` subquery (which itself contained an inner
JOIN to `alt_reads`) was inlined verbatim twice in the same `WITH` block — once for
`locus_counts` and once for the `_filt` inner join. DuckDB's binder failed to reconcile
column types across the two independently expanded copies of the subquery, producing a
spurious type-mismatch assertion.
**Fix:** Materialized `table_expr` into a single leading CTE (`_base`) that is referenced
by name in both `locus_counts` and `_filt`. Also added `CAST(pos AS BIGINT)` in `_base`
to pin the join-key type unambiguously, regardless of what the source expression reports.
**Lesson:** Never inline the same complex subquery expression more than once in a `WITH`
block. Assign it a named CTE so the engine resolves types once and reuses.

---

## Python / Streamlit Explorer

### Depth distribution overrepresented low-depth bins
**Symptom:** The depth distribution histogram in the Coverage Explorer appeared
to have far more low-coverage positions than expected.
**Root cause:** With `--bin-size > 1`, each row in the Parquet represents multiple
genomic positions. Counting rows (`COUNT(*)`) treated a 100bp bin the same as a
1bp position, massively underweighting high-depth bins.
**Fix:** Added a `bin_n` column to `CoverageRecord` tracking positions per bin.
All locus-count queries now use `SUM(bin_n)` instead of `COUNT(*)`.

### Filter defaults silently excluding records at startup
**Symptom:** With no filters active, the "Filtered" Alt records count was
slightly lower than the "Overall" count (e.g. 299,487 vs 299,669). The
discrepancy persisted after clearing all filters.
**Root cause:** Four filter conditions were unconditionally injected into
every query, even when at their neutral defaults:
- `alt_count >= 1` — excluded records with `alt_count = 0` (rare but possible)
- `alt_count * 1.0 / total_depth BETWEEN 0.0 AND 1.0` — excluded records
  where `total_depth = 0` (VAF → inf) or `alt_count > total_depth` (VAF > 1)
- `homopolymer_len BETWEEN 0 AND 20` — excluded records with `homopolymer_len > 20`
- `str_len BETWEEN 0 AND 50` — excluded records with `str_len > 50`

The NULL issue for the repeat columns was identified first and partially fixed
(wrapping with `IS NULL OR`), but the ceiling truncation was still active. The
VAF range condition was the next fix, but the repeat ceiling issue persisted
until a third pass.

**Fix:** Each condition is now only added when the user has actually moved it
away from its default: `min_alt > 1`, `vaf_range != (0.0, 1.0)`,
`homopolymer_range != (0, 20)`, `str_len_range != (0, 50)`. The `where`
clause falls back to `"TRUE"` when no conditions are active.

### VAF distribution charts empty for insertion and deletion
**Symptom:** After fixing filter defaults (removing the always-on VAF BETWEEN
condition), the insertion and deletion VAF distribution charts appeared visually
but showed no bars. SNV worked correctly. Deselecting SNV from the variant type
filter made insertion/deletion charts render correctly.

**Root cause (first attempt — wrong):** Assumed the three `st.altair_chart(...,
on_select="rerun")` calls in the same render pass were colliding because all
used the same Vega-Lite selection name (`"bar_click"`). Fixed by using unique
names per variant type (`bar_click_SNV`, etc.) — did not resolve the issue.

**Root cause (second attempt — wrong):** Assumed Streamlit only properly
initialises the first `st.altair_chart(..., on_select="rerun")` call per pass.
Rewrote to collect all three sub-charts and render as a single `alt.vconcat`
spec. This broke all three charts instead of just two.

**Root cause (actual):** Removing the always-on `alt_count * 1.0 / total_depth
BETWEEN 0.0 AND 1.0` filter exposed records where `total_depth = 0` (VAF → inf)
or `alt_count > total_depth` (VAF > 1.0). These records produce `vaf_bin` values
outside [0, 1] or non-finite. Altair sanitizes `inf` to `null`; Vega-Lite then
renders no bar for that datum, making the chart appear empty. Insertion/deletion
records were more likely to have this edge case than SNV records in the dataset.

**Fix:** Reverted to separate `st.altair_chart` calls; added three guards to the
VAF distribution query:
```sql
AND total_depth > 0
AND alt_count <= total_depth
HAVING vaf_bin IS NOT NULL AND vaf_bin >= 0.0
```

Also fixed `_to_spec96_strat` and `_strat_sbs96_chart` being defined inside
`if not raw.empty:` in the Error Spectrum tab, making them unavailable in the
Reads tab. Moved both definitions above `_trinuc_available` so they are always
defined.

### Top-level Explorer navigation refactor stranded shared helpers inside Summary tab
**Symptom:** After replacing top-level `st.tabs(...)` navigation with a session-state-backed
selector to preserve the active section across reruns, several non-Summary tabs started
throwing `NameError` exceptions:
- `NameError: name 'query_records' is not defined` in the VAF Distribution tab when clicking
  insertion/deletion AF bars
- `NameError: name '_to_spec96_strat' is not defined` in the Reads tab's
  "Family-size stratified Spectrum"

**Root cause:** The refactor changed top-level sections from `with tab_x:` containers to
conditional `if _active_main_tab == ...:` blocks. Two helpers that had historically been
shared across tabs were still defined inside tab-specific blocks:
- `query_records()` and the `_table_cols` setup were defined inside the Summary section
- `_to_spec96_strat()` / `_strat_sbs96_chart()` were defined inside the Error Spectrum section

Under `st.tabs()`, every tab body reran, so those definitions were executed every time and the
scope bug was masked. Once rendering became conditional, opening VAF Distribution or Reads
directly no longer executed Summary or Error Spectrum first, so the helpers were undefined.

**Fix:** Moved `query_records()`, `_table_cols` setup, `_to_spec96_strat()`, and
`_strat_sbs96_chart()` back to shared top-level scope. Tab bodies now only contain rendering
logic; shared query/chart helpers are defined unconditionally before any per-tab branches.

**Lesson:** `st.tabs()` eagerly evaluates every tab body, which can accidentally hide bad
scoping. When converting tab containers into conditional rendering, audit every helper and
shared variable previously defined "inside a tab" and promote cross-tab dependencies to
top-level scope first.

### Re-aggregation mode: COALESCE couldn't distinguish "no reads" from "no reads passing filter"
**Symptom:** In `recompute_vaf=True` mode, loci where every read failed the per-read filter
showed the original (unfiltered) `alt_count` instead of 0.
**Root cause:** The LEFT JOIN to `alt_reads` used `COALESCE(filtered_alt_count, 0)` where
`filtered_alt_count = COUNT(*) FILTER (WHERE <filter>)`. When all reads fail the filter,
`filtered_alt_count` is 0 and the COALESCE correctly returns 0 — but when the locus has
*no rows at all* in `alt_reads` (e.g. indels, which were not yet collected), the LEFT JOIN
produces NULL for `filtered_alt_count` and the COALESCE also returns 0, which is wrong —
the original `alt_count` should be preserved. The two cases were indistinguishable.
**Fix:** Added a sentinel column `TRUE AS has_reads` to the `alt_reads` aggregate subquery.
The outer CASE now branches on `ar_agg.has_reads IS NULL` (no rows → preserve original
`alt_count`) vs `has_reads IS TRUE` (rows exist → use `COALESCE(filtered_alt_count, 0)`).
**Lesson:** A LEFT JOIN that aggregates with `COUNT(*) FILTER` cannot distinguish "no rows
joined" from "rows joined but none passed" via the count alone. A sentinel boolean column
in the aggregate is needed to distinguish the two cases.

### Family-size stratified spectrum silently bypassed per-read filters
**Symptom:** The family-size stratified SBS96 spectrum (singleton vs multi-member) in the
Reads tab classified loci using all reads regardless of the active per-read filter. Setting
a family-size filter (e.g. `family_size >= 2`) did not change the singleton/multi
classification even though it changed all other per-read plots.
**Root cause:** The `locus_fs` CTE queried `alt_reads` with `WHERE family_size IS NOT NULL`
but did not apply `_reads_where`. The per-read filter was threaded through every other
query in the Reads tab but was missed in this CTE.
**Fix:** Appended `AND {_reads_where}` to the `locus_fs` CTE when `_reads_active` is True.
One line change.
**Lesson:** When adding a new CTE that queries `alt_reads`, explicitly check whether the
active `_reads_where` clause should also be applied. Missing it produces silently incorrect
results with no error or warning.

### R1/R2 + Sample combined grouping: `KeyError: 'sample_id'`
**Symptom:** Checking "Show R1/R2" while "Color by = Sample" was selected raised
`KeyError: 'sample_id'` in the normalization step of the Read position bias and
Mean base quality by cycle plots.
**Root cause:** Three separate variables — `_dfe_select_expr` (SQL SELECT), `_dfe_group_expr`
(SQL GROUP BY), and `_dfe_label_col` (Python column name used for `groupby` and Altair
color encoding) — were each derived from independent ternary chains that evaluated the
grouping flags in different priority orders.  With the old "R1/R2 overrides Color by"
logic, the chains were kept consistent by adding `and not _dfe_by_read` to `_dfe_by_sample`
and `_dfe_by_batch`.  When that guard was removed to allow combined grouping, the chains
diverged: `_dfe_label_col` resolved to `"sample_id"` (because `_dfe_by_sample` came first),
but `_dfe_select_expr` and `_dfe_group_expr` resolved to the R1/R2 read expression (because
`_dfe_by_read` came first in those chains).  The SQL query therefore returned a `read`
column, not `sample_id`, and `_dfe_df.groupby("sample_id")` raised a KeyError.
**Fix:** Replaced all three ternary chains with a single explicit `if/elif` block covering
all six grouping combinations (batch+read, sample+read, read-only, batch-only, sample-only,
aggregate).  For the combined cases, a `label` column is built in SQL using string
concatenation (e.g. `ar.sample_id || ' ' || CASE WHEN ar.is_read1 THEN 'R1' ELSE 'R2' END AS label`)
and grouped via DuckDB's alias GROUP BY.  All three derived variables are always assigned
in the same branch, keeping them structurally consistent.
**Lesson:** When multiple derived variables are computed from the same set of boolean flags
via independent ternary chains, any change to one chain's priority order silently diverges
from the others.  Replace parallel ternary chains with a single branching block so each
case assigns all dependent variables together.

### Gene bar chart click-to-drill-down not working
**Symptom:** Clicking a bar in the "Affected loci per gene" chart appeared to
trigger a Streamlit rerun, but the detail table never appeared.
**Root cause (first attempt):** Used `event.selection.point` — but
`event.selection` is an `AttributeDictionary`, not an object with a `.point`
attribute. This raised `AttributeError`.
**Partial fix:** Changed to iterate `event.selection.values()` to find a non-empty
list. The rerun now completes without error but the table still doesn't render
reliably.
**Status:** Unresolved — tracked as a known issue for 0.4.0. The exact structure
of Altair point-selection events in Streamlit's `on_select="rerun"` API is
unclear; needs further investigation.
**See also:** The general `on_select="rerun"` debugging recipe below.

### Sample recurrence filter causes intermittent Position drill-down failures
**Symptom:** After adding a sample recurrence slider (filters loci by how many
samples carry a given alt allele), the Position drill-down table appears
inconsistently — reliably at low recurrence values (e.g. 7–12 samples) but only
about 1-in-5 clicks at high values (e.g. 52–73 samples). The drill-down worked
every time before the recurrence feature was added.

**Root cause 1 — session state out of range:**
When switching between datasets with different cohort sizes, the stored
`sample_recurrence` session state value (e.g. `(1, 73)`) can exceed the new
dataset's `max_value` (e.g. 8 samples). Streamlit raises a silent error rendering
the slider, which prevents the page from executing past that widget — so the
drill-down section is never reached. This explains the "sometimes" breakage
independent of recurrence value.
**Fix:** Clamp the stored session state to `[1, _n_samples_total]` before each
slider render:
```python
_sr = st.session_state["sample_recurrence"]
st.session_state["sample_recurrence"] = (
    max(1, min(_sr[0], _n_samples_total)),
    max(1, min(_sr[1], _n_samples_total)),
)
```

**Root cause 2 — expensive GROUP BY re-executes on every rerun:**
The recurrence condition was a self-referencing subquery embedded directly in the
`WHERE` clause:
```sql
(chrom, pos, ref_allele, alt_allele) IN (
    SELECT chrom, pos, ref_allele, alt_allele FROM alt_bases
    GROUP BY chrom, pos, ref_allele, alt_allele
    HAVING COUNT(DISTINCT sample_id) BETWEEN {lo} AND {hi}
)
```
This full-table GROUP BY runs on every Streamlit rerun — including the rerun
triggered by clicking a row to open the drill-down. At high recurrence values the
query takes several seconds. The user, seeing no immediate response, clicks again.
The second click triggers another rerun that clears the first rerun's row
selection, so the drill-down never appears. This is why it works at low recurrence
(fast query) but fails at high recurrence (slow query).
**Fix:** Wrapped the computation in `@st.cache_data` keyed on `(path, lo, hi)`.
The GROUP BY now runs only when the slider values change; row-click reruns use
the cached result registered as a DuckDB in-memory view:
```python
@st.cache_data
def _compute_recurrence_loci(path, sr_lo, sr_hi): ...
_rec_df = _compute_recurrence_loci(path, _sr_lo, _sr_hi)
con.register("_recurrence_loci", _rec_df)
conditions.append("(chrom, pos, ref_allele, alt_allele) IN "
                  "(SELECT ... FROM _recurrence_loci)")
```

**Root cause 3 — subquery used reads-filtered table_expr:**
The original subquery used the current `table_expr`, which (when per-read filters
are active) is a complex multi-table JOIN subquery. Embedding this inside the
recurrence IN-subquery doubled the join cost and introduced edge cases.
**Fix:** Saved `_base_table_expr = table_expr` before the reads-filter
reassignment and used it in the recurrence subquery. Sample recurrence should count
across the raw data regardless of per-read filter state.

**Root cause 4 — missing `key=` on `st.dataframe` (PRIMARY):**
The main data table used `on_select="rerun"` but had no `key=` parameter:
```python
_tbl_event = st.dataframe(
    df[_table_cols],
    on_select="rerun",
    selection_mode="single-row",
    # no key= !
)
```
Without a stable `key`, Streamlit auto-generates one from the widget's position in
the render tree. Any change above the widget — updated record counts, reworded
captions, new sidebar filters (like the sample recurrence slider) — shifts the
auto-key. On the rerun triggered by clicking a row, Streamlit cannot match the
incoming selection event back to the widget because the key has changed. The widget
returns an empty selection, so the drill-down section sees no selected row and
renders nothing.

This is the same class of bug documented below in the `st.altair_chart` section
(the AB vs BA heatmap fix). The sample recurrence feature introduced new
conditional UI elements above the data table (stats captions, recurrence slider),
which made the auto-key unstable on most reruns — explaining why the drill-down
"worked before recurrence was added" and failed intermittently after.

**Fix:** Added explicit `key=` to every widget using `on_select="rerun"`. A
systematic audit found 7 widgets that needed keys:
```python
st.dataframe(..., key="main_data_table")         # Position drill-down table
st.dataframe(..., key="cohort_data_table")        # Cohort tab table
st.altair_chart(..., key="vaf_chart_{vtype}")     # VAF distribution charts
st.altair_chart(..., key="sbs96_spectrum")        # SBS96 spectrum
st.altair_chart(..., key="sbs96_r1")              # SBS96 R1
st.altair_chart(..., key="sbs96_r2")              # SBS96 R2
st.altair_chart(..., key="snv_error_spectrum")    # SNV error spectrum
st.altair_chart(..., key="strand_bias_scatter")   # Strand bias scatter
```

**Lesson learned:** *Every* widget that uses `on_select="rerun"` must have an
explicit `key=`. This should be treated as a hard rule, not a nice-to-have. The
failure mode is silent (empty selection, no error), making it difficult to diagnose
unless you know to look for it. The symptom worsens as more dynamic content is
added above the widget, which is why it appeared to be caused by the sample
recurrence feature when the underlying issue was pre-existing.

**Current status:** After all four fixes the Position drill-down works reliably.
Root cause 4 was the primary remaining issue — the first three fixes addressed
real problems (stale session state, slow queries, wrong table expression) but the
missing `key=` was responsible for most of the observed failures.

**Possible future approaches:**
- Enforce a lint or code review rule: `on_select="rerun"` → must have `key=`.
- Pre-materialise recurrence loci into a DuckDB temp table at filter-change
  time and use it as a plain table join rather than an IN-clause.
- Use a JOIN instead of IN for better query planning on very large cohorts.

### `st.altair_chart(on_select="rerun")` selection returns `{}` instead of datum
**Symptom:** Clicking a chart cell or bar triggers a Streamlit rerun and the
correct Altair selection name key is present in `event.selection`, but its value
is an empty dict `{}` rather than a list of selected datums. The visual selection
highlight (opacity change) works correctly client-side. No error is raised.

**Root causes and fixes (in order of likelihood):**

1. **Missing `key=` on `st.altair_chart`** *(most common — fixed the AB vs BA heatmap)*
   Without a stable `key`, Streamlit generates an auto-key based on render-tree
   position. Inside conditional blocks or tabs, the auto-key can differ between
   the initial render and the post-click rerun, so Streamlit cannot match the
   incoming selection event to the widget and returns `{}`.
   **Fix:** Add `key="some_unique_string"` to every `st.altair_chart` that uses
   `on_select="rerun"`.

2. **`mark_rect` with `fields=` in `selection_point`** *(obscure Vega-Lite behaviour)*
   `mark_rect` with `fields=["x_field", "y_field"]` sends the event but not the
   field values on some Streamlit/Altair version combinations. Without `fields=`
   (index-based selection) the datum is also not sent.
   **Fix:** Adding `key=` (point 1) resolved this entirely — `mark_rect` +
   `fields=` works correctly once the widget has a stable key.

3. **Multiple charts with the same Vega-Lite selection name**
   Two `st.altair_chart` calls in the same render pass sharing the same Altair
   `selection_point(name=...)` string can collide, causing both to return `{}`.
   **Fix:** Use unique `name=` values per chart (e.g. `"sel_r1_click"`,
   `"sel_r2_click"`).

**Debugging recipe:**
```python
ev = st.altair_chart(chart, on_select="rerun", key="my_chart")
st.warning(f"DEBUG sel: {ev.selection!r}")          # show raw state in UI
import sys; print(ev.selection, file=sys.stderr)    # also to terminal
```
The debug box shows:
- `{}` — `on_select` never fired (no click yet, or `key=` problem → fix point 1)
- `{'my_sel': {}}` — event fired but datum empty → try `fields=` / mark type (point 2)
- `{'my_sel': [{'field': value, ...}]}` — working correctly

**Confirmed working combinations (Streamlit 1.55.0, Altair 6.0.0):**
- `mark_rect` + `fields=["ab_count", "ba_count"]` + `key=` ✓
- `mark_point(shape="square")` + `fields=` + `key=` ✓
- `mark_bar` + `fields=["sbs_label"]` + `key=` ✓

### Drill-down locus changes when toggling "Same alt allele only" checkbox
**Symptom:** Clicking the "Same alt allele only" checkbox in the Position
drill-down section would sometimes change the drill-down to an entirely
different locus (different chrom, pos, and alt allele). The first toggle
usually worked, but subsequent toggles jumped to unrelated positions.

**Root cause 1 — non-deterministic row order:**
`query_records()` had no `ORDER BY` clause. DuckDB returned rows in
arbitrary order that could change between executions. When the checkbox
triggered a Streamlit rerun, `st.dataframe` reported the same row *index*
(e.g. `[2]`) as still selected, but row 2 now pointed to a completely
different locus because the underlying `df` had been re-queried in a
different order. The guard comparing `(chrom, pos, alt_allele)` tuples
correctly detected a "new" locus and dutifully overwrote the persisted
state — with the wrong position.

Debug output confirmed this directly:
```
[DEBUG] persisted _drill_locus=('10', 16232465, 'G')
[DEBUG] row index=2, row locus=('16', 16136558, '-A')
[DEBUG] UPDATING _drill_locus → ('16', 16136558, '-A')
```
Same row index, completely different data.

**Fix:** Added `ORDER BY chrom, pos, alt_allele, sample_id` to
`query_records()` so row indices are stable across reruns.

**Root cause 2 — drill-down gated on transient dataframe selection:**
The entire drill-down block was inside `if _selected_rows:`, which
depended on `st.dataframe`'s `on_select="rerun"` event. When the checkbox
triggered a rerun, the dataframe could lose its selection state, hiding
the drill-down (and the checkbox) entirely — making it appear to uncheck
itself.
**Fix:** Persist the selected `(chrom, pos, alt_allele)` in
`st.session_state["_drill_locus"]` when a row is clicked. The drill-down
renders from the persisted state, not the transient selection event.

**Lesson:** Any query feeding a `st.dataframe` with `on_select="rerun"`
must have a deterministic `ORDER BY`. Without it, the row index reported
by Streamlit is meaningless across reruns — it's a pointer into a
shuffled deck. This is especially dangerous because the failure is silent:
the drill-down confidently renders the wrong locus with no error.

---

## IGV.js Integration

### "Access Unauthorized" on initial BAM load
**Symptom:** IGV.js displayed "Access Unauthorized" immediately when loading a
`gs://` BAM URL.
**Root cause:** The ADC access token was fetched server-side via `google.auth.default`
but not passed to the IGV.js browser instance.
**Fix:** Injected the token three ways in the generated HTML: globally via
`igv.setOauthToken()`, in the browser-level `oauthToken` option, and per-track
`oauthToken` field. Host pattern changed from `"*.googleapis.com"` to
`"storage.googleapis.com"` to match actual GCS request hostnames.

### "Access Unauthorized" on zoom / range requests
**Symptom:** Initial load worked, but zooming in triggered new HTTP range requests
that returned 401.
**Root cause:** Token was set globally but subsequent range requests did not pick
it up. Only the per-track `oauthToken` field reliably covers all requests.
**Fix:** Ensured `oauthToken` is set on every track dict in addition to the global
call.

### IGV.js page freeze on load
**Symptom:** Entire browser tab froze when switching to the IGV tab.
**Root cause (1):** CRAM files require a reference FASTA. Without explicit
`fastaURL`/`indexURL` in the browser config, IGV.js tried to fetch the default
hg38 reference and hung.
**Fix:** Added explicit hg38 reference URLs to the browser options object.
**Root cause (2):** Streamlit rerenders the entire component on every interaction,
reinitializing the IGV browser and triggering repeated heavy network fetches.
**Fix:** Gated the `st.components.v1.html` call behind a "Load IGV" button using
`st.session_state["_igv_loaded"]`, so the browser only initializes once per session.

---

## IGV Desktop Integration

### Sort-by-base on session load silently ignored
**Symptom:** IGV sessions generated by GEAC Explorer were supposed to
sort reads by base at the drill-down locus, but the sort never happened.
No error was raised.

**Root cause 1 — XML session format doesn't support sort:**
The initial implementation added `<RenderOptions sortOption="BASE"
sortByPosition="chr:pos"/>` to each `<Track>` element in the session XML.
IGV Desktop silently ignores this element — confirmed by an IGV maintainer
in [igvteam/igv#224](https://github.com/igvteam/igv/issues/224). There is
no way to specify sort-on-load in IGV's XML session format.

**Root cause 2 — sort sent as HTTP GET to a socket interface:**
After removing the XML approach, the sort was sent as an HTTP GET request:
`http://localhost:60151/sort?option=BASE&locus=chr:pos`. This silently
failed because port 60151 is a **plain-text socket command interface**, not
an HTTP API. Only `/load` and `/goto` have HTTP handlers; all other batch
commands (including `sort`) must be sent as raw text over a TCP socket.

**Fix:** Replaced the HTTP request with a raw socket command:
```python
with socket.create_connection(("localhost", 60151), timeout=5) as sock:
    sock.sendall(f"sort base {locus}\n".encode())
    sock.recv(256)  # read "OK" response
```
A 2-second delay after session load gives IGV time to fetch BAM index data
before the sort command arrives. The sort only works when IGV is reachable
via the socket (auto-launch mode); downloaded session ZIPs cannot include
sort instructions.

**Lesson:** IGV Desktop's port 60151 has two interfaces that look similar
but behave differently: a small HTTP handler for `/load` and `/goto`, and
a text-based socket command interface for everything else (sort, snapshot,
goto, preference, etc.). The HTTP interface returning no error for unknown
paths makes it easy to assume the command was received when it was actually
dropped.

---

## Docker / Deployment

### Podman: collaborator couldn't find database file
**Symptom:** Collaborator entered `cohort.duckdb` in the Explorer UI but got
"database does not exist."
**Root cause:** Podman does not mount the current directory by default. The
`cohort.duckdb` file existed on the host but was not visible inside the container.
**Fix:** Required explicit `-v $(pwd):/data` volume mount and using `/data/cohort.duckdb`
as the path inside the container.

### Docker image was private on ghcr.io
**Symptom:** Collaborator running `podman pull ghcr.io/fleharty/geac` received
"permission denied."
**Root cause:** GitHub Container Registry defaults to private for new packages.
**Fix:** Made the package public via GitHub → Packages → Package Settings →
"Change visibility."

---

## CI / GitHub Actions

### `macos-13` runner no longer available
**Symptom:** The `native-binaries` job for `macos-x86_64` failed with
`"The configuration 'macos-13-us-default' is not supported"`, causing the entire
workflow run to be cancelled.
**Root cause:** GitHub deprecated the `macos-13` runner pool.
**Fix:** Dropped `macos-x86_64` (and all Linux native builds) from CI entirely —
only `macos-arm64` is needed for the Homebrew tap.

### Homebrew install: "No developer tools installed" on headless Mac
**Symptom:** `brew install fleharty/geac/geac` fails immediately with "No developer tools installed" on a collaborator's machine with no GUI access.
**Root cause:** Homebrew requires Xcode Command Line Tools. `xcode-select --install` only works when a GUI is available (it spawns a dialog).
**Fix (headless):**
1. Create a placeholder file to surface CLT in softwareupdate: `touch /tmp/.com.apple.dt.CommandLineTools.installondemand.in-progress`
2. Find the package name: `softwareupdate --list`
3. Install: `softwareupdate --install "Command Line Tools for Xcode-16.4" --agree-to-license`

### Homebrew formula: SHA256 mismatch on source tarball
**Symptom:** `brew install` failed with "Resource reports different checksum".
**Root cause:** CI used `gh api repos/.../tarball/refs/tags/TAG` to download the source tarball for SHA256 computation, but Homebrew fetches from `https://github.com/.../archive/refs/tags/TAG.tar.gz`. GitHub serves slightly different tarballs from these two endpoints, so the SHA256s don't match.
**Fix:** Changed CI to use plain `curl` against the exact `archive/refs/tags` URL that Homebrew uses. Works now that the repo is public (no auth needed).

### Homebrew formula: `#{version}` empty inside `resource` block
**Symptom:** `brew install` attempted to fetch `v.tar.gz` (version missing from URL).
**Root cause:** In Homebrew's DSL, `#{version}` inside a `resource` block refers to the *resource's* own version, not the enclosing formula's version. Since the resource has no version set, it stringifies to empty string.
**Fix:** Changed the resource URL template to use the `FORMULA_VERSION` sed placeholder (`vFORMULA_VERSION.tar.gz`) so CI substitutes the correct version directly, rather than relying on Ruby interpolation at install time.

### Homebrew tap: `HOMEBREW_TAP_TOKEN` secret not set
**Symptom:** CI log shows `TAP_TOKEN: ` (blank); `git push` to `homebrew-geac` would fail silently or with auth error.
**Root cause:** The `HOMEBREW_TAP_TOKEN` secret must be manually created in GitHub repo settings — it is not auto-provisioned like `GITHUB_TOKEN`.
**Fix:** Create a classic PAT with `repo` scope for `fleharty/homebrew-geac`, then add it as a repository secret named `HOMEBREW_TAP_TOKEN` in the GEAC repo settings.

### Homebrew tap: `curl` 404 on source tarball
**Symptom:** `curl: (22) The requested URL returned error: 404` when downloading the GitHub source archive for SHA256 computation.
**Root cause:** Using plain `curl` without authentication on a repo that requires it returns 404 instead of 401.
**Fix:** Replaced `curl -fsSL` with `gh api repos/.../tarball/...`, which uses `GH_TOKEN` automatically.

### `gh release upload` fails with "release not found"
**Symptom:** The `Package and upload binary` step in `native-binaries` exited with
`release not found` when using `gh release upload`.
**Root cause:** `gh release upload` requires the GitHub Release object to already
exist. The `docker` job runs in parallel and does not create a release; nothing
creates the release before the binary job runs.
**Fix:** Reverted to `softprops/action-gh-release@v2`, which creates the release
automatically if it doesn't exist.

### Multi-platform Docker: arm64 runner unavailable → whole run cancelled
**Symptom:** v0.3.7 first release attempt was cancelled mid-run.
**Root cause:** The original Docker workflow used a matrix (amd64 + arm64 via
`ubuntu-22.04-arm`) with a digest-merge approach. Combined with the `macos-13`
failure above, the run was cancelled before the Homebrew tap job ran.
**Fix:** Simplified to a single `docker` job building only `linux/amd64`, removing
the matrix, digest upload/download, and `imagetools create` merge step entirely.

---

### `vec![AccumulatorType::default(); n]` requires `Clone`
**Problem:** During `geac coverage` per-interval accumulation, `IntervalAccumulator`
was initialized with `vec![IntervalAccumulator::default(); n]`. The Rust compiler
rejected this because the repeat-count form of `vec![]` requires the element type to
implement `Clone` (it clones the initial value to fill the remaining slots).
**Root cause:** `IntervalAccumulator` had `#[derive(Default)]` but not `#[derive(Clone)]`.
The compiler error was `the trait Clone is not implemented for IntervalAccumulator`.
**Fix:** Add `Clone` to the derive list: `#[derive(Default, Clone)]`.
**Lesson:** `vec![value; n]` always needs `Clone`. Either derive it or initialise with
`(0..n).map(|_| T::default()).collect()` to avoid the requirement.

---

### DuckDB rejects SQL reserved word `end` in an index definition
**Problem:** Adding a `CREATE INDEX` on `coverage_intervals` with the column list
`(sample_id, chrom, start, end)` caused a DuckDB parser error at runtime.
**Root cause:** `end` is a reserved keyword in SQL/DuckDB. It was accepted as a
column name in `CREATE TABLE` (DuckDB is lenient there) but rejected unquoted inside
an index definition.
**Fix:** Quote the column name: `\"end\"` in the Rust string literal, which renders as
`"end"` in the emitted SQL.
**Lesson:** Avoid reserved words as column names. If you must use them, always
double-quote them everywhere — `CREATE TABLE`, `SELECT`, `JOIN ON`, index definitions.
DuckDB's leniency in `CREATE TABLE` can mask the problem until an index or query hits it.

---

## v0.4.0 Agenda

### Stream `alt_bases`/`alt_reads` to Parquet instead of buffering in memory
**Problem:** `collect_alt_bases` (`src/bam/mod.rs`) returns fully materialized
`Vec<AltBase>` and `Vec<AltRead>` buffers, which are only written to Parquet later
in `src/main.rs` and `src/writer/parquet.rs`. For `--reads-output` mode this is the
dominant scalability risk: a single high-depth sample with many alt loci can produce
millions of `AltRead` rows, all held in RAM simultaneously.
**Proposed fix:** Thread a streaming Parquet writer (following the pattern established
by `CoverageWriter` in `src/writer/parquet_coverage.rs`) through the collect loop,
flushing in fixed-size batches. The `Vec`-return API would be replaced with a
callback or writer-sink pattern.
**Impact:** `alt_bases` output unaffected (Terra re-run not required); `alt_reads`
output format unchanged but memory footprint drops from O(total_alt_reads) to O(batch).

### Gene bar chart click-to-drill-down not working in Low Coverage tab
**Problem:** Clicking a bar in the gene-level bar chart in the Low Coverage Explorer
tab does not trigger the gene drill-down table.
**Root cause:** Same class of bug as the main Explorer's position drill-down:
- `st.altair_chart(..., on_select="rerun")` had no `key=` parameter — Streamlit
  auto-generated a key from render-tree position, which shifted across reruns causing
  the selection event to be dropped.
- `st.dataframe(..., on_select="rerun")` for the low-coverage table was also missing
  `key=`.
- The "deselect" branch in the session-state update unconditionally cleared
  `_low_selected_gene` whenever `event.selection` returned empty `{}` — which is
  exactly what a missing-key widget returns on every non-click rerun, so the drill-down
  was immediately cleared after being set.
**Fix:**
- Added `key="low_coverage_gene_bar"` to `st.altair_chart`.
- Added `key="low_coverage_table"` to the low-coverage `st.dataframe`.
- Tightened the deselect logic: only clear `_low_selected_gene` when
  `event.selection` is non-empty but contains no populated list (explicit
  deselection), not when `event.selection` is falsy (no event yet / key mismatch).
**Lesson:** Every `on_select="rerun"` widget needs an explicit stable `key=`.
The failure mode is always silent (empty selection, no error). See the main
`on_select="rerun"` recipe entry above.

---

### `select_slider` range mode requires explicit `value=` tuple — `key=` alone is insufficient
**Symptom:** The gnomAD AF range slider snapped back to its default position (`("0", "1.0")`)
every time the user released the mouse, even though the drag appeared to work visually.
**Root cause:** Streamlit's `select_slider` determines whether to render as a single-value
or range slider based on the `value=` parameter passed at call time. Using `key=` alone
(relying on session state to supply the value) does not communicate range mode to Streamlit.
On the first render, our normalization code set `st.session_state["gnomad_af_range"] =
("0", "1.0")` (a tuple), but the slider rendered in single-value mode because no `value=`
tuple was passed to the function. After interaction, Streamlit stored a bare string `"0"`
in session state. Our `isinstance(..., tuple)` guard then fired, reset the state to
`("0", "1.0")`, and the slider snapped back.
**Fix:** Removed `key=` from the `select_slider` call. Instead, read the current value
from session state manually, normalize it to a 2-tuple, pass it as `value=` explicitly
(forcing range mode), and then write the returned value back to session state after the
widget renders.
**Lesson:** For `st.select_slider` (and likely other range widgets), range mode is
determined at call time by `value=`, not by what is in session state. If you want a
persistent range slider, manage session state manually rather than relying on `key=`.

---

### `--fill-zeros` filled entire genome when `--region` was set
**Symptom:** Running `geac coverage --fill-zeros --region chr20` would emit zero-depth
records for every position on every chromosome in the BAM header, not just chr20.
**Root cause:** The zero-fill loop iterates over `bam_contigs` (all chromosomes from the
BAM header). The `covered` bitset only has positions marked for the restricted region
(chr20 in this case, since `--region` restricts the pileup fetch). Every position on
every other chromosome passes the `!cov.contains()` check and gets emitted as zero-depth.
**Fix:** Parse `--region` (handling `chrom`, `chrom:start-end`, and `chrom:start` forms)
at the start of the zero-fill phase and skip chromosomes that don't match, restricting
the fill range to the parsed start/end when coordinates are provided.
**Lesson:** Any post-pileup fill pass that iterates over all reference contigs must
respect `--region`. The pileup fetch is the only thing that naturally inherits the
region restriction; anything that follows the pileup loop does not.

---

### Exon shading in Depth Profile — multiple failed approaches

**Feature goal:** Show exon/feature boundaries as visible shaded bands on the Depth
Profile chart in the coverage explorer, to help users see which coverage dips fall
inside vs. outside of exons.

**Data confirmed correct:** Debugging confirmed `_ivl_bounds` contains the right data
(e.g. 6 exon rows with correct coordinates matching the `_prof_df` position range).
The problem is entirely in the visualization layer.

#### Attempt 1 — `mark_rect` layered behind depth chart
Built `mark_rect` bands from `_ivl_bounds` and added them as the first (bottom) layer:
`_prof_chart = _bands + _borders + _labels + _depth_chart`.
**Why it failed:** `mark_area` fills from y=0 up to the depth value across the full
width of each position. It is a solid opaque fill that completely covers any layer
beneath it, regardless of layer order in Altair.

#### Attempt 2 — `mark_rect` layered on top of depth chart
Reversed the order: `_prof_chart = _prof_chart + _bands + _borders + _labels`.
**Why it failed:** Same result. Altair renders the depth `mark_area` over the rect
bands even when the rects are declared last. Altair layer order does not reliably
control z-order in the browser SVG renderer for this combination.

#### Attempt 3 — `mark_rule` border lines instead of filled rects
Replaced `mark_rect` with `mark_rule` vertical lines at exon start/end positions.
**Why it failed:** Lines were present in theory but invisible at gene scale — the
gene's intron span is much wider than any single exon, so the domain auto-scales to
show the entire gene and exon widths compress to sub-pixel size.

#### Attempt 4 — Explicit x-axis domain clamped to exon extents
Added `_x_domain = [_exon_x_min - _x_padding, _exon_x_max + _x_padding]` and
passed it to all chart encodings.
**Why it failed:** The depth `mark_area` chart also needs the matching domain, but
applying it to only the exon layers (with the main depth chart sharing the same
x-axis via layering) caused domain conflicts. Exon numbers appeared as small floating
text but the shaded bands were still not visible.

#### Attempt 5 — `alt.vconcat` genome-browser strip (current approach)
Built a separate 40px exon track below the depth chart using `mark_rect` on a 0–1
y-scale, combined with `alt.vconcat(...).resolve_scale(x="shared")`.
**Status:** Not yet confirmed working. Likely issues to investigate:
- `resolve_scale(x="shared")` may not propagate the explicit `_x_domain` from the
  exon track up to the depth chart's x-axis, so the depth chart may still auto-scale
  to the full gene range while the exon track zooms to exon extents.
- The `_x_domain` is only applied to the exon track's x encoding; the depth chart
  (`_prof_chart`) has no domain restriction and will use the data range of `_prof_df`.
  If `_prof_df` covers the full gene (including introns), the shared x domain will
  expand to match it.
- Fix direction: apply `_x_domain` to the depth chart's x-axis as well, or filter
  `_prof_df` to the exon extent range before plotting.

**Lesson:** Altair's layer z-order does not reliably place `mark_rect` on top of
`mark_area`. When exon annotation must be visible alongside depth area charts, use
`vconcat` with a separate exon strip. Ensure the shared x-domain is applied
consistently to both sub-charts, not just the exon track.

### Depth profile distorted with adaptive depth thresholding
**Symptom:** When a cohort DuckDB was built from runs using `--adaptive-depth-threshold`,
the depth profile showed inflated depth values for some samples (e.g. 120x when expected
30–40x) and the mean line no longer sat within the IQR band.
**Root cause (two issues):**
1. *Wrong `end` coordinate for partial bins (Rust):* When adaptive thresholding flushes
   a bin early (e.g. after 20 of a configured 100 positions), `end` was set to
   `bin_start + bin_size` (100) instead of `bin_start + n` (20). This caused the interval
   to overstate its genomic span.
2. *Mixed-resolution aggregation (Python):* The coverage table contains both 100bp bins
   (covered regions) and 1bp rows (dropout positions). Computing percentiles over all raw
   records let the 1bp dropout rows dominate the IQR, while the weighted mean ignored them
   — the two statistics diverged. An earlier attempt used a two-level CTE to re-bucket
   everything to the largest `bin_n`, but this still aggregated by the left-edge `pos`
   value and didn't correctly handle records whose intervals overlapped a display bin
   boundary.
**Fix:**
1. `src/coverage/mod.rs`: Changed `end: self.bin_start + self.bin_size` to
   `end: self.bin_start + self.n as i64`, so partial bins carry the correct half-open
   interval.
2. `app/coverage_profile.py`: Replaced re-bucketing with interval expansion using DuckDB's
   `UNNEST(range(pos, "end"))`, expanding each coverage row to individual base positions
   before aggregation. Per-sample stats are then computed at base resolution and
   cross-sample statistics follow. Also fixed the WHERE clause for region queries from
   `pos >= start` to `"end" > start` so records spanning the boundary are included.
**Lesson:** When coverage records represent variable-width intervals, always expand to base
positions (or re-bucket using the `end` coordinate, not `pos + bin_size`) before computing
cross-sample statistics. Percentiles computed over mixed-width records are not meaningful.

### Depth profile looked like bars and rendered slowly after interval expansion fix
**Symptom:** After switching to `UNNEST(range(pos, "end"))` for correct interval expansion,
the depth profile plot appeared as a staircase of solid bars rather than smooth area/line
bands, and rendering became noticeably slow (especially on larger genes).
**Root cause:** Expanding each 100bp bin into 100 individual 1bp rows means the aggregated
profile has 100 consecutive positions with identical depth values (one per base within the
original bin). Altair's `mark_area` renders these as a stepped/blocky fill rather than a
smooth curve, and the chart has 100× more data points than the original binned view.
**Fix:** After interval expansion and cross-sample aggregation, re-aggregate to the maximum
bin width in the region (`max_bin_width()` queries `MAX("end" - pos)`). This collapses the
100 identical 1bp rows back to a single 100bp display point — preserving correctness for
mixed-resolution adaptive-threshold data while matching the original visual resolution and
performance. Added `display_step` parameter to `load_expanded_depth_profile` and
`load_expanded_sample_profile` in `app/coverage_profile.py`.
**Lesson:** Interval expansion to 1bp is the right aggregation strategy but not the right
display strategy. Always re-bucket to the native display resolution after aggregation.

### `geac merge` failing when merging two DuckDB files — table existence checks broken
**Symptom:** `geac merge a.duckdb b.duckdb --output cohort.duckdb` failed with:
```
Catalog Error: Table with name tables does not exist!
Did you mean "information_schema.tables"?
```
and then on the second run attempt:
```
failed to insert 'alt_bases' from attached '_src0'
Catalog Error: Table with name alt_bases does not exist!
Did you mean "_src0.alt_bases"?
```
**Root cause (two bugs in `src/merge.rs`):**
1. `src_table_exists` queried `{alias}.information_schema.tables` — not reliably accessible
   for attached databases in DuckDB 1.x.
2. `dst_table_exists` queried `information_schema.tables` (unqualified) — this includes
   tables from all attached databases. So when the source DuckDB was attached, the output
   DB appeared to already have `alt_bases`, causing an `INSERT INTO` on a non-existent
   destination table.
**Fix:** Both functions now use `duckdb_tables()` — `src_table_exists` filters by
`database_name = '{alias}'` and `dst_table_exists` filters by
`database_name = current_database()`.
**Lesson:** In DuckDB, `information_schema.tables` (unqualified) is not scoped to the
current database when other databases are attached. Always use `duckdb_tables()` with an
explicit `database_name` filter when checking table existence in a multi-database context.

### Explorer slowdown after adding pipeline comparison and N-in-reads features
**Symptom:** After the v0.4.5/v0.4.6 additions (pipeline column + pipeline comparison
tab expansion + read-context N-tracking), the Streamlit explorer became noticeably
sluggish. Every sidebar widget interaction — even toggling a filter on a tab that
didn't use any of the new features — took seconds to respond.
**Initial (wrong) diagnosis:** Uncached sidebar `SELECT DISTINCT` queries
(`chroms`, `samples`, `batch`, `pipeline`, `label1/2/3`, `variant_filter`). These
were indeed firing on every rerun and were worth caching in `st.session_state`, but
fixing them did not resolve the perceived slowness.
**Root cause:** Streamlit evaluates the body of **every** `st.tabs()` tab on every
rerun — tabs are pure containers, not lazy. Four heavy queries added by the new
features were running unconditionally inside tab bodies, so interacting with *any*
widget (even on an unrelated tab) forced them to re-execute:
1. `_nctx_df` (Read Context tab): `SELECT ... FROM alt_reads JOIN filtered_loci`
   for N-context fraction metrics. Full `alt_reads` scan.
2. `_nasym_table_df` (Read Context tab): IGV-backing query that builds a huge
   `OR`-chain WHERE clause over the ranked N-asymmetry table and re-scans
   `table_expr`.
3. `_pc_df` (Pipeline Comparison tab): `FULL OUTER JOIN` scanning `table_expr`
   twice (once per pipeline).
4. `_rt_df` (Read Type Comparison tab): another `FULL OUTER JOIN` scanning
   `table_expr` twice (once per read type).
The `alt_reads` table got wider with the new N-tracking columns
(`n_before_alt`, `n_after_alt`, `n_n_before_alt`, `n_n_after_alt`,
`leading_n_run_len`, `trailing_n_run_len`), making query (1) more expensive than
it would have been before that feature landed.
**Fix:** Wrap each query in a `st.session_state` cache guard keyed on the filter
strings that determine its output (`where`, `_r_reads_filter`, `_pc_wa`/`_pc_wb`,
`_rt_wa`/`_rt_wb`). The queries only re-run when those strings change, not on
every widget interaction.
**Lesson:** `st.tabs()` is not lazy — every tab body runs on every rerun. Any
expensive query inside a tab body must be cached on its filter inputs, or it will
bleed performance into every other tab in the app. When investigating "the
explorer got slow after feature X," look for uncached queries in tab bodies
**anywhere** in the app, not just in the tabs that feature X touched.

### Strand bias drill-down crashed with `KeyError: ['pos_display'] not in index`
**Symptom:** Clicking a point on the strand bias scatter triggered a Streamlit rerun that
ended with:
```
KeyError: "['pos_display'] not in index"
```
at the `st.dataframe(_sb_sel_df[_table_cols], ...)` call in the strand bias selection
drill-down block.
**Root cause:** The strand bias selection query was written inline rather than routed
through `query_records()`. When `pos_display` (`pos + 1 AS pos_display`) was introduced
as a 1-based display column and added to `_table_cols`, `query_records()` was updated to
include it, but the handful of inline queries scattered around the file — including the
strand bias one — were not. The mismatch went unnoticed because clicking a strand bias
point is an optional interaction; the rest of the app worked fine.
**Fix:** Add `pos + 1 AS pos_display` to the strand bias inline `SELECT` (line ~3774 in
`geac_explorer.py`), matching what `query_records()` emits.
**Lesson:** Every handwritten `SELECT * FROM {table_expr}` that feeds a `_table_cols`
display is a latent regression point. When adding a column to `_table_cols`, grep for
inline queries that bypass `query_records()` and update them at the same time. Consider
extracting a shared SQL fragment or helper for the computed display columns (`vaf`,
`pos_display`) so there is only one place to update.

### N-context density chart exceeded Streamlit 200 MB message size limit
**Symptom:** Enabling "N-context read distribution" in the Reads tab on a real cohort raised:
```
MessageSizeError: Data of size 2813.1 MB exceeds the message size limit of 200.0 MB.
```
**Root cause:** The original query fetched every row from `alt_reads` — one row per
alt-supporting read — into a Python dataframe, then `_nctx_long` doubled it (one row per
read per side: before/after). The entire dataframe was embedded in the Altair chart spec
as JSON and sent to the browser for client-side `transform_density` computation. At cohort
scale this was millions of rows × multiple columns × JSON overhead.
**Fix:** Replace the full `alt_reads` fetch with two pre-aggregated DuckDB queries:
1. A single-row scalar aggregate for the four metric values (mean frac N before/after,
   fraction with any N before/after).
2. A histogram-bin query (bin width 0.001, range 0–0.05, grouped by `side` and
   `read_group`) that produces at most ~200 rows. The density curve is then drawn from
   these bins in Python instead of using Altair's `transform_density`.
**Lesson:** Never pass a raw `alt_reads` query result to an Altair chart. `alt_reads` is
read-level and grows with cohort size. Any chart driven by it must pre-aggregate in
DuckDB first. `transform_density` and `transform_bin` in Vega-Lite are convenient but
they require the full row-level data to be sent to the browser — always pre-compute these
transformations server-side for tables that scale with cohort size.

### gnomAD AF lookup silently missed every deletion (second bug)
**Symptom:** Even after an earlier fix for deletions, `gnomad_af` was still null for
every deletion in the collect output on v0.4.11.
**Root cause:** `IndelCount.ref_allele` stored **different things for different variant
types**: for insertions it was the single anchor base ("A"), for deletions it was the
deleted bases themselves ("TCG"). The gnomAD call site passed `indel.ref_allele`, and
`GnomadIndex::get` reconstructed the VCF REF as `ref_allele + alt_allele[1..]`. For a
deletion this produced `"TCG" + "TCG" = "TCGTCG"` instead of `"ATCG"`, so no record
ever matched. Separately, insertions also silently returned the wrong allele's AF at
multi-allelic sites because the indel path used `extract_af(record, field, 0)` — the
first alt — rather than finding the specific insertion.
**Fix:** (1) In `src/bam/mod.rs`, pass the reference anchor base (`ref_base.to_string()`)
to `gnomad.get()` for indels, matching the SNV call site and the `AltBase.ref_allele`
contract. (2) Remove the inconsistent `ref_allele` field from `IndelCount` entirely —
it had no other readers. (3) In `src/gnomad.rs`, factor the GEAC→VCF allele translation
into `geac_to_vcf_alleles()` (unit-tested) and use exact alt matching for indels as well
as SNVs, so multi-allelic insertion sites resolve to the right AF.
**Lesson:** A struct field named `ref_allele` with type-dependent semantics is a
landmine. If a field's meaning changes based on another field, extract it into the
discriminated type itself or delete it. Also: the GEAC↔VCF indel translation is the
kind of pure function that belongs in a testable helper — the original inline logic in
`get()` had been wrong for months with no unit test to catch it.

### `input_checksum_sha256` always NULL in merged DuckDB (`geac_inputs` table)
**Symptom:** After running `geac collect --input-checksum-sha256` and then `geac merge`,
querying `SELECT checksum_sha256 FROM geac_inputs` in the merged DuckDB returned NULL for
every row, even though the hash was visibly logged during the collect step.
**Root cause:** The hash is computed and written to the `input_checksum_sha256` column in
the per-sample alt_bases Parquet. But `geac merge`'s `InputProvenance` struct had no
`checksum_sha256` field, `collect_parquet_input_provenance()` never read the column from
the Parquet, and `write_input_provenance()` hardcoded `NULL` in the INSERT SQL. The hash
existed in the Parquet files but was simply never forwarded to the merged database.
**Fix:**
1. Added `checksum_sha256: Option<String>` to the `InputProvenance` struct in `src/merge.rs`.
2. In `collect_parquet_input_provenance()`, after the row/sample-count queries, added a
   query against `alt_bases` Parquets only (other table types don't have this column):
   ```rust
   let checksum_sha256: Option<String> = if spec.table == "alt_bases" {
       conn.query_row(
           &format!("SELECT input_checksum_sha256 FROM read_parquet('{escaped}')
                     WHERE input_checksum_sha256 IS NOT NULL LIMIT 1"),
           [], |row| row.get(0),
       ).optional()?.flatten()
   } else { None };
   ```
3. Added `checksum_sha256: None` to the DuckDB-source `InputProvenance` push (DuckDB
   inputs don't carry the per-BAM hash).
4. In `write_input_provenance()`, replaced the hardcoded `NULL` literal with a proper
   SQL literal derived from `input.checksum_sha256`.
**Lesson:** When adding a new metadata column to `geac collect` output, immediately trace
the full path through `geac merge`: `InputProvenance` struct → `collect_parquet_input_provenance`
→ `write_input_provenance` INSERT. A column that exists in the Parquet but is never read
by the merge step will silently produce NULLs with no error.

---

### `geac fragments` produced empty Parquet (0 records) when no `--region` specified

**Symptom:** Running `geac fragments` on a 1 GB BAM with no `--region` flag produced a
~2 KB Parquet with zero fragment records and completed almost instantly.

**Root cause:** `open_bam` returns a `rust_htslib::bam::IndexedReader`. On an
`IndexedReader`, calling `records()` without a prior `fetch()` silently yields zero
records — htslib requires an explicit fetch to initialise the iterator. When no region
is provided, the code skipped the `bam.fetch()` call, so the inner loop never executed.

**Fix:** Always call `bam.fetch(".")` (all mapped reads) when no region is specified.

**Lesson:** Any code that opens a BAM via `open_bam` (which returns `IndexedReader`)
must call `bam.fetch(...)` before iterating with `records()`, even for whole-file scans.

---

### `end_motif_3p` reported in + strand orientation instead of 5′→3′ of the − strand

**Symptom:** The `end_motif_3p` column in the fragments Parquet contained the raw + strand
sequence at the right cut site rather than the end motif as it would be read in the 5′→3′
direction of the minus strand, making 3′ motifs incomparable to 5′ motifs under the
standard cfDNA end-motif convention (Jiang et al. 2020).

**Root cause:** The code extracted `seq[frag_end-2 : frag_end+2]` directly from the + strand
reference. Because `geac fragments` only processes R1 with positive TLEN (R1 on + strand),
`frag_end` is the right edge of the fragment. The 3′ cut is the 5′ end of the minus strand;
to report its 4-mer in 5′→3′ direction it must be reverse-complemented.

**Fix:** Applied `reverse_complement()` to the extracted sequence for `end_motif_3p`:
```rust
let end_motif_3p = if frag_end as usize >= 2 {
    extract_motif(seq, frag_end as usize - 2, 4).map(|m| reverse_complement(&m))
} else {
    None
};
```
Added a `reverse_complement` helper in `src/fragments/mod.rs`.

**Lesson:** End motifs at the 3′ end of a fragment are on the minus strand. Always
reverse-complement the + strand reference sequence to get the motif in the biologically
meaningful 5′→3′ orientation. The 5′ motif needs no transformation because R1 already
reads in the + strand direction.

---

### `geac fragments` dropped half of all fragments by filtering on `is_first_in_template`

**Symptom:** Fragment counts from `geac fragments` were ~50% lower than expected for
standard FR paired-end libraries — every other proper pair was missing from the output
Parquet.

**Root cause:** The fragment-collection loop required both `record.is_first_in_template()`
(R1) AND `record.insert_size() > 0` to select one read per pair. In FR libraries R1 is
randomly assigned to either strand, so half of all proper pairs have R1 on the −
strand (the rightmost read with negative TLEN). Those pairs satisfied neither side of
the conjunction in a way that mattered: the rightmost R1 had negative TLEN (rejected by
the TLEN check), and the leftmost R2 was not first-in-template (rejected by the R1
check). Net effect: only ~half of pairs (those with R1 on +) were emitted.

**Fix:** Removed the `is_first_in_template()` filter entirely. By the SAM spec, exactly
one read per pair has positive TLEN — the leftmost — and in FR libraries it is on the
+ strand regardless of R1/R2 status. The TLEN-positive check alone is the correct
selector. (`src/fragments/mod.rs`)

**Regression test:** `fragments_emits_record_for_both_r1_strand_orientations` constructs
a BAM with one R1-on-+ pair (flags 99/147) and one R1-on-− pair (flags 83/163) and
asserts both produce a fragment record.

**Lesson:** Don't conflate "R1" with "leftmost / + strand" in paired-end data. R1/R2
designates the order in the FASTQ pair, not the strand. Use TLEN sign or alignment
strand flags (`0x10`, `0x20`) to reason about strand; use `0x40`/`0x80` only when you
need to distinguish the two reads of a pair without regard to position.
