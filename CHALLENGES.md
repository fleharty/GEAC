# Development Challenges

A running log of non-obvious problems encountered during GEAC development, how they
were diagnosed, and what the fix was. Updated as new issues arise.

---

## Rust / Core Tool

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
