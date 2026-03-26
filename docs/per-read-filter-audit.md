# Per-Read Filter Audit

A detailed review of the per-read filtering system: data collection (Rust),
storage (Parquet/DuckDB), and query/visualization (Explorer). Covers correctness,
efficiency, and pitfalls.

---

## Architecture Summary

The per-read filter system spans three layers:

1. **Rust (`src/bam/mod.rs`)** — during `geac collect --reads-output`, one `AltRead`
   record is emitted per alt-supporting read per locus. Written to
   `{stem}.reads.parquet` via `src/writer/parquet_reads.rs`.

2. **DuckDB (`src/merge.rs`)** — `geac merge` loads all `.reads.parquet` files into
   an `alt_reads` table with an index on `(sample_id, chrom, pos, alt_allele)`.

3. **Explorer (`app/geac_explorer.py`)** — sidebar widgets control four filter
   dimensions (family size, cycle number, mapping quality, insert size) with
   include/exclude toggles. Two query modes: **locus-inclusion** (default) and
   **re-aggregation** (`recompute_vaf` checkbox).

### Schema: `alt_reads` table

| Column               | Type  | Nullable | Source                                |
|----------------------|-------|----------|---------------------------------------|
| sample_id            | Utf8  | no       | CLI `--sample-id`                     |
| chrom                | Utf8  | no       | BAM pileup                            |
| pos                  | Int64 | no       | BAM pileup (0-based)                  |
| alt_allele           | Utf8  | no       | observed alt base/indel string        |
| dist_from_read_start | Int32 | no       | `alignment.qpos()`                    |
| dist_from_read_end   | Int32 | no       | `read_len - qpos - 1`                 |
| read_length          | Int32 | no       | `record.seq_len()`                    |
| ab_count             | Int32 | **yes**  | fgbio `aD` tag                        |
| ba_count             | Int32 | **yes**  | fgbio `bD` tag                        |
| family_size          | Int32 | **yes**  | fgbio `cD` tag                        |
| base_qual            | Int32 | no       | `record.qual()[qpos]`                 |
| map_qual             | Int32 | no       | `record.mapq()`                       |
| insert_size          | Int32 | **yes**  | `|TLEN|`; NULL when 0 (unpaired)      |

---

## Confirmed Bugs

### 1. Re-aggregation mode fails to zero-out SNV loci where all reads fail the filter

**Location:** `geac_explorer.py:390–408`

In re-aggregation mode (`recompute_vaf=True`), the query is:
```sql
SELECT
    ab.* EXCLUDE (alt_count),
    COALESCE(ar_agg.filtered_alt_count, ab.alt_count) AS alt_count,
    ...
FROM alt_bases ab
LEFT JOIN (
    SELECT ..., COUNT(*) AS filtered_alt_count
    FROM alt_reads
    WHERE {_reads_where}
    GROUP BY sample_id, chrom, pos, alt_allele
) ar_agg ON ...
```

When a SNV locus has rows in `alt_reads` but **none** pass `_reads_where`, the
subquery returns no row for that locus. The LEFT JOIN produces NULL for
`filtered_alt_count`, and `COALESCE(NULL, ab.alt_count)` falls back to the
**original unfiltered alt_count**.

The intent of this COALESCE was to let **indels** (which have no alt_reads rows
at all) keep their original count. But the current logic conflates "has no reads
in alt_reads" with "has reads but none pass," giving both the original count.

**Impact:** SNV loci that should show alt_count=0 (or be hidden) after filtering
instead appear with their original count — the filter silently has no effect on
them.

**Fix sketch:**
```sql
LEFT JOIN (
    SELECT sample_id, chrom, pos, alt_allele,
           COUNT(*) FILTER (WHERE {_reads_where}) AS filtered_alt_count,
           TRUE AS has_reads
    FROM alt_reads
    GROUP BY sample_id, chrom, pos, alt_allele
) ar_agg ON ...

-- Then:
CASE WHEN ar_agg.has_reads IS NULL THEN ab.alt_count  -- no reads at all (indels)
     ELSE COALESCE(ar_agg.filtered_alt_count, 0)      -- has reads, use filtered
END AS alt_count
```

---

### 2. Warning banner text always claims re-aggregation regardless of mode

**Location:** `geac_explorer.py:820–825`

The warning says:
> "alt_count and VAF are re-aggregated from reads passing the filter.
>  original_vaf shows the unfiltered VAF for comparison."

This is only true when `recompute_vaf=True`. In the default locus-inclusion mode,
alt_count is **not** re-aggregated and the `original_vaf` column does not exist.
The warning text should be conditional on `recompute_vaf`.

---

### 3. Insert size filter missing from the warning banner

**Location:** `geac_explorer.py:809–825`

The warning constructs `_active_parts` for family size, cycle number, and
mapping quality but **not** insert size. When only the insert size filter is
active, the warning reads:

> **Per-read filters active** ().

The empty parentheses are confusing; the insert size filter status is invisible.

---

### 4. Family-size stratified spectrum ignores per-read filters

**Location:** `geac_explorer.py:3197–3220`

The family-size stratified error spectrum CTE reads from `alt_reads` directly:
```sql
WITH locus_fs AS (
    SELECT ..., MEDIAN(family_size) AS median_fs
    FROM alt_reads
    WHERE family_size IS NOT NULL
    GROUP BY sample_id, chrom, pos, alt_allele
)
```

This does **not** apply `_reads_where`. If the user has a family size filter
active (e.g. family_size >= 3), the stratification still classifies loci using
ALL reads' family sizes — including the ones the filter excluded. The locus
set comes from `table_expr` (respects filters) but the singleton/multi
classification comes from unfiltered data.

---

## Semantic / Accuracy Concerns

### 5. Indel `dist_from_read_start` / `dist_from_read_end` / `base_qual` have different semantics than for SNVs

For **SNVs**, these fields refer to the exact position of the mismatched base:
- `dist_from_read_start = qpos` — where the alt base sits in the read
- `base_qual` — quality score of the alt base itself

For **indels**, the anchor base (the reference base immediately before the
indel) is used:
- `dist_from_read_start` = anchor position, not the midpoint of the indel
- `base_qual` = quality of the anchor base, which is a reference base — not
  directly informative about the indel call quality

For **deletions**, `alignment.qpos()` at the pileup anchor position should
return a valid value (the read does have a base here). The `unwrap_or(0)` at
line 702 is a safety fallback that should rarely fire. But if it does, all
affected deletion reads get `dist_from_read_start=0`, concentrating them at
cycle position 0 in plots.

**Recommendations:**
- Document that cycle-number and base-quality metrics are anchor-based for
  indels. Consider adding a note in the Reads tab when indels are included.
- For insertions, `qpos` could reasonably be shifted by `ceil(len/2)` to
  represent the midpoint of the insertion, but this would diverge from the
  anchor convention. Keep the anchor convention but document it.

---

### 6. "Cycle number" label / filter vs visualization column mismatch

The sidebar filter is labeled **"Cycle number"** and filters on
`dist_from_read_end`. The Reads tab visualization is labeled **"Cycle (distance
from read start)"** and plots `dist_from_read_start`.

These are complementary (`start + end + 1 = read_len`), not the same column.
The filter removes reads where the alt base is close to the **end** of the read
(a classic read-end artifact filter). The plot shows position from the
**start** (classic cycle-number plot). This is logically consistent but the
shared "Cycle number" label is misleading — the filter isn't directly controlling
what the x-axis shows.

**Recommendation:** Rename the sidebar filter to **"Min distance from read end"**
or similar. Reserve "cycle number" for the `dist_from_read_start` axis.

---

## Efficiency Observations

### 7. Full `alt_reads` table scan on every filter change

Both query modes (BOOL_OR and COUNT) scan the entire `alt_reads` table, grouping
by locus keys. The composite index on `(sample_id, chrom, pos, alt_allele)`
helps the GROUP BY but DuckDB must still evaluate `_reads_where` against every
row. For a large cohort (millions of reads), this can be slow.

**Mitigation:** DuckDB's columnar engine and Parquet predicate pushdown handle
this reasonably well. The current approach is the correct trade-off — creating
per-filter indexes or materialized views would add complexity for marginal gain.
If performance becomes a problem, consider caching the grouped result in a
DuckDB temp table on the first filter interaction and invalidating on filter
change.

---

### 8. Slider bound queries run on every rerun

`_reads_maxes` (line 211) computes `MAX(family_size)`, `MAX(dist_from_read_end)`,
`MAX(map_qual)`, and `COUNT(insert_size) > 0` from a full table scan of
`alt_reads` on **every** Streamlit rerun. These values never change between
reruns (the database is read-only). The results are stored in `session_state`
but the query still runs.

**Recommendation:** Gate the query behind a session_state check:
```python
if "_cached_fs_max" not in st.session_state:
    _reads_maxes = con.execute(...).fetchone()
    ...
```

---

### 9. `_r_join` re-executes the locus subquery for every plot

The shared join expression `_r_join` (line 975) includes:
```sql
INNER JOIN (
    SELECT DISTINCT sample_id, chrom, pos, alt_allele
    FROM {table_expr}
    WHERE {where}
) _filt ON ...
```

This subquery — which itself contains the LEFT JOIN from `table_expr` — is
textually inlined into every reads-tab query. DuckDB may optimize duplicate
subexpressions, but there's no guarantee. If the reads tab has multiple plots,
this work is repeated.

**Recommendation:** Consider materializing the filtered locus set into a DuckDB
temp table at the top of the render, then referencing it by name.

---

## Pitfalls for Users

### 10. VAF is a lower bound in re-aggregation mode

When `recompute_vaf` is checked, the displayed VAF = `filtered_alt_count /
total_depth`. But `total_depth` still reflects ALL reads at the locus (including
those filtered out). The true filtered VAF would require also re-computing
total_depth from only reads with the matching family size / cycle / MAPQ — which
requires access to ALL reads (not just alt-supporting ones), and is not available
from the current schema.

The Explorer correctly notes this in the help text and sidebar caption, but users
may still be surprised.

---

### 11. Insert size filter silently excludes unpaired reads

`insert_size BETWEEN x AND y` excludes rows where `insert_size IS NULL`
(unpaired reads, mate-unmapped reads). Unlike the family size filter — which has
explicit `IS NULL OR` handling in exclude mode — the insert size filter has no
NULL guard. Activating the insert size filter implicitly removes all unpaired
reads. There is no exclude-mode toggle for insert size.

---

### 12. Locus-level filters interact non-obviously with per-read filters in re-aggregation mode

In re-aggregation mode, the `alt_count` in `table_expr` is the filtered count.
The locus-level VAF range filter then applies to `filtered_alt_count /
total_depth`. A locus with original VAF = 0.05 might have filtered VAF = 0.02
after removing low-family-size reads, causing it to be excluded by a
`vaf_range >= 0.01` filter that would have included it in locus-inclusion mode.

This is correct behavior — the filters compose — but may confuse users who
expect the VAF filter to operate on the original VAF. The `original_vaf` column
is available in the data table for comparison, but the VAF distribution plots
use the filtered VAF.

---

### 13. Per-read filters have no effect on loci without `alt_reads` rows

For databases created before v0.3.8 (or without `--reads-output`), the
`alt_reads` table either doesn't exist (filters are hidden) or contains only
SNV rows (indels pass through unfiltered via the `ar.sample_id IS NULL` path).

Even with v0.3.8+, if the user runs `geac collect` without `--reads-output`,
there are no reads to filter. The sidebar widgets still appear if the table
exists from a previous merge, which could be misleading if the data is stale.

---

## Summary Table

| # | Type       | Severity | Description                                                    |
|---|------------|----------|----------------------------------------------------------------|
| 1 | Bug        | Medium   | Re-aggregation COALESCE gives original count instead of 0      |
| 2 | Bug        | Low      | Warning text wrong for locus-inclusion mode                    |
| 3 | Bug        | Low      | Insert size filter missing from warning banner                 |
| 4 | Bug        | Medium   | Family-size stratified spectrum bypasses per-read filters      |
| 5 | Semantic   | Low      | Indel dist/base_qual use anchor position, not indel position   |
| 6 | Semantic   | Low      | "Cycle number" label mismatch between filter and visualization |
| 7 | Efficiency | Low      | Full alt_reads scan per filter change (acceptable)             |
| 8 | Efficiency | Low      | Slider bound MAX queries re-run on every rerun                 |
| 9 | Efficiency | Low      | `_r_join` re-executes locus subquery per plot                  |
| 10 | Pitfall   | Medium   | VAF is lower bound in re-aggregation mode                     |
| 11 | Pitfall   | Medium   | Insert size filter silently drops unpaired reads              |
| 12 | Pitfall   | Low      | Locus-level VAF filter applies to filtered alt_count          |
| 13 | Pitfall   | Low      | Stale or absent alt_reads → filters have no/partial effect    |
