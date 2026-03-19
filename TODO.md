# GEAC — TODO

## Rust / CLI

- [x] `geac qc` subcommand — per-sample summary of error rates by substitution type, strand bias metrics, and overlap concordance
- [x] `geac cohort` subcommand — per-locus artifact frequencies across samples (flag positions seen in N% of samples); outputs TSV or Parquet
- [x] On-target annotation — `--targets` BED/Picard interval list flag; records `on_target bool?` column in Parquet
- [x] Gene annotation — `--gene-annotations` flag; accepts GFF3, GTF, or UCSC genePred (.txt/.txt.gz); records `gene string?` column
- [x] Locus repetitiveness metrics — `homopolymer_len`, `str_period`, `str_len` columns; `--repeat-window` flag (default 10 bp)
- [x] Trinucleotide context — `trinuc_context` column computed from reference at each SNV locus
- [x] Variant annotation — `--vcf` and `--variants-tsv` flags; annotates `variant_called` / `variant_filter` columns
- [x] Fragment overlap metrics — `overlap_alt_agree`, `overlap_alt_disagree`, `overlap_ref_agree` columns for read-pair concordance
- [x] Integration tests — generate synthetic BAM data and write end-to-end tests
## Per-read detail table (two-table design)

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

### Implementation steps
- [ ] Step 1: Define `AltRead` struct in `src/record.rs` with the columns listed above
- [ ] Step 2: Populate `AltRead` records during BAM pileup in `src/bam/mod.rs`; parse fgbio/DRAGEN family size tags (`cD`/`cE` or `RX`/`MI`) per pipeline
- [ ] Step 3: Add `src/writer/parquet_reads.rs` — write `AltRead` records to `{stem}.reads.parquet`
- [ ] Step 4: Update `geac collect` CLI to emit both files; add `--reads-output` flag (optional; if omitted, reads table is not written)
- [ ] Step 5: Update `geac merge` to accept and merge reads Parquets into a second DuckDB table (`alt_reads`) alongside `alt_bases`
- [ ] Step 6: Update WDL workflows to handle the optional reads Parquet output from `Collect` and pass it to `Merge`
- [ ] Step 7: Explorer — in the position drill-down, JOIN `alt_reads` on `(sample_id, chrom, pos, alt_allele)` to show per-read detail (dist from end, family size, base qual) when reads table is present
- [ ] Step 8: Explorer — add sidebar filters for `min_family_size` and `max_dist_from_end` that apply to the locus table via a subquery against `alt_reads`

## Intra-sample comparison (read-type)

- [ ] Multi-BAM collect — allow `geac collect` to accept multiple input BAMs for the same sample (raw, simplex, duplex) and tag each record with its `read_type`; output a single merged Parquet with a `read_type` column
- [ ] Read-type comparison view in Explorer — given a Parquet or DuckDB with mixed read types, show side-by-side or overlaid metrics (VAF distribution, strand balance, SBS96 spectrum) broken down by `read_type` (raw / simplex / duplex / mixed); goal is to quantify what duplex consensus processing removes vs retains

## Explorer (Streamlit)

- [x] IGV session download — manifest-driven BAM tracks + BED positions zip, capped at 5 samples
- [ ] IGV sample picker — verify that the BED file only contains positions with alt bases in the selected samples; suspected bug where non-selected sample loci still appear in the session
- [x] BED file has some incorrect entries — fixed: deletion loci now span the full deleted region [pos, pos+del_len)
- [x] Position-level drill-down — click a locus and see all samples/alleles at that position
- [x] Export filtered data to CSV — handled by Streamlit's built-in dataframe toolbar download button
- [x] On-target filter — sidebar selectbox "Target bases": All / On target / Off target
- [x] Gene filter — sidebar text input with partial match (ILIKE); depends on `gene` column being populated
- [x] Repeat filter — sidebar range sliders for `homopolymer_len` and `str_len`
- [x] Strand bias plot — dashed y=x diagonal + 95% binomial CI band; gene name in hover tooltip
- [x] Strand bias click drill-down — click/shift-click to select points; shows table of selected loci and IGV session with correct BAMs and BED
- [ ] Strand bias selection: verify `toggle="event.shiftKey"` is valid in the current Altair/Vega-Lite version; may need an alternative approach for shift-click multi-select
- [ ] Strand bias selection: `pos` may be returned as float from Altair selection, causing SQL clause to silently fail; cast to int before building the WHERE clause
- [x] SNV trinucleotide spectrum (SBS96) — 3×2 grid of per-mutation-type panels with shared y-axis; click drill-down
- [ ] Cohort comparison view — side-by-side stats across samples loaded from a DuckDB
  - [x] Step 1: Per-sample summary table — one row per sample_id with n_snv, n_insertion, n_deletion, mean_depth, mean_vaf, strand_balance, overlap_concordance; clicking a row filters all other tabs to that sample
  - [x] Step 2: VAF distribution overlay — all samples on one plot as density curves, colored by sample; highlights shifted VAF distributions
  - [x] Step 3: Strand balance scatter — one dot per sample (x = mean strand balance, y = mean VAF); outliers immediately visible
  - [x] Step 4: SNV count bar chart — n_snv per sample, stacked/colored by SBS6 substitution type breakdown
  - [x] Step 5: SBS96 heatmap — samples as rows, 96 trinucleotide contexts as columns, color = normalized count; reveals samples with unusual mutational profiles
- [x] NMF decomposition — fit the per-sample SBS96 spectrum against COSMIC reference signatures using NNLS; show the largest contributing signatures and their weights

## Coverage Analysis

**Motivation:** `geac collect` only records positions where an alt base was observed.
`geac coverage` fills the denominator — depth at every covered position — enabling true
per-base error rates (`alt_count / total_depth` across all positions, not just alt positions)
and identification of systematically undercovered or low-mappability sites across a cohort.

Coverage is confounded by mappability: a region may appear undercovered simply because
reads that originate there cannot be placed uniquely by the aligner. Without mappability
context, low-coverage sites and low-mappability sites are indistinguishable from each other
or from genuine dropout (e.g. FFPE degradation, GC bias). The schema below captures both.

### CLI

```
geac coverage \
  --input        sample.bam \
  --reference    ref.fa \
  --output       sample.coverage.parquet \
  [--targets     targets.bed]            # BED or Picard interval list; restricts to these positions
  [--region      chr1:1-50000]           # genomic region (alternative to --targets)
  [--mappability mappability.bedgraph]   # pre-computed mappability track (BEDGraph, score col 4)
  [--sample-id   override]
  [--read-type   raw|simplex|duplex]
  [--pipeline    fgbio|dragen|raw]
  [--min-map-qual 20]                    # same as collect; also used to compute frac_low_mapq
  [--min-depth   0]                      # suppress positions with depth below this value
  [--bin-size    1]                      # aggregate N bp into one row (1 = per-position)
  [--threads     1]
```

`--targets` is strongly recommended for targeted panels — it bounds output size and ensures
zero-depth positions (complete dropout) are still recorded. Without `--targets` or `--region`
the whole BAM is scanned; useful for targeted panels but impractical for WGS without `--bin-size`.

### Output schema (`CoverageRecord`)

One row per position (or per bin when `--bin-size > 1`). Positions are 0-based.

```
sample_id:       String

chrom:           String
pos:             i64        # 0-based start of the position or bin
end:             i64        # pos + 1 normally; pos + bin_size with --bin-size > 1

# Depth — reads passing --min-map-qual
total_depth:     i32
fwd_depth:       i32
rev_depth:       i32
overlap_depth:   i32        # number of fragment pairs where both reads cover this position

# BAM-derived mappability signals (computed over ALL reads before any MAPQ filter)
# These are always populated; they cost nothing since we are already doing the pileup.
mean_mapq:       f32        # mean MAPQ of all reads at this position
frac_mapq0:      f32        # fraction of reads with MAPQ = 0 (definitive multi-mappers)
frac_low_mapq:   f32        # fraction with MAPQ < --min-map-qual threshold

# Pre-computed mappability (optional, populated when --mappability is provided)
# Score is 0.0–1.0 where 1.0 = perfectly unique, 0.0 = never uniquely mappable at this k-mer length.
# Common sources: ENCODE GEM mappability tracks (150-mer), genmap, Umap.
mappability:     Option<f32>

# Optional annotations
on_target:       Option<bool>   # populated when --targets is given

# Provenance
read_type:       ReadType
pipeline:        Pipeline
```

**Why two mappability signals?**

`frac_mapq0` is the empirical signal — it reflects what actually happened in this experiment
at this read length and mapper version. `mappability` is the theoretical signal — it reflects
the intrinsic repetitiveness of the reference sequence. Together they support three diagnostic
cases:

| `frac_mapq0` | `mappability` | Interpretation |
|---|---|---|
| High | Low | Classic multi-mapping locus — expected, filter with confidence |
| High | High | Unexpected low MAPQ — structural variant, misassembly, or aligner artifact |
| Low  | Low | Mappability track may be outdated or wrong k-mer length |
| Low  | Low, `total_depth` also low | Genuine biological dropout (GC bias, FFPE, etc.) |

### Pre-computed mappability track (`--mappability`)

Accepted format: **BEDGraph** (chrom, start, end, score). This covers the ENCODE GEM tracks
and output from genmap and Umap. bigWig support can be added later if needed (requires an
additional crate dependency).

Implementation: load into a sorted `Vec<(i64, i64, f32)>` per chromosome, then binary-search
for each pileup position. For targeted panels (tens of thousands of positions) this is fast
and fits comfortably in memory. For WGS the track itself is ~2 GB uncompressed; a streaming
approach (advance through the sorted track in lock-step with the sorted pileup) is more
appropriate and should be used if `--targets` is not provided.

### DuckDB integration

`geac coverage` outputs `{sample}.coverage.parquet`. When this file is passed to `geac merge`,
it is inserted into a `coverage` table in the cohort DuckDB (alongside `alt_bases`). `geac merge`
detects coverage Parquets by schema (presence of `mean_mapq` column; absence of `alt_allele`).

This enables cross-table queries in the Explorer — for example, joining alt bases to their
coverage context:

```sql
-- Alt bases at systematically low-coverage, low-mappability loci
SELECT a.*, c.mean_depth, c.mean_mapq
FROM alt_bases a
JOIN (
    SELECT chrom, pos,
           AVG(total_depth) AS mean_depth,
           AVG(mean_mapq)   AS mean_mapq
    FROM coverage
    GROUP BY chrom, pos
) c ON a.chrom = c.chrom AND a.pos = c.pos
WHERE c.mean_depth < 20 AND c.mean_mapq < 30
ORDER BY c.mean_depth;
```

```sql
-- Systematically undercovered positions across the cohort
SELECT chrom, pos,
       COUNT(DISTINCT sample_id)          AS n_covered,
       AVG(total_depth)                   AS mean_depth,
       AVG(frac_mapq0)                    AS mean_frac_mapq0,
       AVG(mappability)                   AS mean_mappability
FROM coverage
GROUP BY chrom, pos
HAVING AVG(total_depth) < 20
ORDER BY mean_depth;
```

### Implementation steps

- [ ] Step 1: Add `CoverageRecord` struct to `src/record.rs` with the schema above
- [ ] Step 2: Add `CoverageArgs` to `src/cli.rs` and `Command::Coverage` variant
- [ ] Step 3: Add `src/coverage/mod.rs` — pileup loop that records all positions
  (including zero-depth if `--targets` given); compute BAM-derived mappability signals
  by iterating over all reads before MAPQ filter; handle `--bin-size` aggregation
- [ ] Step 4: Add `src/mappability.rs` — `MappabilityTrack` struct; BEDGraph loader;
  `get(chrom, pos) -> Option<f32>` with binary search; streaming mode for WGS
- [ ] Step 5: Add `src/writer/parquet_coverage.rs` — write `CoverageRecord` to Parquet
- [ ] Step 6: Update `src/main.rs` to handle `Command::Coverage`
- [ ] Step 7: Update `geac merge` (`src/merge.rs`) to detect coverage Parquets by schema
  and insert into `coverage` table in DuckDB
- [ ] Step 8: Add `wdl/geac_coverage.wdl` task and scatter/gather workflow
- [ ] Step 9: Integration tests — synthetic BAM; assert row count, depth values,
  `frac_mapq0` correctness for reads with known MAPQ=0
- [ ] Step 10: Explorer — "Coverage" tab (DuckDB mode only):
  - Systematically undercovered positions table (configurable depth and
    fraction-of-samples threshold); click to drill down to per-sample depths
  - Scatter plot of `mean_mapq` vs `mean_depth` per position — low-mappability
    sites cluster in the lower-left; color by `mappability` score if available
  - Per-sample depth distribution histogram across target positions

## WDL / Terra

- [x] WDL task wrapping `geac collect` — single-sample workflow in `wdl/geac_collect.wdl`
- [x] Terra-compatible Docker image — multi-stage `docker/Dockerfile` with htslib + geac binary
- [x] WDL workflow — scatter `geac collect` across a sample list, then gather with `geac merge` (`wdl/geac_cohort.wdl`)
- [x] WDL task wrapping `geac merge` — standalone workflow in `wdl/geac_merge.wdl`
- [ ] Test on Terra with a small cohort
