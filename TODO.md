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
- [ ] Integration tests — generate synthetic BAM data and write end-to-end tests
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

**Motivation:** Identify genomic regions that are systematically undercovered across a cohort —
positions where a meaningful fraction of samples fall below a minimum depth threshold.
This is complementary to artifact detection: low coverage reduces confidence in both
variant calls and clean sites.

### Design

New `geac coverage` subcommand:
- Inputs: BAM/CRAM + `--targets` BED/interval list + `--reference`
- Outputs: `{sample}.coverage.parquet` with one row per target position:
  - `sample_id`, `chrom`, `pos`, `total_depth`
- Covers every position in targets, including positions with zero alt reads
- Merged across samples with `geac merge` (second DuckDB table: `coverage`)

Cohort-level analysis via DuckDB:
```sql
SELECT chrom, pos,
       COUNT(DISTINCT sample_id)  AS n_covered,
       AVG(total_depth)           AS mean_depth,
       MIN(total_depth)           AS min_depth
FROM coverage
WHERE total_depth < 20
GROUP BY chrom, pos
HAVING n_covered > 0.5 * (SELECT COUNT(DISTINCT sample_id) FROM coverage)
ORDER BY n_covered DESC
```

### Implementation steps

- [ ] Step 1: Add `geac coverage` subcommand (`src/coverage.rs` + CLI args):
  `--input`, `--reference`, `--targets` (required), `--output`, `--sample-id`,
  `--min-map-qual`, `--threads`
- [ ] Step 2: Pileup over every position in `--targets`; record `total_depth` even if zero
- [ ] Step 3: Write output as `{stem}.coverage.parquet`
- [ ] Step 4: Update `geac merge` to ingest `.coverage.parquet` files into a `coverage`
  DuckDB table alongside `alt_bases`
- [ ] Step 5: Add WDL task and workflow for `geac coverage` scatter + merge
- [ ] Step 6: Explorer — add "Coverage" tab showing a table of systematically undercovered
  positions (configurable depth threshold and fraction-of-samples threshold);
  click a row to drill down to per-sample depths at that locus

## WDL / Terra

- [x] WDL task wrapping `geac collect` — single-sample workflow in `wdl/geac_collect.wdl`
- [x] Terra-compatible Docker image — multi-stage `docker/Dockerfile` with htslib + geac binary
- [x] WDL workflow — scatter `geac collect` across a sample list, then gather with `geac merge` (`wdl/geac_cohort.wdl`)
- [x] WDL task wrapping `geac merge` — standalone workflow in `wdl/geac_merge.wdl`
- [ ] Test on Terra with a small cohort
