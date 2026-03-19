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
- [x] Strand bias selection: `pos` may be returned as float from Altair selection, causing SQL clause to silently fail; cast to int before building the WHERE clause
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
per-base error rates and identification of systematically undercovered sites across a cohort.

Coverage has three confounders that must be measured, not just controlled for:

1. **Mappability** — a region may appear undercovered because reads cannot be placed uniquely.
   Without a mappability signal, low coverage and multi-mapping are indistinguishable from
   genuine dropout (GC bias, FFPE degradation, probe failure).
2. **Duplicates** — PCR duplicates inflate raw read counts but represent the same original
   molecule. Per-region duplication rates reveal library complexity problems and can differ
   substantially between GC-rich and GC-poor targets.
3. **Fragment overlap** — when paired reads are longer than the insert, both reads cover the
   same bases but provide only one independent observation. High overlap inflates apparent
   depth while providing no additional evidence. The fraction of overlapping fragments is
   itself a useful QC signal (short inserts relative to read length).

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
  [--min-depth     0]                            # suppress positions with total_depth below this
  [--bin-size      1]                            # aggregate N bp into one row (1 = per-position)
  [--threads       1]
```

`--targets` is strongly recommended — it bounds output size and ensures zero-depth positions
are still recorded (complete dropout is important to capture). Without `--targets` or `--region`,
the whole BAM is scanned; fine for targeted panels, impractical for WGS without `--bin-size`.

`--track` can be repeated for multiple annotation tracks (e.g. mappability at two k-mer lengths,
a CpG density track, a GC content track). Each `NAME` becomes a column in the output Parquet.
Common sources: ENCODE GEM tracks (150-mer), genmap, Umap, custom BEDGraph from any tool.

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

This means `total_depth` is directly comparable to the `total_depth` in `alt_bases` from
`geac collect`, which uses the same duplicate-exclusion and overlap-collapsing logic.

### Mappability diagnostic table

| `frac_mapq0` | track score | `total_depth` | Interpretation |
|---|---|---|---|
| High | Low | Low | Classic multi-mapping — expected, filter confidently |
| High | High | Low | Unexpected low MAPQ — SV, misassembly, or aligner artifact |
| Low | Low | Low | Genuine dropout (GC bias, FFPE, probe failure) |
| Low | Low | Normal | Mappability track k-mer length may not match read length |

### Pre-computed annotation tracks (`--track`)

**Format**: BEDGraph (chrom, start, end, score). Covers ENCODE GEM, genmap, Umap, and any
custom track. bigWig support can be added later if needed.

**Multiple tracks** are useful in practice: e.g. mappability at the experiment's read length
alongside a GC-content track lets you disentangle GC bias from repeat-element dropout.

**Implementation**: each track is loaded into a sorted `Vec<(i64, i64, f32)>` per chromosome.
Lookup at each pileup position uses binary search (O(log n)). For targeted panels this fits
comfortably in memory. For WGS without `--targets`, a streaming approach (advance through the
sorted track in lock-step with the sorted pileup) avoids loading the full ~2 GB track.

Because the number of tracks is not known until runtime, `CoverageRecord` cannot be a plain
Rust struct with fixed fields. Instead, the Arrow `Schema` and `RecordBatch` are constructed
dynamically in `src/writer/parquet_coverage.rs` based on the track names provided.

### Gene annotation extensions

`--gene-annotations` reuses `src/gene_annotations.rs` but extends it to also store
`feature_type` (exon / intron / UTR / CDS) and `exon_number` from the GTF/GFF3 attribute
field. This enables Explorer queries like:

```sql
-- Coverage of BRCA1 exon 1 across all samples
SELECT sample_id, pos, total_depth, frac_mapq0, gem150
FROM coverage
WHERE gene = 'BRCA1' AND feature_type = 'exon' AND exon_number = 1
ORDER BY pos;
```

### DuckDB integration

`geac coverage` outputs `{sample}.coverage.parquet`. `geac merge` detects these by schema
(presence of `frac_dup`; absence of `alt_allele`) and inserts into a `coverage` table in
the cohort DuckDB alongside `alt_bases`.

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

### Implementation steps

- [ ] Step 1: Extend `src/gene_annotations.rs` — add `feature_type` and `exon_number` to
  the annotation lookup result; update GTF/GFF3 parser to extract the `exon_number` attribute
- [ ] Step 2: Add `src/track.rs` — `AnnotationTrack` struct; BEDGraph loader; binary-search
  lookup; streaming loader for WGS; `TrackSet` holding multiple named tracks
- [ ] Step 3: Add `CoverageArgs` to `src/cli.rs` with `--track` as a repeatable
  `Vec<String>` parsed as `NAME:path` pairs; add `Command::Coverage` variant
- [ ] Step 4: Add `src/coverage/mod.rs` — pileup loop that records every position in
  `--targets` / `--region` (including zero-depth); three-layer read counting
  (raw → de-dup → de-overlap → total_depth); compute all BAM-derived signals
- [ ] Step 5: Add `src/writer/parquet_coverage.rs` — dynamic Arrow schema built from
  fixed columns + one Float32 column per named track; write to Parquet
- [ ] Step 6: Update `src/main.rs` to handle `Command::Coverage`
- [ ] Step 7: Update `src/merge.rs` — detect coverage Parquets by schema and insert into
  `coverage` table; handle `union_by_name` for variable track columns across samples
- [ ] Step 8: Add `wdl/geac_coverage.wdl` task and scatter/gather workflow
- [ ] Step 9: Integration tests — synthetic BAM with known duplicates (BAM_FDUP reads),
  MAPQ=0 reads, and overlapping pairs; assert `frac_dup`, `frac_mapq0`, `frac_overlap`,
  and `total_depth` correctness
- [ ] Step 10: Explorer — "Coverage" tab (DuckDB mode only):
  - Systematically undercovered positions table with configurable depth threshold and
    minimum fraction-of-samples; columns include gene, exon_number, mean_mappability
  - Scatter plot: `mean_mapq` vs `mean_depth` per position; color by track score if
    available; low-mappability sites cluster visibly in lower-left
  - Per-sample depth distribution histogram; `frac_dup` and `frac_overlap` summary bars
  - Per-exon coverage heatmap when gene annotations are present (genes as rows,
    exons as columns, color = mean depth)

## WDL / Terra

- [x] WDL task wrapping `geac collect` — single-sample workflow in `wdl/geac_collect.wdl`
- [x] Terra-compatible Docker image — multi-stage `docker/Dockerfile` with htslib + geac binary
- [x] WDL workflow — scatter `geac collect` across a sample list, then gather with `geac merge` (`wdl/geac_cohort.wdl`)
- [x] WDL task wrapping `geac merge` — standalone workflow in `wdl/geac_merge.wdl`
- [ ] Test on Terra with a small cohort
