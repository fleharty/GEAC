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
- [ ] Read-end proximity metrics — record distance of each alt base from the nearest read end; enables detection of end-repair artifacts (alt bases enriched near read ends are likely artifactual). New columns: `dist_from_read_end` (int), optionally a boolean `near_read_end` flag with a configurable threshold (e.g. `--end-repair-window`, default 10 bp)
- [ ] Consensus family size metrics — for duplex/simplex reads, record the number of raw fragments that contributed to the consensus. New columns: `ab_count` (top-strand family size), `ba_count` (bottom-strand family size), `family_size` (ab_count + ba_count). Requires parsing fgbio/DRAGEN consensus tags (e.g. `cD`, `cE`, or `RX`/`MI` tags depending on pipeline). Useful for filtering low-confidence consensus calls.

## Intra-sample comparison (read-type)

- [ ] Multi-BAM collect — allow `geac collect` to accept multiple input BAMs for the same sample (raw, simplex, duplex) and tag each record with its `read_type`; output a single merged Parquet with a `read_type` column
- [ ] Read-type comparison view in Explorer — given a Parquet or DuckDB with mixed read types, show side-by-side or overlaid metrics (VAF distribution, strand balance, SBS96 spectrum) broken down by `read_type` (raw / simplex / duplex / mixed); goal is to quantify what duplex consensus processing removes vs retains

## Explorer (Streamlit)

- [x] IGV session download — manifest-driven BAM tracks + BED positions zip, capped at 5 samples
- [x] BED file has some incorrect entries — fixed: deletion loci now span the full deleted region [pos, pos+del_len)
- [x] Position-level drill-down — click a locus and see all samples/alleles at that position
- [x] Export filtered data to CSV — handled by Streamlit's built-in dataframe toolbar download button
- [x] On-target filter — sidebar selectbox "Target bases": All / On target / Off target
- [x] Gene filter — sidebar text input with partial match (ILIKE); depends on `gene` column being populated
- [x] Repeat filter — sidebar range sliders for `homopolymer_len` and `str_len`
- [x] Strand bias plot — dashed y=x diagonal + 95% binomial CI band; gene name in hover tooltip
- [x] SNV trinucleotide spectrum (SBS96) — 3×2 grid of per-mutation-type panels with shared y-axis; click drill-down
- [ ] Cohort comparison view — side-by-side stats across samples loaded from a DuckDB
  - [x] Step 1: Per-sample summary table — one row per sample_id with n_snv, n_insertion, n_deletion, mean_depth, mean_vaf, strand_balance, overlap_concordance; clicking a row filters all other tabs to that sample
  - [x] Step 2: VAF distribution overlay — all samples on one plot as density curves, colored by sample; highlights shifted VAF distributions
  - [ ] Step 3: Strand balance scatter — one dot per sample (x = mean strand balance, y = mean VAF); outliers immediately visible
  - [ ] Step 4: SNV count bar chart — n_snv per sample, stacked/colored by SBS6 substitution type breakdown
  - [ ] Step 5: SBS96 heatmap — samples as rows, 96 trinucleotide contexts as columns, color = normalized count; reveals samples with unusual mutational profiles
- [x] NMF decomposition — fit the per-sample SBS96 spectrum against COSMIC reference signatures using NNLS; show the largest contributing signatures and their weights

## WDL / Terra

- [x] WDL task wrapping `geac collect` — single-sample workflow in `wdl/geac_collect.wdl`
- [x] Terra-compatible Docker image — multi-stage `docker/Dockerfile` with htslib + geac binary
- [x] WDL workflow — scatter `geac collect` across a sample list, then gather with `geac merge` (`wdl/geac_cohort.wdl`)
- [x] WDL task wrapping `geac merge` — standalone workflow in `wdl/geac_merge.wdl`
- [ ] Test on Terra with a small cohort
