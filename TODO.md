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
- [x] NMF decomposition — fit the per-sample SBS96 spectrum against COSMIC reference signatures using NNLS; show the largest contributing signatures and their weights

## WDL / Terra

- [x] WDL task wrapping `geac collect` — single-sample workflow in `wdl/geac_collect.wdl`
- [x] Terra-compatible Docker image — multi-stage `docker/Dockerfile` with htslib + geac binary
- [x] WDL workflow — scatter `geac collect` across a sample list, then gather with `geac merge` (`wdl/geac_cohort.wdl`)
- [x] WDL task wrapping `geac merge` — standalone workflow in `wdl/geac_merge.wdl`
- [ ] Test on Terra with a small cohort
