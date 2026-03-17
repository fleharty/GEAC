# GEAC — TODO

## Rust / CLI

- [ ] `geac qc` subcommand — per-sample summary of error rates by substitution type, strand bias metrics, and overlap concordance
- [ ] Integration tests — generate synthetic BAM data and write end-to-end tests
- [ ] Cohort-level CLI subcommand — per-locus artifact frequencies across samples (e.g. flag positions seen in N% of samples)
- [x] On-target annotation — add `--targets` BED flag to `geac collect`; records `on_target bool?` column in Parquet (null if no BED provided). Locus is on-target if it overlaps any interval in the BED.
- [x] Gene annotation — add `--gene-annotations` flag to `geac collect`; accepts GFF3 or GTF (auto-detected). Records `gene string?` column (null if not provided or intergenic). Explorer shows gene multiselect when column is populated.
- [x] Locus repetitiveness metrics — `homopolymer_len`, `str_period`, `str_len` columns computed from reference window. `--repeat-window` flag (default 10 bp). Explorer sliders for max homopolymer and STR length.

## Explorer (Streamlit)

- [x] IGV session download — manifest-driven BAM tracks + BED positions zip, capped at 5 samples
- [x] BED file has some incorrect entries — fixed: deletion loci now span the full deleted region [pos, pos+del_len) instead of just the anchor base
- [ ] Position-level drill-down — click a locus and see all samples/alleles at that position
- [x] Export filtered data to CSV — handled by Streamlit's built-in dataframe toolbar download button
- [ ] Cohort comparison view — side-by-side stats across samples loaded from a DuckDB
- [x] On-target filter — sidebar selectbox "Target bases": All / On target / Off target; depends on `on_target` column being populated (show warning if all null)
- [x] Gene filter — sidebar multiselect in the chromosome/region section; populated from distinct `gene` values; depends on `gene` column being populated
- [x] Repeat filter — sidebar sliders for max `homopolymer_len` and max `str_len`; show warning if columns are absent. Useful for excluding error-prone repetitive loci from analysis.

## WDL / Terra

- [x] WDL task wrapping `geac collect` — single-sample workflow in `wdl/geac_collect.wdl`
- [x] Terra-compatible Docker image — multi-stage `docker/Dockerfile` with htslib + geac binary
- [ ] WDL workflow — scatter `geac collect` across a sample list, then gather with `geac merge`
- [ ] WDL task wrapping `geac merge` — aggregate per-sample Parquets into a cohort DuckDB
- [ ] Test on Terra with a small cohort
