# GEAC — TODO

## Rust / CLI

- [ ] `geac qc` subcommand — per-sample summary of error rates by substitution type, strand bias metrics, and overlap concordance
- [ ] Integration tests — generate synthetic BAM data and write end-to-end tests
- [ ] Cohort-level CLI subcommand — per-locus artifact frequencies across samples (e.g. flag positions seen in N% of samples)

## Explorer (Streamlit)

- [ ] Position-level drill-down — click a locus and see all samples/alleles at that position
- [ ] Export filtered data to CSV
- [ ] Cohort comparison view — side-by-side stats across samples loaded from a DuckDB

## WDL / Terra

- [ ] WDL task wrapping `geac collect` — process a single sample (BAM/CRAM + reference → Parquet)
- [ ] WDL task wrapping `geac merge` — aggregate per-sample Parquets into a cohort DuckDB
- [ ] WDL workflow — scatter `geac collect` across a sample list, then gather with `geac merge`
- [ ] Terra-compatible Docker image — package the `geac` binary with htslib dependencies
- [ ] Test on Terra with a small cohort
