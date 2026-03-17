# GEAC — TODO

## Rust / CLI

- [ ] `geac qc` subcommand — per-sample summary of error rates by substitution type, strand bias metrics, and overlap concordance
- [ ] Integration tests — generate synthetic BAM data and write end-to-end tests
- [ ] Cohort-level CLI subcommand — per-locus artifact frequencies across samples (e.g. flag positions seen in N% of samples)

## Explorer (Streamlit)

- [x] IGV session download — manifest-driven BAM tracks + BED positions zip, capped at 5 samples
- [ ] Position-level drill-down — click a locus and see all samples/alleles at that position
- [ ] Export filtered data to CSV
- [ ] Cohort comparison view — side-by-side stats across samples loaded from a DuckDB

## WDL / Terra

- [x] WDL task wrapping `geac collect` — single-sample workflow in `wdl/geac_collect.wdl`
- [x] Terra-compatible Docker image — multi-stage `docker/Dockerfile` with htslib + geac binary
- [ ] WDL workflow — scatter `geac collect` across a sample list, then gather with `geac merge`
- [ ] WDL task wrapping `geac merge` — aggregate per-sample Parquets into a cohort DuckDB
- [ ] Test on Terra with a small cohort
