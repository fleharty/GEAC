# GEAC — Roadmap

High-level milestones and release themes. Day-to-day backlog lives in [TODO.md](TODO.md).
Historical context, completed-work archives, and design notes live in
[docs/DEVELOPMENT_LOG.md](docs/DEVELOPMENT_LOG.md).

## Current state

- **Released:** v0.4.29 — Experimental command isolation: fusion / k-mer subcommands (`build-fusion-index`, `fusions`, `extract-gene`, `locate-kmer`) moved under `geac experimental <...>` to make their pre-production status explicit. Stable command surface unchanged. KmerIter position fix; BED output of unique k-mer intervals; `extract-gene` gains `--complement-output`, `--kmer-hits-output`, and redacted `--debug-region` diagnostics; `locate-kmer` gains `--gene-annotations`; `fusions` includes secondary/supplementary/unmapped reads (k-mer matching is alignment-independent) with `--min-mapq` defaulting to 0.
- **Branch policy:** all work goes directly to `main`.
- **Releases are gated on Rust/pipeline changes**, not Explorer UI changes. UI improvements
  ship in whatever version is current.

## Release theme history and forward plan

| Version    | Theme                       | Key deliverables                                                                                                                                                                                                                              |
|------------|-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| v0.3.1     | Bug fix                     | fgbio tag names corrected (`aD`/`bD`/`cD`); `family_size` populated for simplex reads                                                                                                                                                         |
| v0.3.2     | Per-read detail + Explorer  | `--reads-output` validated on Terra; multi-platform Docker (amd64+arm64) via GitHub Actions; Explorer embedded in Docker image; Reads tab with family size, position bias, base quality, mapping quality, and cohort artefact plots; insert size added to reads schema |
| v0.4.0     | Coverage                    | `geac coverage` command; coverage Explorer tab; WDL coverage workflow                                                                                                                                                                         |
| v0.4.23–24 | Sample annotations + refactor | `subject_id` and `sample_type` schema annotations; Explorer tab decomposition into `app/explorer/tabs/`                                                                                                                                   |
| v0.4.25    | MNV detection               | `fragment_id` (FNV-1a hash of qname) in `alt_reads`; MNV candidates tab; Overlap agreement tab                                                                                                                                            |
| v0.4.26    | MNV integration test + IGV  | MNV integration test (`reads_fragment_id_enables_mnv_detection`); IGV session fixes (manifest relative-path resolution, index existence check, always-visible sample multiselect); `make_bed` tested via `igv_helpers`                     |
| v0.4.27    | Pileup-cap fix + indel filters | `--max-pileup-depth` (0 = unlimited) on `collect`, `coverage`, `annotate-normal`; raises the htslib 8000-read default that silently downsamples high-coverage loci and dropped a 42 bp deletion at 14 k× depth. Deletion `on_target` now evaluated by range overlap, not anchor position. `pipeline` made nullable across all output Parquet schemas. Explorer indel-length range slider; MNV inter-substitution min-distance slider |
| v0.4.28    | Embedded resource paths     | `bam_path`, `variants_path`, `gnomad_path` stored in Parquet output and aggregated into DuckDB `samples` table; `--bam-uri` / `--variants-uri` / `--gnomad-uri` flags for Terra/cloud runs (store cloud URI instead of localized path); Explorer pre-populates IGV resources from database with manifest/geac.toml override priority |
| v0.4.29    | Experimental fence + fusion groundwork | `geac experimental <...>` namespace for `build-fusion-index`, `fusions`, `extract-gene`, `locate-kmer` (stable surface unchanged); KmerIter `(kmer, position)` fix; `--bed-output` for build-fusion-index; `--complement-output` / `--kmer-hits-output` / redacted `--debug-region` for extract-gene; `--gene-annotations` for locate-kmer; fusions includes secondary/supplementary/unmapped reads with `--min-mapq` default 0 |
| **v0.5.0** | Analysis (next)             | Remaining customer-facing Coverage Explorer items; Reads tab review; trailing-N investigation                                                                                                                                              |
| v0.6.0     | External beta               | First release shared with external users for feedback; documentation polished; onboarding guide                                                                                                                                               |
| v1.0.0     | Stable                      | Feedback from beta incorporated; API/schema stable; production-ready                                                                                                                                                                          |

## Cross-cutting initiatives

These don't gate any single version, but are tracked alongside feature work:

- **Code health / tech debt** — Python explorer decomposition, dependency pinning,
  lint CI, schema-contract enforcement. See the Code health section of `TODO.md`.
- **Fragmentomics** — long-term direction toward end-motif analysis and nucleosome
  positioning. Design notes in `docs/DEVELOPMENT_LOG.md`.
