# GEAC — Roadmap

High-level milestones and release themes. Day-to-day backlog lives in [TODO.md](TODO.md).
Historical context, completed-work archives, and design notes live in
[docs/DEVELOPMENT_LOG.md](docs/DEVELOPMENT_LOG.md).

## Current state

- **Released:** v0.4.34 — Fusion caller quality overhaul: junction-coherence filter (`--min-coherent-fragments`, `--min-anchor-kmers`); `n_spanning_reads` and `n_coherent_fragments` output columns; seven correctness fixes to `fragment_top_pair`, `assign_gene`, and the secondary-output pipeline (see `docs/DEVELOPMENT_LOG.md`). WDL updated.
- **Released:** v0.4.33 — `geac experimental fusions --breakpoints-output`: per-fusion breakpoint coordinates derived from junction-spanning reads (median position, contributing-read counts, spread), computed in the same second BAM pass as `--kmer-hits-output`. Helper scripts committed: `reconstruct_fusions.sh` (de novo junction reconstruction, experimental) and `make_manifest.py`. Docs refreshed across `EXPERIMENTAL.md`, README WDL table, and `DEVELOPMENT_LOG.md`.
- **Released:** v0.4.32 — `kmer_pos_in_read` column added to `--kmer-hits-output` TSV (0-based k-mer offset within read sequence, enabling fusion breakpoint localization from spanning reads); `reads_bam_index` (.bai) added as a named WDL output in `geac_fusions.wdl`.
- **Released:** v0.4.31 — Fusion WDL workflows (`wdl/experimental/`) and `fusions --kmer-hits-output` memory fix: k-mer hit detail now collected in a targeted second BAM pass over fusion-supporting reads only, eliminating gigabytes of per-read accumulation on deep BAMs.
- **Released:** v0.4.30 — Fusion-index k-mer copy-number quantification and two-tier (panel vs genome) uniqueness. `build-fusion-index` gains `--max-genome-copies` (retain near-unique k-mers), `--gene-stats-output` (per-gene uniqueness + copy buckets), `--copy-histogram-output`, and `--bed-output-by-copies` (per-tier BEDs with a full-GTF gene-annotation column on tiers ≥2); the index now stores per-k-mer `genome_copies`. `fusions` gains `--max-kmer-copies` (re-tighten/relax uniqueness at call time) and tags `--reads-output` records with `FX:Z:GENEA::GENEB`. Experimental commands documented in [docs/EXPERIMENTAL.md](docs/EXPERIMENTAL.md).
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
| v0.4.30    | Fusion k-mer uniqueness tiers | Copy-number quantification + panel/genome two-tier uniqueness: `--max-genome-copies`, `--gene-stats-output`, `--copy-histogram-output`, `--bed-output-by-copies` (gene-annotated per-tier BEDs), per-k-mer `genome_copies` in the index; `fusions --max-kmer-copies` and `FX:Z:` fusion tag on `--reads-output`; experimental commands documented in `docs/EXPERIMENTAL.md` |
| v0.4.31    | Fusion WDLs + memory fix    | `wdl/experimental/geac_build_fusion_index.wdl` and `wdl/experimental/geac_fusions.wdl` for Terra; `fusions --kmer-hits-output` memory fix: k-mer detail collected in a second BAM pass over fusion-supporting reads only (was: all gene-assigned reads, causing gigabytes of peak RSS on deep BAMs) |
| v0.4.32    | Breakpoint localization     | `kmer_pos_in_read` (0-based k-mer offset within read sequence) added to `--kmer-hits-output` TSV; `reads_bam_index` (.bai) exposed as named output in `geac_fusions.wdl` |
| v0.4.34    | Fusion quality overhaul     | Junction-coherence filter (`--min-coherent-fragments`, `--min-anchor-kmers`); `n_spanning_reads` / `n_coherent_fragments` output columns; 7 correctness fixes to `fragment_top_pair` / `assign_gene` / secondary-output pipeline; WDL updated |
| v0.4.33    | Breakpoint output + scripts | `fusions --breakpoints-output` (per-fusion coordinates from junction-spanning reads, shares the `--kmer-hits-output` second BAM pass); `reconstruct_fusions.sh` and `make_manifest.py` committed; experimental fusion docs refreshed |
| **v0.5.0** | Analysis (next)             | Remaining customer-facing Coverage Explorer items; Reads tab review; trailing-N investigation                                                                                                                                              |
| v0.6.0     | External beta               | First release shared with external users for feedback; documentation polished; onboarding guide                                                                                                                                               |
| v1.0.0     | Stable                      | Feedback from beta incorporated; API/schema stable; production-ready                                                                                                                                                                          |

## Cross-cutting initiatives

These don't gate any single version, but are tracked alongside feature work:

- **Code health / tech debt** — Python explorer decomposition, dependency pinning,
  lint CI, schema-contract enforcement. See the Code health section of `TODO.md`.
- **Fragmentomics** — long-term direction toward end-motif analysis and nucleosome
  positioning. Design notes in `docs/DEVELOPMENT_LOG.md`.
