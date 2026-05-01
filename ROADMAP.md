# GEAC — Roadmap

High-level milestones and release themes. Day-to-day backlog lives in [TODO.md](TODO.md).
Historical context, completed-work archives, and design notes live in
[docs/DEVELOPMENT_LOG.md](docs/DEVELOPMENT_LOG.md).

## Current state

- **Released:** v0.5.0 — `subject_id` and `sample_type` as schema-level sample annotations.
- **Branch policy:** all work goes directly to `main`.
- **Releases are gated on Rust/pipeline changes**, not Explorer UI changes. UI improvements
  ship in whatever version is current.

## Release theme history and forward plan

| Version    | Theme                       | Key deliverables                                                                                                                                                                                                                              |
|------------|-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| v0.3.1     | Bug fix                     | fgbio tag names corrected (`aD`/`bD`/`cD`); `family_size` populated for simplex reads                                                                                                                                                         |
| v0.3.2     | Per-read detail + Explorer  | `--reads-output` validated on Terra; multi-platform Docker (amd64+arm64) via GitHub Actions; Explorer embedded in Docker image; Reads tab with family size, position bias, base quality, mapping quality, and cohort artefact plots; insert size added to reads schema |
| v0.4.0     | Coverage                    | `geac coverage` command; coverage Explorer tab; WDL coverage workflow                                                                                                                                                                         |
| v0.5.0     | Sample annotations          | `subject_id` and `sample_type` as schema-level sample annotations                                                                                                                                                                             |
| **v0.6.0** | Analysis (next)             | MNV detection; remaining customer-facing Coverage Explorer items; Reads tab review; trailing-N investigation                                                                                                                                  |
| v0.7.0     | External beta               | First release shared with external users for feedback; documentation polished; onboarding guide                                                                                                                                               |
| v1.0.0     | Stable                      | Feedback from beta incorporated; API/schema stable; production-ready                                                                                                                                                                          |

## Cross-cutting initiatives

These don't gate any single version, but are tracked alongside feature work:

- **Code health / tech debt** — Python explorer decomposition, dependency pinning,
  lint CI, schema-contract enforcement. See the Code health section of `TODO.md`.
- **Fragmentomics** — long-term direction toward end-motif analysis and nucleosome
  positioning. Design notes in `docs/DEVELOPMENT_LOG.md`.
