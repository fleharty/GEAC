# CLAUDE.md — GEAC project guide for Claude Code

## Architecture

GEAC has two independent layers:

1. **Rust CLI (`src/`)** — pileup-based BAM/CRAM processing. Produces Parquet files
   (one per sample) and DuckDB cohort databases.
2. **Python Streamlit apps (`app/`)** — interactive explorers that consume those files.

The schema contract between them lives in **`schema/geac_schema.json`**. Any new column
added to Rust output must be reflected there and checked against `app/explorer/schema.py`.

### Two separate explorers

| App | Launcher | Data |
|-----|----------|------|
| `app/geac_explorer.py` | `geac-cohort` | Alt-base locus/reads tables |
| `app/geac_coverage_explorer.py` | `geac-coverage-explorer` | Per-position coverage tables |

Shared code lives under `app/explorer/`. The two apps are otherwise independent.

## Build & test

```bash
# Rust
cargo build --release
cargo test                      # unit + integration tests

# Python
pytest app/tests/               # pure-Python helpers (no Streamlit runtime needed)

# Run explorers locally
streamlit run app/geac_explorer.py
streamlit run app/geac_coverage_explorer.py
```

## Key files

| Path | Purpose |
|------|---------|
| `schema/geac_schema.json` | Source of truth for all table column definitions |
| `app/explorer/schema.py` | Python mirror of the schema; `GEAC_VERSION` lives here |
| `app/explorer/filter_state.py` | `MAIN_FILTER_STATE` and `COVERAGE_FILTER_STATE` — all sidebar filter keys and defaults |
| `src/cli.rs` | All CLI flag definitions |
| `wdl/geac_collect.wdl` | Terra WDL for `geac collect` (mirrors CLI flags) |
| `wdl/geac_cohort.wdl` | Terra WDL for full cohort pipeline |
| `Cargo.toml` | Rust version — bump here when releasing |
| `CHALLENGES.md` | Log of non-obvious bugs and multi-attempt fixes |

## Rules of thumb

### When adding a CLI flag to Rust
1. Add it to `src/cli.rs`
2. Expose it in `wdl/geac_collect.wdl` and/or `wdl/geac_cohort.wdl`
3. If it produces a new output column, add it to `schema/geac_schema.json`
4. Update `README.md`

### When changing an Explorer feature
1. Check `README.md` — the Explorer UI section often needs updating
2. Check `docs/` for relevant doc files (e.g. `docs/provenance.md`)
3. If adding a new Streamlit widget with a `key=`, initialize its default in
   `MAIN_FILTER_STATE.defaults` (or `COVERAGE_FILTER_STATE`) so the session-state
   conflict warning doesn't fire after "Clear all"

### When releasing
Releases are gated on Rust/pipeline changes, not Explorer UI changes. UI improvements
ship in whatever version is current. To cut a release:
1. Bump `version` in `Cargo.toml`
2. Bump `GEAC_VERSION` in `app/explorer/schema.py`
3. Tag the commit (`git tag v0.X.Y && git push --tags`)

## Sensitive data

This is a **public** repository. Flag any content that could be patient-identifiable
(real sample IDs, file paths containing patient identifiers, genomic coordinates from
real clinical data) before committing.

## CHALLENGES.md

Append non-obvious bugs and fixes that required multiple attempts to `CHALLENGES.md`
as they occur. This prevents re-discovering the same pitfalls.

## Branch workflow

All work goes directly to `main`.
