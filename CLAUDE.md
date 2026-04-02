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

This is a **public** repository. Before committing anything, verify that **none** of
the following are present:

- **Real genomic data** — no BAM/CRAM/FASTQ/VCF files or their derivatives; no real
  coverage Parquet or DuckDB files; no actual sequencing reads or variant calls.
- **External file paths** — no absolute paths or filenames that reference locations
  outside this repository (e.g. `/home/user/...`, `/data/projects/...`, GCS bucket
  URIs from real runs). These can reveal infrastructure details or project names.
- **Real sample identifiers** — the only sample names permitted in this repo are
  public reference standards: **HG001 / NA12878** and **HG002 / NA24385**.
  Any other sample name must not appear in code, tests, documentation, or data files.
- **Synthetic data** — even synthetic/simulated data should be reviewed before
  committing. If in doubt, flag it for review rather than committing immediately.

When in doubt, **stop and ask** rather than committing. It is much easier to add
something later than to scrub it from git history after the fact.

## CHALLENGES.md

Append non-obvious bugs and fixes that required multiple attempts to `CHALLENGES.md`
as they occur. This prevents re-discovering the same pitfalls.

## Branch workflow

All work goes directly to `main`.
