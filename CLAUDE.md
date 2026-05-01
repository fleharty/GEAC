# CLAUDE.md — GEAC project guide for Claude Code

## Architecture

GEAC has three layers:

1. **Rust CLI (`src/`)** — pileup-based BAM/CRAM processing. Produces Parquet files
   (one per sample) and DuckDB cohort databases.
2. **WDL workflows (`wdl/`)** — the primary way the Rust CLI runs in production on Terra
   or other cloud compute platforms. Wraps the CLI flags as WDL inputs and handles
   scatter-gather across cohorts.
3. **Python Streamlit apps (`app/`)** — interactive explorers that consume the Parquet
   and DuckDB files produced by the CLI.

The schema contract between layers 1 and 3 lives in **`schema/geac_schema.json`**. Any
new column added to Rust output must be reflected there and checked against
`app/explorer/schema.py`. When adding a CLI flag, also expose it in the relevant WDL.

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
| `TODO.md` | Active backlog — only unchecked, actionable items |
| `ROADMAP.md` | Release theme history and forward milestones |
| `docs/DEVELOPMENT_LOG.md` | Archive of completed work and design notes |

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
3. Bump the version in `VERSION` (used by `scripts/build_docker.sh` and the `run_*.sh` explorer launch scripts to tag/pull Docker images)
4. Update `ROADMAP.md` — mark the new version's row as released and update "Current state" if the theme has shifted
5. Tag the commit (`git tag v0.X.Y && git push --tags`)

### When completing a TODO item
1. Remove the item from `TODO.md` (don't leave a `[x]` — `TODO.md` is for unchecked work only)
2. If the work is non-trivial, has a useful design rationale, or the *why* is
   non-obvious, append a short entry to `docs/DEVELOPMENT_LOG.md` under the
   appropriate section. Trivial fixes don't need an entry — the commit message
   covers them.
3. If the item was a multi-step design (like the per-read detail table or coverage
   work), move the design context block into `docs/DEVELOPMENT_LOG.md` rather than
   discarding it.

### When starting forward work
1. Check `ROADMAP.md` to see which milestone the work belongs to.
2. Check `docs/DEVELOPMENT_LOG.md` for prior context — many features have design
   notes that explain constraints or rejected alternatives.
3. Add an actionable entry to `TODO.md` under the appropriate area heading.

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

## Project memory files

Four append-mostly markdown files capture different slices of project history.
Pick the right one when recording something:

- **`CHALLENGES.md`** — non-obvious bugs and fixes that required multiple attempts.
  Append as they occur to prevent re-discovering the same pitfalls.
- **`docs/DEVELOPMENT_LOG.md`** — completed features, shipped design decisions,
  and rejected alternatives ("won't do, because..."). Append when finishing a
  non-trivial item from `TODO.md`.
- **`TODO.md`** — only *unchecked*, actionable backlog items. When an item is
  done, move its useful context to `docs/DEVELOPMENT_LOG.md` and remove it from
  `TODO.md`.
- **`ROADMAP.md`** — release theme history and forward milestones. Update when
  cutting a release or when the next-version theme changes.

## Branch workflow

All work goes directly to `main`.
