# GEAC — TODO

Active backlog. Only unchecked, actionable items live here. For shipped work, design
notes, and historical context, see [docs/DEVELOPMENT_LOG.md](docs/DEVELOPMENT_LOG.md).
For high-level milestones and release themes, see [ROADMAP.md](ROADMAP.md).

Conventions:
- Items are grouped by area, not chronology.
- Each top-level bullet should be self-contained enough to start work on. If it
  needs deep context, link to the relevant section of `docs/DEVELOPMENT_LOG.md`.
- When an item is completed, move it (with any useful context) into
  `docs/DEVELOPMENT_LOG.md` rather than leaving a `[x]` here.

---

## Rust / CLI

- [ ] **MNV detection — Step 5: Integration test** — synthetic BAM with reads carrying
  substitutions at two adjacent positions vs reads carrying only one; verify `fragment_id`
  matches and the MNV candidates query returns correct `co_count` / `frac_cooccurring`.
  Steps 1–4 are complete (see `docs/DEVELOPMENT_LOG.md`).

---

## Per-read filter validation

- [ ] **Single-site read-level filter test** — run `geac collect --reads-output
  --region chr:start-end` on one sample restricted to the region of a known het
  variant. Single locus record with a known number of alt-supporting reads. Use
  the Explorer to manipulate family size and dist-from-read-end sliders; verify
  `alt_count` changes as expected. Confirm include vs exclude toggle behaviour
  matches intuition.

---

## Reads tab (Explorer)

Requires `alt_reads` table. Plots already shipped are listed in
`docs/DEVELOPMENT_LOG.md`.

- [ ] **Reads tab review** — work through all plots with real cohort data; assess
  usefulness, make changes, remove plots that don't add value.
- [ ] **Investigate trailing-N pattern after alt-supporting bases** — in some
  reads, an alt base is followed by a run of `N` bases. Determine whether this
  is primarily a property of the input BAM/consensus pipeline, an alignment/context
  effect (e.g. indels, read-end failure, masking), or something GEAC should
  explicitly annotate or visualize. Start with real examples in IGV and the Reads
  tab; collect compact read-context metrics in `alt_reads` (`n_before_alt`,
  `n_after_alt`, `n_n_before_alt`, `n_n_after_alt`, `leading_n_run_len`,
  `trailing_n_run_len`) and compare the before/after N burden distributions.
  Compare by pipeline/read_type where possible; assess whether these loci are
  enriched for low-confidence or pipeline-unique calls.
- [ ] **Family size vs VAF click-through** — add click/shift-click selection with
  drill-down table and IGV buttons, same as the strand bias plot. Currently
  blocked: `selection_point` with `on_select="rerun"` returns `{"fsvaf_select":
  {}}` regardless of what is clicked; strand bias works identically so the root
  cause is unknown. Investigate when Altair/Streamlit version context is clearer.

---

## Coverage Explorer (Customer-facing)

- [ ] Surface `feature_type` and `exon_number` in the intervals heatmap — available
  from v0.4.0; use `exon_number` as the x-axis instead of `interval_name` when
  present, so exons line up consistently across genes regardless of naming
  conventions.
- [ ] Surface BEDGraph track columns in the intervals GC scatter — detect any
  `Float32` columns beyond the fixed schema in `coverage_intervals` and offer them
  as a color-by selector.
- [ ] `config/panel_example.toml` — example config file for operators deploying
  the Coverage Explorer for a specific panel; documents all supported `geac.toml`
  keys with comments.
- [ ] **Longitudinal tracking** — add `run_id` and `run_date` provenance fields
  to `CoverageArgs` and `CoverageRecord`; surface a run-over-run depth trend line
  in the Summary tab.

---

## Coverage Analysis — Explorer integration

Steps 1–10 of `geac coverage` shipped in v0.4.0 (see `docs/DEVELOPMENT_LOG.md`).
Remaining:

- [ ] **Step 11: Explorer — "Coverage" tab (DuckDB mode only)**
  - Systematically undercovered intervals table (from `coverage_intervals`):
    configurable depth threshold and fraction-of-samples; columns include gene,
    exon_number, mean_gc, mean_mappability; explains whether dropout is GC,
    mappability, or other.
  - Scatter plot: `mean_gc_content` vs `mean_depth` per interval — GC bias
    visible as a U-shaped curve; color by `mean_frac_mapq0` to overlay
    mappability.
  - Per-sample depth distribution histogram; `frac_dup`, `frac_overlap`,
    `frac_soft_clipped` summary bars for QC overview.
  - Per-exon coverage heatmap (genes as rows, exons as columns, color = mean
    depth or `frac_at_30x`) — the primary customer-facing view.

---

## Fragmentomics

Long-term: end-motif table comparing alt-supporting vs reference-supporting reads
as a bait-bias signal. Design context in `docs/DEVELOPMENT_LOG.md`.

- [ ] Design the new table schema (per-read, not per-alt-read).
- [ ] Implement collection in `geac collect` (opt-in flag, like `--reads-output`).
- [ ] Add Explorer plot comparing end-motif frequency for alt vs ref reads.

---

## Code health / tech debt

Captured 2026-05-01 from a code review. Items here are not blocking any release;
they are organizational/structural improvements to take on between feature
milestones. Each "what to do" item below maps to a "why" weakness for context.

### Weaknesses observed

- **Python explorer is monolithic** — `app/geac_explorer.py` is ~7,853 LOC in a
  single file. Tab modules under `app/explorer/tabs/` are mostly stubs; the real
  logic is inline in the main file.
- **No app-level Python tests** — `app/tests/` covers ~400 LOC of helper
  functions, but neither `geac_explorer.py` nor `geac_coverage_explorer.py` has
  integration tests.
- **No reproducibility pinning** — there is no `rust-toolchain.toml`, no
  `=`-pinned versions in `Cargo.toml` (e.g. `rust-htslib`, `parquet`, `duckdb`),
  and no `pyproject.toml`/`ruff`/`mypy` config.
- **CLI metadata-flag boilerplate** — `src/cli.rs` has ~200 lines of copy-paste
  for the optional `--label1` / `--label2` / `--label3` style metadata flags.
- **Sparse module-level docs in Rust** — most `.rs` files have function-level
  doc comments but no `//!` describing module intent, invariants, or how the
  module fits into the larger pipeline.
- **No lint enforcement in CI** — no `clippy`, `rustfmt`, or `ruff` step gates
  merges. Style and lint regressions can land silently.
- **Schema contract is one-way at runtime** — `app/explorer/schema.py` loads
  `schema/geac_schema.json` at startup, but nothing fails the build if a column
  added to the schema is never referenced by the app, or vice versa.

### Proposed work items

- [ ] **Decompose `app/geac_explorer.py`** — introduce a `BaseTab` abstract class
  with a `render()` method; move each tab's logic into its own file under
  `app/explorer/tabs/` (the directory already exists). Goal: no single Python
  file >1,500 LOC. Unlocks unit tests of tab logic without a Streamlit runtime.
- [ ] **Pin Rust toolchain and critical dependencies** — add `rust-toolchain.toml`
  with the channel currently used in CI; switch `rust-htslib`, `parquet`, `arrow`,
  and `duckdb` in `Cargo.toml` to `=`-pinned versions. Document the bump cadence
  in `CLAUDE.md`.
- [ ] **Generate schema bindings at build time** — add a `build.rs` that reads
  `schema/geac_schema.json` and emits a `const SCHEMA_COLUMNS: &[&str]` per
  table; add a Python test that fails if any schema column is unreferenced in
  the explorer code.
- [ ] **Add module-level `//!` docs to all public Rust modules** — at minimum:
  `src/coverage/mod.rs`, `src/bam/pileup.rs`, `src/merge.rs`, `src/track.rs`,
  `src/writer/*`. Each should describe module purpose, key invariants, and
  dependencies.
- [ ] **Add lint/format CI** — new GitHub Actions job that runs `cargo fmt --check`,
  `cargo clippy -- -D warnings`, and `ruff check app/` on every PR. Fix existing
  warnings before turning the gate on.
- [ ] **Collapse CLI metadata flag boilerplate** — replace `--label1`/`--label2`/
  `--label3` (and any parallel groups) in `src/cli.rs` with a single repeatable
  `--label KEY=VALUE` flag backed by `Vec<(String, String)>`. Update WDLs and
  README.
- [ ] **Add Streamlit integration smoke tests** — at minimum, a test that imports
  each tab module and instantiates its render function against a tiny synthetic
  Parquet fixture. Catches import-time and obvious-render regressions without
  requiring a full browser.
- [ ] **Migrate this backlog into GitHub Issues** — once the team grows beyond a
  single contributor. Keep `ROADMAP.md` for milestones; let issues hold
  individual work items with priority/milestone labels.
