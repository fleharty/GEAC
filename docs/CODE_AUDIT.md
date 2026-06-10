# GEAC Code Audit Log

Ongoing notes from targeted code audits. These are potential bugs, design risks,
missing safeguards, or questionable feature behavior found during review. Items
are not necessarily release blockers; promote them into `TODO.md` or GitHub
Issues when they become scheduled work.

Conventions:
- Newest audits go at the top.
- Keep entries scoped to the code area reviewed.
- Prefer concrete file references and user-visible impact over broad refactor
  suggestions.
- Include verification performed, even when no code changes were made.

---

## 2026-06-10 — whole-codebase pass: explorer layer, cross-layer contracts, structural debt

Scope:
- Whole-repo evaluation requested ("evaluate the code as a whole, not just recent
  features"), with emphasis on the areas the 2026-06-09 cross-cutting entry under-covered:
  the ~15 K-line Python explorer layer (`app/`), the cross-layer contracts
  (`schema/geac_schema.json` ↔ Rust writers ↔ `app/explorer/schema.py`), the test strategy,
  and documentation. Deliberately deduped against the 2026-06-09 entry and the `TODO.md`
  "Code health / tech debt" section — items already logged there are referenced, not
  repeated.

Verification:
- Test inventory: 101 unit + 45 integration (Rust) + 99 `pytest app/tests` = 245 total.
- Structural checks via grep/read over `src/`, `app/`, `wdl/`, `schema/`, `docs/`.
- No code changes; findings only. Magnitude note: a full count of
  `unwrap/expect/panic/unreachable` in `src/` is 60 total, the majority inside
  `#[cfg(test)]` modules — the project is not panic-heavy.

Findings — cross-layer contract:

- [ ] **Schema-contract test is one-directional.** `tests/integration.rs:30`
  (`assert_schema_columns_present`) verifies schema-declared columns are present in emitted
  tables, but only for explicitly named tables, and never the reverse: nothing asserts that
  every table `src/merge.rs` creates appears in `schema/geac_schema.json`. That gap is the
  root cause of the drift the 2026-06-09 entry flagged — `coverage_intervals` and `fusions`
  are created by `src/merge.rs`/`src/main.rs` but absent from both the `tables` map and the
  `feature_tables` array in the JSON. Because `app/explorer/schema.py` reads
  `raw["feature_tables"]` at runtime, `DataSource.has_optional_table("fusions")` /
  `("coverage_intervals")` can never return True. Fix: add both tables to the JSON and add a
  test enumerating merge `TableSpec`s that asserts each is in the schema.
- [ ] **`DUCKDB_SCHEMA_VERSION` is stamped but never validated.** Defined at
  `src/merge.rs:12` (`"duckdb-v4"`) and written into the metadata table (`:842`), but merge
  never checks the schema version of incoming data (no compare/bail anywhere). Combined with
  the open forward item to stamp version into per-Parquet metadata, merging inputs produced
  by an incompatible GEAC version can silently truncate columns rather than failing.

Findings — refactor / duplication:

- [ ] **Parquet-writer boilerplate across 8 modules.** Every `src/writer/*.rs`
  (`parquet.rs`, `parquet_coverage.rs`, `parquet_coverage_intervals.rs`, `parquet_reads.rs`,
  `parquet_fragments.rs`, `parquet_normal.rs`, `parquet_pon.rs`, `parquet_sample_metrics.rs`)
  repeats the same `File::create` → `WriterProperties::builder().set_compression(SNAPPY)` →
  `ArrowWriter::try_new` → `write` → `close` shell, varying only the schema and
  record-to-batch functions. Extract one `write_parquet_generic<T>(records, path, schema_fn,
  batch_fn)` to remove ~150-200 LOC and centralize error/compression policy.
- [ ] **`reads.py` is one 1345-line `render()` with nested closures.**
  `app/explorer/tabs/reads.py:14` is the only top-level `def`; all helpers
  (`_nctx_cache_key`, `_nasym_cache_key`, `_cache_get`, …) are closures that capture
  `where`/`_r_reads_filter`, and cache keys are built from those captures. Works today, but
  fragile under refactor and not unit-testable in isolation. The existing "decompose
  `geac_explorer.py`" TODO targets a different file and a >1500-LOC bar this file slips
  under, so call it out separately.
- [ ] **Two SQL-escaping idioms across the two explorers.** `geac_explorer.py` uses the
  `_sql_str()` helper; `geac_coverage_explorer.py` (~line 194) uses inline
  `replace(chr(39), chr(39)*2)`; many tabs interpolate column names/values directly. The
  2026-06-09 entry and `TODO.md` already flag the f-string SQL pattern — the point here is
  that the right fix is one shared escaping/parameterization helper used by both apps, not
  per-tab patches.
- [ ] **Cross-explorer sidebar duplication.** Region input, sample selection, optional-label
  multiselects, and manifest wiring are duplicated between `geac_explorer.py` and
  `geac_coverage_explorer.py`. Extract `app/explorer/sidebar.py` shared by both entrypoints
  (distinct from the tab-decomposition item already in `TODO.md`).

Findings — missing safeguards:

- [ ] **No end-to-end Rust → DuckDB → Python test.** Rust integration tests verify output
  against the schema, and Python tests read a DuckDB, but nothing chains
  `collect → merge → open the DuckDB from Python → assert a known locus's counts`. The
  cross-layer contract is verified only piecewise.
- [ ] **No bounds validation on `build-fusion-index` numeric flags.** `min_k`/`max_k` are
  not cross-checked (e.g. `min_k=30, max_k=20`), `min_gene_kmers` has no upper bound, and
  `max_genome_copies` is only meaningful with `--check-genome-uniqueness` but is not gated on
  it. Same class as the open `--min-sample-fraction` range item.
- [ ] **A few non-test `unwrap()`s on possibly-empty vectors / post-check Options.**
  `src/fusions.rs:1097-1100` (`max()/min().unwrap()` on `a_positions`/`b_positions` in
  spanning-read logic), `:1122`/`:1127` (`as_mut().unwrap()`, `get_mut(...).unwrap()` after a
  guard), and `src/coverage/mod.rs:810-811,1144-1145,1156-1157` (`min()/max().unwrap()` in
  stats helpers). Replace with `ok_or_else`/empty-guards so edge-case data yields a clean
  error rather than a panic.

Findings — polish:

- [ ] **WDL resource doc/default drift.** `wdl/geac_collect.wdl:50-52` documents memory
  default 8 / disk 100 / preemptible 2, but the actual input defaults are 32 / 100 / 2
  (`:107-109`). Resources are correctly exposed as overridable inputs; only the header
  comment is wrong. Reconcile across sibling WDLs.
- [ ] **Doc duplication.** `README.md:148` "Data Model" suffix→table table duplicates
  `docs/schema.md:6` "Table Routing"; keep one canonical copy and link the other.
  `docs/EXPERIMENTAL.md` (514 L) and `docs/FUSION_DEVELOPMENT.md` (152 L) overlap on fusion
  scope — add a one-line scope header to each pointing at the other.
- [ ] **Memory-file growth.** `CHALLENGES.md` (~72 KB) and `docs/DEVELOPMENT_LOG.md`
  (~61 KB) remain valuable and well-organized; no action now, but if either passes ~100 KB
  consider archiving older sections under `docs/archive/`.

Good (recorded so it is preserved, not changed):
- Streaming writers, the `experimental` namespace isolation (extraction-ready), trait-based
  annotation (`VariantAnnotator`/`TargetIntervals`), `anyhow`+`.context()` error
  propagation, the explorer's `DataSource`/`FilterState`/`TabContext` abstractions (and
  `error_spectrum.py`'s 13-helper decomposition as the model for other tabs), the partial
  schema-contract test that already exists, and the `CHALLENGES.md`/`DEVELOPMENT_LOG.md`
  memory discipline.

Promoted to `TODO.md` (this pass): test-running CI job; schema-JSON drift fix + reverse
contract test; shared Parquet-writer helper; end-to-end Rust→DuckDB→Python test; WDL
resource doc/default reconciliation.

---

## 2026-06-09 — cross-cutting audit: CI, scope, collect performance, schema contract

Scope:
- Whole-repo pass rather than a single code area: pileup core (`src/bam/pileup.rs`,
  `src/bam/mod.rs`), `src/merge.rs`, `src/record.rs`, `src/inspect.rs`, the full CLI
  surface, the explorer data layer (`app/explorer/data_source.py`),
  `schema/geac_schema.json`, `.github/workflows/`, and the Homebrew formula.
- Focus on architecture, build/release enforcement, schema-contract integrity,
  and feature-level scope.

Verification:
- `cargo test` against the existing target: 101 unit + 45 integration passing.
- `pytest app/tests`: 99 passing.
- No code changes made; findings only.

Findings — build & release enforcement:

- [ ] **No CI runs the test suite.** `.github/workflows/release.yml` is the only
  workflow and only builds/pushes the Docker image on `v*` tags. Nothing runs
  `cargo test`, `pytest app/tests`, `cargo clippy`, or `cargo fmt --check` on
  push/PR. The project has 200+ tests but no automated proof the tree is green. A
  minimal `ci.yml` running `cargo test` + `pytest app/tests` on push is the
  highest-leverage, lowest-effort fix in the repo and a prerequisite before the
  lint gate already noted in `TODO.md`.

Findings — scope / strategic:

- [ ] **Experimental fusion subsystem is >20% of the Rust code and gates the
  release cadence.** `fusions.rs` (1,326), `build_fusion_index.rs` (1,011),
  `extract_gene.rs` (458), plus `locate_kmer`/`lookup_kmer`/`scan_read`/
  `compute_uniqueness_map`/`build_fusion_kmer_blacklist` total ~3,400 LOC, all
  under the `geac experimental` "not for production" banner. `ROADMAP.md` shows
  v0.4.29–v0.4.37 were almost entirely fusion work while customer-facing
  coverage/analysis (v0.5/v0.6 beta) waited. Decide deliberately: (a) freeze and
  extract into a separate `geac-fusion` crate/repo so the core stays lean and
  fusion work matures on its own cadence (the clean `experimental` namespace makes
  this low-friction), or (b) commit to it as a product bet and build the
  benchmarking harness — already a TODO — first, since it competes with validated
  tools (Arriba, STAR-Fusion, FusionCatcher). Either way, stop letting
  pre-production work set release timing.
- [ ] **Demote fusion dev/debug subcommands from the user surface.** `locate-kmer`,
  `lookup-kmer`, and `scan-read` are developer diagnostics for the fusion index,
  not end-user features. If the fusion code stays, fold them behind
  `geac experimental debug <...>` (or a debug build) to shrink the advertised CLI.

Findings — `collect` performance & memory:

- [ ] **`collect` buffers all output in RAM.** `collect_alt_bases` returns fully
  materialized `Vec<AltBase>`/`Vec<AltRead>` and `src/main.rs:144` writes them only
  after the whole BAM is processed, whereas `coverage` and `fragments` stream via
  `CoverageWriter`/`FragmentsWriter`. On a deep WGS sample with `--reads-output`
  the per-read vector alone can be many GB of peak RSS. Refactor `collect` to a
  streaming writer; this also unblocks the parallelism item below.
- [ ] **No single-sample parallelism.** `rayon` is a dependency but used only in
  `src/sample_identity.rs`. The `collect` pileup loop is strictly sequential, one
  contig/region at a time. Cohort scale is handled by WDL scatter across *samples*,
  but a single deep BAM cannot use more than one core. Scatter the pileup across
  contigs/regions and merge per-region record vectors — the largest untapped
  single-sample throughput lever.

Findings — correctness & safety guards:

- [ ] **No reference/BAM concordance check.** `collect` opens the reference FASTA
  and the BAM independently (`src/bam/mod.rs:51`). If they disagree (wrong build,
  contig-naming mismatch, shifted coordinates) every base silently reads as alt and
  the output is quietly, catastrophically wrong. Compare `@SQ` contig names/lengths
  (and `M5` where present) against the reference `.fai` at startup and fail fast.
  Cheap, and prevents an entire class of silently-wrong runs.
- [ ] **Per-locus output is not deterministic.** Alt-base tallies are built in a
  `HashMap<char, BaseTally>` (`src/bam/pileup.rs:327`) and iterated in hash order at
  `src/bam/mod.rs:196`, so two runs can emit a locus's alt alleles in different row
  order. Given the project's investment in checksums/provenance, sort alt alleles
  before emit to guarantee byte-identical output across runs.

Findings — schema contract:

- [ ] **`schema/geac_schema.json` has drifted from the actual tables.** `CLAUDE.md`
  calls it the source of truth, but it omits the `coverage_intervals` and `fusions`
  tables entirely, even though `src/merge.rs` and `src/inspect.rs` create/expect
  them. Concrete downstream bug: the Python `feature_tables` list (mirrored from
  the JSON) omits both, so `DataSource.has_optional_table("fusions")` and
  `("coverage_intervals")` can never return True regardless of whether the table
  exists. Add both tables to the schema, or document why they are intentionally
  outside the contract. This is the exact failure the `TODO.md` "schema contract is
  one-way" item predicted.

Findings — explorer:

- [ ] **Systemic raw f-string SQL.** Nearly every tab builds SQL with f-strings.
  The Rust layer escapes literals; the Python layer escapes *paths* (`_sql_str`)
  but interpolates column names and selected values directly. Inputs are mostly
  trusted (sample IDs from the data itself), so this is defense-in-depth, but
  `TODO.md` already flags three tabs (fragmentomics, sample_identity,
  pipeline_comparison) with the pattern — better solved with one shared
  escaping/parameterization helper than three point fixes.

Findings — forward features (not yet in `TODO.md`/`ROADMAP.md`):

- [ ] **Add a cohort background-error / recurrence-statistics model — the actual
  "Atlas."** `src/cohort.rs` only counts recurrence (`COUNT(DISTINCT sample_id)`).
  The differentiating feature latent in the collected data is a position-specific
  (and per-trinucleotide-context) background error model: fit a per-site noise rate
  across the cohort (beta-binomial or empirical), then score each observation as
  consistent-with-background vs. outlier. This separates recurrent artifacts from
  recurrent biology and turns raw alt-base evidence into a confidence signal that
  per-sample callers don't provide. Highest-value forward feature.
- [ ] **Stamp `geac_version`/`schema_version` into per-Parquet file metadata.** The
  merged DuckDB carries `geac_metadata`, but individual `.parquet` files don't embed
  version/schema info in Parquet key-value file metadata. Merging an old Parquet
  later leaves no in-file record of which binary/schema produced it. Stamp at write
  time to complete the provenance story.
- [ ] **Add a computed strand-bias statistic column.** `fwd_alt_count`/
  `rev_alt_count` are collected and there is a strand-bias tab, but no per-locus
  statistic (Fisher exact p or SOR) exists in the schema. Compute it once in Rust so
  the explorer can filter/sort on a number instead of re-deriving balance visually.
- [ ] **Parser fuzz/property tests.** Prior audit entries list many BED/Picard/GTF/
  VCF/TSV robustness issues (silent row-skips, `i64→u32` casts, `end<=start`). A
  small `cargo-fuzz`/`proptest` harness over the parsers would surface this class in
  bulk and guard against regressions instead of fixing them one at a time.
- [ ] **`geac doctor` input preflight.** `inspect` validates *outputs*; there is no
  input-side analogue. A preflight that checks BAM index present, reference
  concordance (above), target-file parse, and required consensus tags present for
  the chosen `--pipeline` would save cloud time on long runs.

Findings — smaller items:

- [ ] **Duplicate inner attribute.** `#![deny(unsafe_code)]` appears twice
  (`src/main.rs:1` and `:2`). Harmless, but exactly what a `clippy`/`fmt` gate would
  catch.
- [ ] **Phantom `geac compare` command.** `src/cli.rs:123` and `src/record.rs:100`
  document `subject_id` as "Used by `geac compare`," but no `Compare` variant exists
  in the `Command` enum and the README lists no such command. Build it or strike the
  references.
- [ ] **Stray root file.** `cic_dux4_genes.txt` (3-line CIC::DUX4 gene list) is
  tracked at repo root; example/panel files belong in `examples/` or `config/`.
- [ ] **Homebrew formula ships placeholder hashes.** `homebrew/Formula/geac.rb`
  contains literal `sha256 "SHA256_MACOS_ARM64"` / `"SHA256_SOURCE"`. If
  `brew install fleharty/geac/geac` is meant to work today, verify the published tap
  has real SHAs. The formula is arm64-mac only (Intel mac/Linux fall back to
  Docker/source).

---

## 2026-06-08 — variant annotation from VCF and TSV

Scope:
- `src/vcf.rs`
- `src/variants_tsv.rs`
- `src/bam/mod.rs` `vcf_annotation()` consumer
- `src/cli.rs` `--vcf` / `--variants-tsv` arguments

Verification:
- `cargo test -q`
- Result: full Rust test suite passed: 101 unit tests and 44 integration tests.

Findings:

- [ ] **Remove SNV position fallback or make it explicit** — both VCF and TSV
  annotators first try exact `(chrom, pos, alt)` matching for SNVs, then fall
  back to `(chrom, pos)` when the allele is absent. That means an observed `G`
  at a site where the caller reported only `A` is marked `variant_called=true`
  and inherits the caller's filter. This can inflate called status and confuse
  pipeline comparison/error-spectrum views. Keep position fallback for indels if
  necessary, but SNVs should probably require exact allele matching.

- [ ] **Make indel annotation allele-aware instead of position-only** —
  indels are annotated by position alone because GEAC `+seq`/`-seq` notation
  differs from VCF/TSV representation. A site with multiple different indels, or
  an SNV and indel at the same anchor, can mark the wrong GEAC indel as called
  and assign the wrong filter. Reuse the GEAC-to-VCF indel conversion logic from
  gnomAD annotation, or add a canonical indel normalization step for both VCF
  and TSV sources.

- [ ] **Validate and report malformed `variants_tsv` rows** —
  `src/variants_tsv.rs` silently skips rows with fewer than five fields and rows
  whose `pos_start` does not parse. A truncated or malformed variant list can
  annotate only part of a sample without any warning. Header/comment rows can be
  skipped deliberately, but malformed data rows should fail with line number and
  file context.

- [ ] **Use `pos_end` and `ref` from `variants_tsv` for consistency checks** —
  the TSV loader documents columns `chrom pos_start pos_end ref var`, but it only
  uses `chrom`, `pos_start`, and `var`. It does not validate coordinate span,
  reference allele, or whether a row is SNV vs indel. This makes off-by-one or
  incompatible TSVs difficult to detect and contributes to false-positive
  position fallback. Validate `pos_end > pos_start`, compare `ref` where
  possible, and route indel rows through explicit allele normalization.

- [ ] **Surface duplicate annotation records deterministically** — `by_position`
  and `by_allele` use `or_insert_with`, so the first record encountered wins for
  duplicate VCF/TSV records at the same allele or position. File order can decide
  whether a locus gets `PASS` or a filtered status. If duplicate records are
  expected, define a precedence rule (for example non-PASS over PASS, or all
  filters joined); otherwise detect duplicates and warn/fail.

- [ ] **Add direct tests for annotation loaders** — there are no focused tests for
  `VcfIndex` or `VariantsTsv` behavior. Add fixtures covering exact SNV match,
  same-position different-allele non-match, indel matching, duplicate records,
  malformed TSV rows, and `post_filter_flag` parsing (`NoRules` -> `PASS`).

---

## 2026-06-08 — target interval parsing and on-target lookup

Scope:
- `src/targets.rs`
- `src/bam/mod.rs` on-target calls for SNVs/indels
- `src/coverage/mod.rs` target zero-fill and interval summaries
- CLI target documentation

Verification:
- `cargo test -q targets`
- Result: 18 target unit tests and 3 target-related integration tests passed.

Findings:

- [ ] **Reject invalid or empty target intervals at parse time** —
  `src/targets.rs` accepts intervals where `end <= start` for BED, and Picard
  interval rows where the converted 0-based start is greater than or equal to
  the end. Zero-length intervals silently disappear from position iteration but
  remain in `named_intervals()`, while reversed intervals can underflow in
  `total_bases()` (`e - s` on `u32`) or create nonsensical interval summaries.
  Validate `start < end` after coordinate conversion and fail with line context.

- [ ] **Avoid casting negative query positions to `u32`** —
  `contains()`, `overlaps()`, and `interval_index()` cast `i64` positions to
  `u32`. A negative position becomes a very large coordinate, so callers with a
  bad upstream coordinate can get false positives against high-coordinate
  intervals instead of a clear false/validation failure. Return false for
  negative `pos`, `start`, or `end_exclusive`, and handle empty query intervals
  explicitly.

- [ ] **Clarify or fix overlapping named interval assignment** —
  merged intervals are used for fast `contains()`/`overlaps()`, but coverage
  interval summaries use original named intervals and `interval_index()`.
  For overlapping named intervals, `interval_index()` returns only the last
  interval whose start is <= the position, so bases in earlier overlapping
  intervals are not assigned to every interval they belong to. Either reject
  overlapping named targets when interval summaries are requested, or update
  interval summary accumulation to assign a position to all overlapping named
  intervals.

- [ ] **Make target format detection stricter** — `TargetIntervals::parse()`
  treats the entire file as Picard interval-list format if any line starts with
  `@`; otherwise every data line is parsed as BED. A BED file with an accidental
  header/comment line beginning with `@`, or a malformed mixed file, can be
  interpreted with Picard coordinates and name columns. Prefer explicit format
  detection from valid Picard headers (`@HD`/`@SQ`) plus row shape validation, or
  add a `--targets-format` override for ambiguous files.

- [ ] **Do not silently skip malformed short target rows** — rows with fewer than
  three tab-delimited fields are ignored. That is convenient for comments, but
  it also hides truncated target rows and can produce an incomplete target set
  without warning. Skip only blank/comment/header lines; otherwise fail on
  malformed rows with the line number.

---

## 2026-06-08 — gnomAD allele-frequency annotation

Scope:
- `src/gnomad.rs`
- `src/bam/mod.rs` calls into `GnomadIndex::get()`
- `src/cli.rs` gnomAD-related collect arguments
- README gnomAD annotation documentation

Verification:
- `cargo test -q gnomad`
- Result: 4 unit tests passed.

Findings:

- [ ] **Add integration tests with a real indexed VCF/BCF fixture** —
  `src/gnomad.rs` currently has unit tests for GEAC-to-VCF allele string
  conversion, but no test exercises `GnomadIndex::open()`, chromosome RID
  resolution, tabix/CSI fetch behavior, INFO extraction, multi-allelic records,
  or custom `--gnomad-af-field`. A tiny bgzip+tabix fixture should cover SNV,
  insertion, deletion, chr-prefix mismatch, multi-allelic AF indexing, and a
  missing/custom AF field.

- [ ] **Warn or fail when the requested AF INFO field is absent everywhere** —
  `GnomadIndex::extract_af()` returns `None` when the configured INFO field is
  missing or not a float. That is correct for an allele with no AF, but a typo
  such as `--gnomad-af-field AF_jonit` silently makes every `gnomad_af` null.
  Track whether the field is present in the header and/or observed during lookup,
  then emit a clear warning or fail at open time when the requested field is not
  available.

- [ ] **Do not silently treat all gnomAD fetch errors as allele absence** —
  `GnomadIndex::get()` maps any `IndexedReader::fetch()` error to `Ok(None)`.
  That hides index/reference problems as "allele absent from gnomAD" and can
  produce a cohort with systematically null `gnomad_af` values. Missing contigs
  can remain non-fatal, but malformed indices, invalid regions, or reader errors
  should be counted and reported, or promoted to an error.

- [ ] **Document and/or handle normalized-equivalent indel misses** —
  gnomAD lookup requires exact same `record.pos()`, REF, and ALT after converting
  GEAC `+seq`/`-seq` notation to VCF form. Equivalent indels can be represented
  at different left-aligned positions in repetitive sequence, so an allele that
  is biologically present in gnomAD may still return null. Either document this
  limitation explicitly for `gnomad_af`, or add normalization-aware lookup using
  reference context and a wider fetch window.

- [ ] **Validate custom AF field cardinality expectations** — the code assumes the
  selected AF field returns one float per ALT allele and indexes by `alt_idx`.
  That matches standard `Number=A` AF fields, but the CLI allows "any other INFO
  key"; scalar or differently-cardinalized fields can return misleading nulls
  for non-first ALT alleles. Restrict accepted fields to allele-indexed float
  INFO definitions when possible, or document the requirement and surface a
  warning for incompatible header metadata.

---

## 2026-06-08 — `geac cohort` and `geac qc`

Scope:
- `src/cohort.rs`
- `src/qc.rs`
- `src/cli.rs` argument definitions for `CohortArgs` and `QcArgs`
- `tests/integration.rs` cohort/QC happy-path coverage

Verification:
- `cargo test -q cohort_finds_recurrent_locus`
- `cargo test -q qc_exits_successfully`
- Result: both targeted integration tests passed.
- Note: an initial combined invocation
  `cargo test -q cohort_finds_recurrent_locus qc_exits_successfully` was rejected
  by Cargo before running tests because `cargo test` accepts a single test-name
  filter.

Findings:

- [ ] **`geac cohort` can inflate sample fractions after filtering** —
  `src/cohort.rs` computes `total_samples` from rows that survive the `data`
  view. With `--on-target-only`, any sample with no on-target alt rows drops out
  of the denominator, inflating `sample_fraction`. The same denominator issue
  exists for samples with zero alt-base rows in general. For recurrence
  reporting, the denominator should probably come from the input sample set or a
  manifest/sample table rather than from samples visible after locus filtering.

- [ ] **Duplicate `sample_id` inputs distort cohort aggregate metrics** —
  `src/cohort.rs` counts recurrence with `COUNT(DISTINCT sample_id)`, but VAF,
  depth, alt-count, and strand metrics average all rows for a locus. If the same
  `sample_id` appears in multiple files, reruns, pipelines, or read types, it
  counts as one sample while duplicated rows still affect aggregate metrics.
  Validate unique `sample_id` across inputs, or make the grouping key include
  pipeline/read_type/run metadata when those dimensions are intended to coexist.

- [ ] **Validate `--min-sample-fraction` range before SQL interpolation** —
  `src/cli.rs` documents the value as `0.0-1.0`, but no parser or runtime check
  enforces that. Values such as `-1`, `2`, `NaN`, or infinity can reach the SQL
  in `src/cohort.rs`, producing misleading filters or confusing DuckDB errors.
  Add a Clap value parser or explicit range check.

- [ ] **Improve `--on-target-only` errors when `on_target` is absent** —
  both `src/cohort.rs` and `src/qc.rs` append `WHERE on_target = true` without
  first checking that the column exists. The CLI says the option requires the
  column, but older/no-target Parquets fail with a generic DuckDB binder/create
  view error. Detect the missing column and fail with a clear message such as
  "rerun collect with --targets or omit --on-target-only."

---

## 2026-06-08 — Tumor/Normal Explorer and `geac annotate-normal`

Scope:
- `src/normal.rs`
- `app/explorer/tabs/tumor_normal.py`
- `src/writer/parquet_normal.rs`
- `tests/integration.rs` annotate-normal coverage

Verification:
- `cargo test -q annotate_normal`
- Result: 3 annotate-normal integration tests passed.

Findings:

- [ ] **Tumor/Normal tab can duplicate or mix normal evidence rows** —
  `app/explorer/tabs/tumor_normal.py` builds `ne_anchor` directly from
  `normal_evidence` anchor rows and joins it to tumor loci without enforcing one
  row per `(tumor_sample_id, chrom, pos, tumor_alt_allele)`. If duplicate
  normal-evidence files, reruns, or multiple normals for one tumor are merged,
  tumor rows can be duplicated. `ne_match` also sums matching alt counts without
  carrying `normal_sample_id`, while depth comes from anchor rows separately,
  which can create mixed normal VAFs. Consider enforcing uniqueness during
  merge/annotation, or aggregate by
  `(tumor_sample_id, chrom, pos, tumor_alt_allele, normal_sample_id)` and expose
  normal-sample selection in the tab.

- [ ] **`annotate-normal` should reject multi-sample tumor Parquet inputs** —
  `src/normal.rs` reads `ANY_VALUE(sample_id)` from the tumor Parquet and writes
  that same `tumor_sample_id` onto every output row. If a user accidentally
  passes a multi-sample Parquet/cohort output instead of a single-sample
  `geac collect` Parquet, normal evidence for all tumor loci is mislabeled as
  one arbitrary sample. Validate `COUNT(DISTINCT sample_id) = 1` and fail with a
  clear message when the input is not a single tumor sample.

- [ ] **`annotate-normal` treats all `bam.fetch()` errors as zero coverage** —
  `src/normal.rs` maps any fetch error to `found = false`, which emits
  `normal_depth = 0`. That is reasonable for a missing contig, but can hide
  reference/index/region problems as "No normal coverage" in downstream review.
  Distinguish missing-contig errors from other fetch failures, and report a
  warning/count of skipped loci if zero-depth fallback is retained.

- [ ] **Normal pileup de-duplication is iteration-order dependent** —
  `src/normal.rs::tally_position()` groups by query name and keeps the first
  observed base. For overlapping read pairs with discordant bases or different
  base qualities, normal alt support depends on pileup iteration order. The main
  collect path tracks overlap agreement/disagreement more explicitly; normal
  evidence should choose a deterministic higher-quality base or expose pair
  agreement state so low-VAF normal calls are not arbitrary.
