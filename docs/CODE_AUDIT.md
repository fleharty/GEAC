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
