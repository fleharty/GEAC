# Development Challenges

A running log of non-obvious problems encountered during GEAC development, how they
were diagnosed, and what the fix was. Updated as new issues arise.

---

## Rust / Core Tool

### Coverage memory growth on large genomes
**Symptom:** `geac coverage` on a whole-genome BAM would accumulate all
`CoverageRecord` objects in a `Vec` before writing, consuming multiple GB of RAM.
**Root cause:** The original design collected all records, sorted them, then wrote
the Parquet file in one shot.
**Fix:** Replaced with a streaming `CoverageWriter` that buffers 100k records at a
time and flushes to disk as it goes (`src/writer/parquet_coverage.rs`). The sort was
also dropped — pileup order is already genomic order.
**Lesson:** For genome-scale outputs, never buffer the full result set in memory.

### Rust import path: `ReadType` not in scope in submodule
**Symptom:** Compiler error when `BinAccumulator` in `src/coverage/mod.rs` tried to
reference `crate::cli::ReadType`.
**Root cause:** `ReadType` is defined in `crate::record`, not `crate::cli`. The CLI
module re-exports it for clap parsing but the canonical definition lives in `record`.
**Fix:** Changed import to `use crate::record::ReadType`.

### Edit tool conflict: ambiguous surrounding context
**Symptom:** An edit to add `adaptive_depth_threshold` to `CoverageArgs` in
`src/cli.rs` failed because two structs had identical surrounding lines, making
the match non-unique.
**Fix:** Extended the `old_string` context to include a nearby unique field
(`min_depth`) so the replacement target was unambiguous.

---

## Python / Streamlit Explorer

### Depth distribution overrepresented low-depth bins
**Symptom:** The depth distribution histogram in the Coverage Explorer appeared
to have far more low-coverage positions than expected.
**Root cause:** With `--bin-size > 1`, each row in the Parquet represents multiple
genomic positions. Counting rows (`COUNT(*)`) treated a 100bp bin the same as a
1bp position, massively underweighting high-depth bins.
**Fix:** Added a `bin_n` column to `CoverageRecord` tracking positions per bin.
All locus-count queries now use `SUM(bin_n)` instead of `COUNT(*)`.

### Gene bar chart click-to-drill-down not working
**Symptom:** Clicking a bar in the "Affected loci per gene" chart appeared to
trigger a Streamlit rerun, but the detail table never appeared.
**Root cause (first attempt):** Used `event.selection.point` — but
`event.selection` is an `AttributeDictionary`, not an object with a `.point`
attribute. This raised `AttributeError`.
**Partial fix:** Changed to iterate `event.selection.values()` to find a non-empty
list. The rerun now completes without error but the table still doesn't render
reliably.
**Status:** Unresolved — tracked as a known issue for 0.4.0. The exact structure
of Altair point-selection events in Streamlit's `on_select="rerun"` API is
unclear; needs further investigation.

---

## IGV.js Integration

### "Access Unauthorized" on initial BAM load
**Symptom:** IGV.js displayed "Access Unauthorized" immediately when loading a
`gs://` BAM URL.
**Root cause:** The ADC access token was fetched server-side via `google.auth.default`
but not passed to the IGV.js browser instance.
**Fix:** Injected the token three ways in the generated HTML: globally via
`igv.setOauthToken()`, in the browser-level `oauthToken` option, and per-track
`oauthToken` field. Host pattern changed from `"*.googleapis.com"` to
`"storage.googleapis.com"` to match actual GCS request hostnames.

### "Access Unauthorized" on zoom / range requests
**Symptom:** Initial load worked, but zooming in triggered new HTTP range requests
that returned 401.
**Root cause:** Token was set globally but subsequent range requests did not pick
it up. Only the per-track `oauthToken` field reliably covers all requests.
**Fix:** Ensured `oauthToken` is set on every track dict in addition to the global
call.

### IGV.js page freeze on load
**Symptom:** Entire browser tab froze when switching to the IGV tab.
**Root cause (1):** CRAM files require a reference FASTA. Without explicit
`fastaURL`/`indexURL` in the browser config, IGV.js tried to fetch the default
hg38 reference and hung.
**Fix:** Added explicit hg38 reference URLs to the browser options object.
**Root cause (2):** Streamlit rerenders the entire component on every interaction,
reinitializing the IGV browser and triggering repeated heavy network fetches.
**Fix:** Gated the `st.components.v1.html` call behind a "Load IGV" button using
`st.session_state["_igv_loaded"]`, so the browser only initializes once per session.

---

## Docker / Deployment

### Podman: collaborator couldn't find database file
**Symptom:** Collaborator entered `cohort.duckdb` in the Explorer UI but got
"database does not exist."
**Root cause:** Podman does not mount the current directory by default. The
`cohort.duckdb` file existed on the host but was not visible inside the container.
**Fix:** Required explicit `-v $(pwd):/data` volume mount and using `/data/cohort.duckdb`
as the path inside the container.

### Docker image was private on ghcr.io
**Symptom:** Collaborator running `podman pull ghcr.io/fleharty/geac` received
"permission denied."
**Root cause:** GitHub Container Registry defaults to private for new packages.
**Fix:** Made the package public via GitHub → Packages → Package Settings →
"Change visibility."

---

## CI / GitHub Actions

### `macos-13` runner no longer available
**Symptom:** The `native-binaries` job for `macos-x86_64` failed with
`"The configuration 'macos-13-us-default' is not supported"`, causing the entire
workflow run to be cancelled.
**Root cause:** GitHub deprecated the `macos-13` runner pool.
**Fix:** Dropped `macos-x86_64` (and all Linux native builds) from CI entirely —
only `macos-arm64` is needed for the Homebrew tap.

### Homebrew formula: `#{version}` empty inside `resource` block
**Symptom:** `brew install` attempted to fetch `v.tar.gz` (version missing from URL).
**Root cause:** In Homebrew's DSL, `#{version}` inside a `resource` block refers to the *resource's* own version, not the enclosing formula's version. Since the resource has no version set, it stringifies to empty string.
**Fix:** Changed the resource URL template to use the `FORMULA_VERSION` sed placeholder (`vFORMULA_VERSION.tar.gz`) so CI substitutes the correct version directly, rather than relying on Ruby interpolation at install time.

### Homebrew tap: `HOMEBREW_TAP_TOKEN` secret not set
**Symptom:** CI log shows `TAP_TOKEN: ` (blank); `git push` to `homebrew-geac` would fail silently or with auth error.
**Root cause:** The `HOMEBREW_TAP_TOKEN` secret must be manually created in GitHub repo settings — it is not auto-provisioned like `GITHUB_TOKEN`.
**Fix:** Create a classic PAT with `repo` scope for `fleharty/homebrew-geac`, then add it as a repository secret named `HOMEBREW_TAP_TOKEN` in the GEAC repo settings.

### Homebrew tap: `curl` 404 on source tarball
**Symptom:** `curl: (22) The requested URL returned error: 404` when downloading the GitHub source archive for SHA256 computation.
**Root cause:** Using plain `curl` without authentication on a repo that requires it returns 404 instead of 401.
**Fix:** Replaced `curl -fsSL` with `gh api repos/.../tarball/...`, which uses `GH_TOKEN` automatically.

### `gh release upload` fails with "release not found"
**Symptom:** The `Package and upload binary` step in `native-binaries` exited with
`release not found` when using `gh release upload`.
**Root cause:** `gh release upload` requires the GitHub Release object to already
exist. The `docker` job runs in parallel and does not create a release; nothing
creates the release before the binary job runs.
**Fix:** Reverted to `softprops/action-gh-release@v2`, which creates the release
automatically if it doesn't exist.

### Multi-platform Docker: arm64 runner unavailable → whole run cancelled
**Symptom:** v0.3.7 first release attempt was cancelled mid-run.
**Root cause:** The original Docker workflow used a matrix (amd64 + arm64 via
`ubuntu-22.04-arm`) with a digest-merge approach. Combined with the `macos-13`
failure above, the run was cancelled before the Homebrew tap job ran.
**Fix:** Simplified to a single `docker` job building only `linux/amd64`, removing
the matrix, digest upload/download, and `imagetools create` merge step entirely.
