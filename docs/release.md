# Release Checklist

GEAC releases are tag-driven. Pushing `v*.*.*` starts the GitHub Actions release
workflow.

## Before Tagging

1. Pick the next semantic version, for example `0.4.40`.
2. Update the version in:
   - `Cargo.toml`
   - `Cargo.lock`
   - `VERSION`
   - `app/explorer/schema.py`
3. Run:

```bash
cargo check -q
cargo test -q
cargo build --release
git diff --check
target/release/geac --version
target/release/geac --help
```

`cargo fmt --check` is useful, but the current release workflow does not enforce it.
If the repo has broad formatting drift, handle that in a separate formatting-only
change rather than mixing it into a release commit.

## Commit And Tag

```bash
git add Cargo.toml Cargo.lock VERSION app/explorer/schema.py
git commit -m "release 0.X.Y"
git push origin main
git tag v0.X.Y
git push origin v0.X.Y
```

## Automated Release Workflow

The tag workflow:

1. Builds and pushes the `linux/amd64` Docker image to GHCR.
2. Builds a native `macos-arm64` binary.
3. Attaches the macOS binary archive to the GitHub release.
4. Updates the Homebrew tap formula.

Published artifacts:

```text
ghcr.io/fleharty/geac:<version>
ghcr.io/fleharty/geac:latest
geac-macos-arm64.tar.gz
```
