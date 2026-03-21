#!/bin/bash
# Launch the GEAC Explorer in a browser.
#
# Usage:
#   ./run_explorer.sh [/path/to/data/directory]
#
# The data directory is mounted at /data inside the container.
# Place your cohort.duckdb and geac.toml there.
# Defaults to the current directory if no argument is given.
#
# Requires Docker or Podman.

set -euo pipefail

DATA_DIR="${1:-.}"
DATA_DIR="$(cd "${DATA_DIR}" && pwd)"   # make absolute

VERSION="$(cat "$(dirname "$0")/VERSION")"
IMAGE="ghcr.io/fleharty/geac:${VERSION}"

if ! command -v podman &>/dev/null; then
    echo "Error: podman not found on PATH." >&2
    exit 1
fi

echo "Starting GEAC Explorer v${VERSION}..."
echo "Data directory: ${DATA_DIR}"
echo "Open http://localhost:8501 in your browser."
echo "Press Ctrl-C to stop."
echo

exec podman run --rm \
    -p 8501:8501 \
    -v "${DATA_DIR}:/data" \
    "${IMAGE}" \
    explorer -- --server.address 0.0.0.0 --server.headless true
