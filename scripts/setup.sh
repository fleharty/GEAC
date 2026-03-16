#!/usr/bin/env bash
# Setup script for GEAC development environment.
# Run once after cloning the repository.
#
# Usage:
#   bash scripts/setup.sh

set -euo pipefail

# ── 1. Conda environment ──────────────────────────────────────────────────────
if ! command -v conda &>/dev/null; then
    echo "ERROR: conda not found. Install Miniforge or Anaconda first." >&2
    exit 1
fi

echo "==> Creating conda environment 'geac'..."
conda env create -f environment.yml || conda env update -f environment.yml

# Resolve the conda prefix for the geac env
CONDA_ENV_PREFIX=$(conda run -n geac conda info --base 2>/dev/null || true)
CONDA_ENV_PREFIX=$(conda env list | awk '/^geac /{print $NF}')

if [[ -z "$CONDA_ENV_PREFIX" ]]; then
    echo "ERROR: could not determine conda prefix for 'geac' environment." >&2
    exit 1
fi

echo "==> Conda prefix: $CONDA_ENV_PREFIX"

# ── 2. activate.d hook ────────────────────────────────────────────────────────
# This script runs automatically on `conda activate geac` and sets the
# environment variables that rust-htslib needs to find htslib.
ACTIVATE_D="$CONDA_ENV_PREFIX/etc/conda/activate.d"
DEACTIVATE_D="$CONDA_ENV_PREFIX/etc/conda/deactivate.d"
mkdir -p "$ACTIVATE_D" "$DEACTIVATE_D"

cat > "$ACTIVATE_D/geac_env.sh" <<'EOF'
#!/usr/bin/env bash
# Set build environment for rust-htslib -> htslib linkage
export LIBRARY_PATH="${CONDA_PREFIX}/lib${LIBRARY_PATH:+:${LIBRARY_PATH}}"
export CPATH="${CONDA_PREFIX}/include${CPATH:+:${CPATH}}"
export PKG_CONFIG_PATH="${CONDA_PREFIX}/lib/pkgconfig${PKG_CONFIG_PATH:+:${PKG_CONFIG_PATH}}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
# Bake conda lib rpath into binaries at link time (required on macOS)
export RUSTFLAGS="-C link-arg=-Wl,-rpath,${CONDA_PREFIX}/lib${RUSTFLAGS:+ ${RUSTFLAGS}}"
EOF

cat > "$DEACTIVATE_D/geac_env.sh" <<'EOF'
#!/usr/bin/env bash
# Restore environment variables on deactivate (best-effort)
unset LIBRARY_PATH
unset CPATH
unset PKG_CONFIG_PATH
unset LD_LIBRARY_PATH
unset RUSTFLAGS
EOF

echo "==> Installed conda activate hooks."

# ── 3. Rust via rustup ───────────────────────────────────────────────────────
if command -v rustup &>/dev/null; then
    echo "==> rustup already installed, updating toolchain..."
    rustup update stable
else
    echo "==> Installing Rust via rustup..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
    # shellcheck source=/dev/null
    source "$HOME/.cargo/env"
fi

echo ""
echo "Setup complete."
echo ""
echo "Next steps:"
echo "  conda activate geac"
echo "  cargo build --release"
echo "  ./target/release/geac --help"
