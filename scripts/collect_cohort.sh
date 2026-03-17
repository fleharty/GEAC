#!/usr/bin/env bash
# collect_cohort.sh — run geac collect on every *.bqsr.bam in a directory,
# writing one Parquet file per sample to the same directory.
#
# Usage:
#   bash scripts/collect_cohort.sh -b /path/to/bams -r /path/to/ref.fa
#
# Options:
#   -b BAM_DIR    Directory containing {sample}.bqsr.bam files (required)
#   -r REFERENCE  Path to reference FASTA (required)
#   -j JOBS       Number of samples to process in parallel (default: 4)
#   -t THREADS    CPU threads per geac collect job (default: 4)
#   -f            Force re-run even if the output Parquet already exists

set -euo pipefail

usage() {
    sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
}

BAM_DIR=""
REFERENCE=""
JOBS=4
THREADS=4
FORCE=0

while getopts "b:r:j:t:fh" opt; do
    case $opt in
        b) BAM_DIR="$OPTARG"    ;;
        r) REFERENCE="$OPTARG" ;;
        j) JOBS="$OPTARG"      ;;
        t) THREADS="$OPTARG"   ;;
        f) FORCE=1             ;;
        h) usage               ;;
        *) usage               ;;
    esac
done

[[ -n "$BAM_DIR"   ]] || { echo "Error: -b BAM_DIR is required";   usage; }
[[ -n "$REFERENCE" ]] || { echo "Error: -r REFERENCE is required"; usage; }
[[ -d "$BAM_DIR"   ]] || { echo "Error: directory not found: $BAM_DIR";   exit 1; }
[[ -f "$REFERENCE" ]] || { echo "Error: reference not found: $REFERENCE"; exit 1; }

# Locate the geac binary: prefer the release build in this repo, then PATH.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
if [[ -x "$REPO_ROOT/target/release/geac" ]]; then
    GEAC="$REPO_ROOT/target/release/geac"
elif command -v geac &>/dev/null; then
    GEAC="geac"
else
    echo "Error: geac binary not found. Run 'cargo build --release' or add geac to PATH."
    exit 1
fi

shopt -s nullglob
bams=("$BAM_DIR"/*.bqsr.bam)
shopt -u nullglob

if (( ${#bams[@]} == 0 )); then
    echo "Error: no *.bqsr.bam files found in $BAM_DIR"
    exit 1
fi

echo "Found ${#bams[@]} BAM(s) in $BAM_DIR"
echo "Geac binary : $GEAC"
echo "Reference   : $REFERENCE"
echo "Parallel    : $JOBS jobs x $THREADS threads"
echo ""

pids=()
statuses=()

run_sample() {
    local bam="$1"
    local sample
    sample="$(basename "$bam" .bqsr.bam)"
    local output="$BAM_DIR/${sample}.parquet"
    local variants_tsv="$BAM_DIR/${sample}.bqsr.annotated_variants.txt"

    if [[ -f "$output" && "$FORCE" -eq 0 ]]; then
        echo "[skip] $sample — output already exists ($output)"
        return 0
    fi

    local cmd=("$GEAC" collect
        --input     "$bam"
        --reference "$REFERENCE"
        --output    "$output"
        --read-type duplex
        --pipeline  fgbio
        --threads   "$THREADS"
    )

    if [[ -f "$variants_tsv" ]]; then
        cmd+=(--variants-tsv "$variants_tsv")
    else
        echo "[warn] $sample — no variants TSV found, skipping --variants-tsv"
    fi

    echo "[run]  $sample"
    if "${cmd[@]}"; then
        echo "[done] $sample"
    else
        echo "[fail] $sample (exit $?)" >&2
        return 1
    fi
}

export -f run_sample
export BAM_DIR REFERENCE THREADS FORCE GEAC

failed=0

for bam in "${bams[@]}"; do
    run_sample "$bam" &
    pids+=($!)

    if (( ${#pids[@]} >= JOBS )); then
        if ! wait "${pids[0]}"; then
            (( failed++ )) || true
        fi
        pids=("${pids[@]:1}")
    fi
done

for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
        (( failed++ )) || true
    fi
done

echo ""
if (( failed > 0 )); then
    echo "Finished with $failed failure(s)."
    exit 1
else
    echo "All samples completed successfully."
fi
