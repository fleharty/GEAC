#!/usr/bin/env bash
# make_manifest.sh — generate a GEAC manifest TSV from a directory of BAM files.
#
# Reads the SM tag from each BAM's @RG header line and writes a manifest with
# columns: sample_id, bam_path, bai_path
#
# Usage:
#   bash scripts/make_manifest.sh -b /path/to/bams [-o manifest.tsv] [-p PATTERN]
#
# Options:
#   -b BAM_DIR   Directory containing BAM files (required)
#   -o OUTPUT    Output manifest path (default: BAM_DIR/manifest.tsv)
#   -p PATTERN   Glob pattern for BAMs (default: *.bam)

set -euo pipefail

usage() {
    sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
}

BAM_DIR=""
OUTPUT=""
PATTERN="*.bam"

while getopts "b:o:p:h" opt; do
    case $opt in
        b) BAM_DIR="$OPTARG"  ;;
        o) OUTPUT="$OPTARG"   ;;
        p) PATTERN="$OPTARG"  ;;
        h) usage              ;;
        *) usage              ;;
    esac
done

[[ -n "$BAM_DIR" ]] || { echo "Error: -b BAM_DIR is required"; usage; }
[[ -d "$BAM_DIR" ]] || { echo "Error: directory not found: $BAM_DIR"; exit 1; }

if [[ -z "$OUTPUT" ]]; then
    OUTPUT="$BAM_DIR/manifest.tsv"
fi

# Require samtools
if ! command -v samtools &>/dev/null; then
    echo "Error: samtools not found. Run 'conda activate geac' or install samtools."
    exit 1
fi

shopt -s nullglob
bams=("$BAM_DIR"/$PATTERN)
shopt -u nullglob

if (( ${#bams[@]} == 0 )); then
    echo "Error: no files matching '$PATTERN' found in $BAM_DIR"
    exit 1
fi

echo "Found ${#bams[@]} BAM(s)"
echo "Writing manifest to $OUTPUT"
echo ""

printf "sample_id\tbam_path\tbai_path\n" > "$OUTPUT"

failed=0
for bam in "${bams[@]}"; do
    bam="$(cd "$(dirname "$bam")" && pwd)/$(basename "$bam")"

    # Extract SM tag from the first @RG line that has one
    sample_id=$(samtools view -H "$bam" \
        | awk '/^@RG/ { for(i=1;i<=NF;i++) if ($i ~ /^SM:/) { sub(/^SM:/,"",$i); print $i; exit } }')

    if [[ -z "$sample_id" ]]; then
        echo "[warn] $(basename "$bam") — no SM tag found, skipping"
        (( failed++ )) || true
        continue
    fi

    # Look for index alongside the BAM (.bai or .bam.bai)
    bai=""
    if   [[ -f "${bam}.bai"                           ]]; then bai="${bam}.bai"
    elif [[ -f "${bam%.bam}.bai"                      ]]; then bai="${bam%.bam}.bai"
    fi

    printf "%s\t%s\t%s\n" "$sample_id" "$bam" "$bai" >> "$OUTPUT"
    echo "[ok]  $sample_id  →  $(basename "$bam")"
done

echo ""
echo "Manifest written to $OUTPUT"
if (( failed > 0 )); then
    echo "Warning: $failed BAM(s) skipped (no SM tag)."
fi
