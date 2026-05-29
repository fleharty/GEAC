#!/usr/bin/env python3
"""Generate a GEAC manifest TSV by scanning pipeline subdirectories.

Usage:
    python scripts/make_manifest.py <root_dir> [--output manifest.tsv]

Expects <root_dir> to contain subdirectories (e.g. fgbio/, dragen/), each
holding files named:
    {sample_name}.bam
    {sample_name}.bai
    {sample_name}.annotated_filtered_variants.txt

The subdirectory name becomes the `pipeline` column value.
"""

import argparse
import csv
import sys
from pathlib import Path


def scan_pipeline_dir(pipeline_dir: Path) -> list[dict]:
    rows = []
    bams = sorted(pipeline_dir.glob("*.bam"))
    for bam in bams:
        sample = bam.stem
        bai = pipeline_dir / f"{sample}.bai"
        variants = pipeline_dir / f"{sample}.annotated_filtered_variants.txt"
        rows.append({
            "collaborator_sample_id": sample,
            "duplex_output_bam": str(bam),
            "duplex_output_bam_index": str(bai) if bai.exists() else "",
            "final_annotated_variants": str(variants) if variants.exists() else "",
            "pipeline": pipeline_dir.name,
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description="Generate GEAC manifest TSV")
    parser.add_argument("root_dir", help="Directory containing pipeline subdirectories")
    parser.add_argument("--output", default="manifest.tsv", help="Output TSV path (default: manifest.tsv)")
    parser.add_argument("--pipelines", nargs="+", default=["fgbio", "dragen"],
                        help="Subdirectory names to scan (default: fgbio dragen)")
    args = parser.parse_args()

    root = Path(args.root_dir).resolve()
    rows = []
    for pipeline_name in args.pipelines:
        pipeline_dir = root / pipeline_name
        if not pipeline_dir.is_dir():
            print(f"Warning: {pipeline_dir} not found, skipping", file=sys.stderr)
            continue
        found = scan_pipeline_dir(pipeline_dir)
        print(f"  {pipeline_name}: {len(found)} sample(s)", file=sys.stderr)
        rows.extend(found)

    if not rows:
        print("No samples found. Check your root directory and pipeline names.", file=sys.stderr)
        sys.exit(1)

    fieldnames = [
        "collaborator_sample_id",
        "duplex_output_bam",
        "duplex_output_bam_index",
        "final_annotated_variants",
        "pipeline",
    ]
    output_path = Path(args.output)
    with output_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {output_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
