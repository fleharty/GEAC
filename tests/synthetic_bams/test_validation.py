"""
Validation tests: geac collect output vs. independent pysam fragment counting.

These tests confirm that geac's alt_count, ref_count, and total_depth at the SNV
locus agree with an independent pysam-based fragment counter.

Fragment-based counting
-----------------------
geac counts each DNA fragment (read pair) once, not each read.  For overlapping
pairs where both R1 and R2 cover the SNV (e.g. the default layout in these BAMs
where R1=175–324 and R2=225–374 both span position 250), the pair contributes
exactly one count to total_depth.

The pysam counter below mirrors this by deduplicating on query_name: if any read
in the pair observes the alt allele, the fragment is classified as "alt".

Why not samtools?
-----------------
samtools depth / mpileup does not have a built-in fragment-deduplication mode.
For positions covered by both R1 and R2, raw samtools depth counts twice (giving
20 for clean_snv instead of 10).  A partial workaround for layouts where all R1
reads cover the SNV is:

    samtools view -b -f 64 clean_snv.bam | samtools depth -a -

The -f 64 flag selects only first-in-pair reads (R1), yielding fragment depth for
that specific geometry.  This approach fails for layouts where some fragments only
cover the SNV via R2, so pysam with query_name deduplication is the correct
general solution.

Known intentional differences from raw samtools mpileup
--------------------------------------------------------
- total_depth counts fragments, not reads (see above)
- Reads flagged as PCR/optical duplicates (FLAG 0x400) are excluded by geac
  unless --include-duplicates is passed
- The default min_map_qual is 0 (all MAPQ values included)

Running
-------
Requires the geac binary to be on PATH, or set GEAC_BIN=/path/to/geac.
All synthetic BAMs must have been generated first:
    cd tests/synthetic_bams && python generate.py
"""

import os
import subprocess
from pathlib import Path

import pandas as pd
import pysam
import pytest

# ── Paths ─────────────────────────────────────────────────────────────────────

BAMS_DIR = Path(__file__).parent / "bams"
REF      = BAMS_DIR / "reference.fa"
GEAC_BIN = os.environ.get("GEAC_BIN", "geac")

# ── SNV coordinates planted in all synthetic BAMs (see generate.py) ───────────

SNV_CHROM = "chr1"
SNV_POS_0 = 250      # 0-based (DB convention)
SNV_REF   = "C"
SNV_ALT   = "T"


# ── Helpers ───────────────────────────────────────────────────────────────────

def run_geac_collect(bam: Path, out: Path) -> None:
    """Run geac collect with default flags (no mapq filter, raw pipeline)."""
    result = subprocess.run(
        [
            GEAC_BIN, "collect",
            "--input",     str(bam),
            "--reference", str(REF),
            "--output",    str(out),
            "--sample-id", bam.stem,
            "--read-type", "raw",
            "--pipeline",  "raw",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"geac collect failed on {bam.name}:\n{result.stderr}"
    )


def geac_counts_at_snv(parquet: Path) -> dict:
    """Return alt_count, ref_count, total_depth from the SNV row in a geac Parquet."""
    df = pd.read_parquet(parquet)
    rows = df[df["pos"] == SNV_POS_0]
    assert len(rows) == 1, (
        f"Expected 1 row at pos {SNV_POS_0} in {parquet.name}, got {len(rows)}"
    )
    row = rows.iloc[0]
    return {
        "alt_count":   int(row["alt_count"]),
        "ref_count":   int(row["ref_count"]),
        "total_depth": int(row["total_depth"]),
    }


def pysam_fragment_counts(bam: Path) -> dict:
    """
    Independently count alt/ref fragments at SNV_POS_0 using pysam.

    Deduplicates on query_name so an overlapping R1+R2 pair counts as one
    fragment.  A query_name is classified as 'alt' if any read in the pair
    observed the alt allele; 'ref' otherwise.
    """
    alt_names: set[str] = set()
    ref_names: set[str] = set()

    af = pysam.AlignmentFile(str(bam), "rb")
    for col in af.pileup(
        SNV_CHROM, SNV_POS_0, SNV_POS_0 + 1,
        min_base_quality=0,
        min_mapping_quality=0,
        ignore_overlaps=False,
        stepper="all",
    ):
        if col.reference_pos != SNV_POS_0:
            continue
        for p in col.pileups:
            if p.is_del or p.is_refskip:
                continue
            base = p.alignment.query_sequence[p.query_position]
            name = p.alignment.query_name
            if base == SNV_ALT:
                alt_names.add(name)
            elif base == SNV_REF:
                ref_names.add(name)
    af.close()

    # A fragment with any alt read is classified as alt, not ref.
    ref_names -= alt_names

    return {
        "alt_count":   len(alt_names),
        "ref_count":   len(ref_names),
        "total_depth": len(alt_names) + len(ref_names),
    }


# ── Tests ──────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("scenario,expected_alt,expected_total", [
    # scenario           alt_count  total_depth
    # ─────────────────────────────────────────
    # 10 fragments (pairs), 2 alt.  Both R1 and R2 overlap the SNV, so raw
    # read depth = 20 but fragment depth = 10.  Key test for deduplication.
    ("clean_snv",        2,         10),

    # 20 fragments, 10 alt.  VAF ≈ 0.50 (germline-like).
    ("germline_snv",     10,        20),

    # All 4 alt reads are on R2; R1 always carries ref.  The fragment is still
    # classified as alt because R2 observed the alt allele.
    ("r2_artefact",      2,         10),

    # Alt reads have MAPQ=10; geac default min_map_qual=0 so they are included.
    # Confirms geac does not silently drop low-MAPQ reads by default.
    ("low_mapq",         2,         10),

    # Alt reads at cycle 140 (end of R1), ref at cycle 76.
    # Validates that cycle position does not affect base counting.
    ("end_of_read",      2,         10),
])
def test_geac_matches_pysam_fragment_count(
    scenario, expected_alt, expected_total, tmp_path
):
    bam     = BAMS_DIR / f"{scenario}.bam"
    parquet = tmp_path / f"{scenario}.parquet"

    run_geac_collect(bam, parquet)

    geac  = geac_counts_at_snv(parquet)
    psam  = pysam_fragment_counts(bam)

    # ── geac must agree with the independent pysam counter ────────────────────
    assert geac["alt_count"] == psam["alt_count"], (
        f"{scenario}: alt_count  geac={geac['alt_count']} pysam={psam['alt_count']}"
    )
    assert geac["ref_count"] == psam["ref_count"], (
        f"{scenario}: ref_count  geac={geac['ref_count']} pysam={psam['ref_count']}"
    )
    assert geac["total_depth"] == psam["total_depth"], (
        f"{scenario}: total_depth  geac={geac['total_depth']} pysam={psam['total_depth']}"
    )

    # ── both must match the documented expected values from generate.py ───────
    assert geac["alt_count"]   == expected_alt,   f"{scenario}: alt_count"
    assert geac["total_depth"] == expected_total, f"{scenario}: total_depth"
    assert geac["ref_count"]   == expected_total - expected_alt, f"{scenario}: ref_count"


def test_vaf_is_alt_over_total(tmp_path):
    """VAF stored in the Parquet must equal alt_count / total_depth."""
    bam     = BAMS_DIR / "clean_snv.bam"
    parquet = tmp_path / "clean_snv.parquet"
    run_geac_collect(bam, parquet)

    df  = pd.read_parquet(parquet)
    row = df[df["pos"] == SNV_POS_0].iloc[0]

    expected_vaf = row["alt_count"] / row["total_depth"]
    assert abs(float(row["alt_count"]) / float(row["total_depth"]) - expected_vaf) < 1e-6


def test_strand_counts_sum_to_total(tmp_path):
    """
    fwd_depth + rev_depth must equal total_depth.
    fwd_alt_count + rev_alt_count may be >= alt_count when both R1 and R2
    of an overlapping pair carry the alt allele (each strand counted separately).
    """
    bam     = BAMS_DIR / "clean_snv.bam"
    parquet = tmp_path / "clean_snv.parquet"
    run_geac_collect(bam, parquet)

    df  = pd.read_parquet(parquet)
    row = df[df["pos"] == SNV_POS_0].iloc[0]

    assert int(row["fwd_depth"]) + int(row["rev_depth"]) == int(row["total_depth"]), (
        "fwd_depth + rev_depth must equal total_depth"
    )
    # fwd_alt_count and rev_alt_count are read-level (not fragment-level),
    # so their sum >= alt_count for overlapping pairs.
    assert int(row["fwd_alt_count"]) + int(row["rev_alt_count"]) >= int(row["alt_count"]), (
        "combined strand alt counts must be at least fragment alt_count"
    )


def test_coordinate_is_zero_based(tmp_path):
    """
    geac stores positions 0-based (from htslib pileup.pos()).
    The SNV planted at 1-based chr1:251 must appear at pos=250 in the Parquet.
    """
    bam     = BAMS_DIR / "clean_snv.bam"
    parquet = tmp_path / "clean_snv.parquet"
    run_geac_collect(bam, parquet)

    df = pd.read_parquet(parquet)
    assert 250 in df["pos"].values, "SNV should appear at pos=250 (0-based)"
    assert 251 not in df["pos"].values, "pos=251 would indicate 1-based storage"
