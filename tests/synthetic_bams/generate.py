#!/usr/bin/env python3
"""
Generate synthetic BAMs for GEAC end-to-end testing.

Creates a minimal synthetic reference (500 bp chr1) and a set of BAM files,
each representing a specific property at a known locus.  All BAMs can be
loaded in IGV alongside reference.fa for visual inspection.

SNV locus: chr1:251 (1-based) = position 250 (0-based)
  REF = C   ALT = T

Read layout (default scenario):
  Read length : 150 bp
  Insert size : 200 bp
  R1 start    : 175  → SNV at R1 qpos 75  → cycle 76
  R2 start    : 225  → SNV at R2 qpos 25  → cycle 26

Usage:
    python generate.py [--outdir <dir>]

Outputs:
    reference.fa / reference.fa.fai   — synthetic 500 bp reference
    clean_snv.bam                     — 20 reads, VAF=0.20, no artefacts
    r2_artefact.bam                   — all 4 alt reads on R2 only
    end_of_read.bam                   — alt reads at cycle 140+ (near end of R1)
    low_mapq.bam                      — alt reads have MAPQ=10; ref reads MAPQ=60
    short_insert.bam                  — alt reads have insert=80 bp; ref reads insert=200 bp
    singleton_family.bam              — alt reads have cD=1; ref reads have cD=3
    germline_snv.bam                  — VAF≈0.50, 40 reads, no artefact properties
"""

import argparse
import os
import random
import struct
from pathlib import Path

import pysam

# ── Reference ─────────────────────────────────────────────────────────────────

CHROM       = "chr1"
REF_LEN     = 500
SNV_POS     = 250   # 0-based; 1-based = 251
REF_BASE    = "C"
ALT_BASE    = "T"

def _make_ref_seq() -> str:
    """Deterministic, repeat-free 500 bp sequence with C at position 250."""
    # Use a simple LCG so the sequence is reproducible without numpy/random dependency
    bases = "ACGT"
    h = 0xDEADBEEF
    seq = []
    for i in range(REF_LEN):
        h = (h * 1664525 + 1013904223) & 0xFFFFFFFF
        seq.append(bases[h % 4])
    seq[SNV_POS] = REF_BASE
    return "".join(seq)

REF_SEQ = _make_ref_seq()


def write_reference(outdir: Path) -> Path:
    fa_path = outdir / "reference.fa"
    with open(fa_path, "w") as fh:
        fh.write(f">{CHROM}\n")
        # 60 bp per FASTA line
        for i in range(0, REF_LEN, 60):
            fh.write(REF_SEQ[i:i+60] + "\n")
    pysam.faidx(str(fa_path))
    return fa_path


# ── BAM helpers ───────────────────────────────────────────────────────────────

READ_LEN = 150

def _bam_header() -> pysam.AlignmentHeader:
    return pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": REF_LEN}],
    })


def _make_read(
    *,
    name:        str,
    r1:          bool,          # True = R1, False = R2
    ref_start:   int,           # 0-based alignment start
    seq:         str,           # full read sequence (length = READ_LEN)
    mapq:        int,
    insert_size: int,           # abs(TLEN); sign is applied per SAM convention
    tags:        list,          # list of (tag, value) tuples
    header:      pysam.AlignmentHeader,
) -> pysam.AlignedSegment:
    seg = pysam.AlignedSegment(header)
    seg.query_name = name
    seg.query_sequence = seq
    seg.flag = _flags(r1, insert_size)
    seg.reference_id = 0
    seg.reference_start = ref_start
    seg.mapping_quality = mapq
    seg.cigar = [(0, READ_LEN)]   # 150M
    seg.next_reference_id = 0
    # Mate start
    if r1:
        seg.next_reference_start = ref_start + insert_size - READ_LEN
        seg.template_length = insert_size
    else:
        seg.next_reference_start = ref_start - insert_size + READ_LEN
        seg.template_length = -insert_size
    seg.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)  # BQ=40
    for tag, val in tags:
        seg.set_tag(tag, val)
    return seg


def _flags(r1: bool, insert_size: int) -> int:
    """Return SAM flags for a properly paired read.

    R1 (forward): 0x1|0x2|0x40|0x20  paired, proper, first-in-pair, mate-reverse
    R2 (reverse): 0x1|0x2|0x80|0x10  paired, proper, second-in-pair, read-reverse
    """
    f = 0x1   # paired
    if insert_size > 0:
        f |= 0x2  # properly paired
    f |= (0x40 if r1 else 0x80)
    if r1:
        f |= 0x20   # mate (R2) is on reverse strand
    else:
        f |= 0x10   # R2 is on reverse strand (mate=R1 is forward, so 0x20 not set)
    return f


def _read_seq(is_alt: bool, r1: bool, r1_start: int, insert_size: int) -> str:
    """
    Build the read sequence from REF_SEQ, inserting ALT_BASE at SNV_POS if is_alt.
    R2 is reverse-complemented relative to the forward strand, but we store
    the original-strand sequence (pysam handles orientation via the flag).
    R2 start = r1_start + insert_size - READ_LEN.
    """
    if r1:
        start = r1_start
    else:
        start = r1_start + insert_size - READ_LEN

    region = list(REF_SEQ[start : start + READ_LEN])
    # Place the alt allele if this read covers the SNV position
    snv_qpos = SNV_POS - start
    if 0 <= snv_qpos < READ_LEN:
        region[snv_qpos] = ALT_BASE if is_alt else REF_BASE

    # BAM convention: query_sequence is stored left-to-right at the alignment
    # position regardless of strand.  For reverse-strand R2 (FLAG 0x10), the
    # sequence is the reverse complement of the original read — but since we
    # build it directly from the forward-strand reference slice, we do NOT
    # revcomp here.  pysam pileup reads query_sequence[qpos] directly.
    return "".join(region)


def _cycle_for(r1: bool, r1_start: int, insert_size: int) -> int:
    """1-based cycle at which the SNV appears in this read."""
    start = r1_start if r1 else r1_start + insert_size - READ_LEN
    snv_qpos = SNV_POS - start
    return snv_qpos + 1  # 1-based


# ── Fragment factory ───────────────────────────────────────────────────────────

def make_fragment(
    *,
    frag_id:     int,
    is_alt:      bool,
    r1_start:    int,
    insert_size: int,
    mapq:        int,
    tags:        list,
    header:      pysam.AlignmentHeader,
) -> tuple:
    """Return (r1_seg, r2_seg) sorted by position (r1 always left of r2 here)."""
    name = f"frag_{frag_id:04d}"
    r1 = _make_read(
        name=name, r1=True, ref_start=r1_start,
        seq=_read_seq(is_alt, True, r1_start, insert_size),
        mapq=mapq, insert_size=insert_size, tags=tags, header=header,
    )
    r2_start = r1_start + insert_size - READ_LEN
    r2 = _make_read(
        name=name, r1=False, ref_start=r2_start,
        seq=_read_seq(is_alt, False, r1_start, insert_size),
        mapq=mapq, insert_size=insert_size, tags=tags, header=header,
    )
    return r1, r2


# ── Write a BAM ───────────────────────────────────────────────────────────────

def write_bam(path: Path, fragments: list, header: pysam.AlignmentHeader) -> None:
    """Write fragments (list of (r1, r2) tuples) to a sorted, indexed BAM."""
    tmp = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=header) as bam:
        for r1, r2 in fragments:
            bam.write(r1)
            bam.write(r2)
    pysam.sort("-o", str(path), tmp)
    os.remove(tmp)
    pysam.index(str(path))


# ── Scenario builders ──────────────────────────────────────────────────────────
#
# Default layout: R1 starts at 175 → SNV at R1 cycle 76
#                 R2 starts at 225 → SNV at R2 cycle 26
# 20 reads (10 pairs), 4 alt reads (VAF = 0.2)
#
# Tags: no fgbio tags unless scenario requires them

DEFAULT_R1_START  = SNV_POS - 75   # 175; SNV at R1 qpos 75 → cycle 76
DEFAULT_INSERT    = 200
N_TOTAL_PAIRS     = 10             # 20 reads per scenario
N_ALT_PAIRS       = 2              # 4 alt reads (VAF ≈ 0.2)


def scenario_clean_snv(outdir: Path, header: pysam.AlignmentHeader) -> None:
    """
    Clean SNV with no artefact properties.
    20 reads (10 pairs), VAF=0.20, MAPQ=60, equal R1/R2, insert=200.
    SNV at R1 cycle 76 (qpos 75), R2 cycle 26 (qpos 25).
    Expected geac collect output: 1 locus at chr1:251, alt_count=2, total_depth=10.
    (geac counts one read per fragment; 10 fragments, 2 alt fragments → VAF=0.20)
    """
    frags = []
    for i in range(N_TOTAL_PAIRS):
        is_alt = i < N_ALT_PAIRS
        frags.append(make_fragment(
            frag_id=i, is_alt=is_alt, r1_start=DEFAULT_R1_START,
            insert_size=DEFAULT_INSERT, mapq=60, tags=[], header=header,
        ))
    write_bam(outdir / "clean_snv.bam", frags, header)
    _print_scenario("clean_snv", frags)


def scenario_r2_artefact(outdir: Path, header: pysam.AlignmentHeader) -> None:
    """
    All 4 alt reads are R2; R1 at the same locus carries the reference allele.
    Filter: R1 only → VAF drops to 0.
    Filter: R2 only → all 4 reads remain.
    Expected: is_read1=False for all 4 alt reads in alt_reads table.
    """
    frags = []
    for i in range(N_TOTAL_PAIRS):
        # Alt reads: force R2 to carry alt, R1 to carry ref
        # We achieve this by making the 'alt' pair have is_alt=False for R1
        # and then patching R2's sequence manually.
        is_alt = i < N_ALT_PAIRS
        frags.append(make_fragment(
            frag_id=i, is_alt=False, r1_start=DEFAULT_R1_START,
            insert_size=DEFAULT_INSERT, mapq=60, tags=[], header=header,
        ))
        if is_alt:
            # Replace R2 sequence with alt allele
            r1, r2 = frags[-1]
            r2_start = DEFAULT_R1_START + DEFAULT_INSERT - READ_LEN
            r2.query_sequence = _read_seq(True, False, DEFAULT_R1_START, DEFAULT_INSERT)
            r2.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            frags[-1] = (r1, r2)
    write_bam(outdir / "r2_artefact.bam", frags, header)
    _print_scenario("r2_artefact", frags)


def scenario_end_of_read(outdir: Path, header: pysam.AlignmentHeader) -> None:
    """
    Alt reads have the SNV at cycle 140 (near end of R1, qpos 139).
    Ref reads have the SNV at cycle 76 (middle of R1).
    Achieved by shifting R1 start so SNV_POS = r1_start + 139.
    Filter: cycle range 1–120 → alt reads excluded.
    Expected: alt reads have cycle=140 in alt_reads; ref reads have cycle=76.
    """
    alt_r1_start = SNV_POS - 139   # r1 starts at 111 → SNV at qpos 139 → cycle 140
    frags = []
    for i in range(N_TOTAL_PAIRS):
        is_alt = i < N_ALT_PAIRS
        r1_start = alt_r1_start if is_alt else DEFAULT_R1_START
        frags.append(make_fragment(
            frag_id=i, is_alt=is_alt, r1_start=r1_start,
            insert_size=DEFAULT_INSERT, mapq=60, tags=[], header=header,
        ))
    write_bam(outdir / "end_of_read.bam", frags, header)
    _print_scenario("end_of_read", frags)


def scenario_low_mapq(outdir: Path, header: pysam.AlignmentHeader) -> None:
    """
    Alt reads have MAPQ=10; ref reads have MAPQ=60.
    Filter: MAPQ range 20–60 → alt reads excluded.
    Expected: alt reads have map_qual=10 in alt_reads.
    """
    frags = []
    for i in range(N_TOTAL_PAIRS):
        is_alt = i < N_ALT_PAIRS
        mapq = 10 if is_alt else 60
        frags.append(make_fragment(
            frag_id=i, is_alt=is_alt, r1_start=DEFAULT_R1_START,
            insert_size=DEFAULT_INSERT, mapq=mapq, tags=[], header=header,
        ))
    write_bam(outdir / "low_mapq.bam", frags, header)
    _print_scenario("low_mapq", frags)


def scenario_short_insert(outdir: Path, header: pysam.AlignmentHeader) -> None:
    """
    Alt reads have insert_size=80 bp; ref reads have insert_size=200 bp.
    Both sets of reads cover the SNV (R1 starts at SNV-39=211 for short insert).
    Filter: insert size range 100–300 → alt reads excluded.
    Expected: alt reads have insert_size=80 in alt_reads.
    """
    alt_insert  = 80
    alt_r1_start = SNV_POS - 39   # 211; SNV at R1 qpos 39 → cycle 40
    frags = []
    for i in range(N_TOTAL_PAIRS):
        is_alt = i < N_ALT_PAIRS
        r1_start  = alt_r1_start    if is_alt else DEFAULT_R1_START
        insert    = alt_insert      if is_alt else DEFAULT_INSERT
        frags.append(make_fragment(
            frag_id=i, is_alt=is_alt, r1_start=r1_start,
            insert_size=insert, mapq=60, tags=[], header=header,
        ))
    write_bam(outdir / "short_insert.bam", frags, header)
    _print_scenario("short_insert", frags)


def scenario_singleton_family(outdir: Path, header: pysam.AlignmentHeader) -> None:
    """
    Alt reads have fgbio cD tag = 1 (singleton molecule).
    Ref reads have cD = 3 (supported by 3 raw reads each).
    Also sets aD=1/bD=0 for alt reads, aD=2/bD=1 for ref reads.
    Filter: family size range 2–∞ (exclude mode off) → alt reads excluded.
    Expected: alt reads have family_size=1 in alt_reads.
    Use with: geac collect --pipeline fgbio --read-type simplex
    """
    frags = []
    for i in range(N_TOTAL_PAIRS):
        is_alt = i < N_ALT_PAIRS
        if is_alt:
            tags = [("cD", 1), ("aD", 1), ("bD", 0)]
        else:
            tags = [("cD", 3), ("aD", 2), ("bD", 1)]
        frags.append(make_fragment(
            frag_id=i, is_alt=is_alt, r1_start=DEFAULT_R1_START,
            insert_size=DEFAULT_INSERT, mapq=60, tags=tags, header=header,
        ))
    write_bam(outdir / "singleton_family.bam", frags, header)
    _print_scenario("singleton_family", frags)


def scenario_germline_snv(outdir: Path, header: pysam.AlignmentHeader) -> None:
    """
    VAF ≈ 0.50: 40 reads (20 pairs), 10 alt pairs (20 alt reads).
    Resembles a germline heterozygous SNV. No artefact properties.
    Expected: alt_count≈20, total_depth≈40, VAF≈0.50.
    Should NOT be filtered by any reasonable per-read threshold.
    """
    n_pairs     = 20
    n_alt_pairs = 10
    frags = []
    for i in range(n_pairs):
        is_alt = i < n_alt_pairs
        frags.append(make_fragment(
            frag_id=i, is_alt=is_alt, r1_start=DEFAULT_R1_START,
            insert_size=DEFAULT_INSERT, mapq=60, tags=[], header=header,
        ))
    write_bam(outdir / "germline_snv.bam", frags, header)
    _print_scenario("germline_snv", frags)


# ── Summary helper ─────────────────────────────────────────────────────────────

def _print_scenario(name: str, frags: list) -> None:
    n_total = len(frags) * 2
    n_alt   = 0
    for r1, r2 in frags:
        for seg in (r1, r2):
            seq = seg.query_sequence or ""
            qpos = SNV_POS - seg.reference_start
            if 0 <= qpos < len(seq) and seq[qpos] == ALT_BASE:
                n_alt += 1
    print(f"  {name:<22}  {len(frags):>3} pairs  {n_total:>3} reads  {n_alt:>3} alt reads  VAF≈{n_alt/n_total:.2f}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=str(Path(__file__).parent / "bams"),
                        help="Output directory (default: bams/ next to this script)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Writing to: {outdir}")
    print(f"Reference:  {CHROM}:{SNV_POS+1} (1-based)  REF={REF_BASE}  ALT={ALT_BASE}")
    print(f"Read len:   {READ_LEN} bp  |  Default insert: {DEFAULT_INSERT} bp")
    print()

    ref_path = write_reference(outdir)
    print(f"  reference.fa written ({REF_LEN} bp)")
    print()

    header = _bam_header()
    print("Scenario              Pairs  Reads  Alt reads  VAF")
    print("─" * 60)
    scenario_clean_snv(outdir, header)
    scenario_r2_artefact(outdir, header)
    scenario_end_of_read(outdir, header)
    scenario_low_mapq(outdir, header)
    scenario_short_insert(outdir, header)
    scenario_singleton_family(outdir, header)
    scenario_germline_snv(outdir, header)

    print()
    print("Done.  Load reference.fa + any .bam in IGV to inspect.")
    print()
    print("Example geac collect commands:")
    print(f"  geac collect --bam bams/clean_snv.bam          --reference bams/reference.fa --sample-id clean_snv          --output clean_snv.parquet")
    print(f"  geac collect --bam bams/singleton_family.bam   --reference bams/reference.fa --sample-id singleton_family  --pipeline fgbio --read-type simplex --output singleton_family.parquet")


if __name__ == "__main__":
    main()
