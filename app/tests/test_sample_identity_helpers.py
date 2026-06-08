"""Tests for the Sample Identity / Duplicates SQL helpers.

Builds a tiny synthetic cohort DuckDB with a known duplicate pair (same genotypes,
different subject_id) and a known swap (same subject_id, divergent genotypes), then
asserts the pairwise fingerprint metrics and flags come out as expected.
"""
import os
import sys

import duckdb
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from explorer.sample_identity_helpers import (
    IdentityParams,
    classify_flags,
    compute_all_loci_jaccard,
    compute_contamination,
    compute_pairwise_identity,
    subject_map,
)


def _rows(sample_id, subject_id, positions, vaf):
    alt_count = int(round(vaf * 100))
    for pos in positions:
        yield {
            "sample_id": sample_id,
            "subject_id": subject_id,
            "chrom": "chr1",
            "pos": pos,
            "ref_allele": "A",
            "alt_allele": "G",
            "variant_type": "SNV",
            "total_depth": 100,
            "alt_count": alt_count,
        }


def _make_db() -> duckdb.DuckDBPyConnection:
    common = list(range(1, 11))        # 10 markers present in all samples
    dup_private = list(range(101, 141))  # 40 markers shared only by S1 & S2

    records = []
    # S1 / S2: identical germline het genotypes, DIFFERENT subjects -> duplicate
    records += list(_rows("S1", "P1a", common + dup_private, 0.50))
    records += list(_rows("S2", "P1b", common + dup_private, 0.50))
    # S3 / S4: SAME subject P2 but divergent genotypes at the common panel -> swap
    records += list(_rows("S3", "P2", common, 0.50))   # het
    records += list(_rows("S4", "P2", common, 0.95))   # hom-alt

    con = duckdb.connect()
    con.execute(
        """
        CREATE TABLE alt_bases (
            sample_id VARCHAR, subject_id VARCHAR, chrom VARCHAR, pos BIGINT,
            ref_allele VARCHAR, alt_allele VARCHAR, variant_type VARCHAR,
            total_depth INTEGER, alt_count INTEGER
        )
        """
    )
    con.register("rec_df", pd.DataFrame(records))
    con.execute("INSERT INTO alt_bases SELECT sample_id, subject_id, chrom, pos, ref_allele, alt_allele, variant_type, total_depth, alt_count FROM rec_df")
    return con


def _pair(df: pd.DataFrame, a: str, b: str) -> pd.Series:
    row = df[(df.sample_a == a) & (df.sample_b == b)]
    assert len(row) == 1, f"expected exactly one row for {a},{b}, got {len(row)}"
    return row.iloc[0]


def test_duplicate_pair_flagged_unknown_duplicate():
    con = _make_db()
    p = IdentityParams()
    pairs = compute_pairwise_identity(con, p)
    flagged = classify_flags(pairs, subject_map(con), p)

    s1s2 = _pair(flagged, "S1", "S2")
    assert s1s2.jaccard == 1.0
    assert s1s2.concordance == 1.0
    assert s1s2.flag == "UNKNOWN_DUPLICATE"


def test_same_subject_divergent_flagged_swap():
    con = _make_db()
    p = IdentityParams()
    pairs = compute_pairwise_identity(con, p)
    flagged = classify_flags(pairs, subject_map(con), p)

    s3s4 = _pair(flagged, "S3", "S4")
    assert s3s4.concordance == 0.0      # het vs hom-alt at every shared marker
    assert s3s4.flag == "POSSIBLE_SWAP"


def test_unrelated_pair_not_flagged():
    con = _make_db()
    p = IdentityParams()
    pairs = compute_pairwise_identity(con, p)
    flagged = classify_flags(pairs, subject_map(con), p)

    s1s3 = _pair(flagged, "S1", "S3")
    assert s1s3.jaccard < p.t_high       # only the common panel overlaps
    assert s1s3.flag == ""


def test_all_loci_jaccard_identical_for_duplicate():
    con = _make_db()
    p = IdentityParams()
    res = compute_all_loci_jaccard(con, ["S1", "S2"], p)
    assert _pair(res, "S1", "S2").all_loci_jaccard == 1.0


def test_contamination_runs_and_covers_samples():
    con = _make_db()
    p = IdentityParams()
    cont = compute_contamination(con, p)
    # S1, S2, S3 carry het markers; S4 is hom-alt only (no het markers).
    assert set(cont.sample_id) >= {"S1", "S2", "S3"}
    # Clean het samples have ~0 dispersion (VAF exactly 0.5).
    s1 = cont[cont.sample_id == "S1"].iloc[0]
    assert abs(s1.het_vaf_dispersion) < 1e-9
