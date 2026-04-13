import math

import duckdb

from bait_bias_helpers import (
    compute_bait_bias_candidates,
    compute_bait_bias_locus_detail,
)


def test_compute_bait_bias_candidates_uses_sample_normalized_expected_depth():
    con = duckdb.connect()
    con.execute("""
        CREATE TABLE alt_bases AS
        SELECT * FROM (VALUES
            ('C1', 'chr1', 100, 'A', 'insertion', 'GENE1', 80, 80, TRUE),
            ('C2', 'chr1', 100, 'A', 'insertion', 'GENE1', 160, 160, TRUE),
            ('N1', 'chr2', 300, 'T', 'SNV', 'GENE2', 95, 100, TRUE)
        ) t(sample_id, chrom, pos, alt_allele, variant_type, gene, alt_count, total_depth, on_target)
    """)
    con.execute("""
        CREATE TABLE ref_bases AS
        SELECT * FROM (VALUES
            ('C1', 'chr1', 200, 100, TRUE),
            ('C1', 'chr1', 300, 100, TRUE),
            ('C2', 'chr1', 200, 200, TRUE),
            ('C2', 'chr1', 300, 200, TRUE),
            ('N1', 'chr1', 100, 150, TRUE),
            ('N1', 'chr1', 200, 150, TRUE),
            ('N1', 'chr1', 300, 150, TRUE)
        ) t(sample_id, chrom, pos, total_depth, on_target)
    """)

    out = compute_bait_bias_candidates(
        con,
        "alt_bases",
        "TRUE",
        ab_on_tgt="AND on_target = TRUE",
        rb_on_tgt="AND rb.on_target = TRUE",
    )

    assert len(out) == 2
    row = out[(out["chrom"] == "chr1") & (out["pos"] == 100)].iloc[0]
    assert row["n_carrier_samples"] == 2
    assert row["n_noncarrier_samples"] == 1
    assert row["median_observed_depth"] == 120.0
    assert row["median_expected_depth"] == 150.0
    assert math.isclose(row["median_depth_retained"], 0.8, rel_tol=1e-9)
    assert math.isclose(row["locus_relative_depth"], 1.0, rel_tol=1e-9)


def test_compute_bait_bias_locus_detail_reports_per_carrier_and_noncarrier_rows():
    con = duckdb.connect()
    con.execute("""
        CREATE TABLE alt_bases AS
        SELECT * FROM (VALUES
            ('C1', 'chr1', 100, 'A', 'insertion', 'GENE1', 80, 80, TRUE),
            ('C2', 'chr1', 100, 'A', 'insertion', 'GENE1', 160, 160, TRUE)
        ) t(sample_id, chrom, pos, alt_allele, variant_type, gene, alt_count, total_depth, on_target)
    """)
    con.execute("""
        CREATE TABLE ref_bases AS
        SELECT * FROM (VALUES
            ('C1', 'chr1', 200, 100, TRUE),
            ('C1', 'chr1', 300, 100, TRUE),
            ('C2', 'chr1', 200, 200, TRUE),
            ('C2', 'chr1', 300, 200, TRUE),
            ('N1', 'chr1', 100, 150, TRUE),
            ('N1', 'chr1', 200, 150, TRUE),
            ('N1', 'chr1', 300, 150, TRUE)
        ) t(sample_id, chrom, pos, total_depth, on_target)
    """)

    carriers, noncarriers = compute_bait_bias_locus_detail(
        con,
        "alt_bases",
        "TRUE",
        chrom="chr1",
        pos=100,
        alt_allele="A",
        ab_on_tgt="AND on_target = TRUE",
        rb_on_tgt="AND rb.on_target = TRUE",
    )

    assert list(carriers["sample_id"]) == ["C1", "C2"]
    assert math.isclose(carriers.iloc[0]["sample_baseline_depth"], 100.0, rel_tol=1e-9)
    assert math.isclose(carriers.iloc[1]["sample_baseline_depth"], 200.0, rel_tol=1e-9)
    assert math.isclose(carriers.iloc[0]["expected_depth"], 100.0, rel_tol=1e-9)
    assert math.isclose(carriers.iloc[1]["expected_depth"], 200.0, rel_tol=1e-9)
    assert math.isclose(carriers.iloc[0]["depth_retained"], 0.8, rel_tol=1e-9)
    assert math.isclose(carriers.iloc[1]["depth_retained"], 0.8, rel_tol=1e-9)

    assert len(noncarriers) == 1
    assert noncarriers.iloc[0]["sample_id"] == "N1"
    assert math.isclose(noncarriers.iloc[0]["sample_baseline_depth"], 150.0, rel_tol=1e-9)
    assert math.isclose(noncarriers.iloc[0]["relative_depth"], 1.0, rel_tol=1e-9)
