"""Tests for bin-weighted coverage profile aggregation helpers."""

from __future__ import annotations

import os
import sys

import duckdb

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from coverage_profile import (
    expanded_profile_position_count,
    load_expanded_depth_profile,
    load_expanded_sample_profile,
    max_bin_width,
)


def test_depth_profile_bin_weighted_aggregation():
    """Bins are plotted at their own pos; no expansion across bin width.

    sample_a has a single 4bp bin [100,104) with depth=40.
    sample_b has four 1bp bins [100,101), [101,102), [102,103), [103,104)
    with depths 10, 20, 30, 40.

    Without expansion, sample_a only contributes at pos=100 (not 101–103).
    Cross-sample stats therefore differ at each position:
      pos=100: both samples → mean=(40+10)/2=25, n=2
      pos=101–103: sample_b only → mean=20/30/40, n=1
    """
    con = duckdb.connect()
    con.execute(
        """
        CREATE TABLE coverage (
            sample_id VARCHAR,
            pos BIGINT,
            "end" BIGINT,
            total_depth INTEGER,
            mean_mapq DOUBLE,
            frac_mapq0 DOUBLE,
            frac_low_mapq DOUBLE,
            gc_content DOUBLE
        )
        """
    )
    con.executemany(
        """
        INSERT INTO coverage
            (sample_id, pos, "end", total_depth, mean_mapq, frac_mapq0, frac_low_mapq, gc_content)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            ("sample_a", 100, 104, 40, 60, 0.0, 0.0, 0.40),
            ("sample_b", 100, 101, 10, 50, 0.1, 0.2, 0.10),
            ("sample_b", 101, 102, 20, 51, 0.2, 0.3, 0.20),
            ("sample_b", 102, 103, 30, 52, 0.3, 0.4, 0.30),
            ("sample_b", 103, 104, 40, 53, 0.4, 0.5, 0.40),
        ],
    )

    # Genomic span: MAX(end=104) - MIN(pos=100) = 4
    assert expanded_profile_position_count(con, "coverage", "") == 4

    profile = load_expanded_depth_profile(con, "coverage", "")
    # Four rows: one per distinct pos value present in the data
    assert profile["pos"].tolist() == [100, 101, 102, 103]
    # pos=100: mean of sample_a(40) and sample_b(10) = 25
    # pos=101–103: only sample_b → 20, 30, 40
    assert profile["mean_depth"].round(2).tolist() == [25.0, 20.0, 30.0, 40.0]
    assert profile["min_depth"].tolist() == [10.0, 20.0, 30.0, 40.0]
    assert profile["max_depth"].tolist() == [40.0, 20.0, 30.0, 40.0]
    assert profile["n_samples"].tolist() == [2, 1, 1, 1]

    sample_profile = load_expanded_sample_profile(con, "coverage", "")
    # sample_a only appears at pos=100; sample_b appears at all four positions
    assert sample_profile["pos"].tolist() == [100, 100, 101, 102, 103]
    assert sample_profile["sample_id"].tolist() == [
        "sample_a",
        "sample_b",
        "sample_b",
        "sample_b",
        "sample_b",
    ]
    assert sample_profile["depth"].tolist() == [40, 10, 20, 30, 40]


def test_display_step_rebuckets_to_bin_resolution():
    """With display_step=4, four 1bp positions should collapse to one display bin."""
    con = duckdb.connect()
    con.execute(
        """
        CREATE TABLE cov2 (
            sample_id VARCHAR,
            pos BIGINT,
            "end" BIGINT,
            total_depth INTEGER,
            mean_mapq DOUBLE,
            frac_mapq0 DOUBLE,
            frac_low_mapq DOUBLE,
            gc_content DOUBLE
        )
        """
    )
    # Two samples, each with a single 4bp bin — mirrors 100bp bin behaviour.
    con.executemany(
        """
        INSERT INTO cov2
            (sample_id, pos, "end", total_depth, mean_mapq, frac_mapq0, frac_low_mapq, gc_content)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            ("sample_a", 0, 4, 20, 60, 0.0, 0.0, 0.50),
            ("sample_b", 0, 4, 40, 55, 0.0, 0.0, 0.50),
        ],
    )

    assert max_bin_width(con, "cov2", "") == 4

    profile = load_expanded_depth_profile(con, "cov2", "", display_step=4)
    assert len(profile) == 1, "4bp bins with display_step=4 should yield 1 display point"
    assert profile["pos"].tolist() == [0]
    assert profile["mean_depth"].tolist() == [30.0]

    sample_profile = load_expanded_sample_profile(con, "cov2", "", display_step=4)
    assert len(sample_profile) == 2
    assert set(sample_profile["sample_id"]) == {"sample_a", "sample_b"}
