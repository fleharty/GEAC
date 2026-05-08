"""Tests for bin-weighted coverage profile aggregation helpers."""

from __future__ import annotations

import os
import sys

import duckdb

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from coverage_profile import (
    profile_position_count,
    load_expanded_depth_profile,
    load_expanded_sample_profile,
)


def test_depth_profile_bin_weighted_aggregation():
    """Bin-width weighted aggregation produces correct cross-sample stats.

    Two samples, each with uniform 100bp bins over two positions:
      sample_a: [100,200) depth=40, [200,300) depth=20
      sample_b: [100,200) depth=10, [200,300) depth=40

    Cross-sample mean:
      pos=100: (40+10)/2 = 25
      pos=200: (20+40)/2 = 30
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
            ("sample_a", 100, 200, 40, 60, 0.0, 0.0, 0.40),
            ("sample_a", 200, 300, 20, 58, 0.0, 0.0, 0.30),
            ("sample_b", 100, 200, 10, 50, 0.1, 0.2, 0.10),
            ("sample_b", 200, 300, 40, 53, 0.0, 0.1, 0.40),
        ],
    )

    # Genomic span: MAX(end=300) - MIN(pos=100) = 200
    assert profile_position_count(con, "coverage", "") == 200

    profile = load_expanded_depth_profile(con, "coverage", "")
    assert profile["pos"].tolist() == [100, 200]
    assert profile["mean_depth"].round(2).tolist() == [25.0, 30.0]
    assert profile["min_depth"].tolist() == [10.0, 20.0]
    assert profile["max_depth"].tolist() == [40.0, 40.0]
    assert profile["n_samples"].tolist() == [2, 2]

    sample_profile = load_expanded_sample_profile(con, "coverage", "")
    assert sample_profile["pos"].tolist() == [100, 100, 200, 200]
    assert sample_profile["sample_id"].tolist() == [
        "sample_a", "sample_b", "sample_a", "sample_b"
    ]
    assert sample_profile["depth"].tolist() == [40.0, 10.0, 20.0, 40.0]
