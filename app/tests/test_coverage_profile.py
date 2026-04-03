"""Tests for interval-aware coverage profile expansion helpers."""

from __future__ import annotations

import os
import sys

import duckdb

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from coverage_profile import (
    expanded_profile_position_count,
    load_expanded_depth_profile,
    load_expanded_sample_profile,
)


def test_expanded_depth_profile_uses_interval_spans():
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

    assert expanded_profile_position_count(con, "coverage", "") == 4

    profile = load_expanded_depth_profile(con, "coverage", "")
    assert profile["pos"].tolist() == [100, 101, 102, 103]
    assert profile["mean_depth"].round(2).tolist() == [25.0, 30.0, 35.0, 40.0]
    assert profile["min_depth"].tolist() == [10.0, 20.0, 30.0, 40.0]
    assert profile["max_depth"].tolist() == [40.0, 40.0, 40.0, 40.0]
    assert profile["n_samples"].tolist() == [2, 2, 2, 2]

    sample_profile = load_expanded_sample_profile(con, "coverage", "")
    assert sample_profile["pos"].tolist() == [100, 100, 101, 101, 102, 102, 103, 103]
    assert sample_profile["sample_id"].tolist() == [
        "sample_a",
        "sample_b",
        "sample_a",
        "sample_b",
        "sample_a",
        "sample_b",
        "sample_a",
        "sample_b",
    ]
    assert sample_profile["depth"].tolist() == [40, 10, 40, 20, 40, 30, 40, 40]
