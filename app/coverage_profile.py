"""Helpers for interval-aware coverage profile aggregation."""

from __future__ import annotations

import duckdb
import pandas as pd


_PROFILE_BASE_COLUMNS = """
    sample_id,
    UNNEST(range(pos, "end")) AS prof_pos,
    total_depth,
    mean_mapq,
    frac_mapq0,
    frac_low_mapq,
    gc_content
"""


def build_expanded_profile_expr(table_expr: str, where_clause: str) -> str:
    """Return a subquery that expands coverage rows to base-resolution positions."""
    return f"""
        (
            SELECT
                {_PROFILE_BASE_COLUMNS}
            FROM {table_expr}
            {where_clause}
        )
    """


def expanded_profile_position_count(
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    where_clause: str,
) -> int:
    """Count distinct genomic positions after interval expansion."""
    base_expr = build_expanded_profile_expr(table_expr, where_clause)
    return int(
        con.execute(f"SELECT COUNT(DISTINCT prof_pos) FROM {base_expr}").fetchone()[0] or 0
    )


def load_expanded_depth_profile(
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    where_clause: str,
) -> pd.DataFrame:
    """Load cross-sample depth profile statistics at base resolution."""
    base_expr = build_expanded_profile_expr(table_expr, where_clause)
    return con.execute(
        f"""
        SELECT
            prof_pos AS pos,
            AVG(total_depth)   AS mean_depth,
            MIN(total_depth)   AS min_depth,
            MAX(total_depth)   AS max_depth,
            PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY total_depth) AS p25_depth,
            PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY total_depth) AS p75_depth,
            COUNT(DISTINCT sample_id) AS n_samples,
            AVG(mean_mapq)     AS mean_mapq,
            AVG(frac_mapq0)    AS mean_frac_mapq0,
            AVG(frac_low_mapq) AS mean_frac_low_mapq,
            AVG(gc_content)    AS mean_gc_content
        FROM {base_expr}
        GROUP BY prof_pos
        ORDER BY prof_pos
        """
    ).df()


def load_expanded_sample_profile(
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    where_clause: str,
) -> pd.DataFrame:
    """Load per-sample depth values at base resolution for profile overlays."""
    base_expr = build_expanded_profile_expr(table_expr, where_clause)
    return con.execute(
        f"""
        SELECT
            prof_pos AS pos,
            sample_id,
            total_depth AS depth
        FROM {base_expr}
        ORDER BY prof_pos, sample_id
        """
    ).df()
