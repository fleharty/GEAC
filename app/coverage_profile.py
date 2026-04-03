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


def max_bin_width(
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    where_clause: str,
) -> int:
    """Return the largest (end - pos) value in the region, used as display step."""
    result = con.execute(
        f'SELECT MAX("end" - pos) FROM {table_expr} {where_clause}'
    ).fetchone()[0]
    return int(result) if result and result > 0 else 1


def load_expanded_depth_profile(
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    where_clause: str,
    display_step: int = 1,
) -> pd.DataFrame:
    """Load cross-sample depth profile statistics.

    Expands each coverage interval to individual base positions for correct
    cross-sample aggregation, then re-buckets to ``display_step`` base-pair
    bins so the plot matches the original data resolution and renders quickly.
    """
    base_expr = build_expanded_profile_expr(table_expr, where_clause)
    if display_step > 1:
        pos_expr = f"FLOOR(prof_pos / {display_step}) * {display_step}"
    else:
        pos_expr = "prof_pos"
    return con.execute(
        f"""
        SELECT
            {pos_expr}                AS pos,
            AVG(total_depth)          AS mean_depth,
            MIN(total_depth)          AS min_depth,
            MAX(total_depth)          AS max_depth,
            PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY total_depth) AS p25_depth,
            PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY total_depth) AS p75_depth,
            COUNT(DISTINCT sample_id) AS n_samples,
            AVG(mean_mapq)            AS mean_mapq,
            AVG(frac_mapq0)           AS mean_frac_mapq0,
            AVG(frac_low_mapq)        AS mean_frac_low_mapq,
            AVG(gc_content)           AS mean_gc_content
        FROM {base_expr}
        GROUP BY {pos_expr}
        ORDER BY {pos_expr}
        """
    ).df()


def load_expanded_sample_profile(
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    where_clause: str,
    display_step: int = 1,
) -> pd.DataFrame:
    """Load per-sample depth values for profile overlays.

    Re-buckets to ``display_step`` bins matching the display resolution used
    for the aggregate profile.
    """
    base_expr = build_expanded_profile_expr(table_expr, where_clause)
    if display_step > 1:
        pos_expr = f"FLOOR(prof_pos / {display_step}) * {display_step}"
    else:
        pos_expr = "prof_pos"
    return con.execute(
        f"""
        SELECT
            {pos_expr}           AS pos,
            sample_id,
            AVG(total_depth)     AS depth
        FROM {base_expr}
        GROUP BY {pos_expr}, sample_id
        ORDER BY {pos_expr}, sample_id
        """
    ).df()
