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

    Uses two-level aggregation:
      1. Expand each coverage interval to base positions and collapse to one
         value per (display_bin, sample_id) — so each sample contributes
         equally regardless of how many 1bp dropout records it has.
      2. Compute cross-sample statistics (mean, IQR, min/max) over the
         per-sample values.

    This ensures the IQR band reflects the spread across samples, not the
    spread of individual base positions across all samples mixed together.
    """
    base_expr = build_expanded_profile_expr(table_expr, where_clause)
    if display_step > 1:
        pos_expr = f"FLOOR(prof_pos / {display_step}) * {display_step}"
    else:
        pos_expr = "prof_pos"
    return con.execute(
        f"""
        WITH per_sample AS (
            SELECT
                {pos_expr}       AS pos,
                sample_id,
                AVG(total_depth) AS depth,
                AVG(mean_mapq)   AS mean_mapq,
                AVG(frac_mapq0)  AS frac_mapq0,
                AVG(frac_low_mapq) AS frac_low_mapq,
                AVG(gc_content)  AS gc_content
            FROM {base_expr}
            GROUP BY {pos_expr}, sample_id
        )
        SELECT
            pos,
            AVG(depth)   AS mean_depth,
            MIN(depth)   AS min_depth,
            MAX(depth)   AS max_depth,
            PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY depth) AS p25_depth,
            PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY depth) AS median_depth,
            PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY depth) AS p75_depth,
            COUNT(DISTINCT sample_id) AS n_samples,
            AVG(mean_mapq)       AS mean_mapq,
            AVG(frac_mapq0)      AS mean_frac_mapq0,
            AVG(frac_low_mapq)   AS mean_frac_low_mapq,
            AVG(gc_content)      AS mean_gc_content
        FROM per_sample
        GROUP BY pos
        ORDER BY pos
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
            {pos_expr}       AS pos,
            sample_id,
            AVG(total_depth) AS depth
        FROM {base_expr}
        GROUP BY {pos_expr}, sample_id
        ORDER BY {pos_expr}, sample_id
        """
    ).df()
