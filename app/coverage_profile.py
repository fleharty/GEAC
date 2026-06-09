"""Helpers for bin-weighted coverage profile aggregation."""

from __future__ import annotations

import duckdb
import pandas as pd


def _bin_width() -> str:
    """SQL expression for the width of a coverage record in base positions."""
    return '("end" - pos)'


def profile_position_count(
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    where_clause: str,
) -> int:
    """Return the genomic span (MAX(end) - MIN(pos)) of the queried records.

    Used as a threshold check before rendering the depth profile — if the span
    is very large the user should narrow their region or increase --bin-size.
    """
    result = con.execute(
        f'SELECT MAX("end") - MIN(pos) FROM {table_expr} {where_clause}'
    ).fetchone()[0]
    return int(result) if result else 0


def profile_bin_count(
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    where_clause: str,
) -> int:
    """Return the number of distinct plotted bin starts in the queried records."""
    result = con.execute(
        f"SELECT COUNT(DISTINCT pos) FROM {table_expr} {where_clause}"
    ).fetchone()[0]
    return int(result) if result else 0


def load_expanded_depth_profile(
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    where_clause: str,
) -> pd.DataFrame:
    """Load cross-sample depth profile statistics.

    Bins are aggregated using bin_n-weighted averages so wider bins contribute
    proportionally to the display. Two-level aggregation:
      1. Collapse to one weighted value per (pos, sample_id).
      2. Compute cross-sample statistics (mean, IQR, min/max) over per-sample values.
    """
    w = _bin_width()
    return con.execute(
        f"""
        WITH per_sample AS (
            SELECT
                pos,
                sample_id,
                SUM(total_depth * {w}) / SUM({w})       AS depth,
                SUM(mean_mapq   * {w}) / SUM({w})       AS mean_mapq,
                SUM(frac_mapq0  * {w}) / SUM({w})       AS frac_mapq0,
                SUM(frac_low_mapq * {w}) / SUM({w})     AS frac_low_mapq,
                SUM(gc_content  * {w}) / SUM({w})       AS gc_content
            FROM {table_expr}
            {where_clause}
            GROUP BY pos, sample_id
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
) -> pd.DataFrame:
    """Load per-sample depth values for profile overlays."""
    w = _bin_width()
    return con.execute(
        f"""
        SELECT
            pos,
            sample_id,
            SUM(total_depth * {w}) / SUM({w})       AS depth
        FROM {table_expr}
        {where_clause}
        GROUP BY pos, sample_id
        ORDER BY pos, sample_id
        """
    ).df()
