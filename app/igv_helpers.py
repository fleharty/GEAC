"""Pure-Python / DuckDB helpers for IGV session generation.

Kept separate from geac_explorer.py so they can be unit-tested without
importing the Streamlit runtime or any visualisation libraries.
"""

from __future__ import annotations

import duckdb


def query_distinct_samples(
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    where: str,
) -> list[str]:
    """Return sorted list of distinct sample_ids matching *where* in *table_expr*.

    Must query the database directly rather than inspecting a display-limited
    DataFrame: the IGV cap warning must reflect the full dataset regardless of
    how many rows are shown in the UI.
    """
    return (
        con.execute(
            f"SELECT DISTINCT sample_id FROM {table_expr} WHERE {where} ORDER BY sample_id"
        )
        .df()["sample_id"]
        .tolist()
    )
