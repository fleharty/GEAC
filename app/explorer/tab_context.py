"""Shared context passed to each tab's `render(ctx)` function.

`TabContext` is built once in `geac_explorer.py` after sidebar and main-query
setup, then handed to whichever tab the user has selected. Tabs read from
`ctx` rather than reaching back into module-level globals — this is what
makes them extractable into their own files.

`PerReadFilters` is a snapshot of the per-read filter slider state. The
`*_has_data` flags are False when the underlying column is absent or empty;
in that case the slider was not rendered and lo/hi/max should not be read.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import duckdb
import pandas as pd


@dataclass(frozen=True)
class PerReadFilters:
    active: bool
    recompute_vaf: bool
    read_strand_sel: str          # "All" / "R1 only" / "R2 only"

    fs_has_data: bool
    fs_lo: int
    fs_hi: int
    fs_max: int

    cycle_lo: int
    cycle_hi: int
    cycle_max: int

    mq_lo: int
    mq_hi: int
    mq_max: int

    n_total_has_data: bool
    n_total_lo: int
    n_total_hi: int
    n_total_max: int

    is_has_data: bool
    is_lo: int
    is_hi: int
    is_min: int
    is_max: int

    gc_has_data: bool
    gc_lo: float
    gc_hi: float


@dataclass(frozen=True)
class TabContext:
    # Database / query primitives
    con: duckdb.DuckDBPyConnection
    table_expr: str               # may be a CTE wrapping per-read joins
    where: str                    # WHERE clause for the active filter set
    schema_cols: set[str]
    has_alt_reads: bool
    path: str                     # source DB or parquet path (e.g. for ".duckdb" suffix gating)
    r_join: str                   # SQL join expression connecting alt_reads to the filtered locus set
    reads_where: str              # SQL WHERE fragment for active per-read filters; "" if reads.active is False

    # Pre-computed summaries
    stats: pd.DataFrame
    pct_called: str
    total_count: int
    table_cols: list[str]

    # Per-read filter state
    reads: PerReadFilters

    # Helpers — passed as callables so tabs don't import from geac_explorer.py
    timed: Callable[[str], Any]
    query_records: Callable[..., pd.DataFrame]
    igv_buttons: Callable[..., Any]
    sql_str: Callable[[str], str]
    has_data: Callable[[str], bool]      # True iff column exists in schema AND has non-null values
    alt_reads_has: Callable[[str], bool]  # True iff column exists in alt_reads table
