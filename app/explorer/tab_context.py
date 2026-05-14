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

from explorer.data_source import DataSource


@dataclass(frozen=True)
class LociThresholds:
    """Numeric values of threshold-based locus filters (alt count, VAF, depth, etc.).

    These filters compare per-pipeline count columns and must be applied *after* the
    pipeline-comparison FULL OUTER JOIN so that a locus with sub-threshold evidence in
    one pipeline is not misclassified as absent. Categorical filters (chrom, sample,
    variant type, on_target, …) are unaffected and stay in the CTE WHERE clauses.
    """
    min_alt: int = 1       # default 1 = inactive (condition only added when > 1)
    max_alt: int = 0       # 0 = no cap
    vaf_lo: float = 0.0
    vaf_hi: float = 1.0
    min_depth: int = 0
    max_depth: int = 0     # 0 = no cap
    min_fwd_alt: int = 0
    min_rev_alt: int = 0
    min_overlap_agree: int = 0
    min_overlap_disagree: int = 0

    @property
    def any_active(self) -> bool:
        return (
            self.min_alt > 1
            or self.max_alt > 0
            or self.vaf_lo > 0.0
            or self.vaf_hi < 1.0
            or self.min_depth > 0
            or self.max_depth > 0
            or self.min_fwd_alt > 0
            or self.min_rev_alt > 0
            or self.min_overlap_agree > 0
            or self.min_overlap_disagree > 0
        )


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
    data_source: DataSource       # underlying source (Parquet or DuckDB) — exposes is_duckdb, distinct_values, etc.
    table_expr: str               # may be a CTE wrapping per-read joins
    where: str                    # WHERE clause for the active filter set
    conditions: list[str]         # individual WHERE-clause fragments before AND-joining; needed by comparison tabs that re-build the WHERE
    schema_cols: set[str]
    has_alt_reads: bool
    has_normal_evidence: bool     # True iff normal_evidence table is present
    has_pon_evidence: bool        # True iff pon_evidence table is present
    path: str                     # source DB or parquet path (e.g. for ".duckdb" suffix gating)
    r_join: str                   # SQL join expression connecting alt_reads to the filtered locus set
    reads_where: str              # SQL WHERE fragment for active per-read filters; "" if reads.active is False
    has_fragments: bool           # True iff fragments table/view is present
    fragments_where_parts: list[str]  # pre-built sidebar-filter SQL fragments against fragments table (chrom, midpoint, batch, label1/2/3, timepoint)

    # Threshold filter values — used by pipeline comparison to apply count/VAF/depth
    # filters post-join rather than pre-join, so sub-threshold evidence is visible
    loci_thresholds: LociThresholds
    threshold_conditions: list[str]   # subset of conditions[] that are threshold-based

    # Pre-computed summaries
    stats: pd.DataFrame
    pct_called: str
    total_count: int
    table_cols: list[str]
    alt_reads_cols: set[str]      # column names of the alt_reads table (empty set if absent)
    cfg: dict                     # parsed config.toml (e.g. anthropic_api_key for AI Plot Builder)

    # Per-read filter state
    reads: PerReadFilters

    # Helpers — passed as callables so tabs don't import from geac_explorer.py
    timed: Callable[[str], Any]
    query_records: Callable[..., pd.DataFrame]
    igv_buttons: Callable[..., Any]
    sql_str: Callable[[str], str]
    has_data: Callable[[str], bool]      # True iff column exists in schema AND has non-null values
    alt_reads_has: Callable[[str], bool]  # True iff column exists in alt_reads table
    build_provenance: Callable[..., pd.DataFrame]  # captures sidebar-filter state for signature-discovery downloads
