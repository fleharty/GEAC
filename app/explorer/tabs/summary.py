"""Summary tab — top-line metrics, the main records table, and position drill-down."""
from __future__ import annotations

import streamlit as st

from explorer.main_table import (
    render_position_drilldown,
    render_records_table,
    render_summary_metrics,
)
from explorer.tab_context import TabContext
from igv_helpers import per_read_warning_note

LABEL = "📋 Summary"


def render(ctx: TabContext) -> None:
    fstats = ctx.con.execute(f"""
        SELECT
            COUNT(*)                                            AS n_records,
            COUNT(DISTINCT sample_id)                           AS n_samples,
            COALESCE(SUM(alt_count), 0)                         AS total_alt_bases,
            COALESCE(ROUND(AVG(alt_count * 1.0 / total_depth), 4), 0) AS mean_vaf,
            COALESCE(ROUND(AVG(total_depth), 1), 0)             AS mean_depth,
            COUNT(*) FILTER (WHERE variant_called IS NOT NULL)  AS n_annotated,
            COUNT(*) FILTER (WHERE variant_called = true)       AS n_called
        FROM {ctx.table_expr}
        WHERE {ctx.where}
    """).df()

    fn_annotated = int(fstats["n_annotated"][0])
    fn_called    = int(fstats["n_called"][0])
    fpct_called  = f"{100 * fn_called / fn_annotated:.1f}%" if fn_annotated > 0 else "N/A"

    render_summary_metrics(ctx.stats, fstats, ctx.pct_called, fpct_called)

    st.info(f"**{ctx.total_count:,}** records match the current filters.", icon="✅")

    _reads_banner = st.empty()
    r = ctx.reads
    if r.active:
        _active_parts = []
        if r.fs_has_data and (r.fs_lo > 0 or r.fs_hi < r.fs_max):
            _active_parts.append(f"family size: {r.fs_lo}–{r.fs_hi}")
        if r.cycle_lo > 1 or r.cycle_hi < r.cycle_max:
            _active_parts.append(f"cycle number: {r.cycle_lo}–{r.cycle_hi}")
        if r.mq_lo > 0 or r.mq_hi < r.mq_max:
            _active_parts.append(f"map qual: {r.mq_lo}–{r.mq_hi}")
        if r.n_total_has_data and (r.n_total_lo > 0 or r.n_total_hi < r.n_total_max):
            _active_parts.append(f"N in read: {r.n_total_lo}–{r.n_total_hi}")
        if r.is_has_data and (r.is_lo > r.is_min or r.is_hi < r.is_max):
            _active_parts.append(f"insert size: {r.is_lo}–{r.is_hi}")
        if r.gc_has_data and (r.gc_lo > 0.0 or r.gc_hi < 1.0):
            _active_parts.append(f"frag GC: {r.gc_lo:.2f}–{r.gc_hi:.2f}")
        if r.read_strand_sel != "All":
            _active_parts.append(r.read_strand_sel.lower())
        _reads_banner.warning(
            f"**Per-read filters active** ({'; '.join(_active_parts)}). "
            f"{per_read_warning_note(r.recompute_vaf)}"
        )

    # ── Data table ────────────────────────────────────────────────────────────────
    _tbl_limit_options = [100, 500, 1000, 5000, 10000, 50000, "All"]
    _tbl_limit_sel = st.selectbox(
        "Table row limit",
        _tbl_limit_options,
        index=1,
        key="table_limit_sel",
        help="Limits rows shown in the data table below. Has no effect on plots — they always use the full filtered dataset.",
    )
    _tbl_limit = None if _tbl_limit_sel == "All" else int(_tbl_limit_sel)

    with ctx.timed("main_table_query"):
        df = ctx.query_records(limit=_tbl_limit)

    _tbl_event = render_records_table(
        df,
        ctx.total_count,
        ctx.table_cols,
        igv_buttons=ctx.igv_buttons,
        key="main",
    )

    # ── Position-level drill-down ──────────────────────────────────────────────────
    _selected_rows = (_tbl_event.selection or {}).get("rows", [])
    if _selected_rows:
        _row = df.iloc[_selected_rows[0]]
        _new_locus = (str(_row["chrom"]), int(_row["pos"]), str(_row["alt_allele"]))
        if st.session_state.get("_drill_locus") != _new_locus:
            st.session_state["_drill_locus"] = _new_locus

    if "_drill_locus" in st.session_state:
        _chrom, _pos, _selected_alt = st.session_state["_drill_locus"]
        render_position_drilldown(
            con=ctx.con,
            table_expr=ctx.table_expr,
            schema_cols=ctx.schema_cols,
            chrom=_chrom,
            pos=_pos,
            selected_alt=_selected_alt,
            has_alt_reads=ctx.has_alt_reads,
            sql_str=ctx.sql_str,
            igv_buttons=ctx.igv_buttons,
        )
