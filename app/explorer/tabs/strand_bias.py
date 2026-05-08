"""Strand bias tab — fwd/rev alt-read scatter with 95% Binomial CI band and click-through drill-down."""
from __future__ import annotations

import altair as alt
import hashlib
import numpy as np
import pandas as pd
import streamlit as st

from explorer.tab_context import TabContext

LABEL = "↕️ Strand bias"


def render(ctx: TabContext) -> None:
    _genes_available = ctx.has_data("gene")

    _sb_col1, _sb_col2 = st.columns(2)
    _sb_scale = _sb_col1.radio(
        "Axis scale", ["Linear", "log1p"], horizontal=True, key="sb_scale",
        help="log1p = log(1 + x), compresses large counts while keeping zeros visible.",
    )
    _use_log1p = _sb_scale == "log1p"

    _color_options = ["Variant type", "Sample"]
    if ctx.has_data("batch"):
        _color_options.append("Batch")
    if ctx.has_data("label1"):
        _color_options.append("Label 1")
    if ctx.has_data("label2"):
        _color_options.append("Label 2")
    if ctx.has_data("label3"):
        _color_options.append("Label 3")
    if ctx.has_data("timepoint"):
        _color_options.append("Timepoint")
    if ctx.has_data("on_target"):
        _color_options.append("On target")
    if ctx.has_data("variant_called"):
        _color_options.append("Called variant")
    _color_by = _sb_col2.radio(
        "Color by", _color_options, horizontal=True, key="sb_color_by",
    )

    _sb_opt_cols = ", ".join(filter(None, [
        "batch"          if ctx.has_data("batch")          else None,
        "label1"         if ctx.has_data("label1")         else None,
        "label2"         if ctx.has_data("label2")         else None,
        "label3"         if ctx.has_data("label3")         else None,
        "timepoint"      if ctx.has_data("timepoint")      else None,
        "on_target"      if ctx.has_data("on_target")      else None,
        "variant_called" if ctx.has_data("variant_called") else None,
        "gene"           if _genes_available             else None,
    ]))
    _SB_SAMPLE_THRESHOLD = 5_000
    _sb_needs_sample = ctx.total_count > _SB_SAMPLE_THRESHOLD
    _sb_show_all = False
    if _sb_needs_sample:
        st.warning(
            f"**{ctx.total_count:,} loci** exceed the scatter plot threshold of "
            f"{_SB_SAMPLE_THRESHOLD:,}. Showing a random sample of "
            f"{_SB_SAMPLE_THRESHOLD:,} points.",
            icon="⚠️",
        )
        _sb_show_all = st.checkbox(
            "Show all loci (may be slow)",
            value=False,
            key="sb_show_all",
        )

    _sb_sql_base = f"""
        SELECT sample_id, chrom, pos, ref_allele, alt_allele, variant_type,
               fwd_alt_count, rev_alt_count,
               ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
               {(', ' + _sb_opt_cols) if _sb_opt_cols else ''}
        FROM {ctx.table_expr}
        WHERE {ctx.where}
    """
    if _sb_needs_sample and not _sb_show_all:
        sample_df = ctx.con.execute(
            f"{_sb_sql_base} USING SAMPLE reservoir({_SB_SAMPLE_THRESHOLD}) REPEATABLE(42)"
        ).df()
    else:
        if _sb_show_all:
            st.warning("Loading all loci — this may take a moment and slow down your browser.", icon="🐢")
        sample_df = ctx.con.execute(_sb_sql_base).df()

    _log1p_ticks_linear = [1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000]

    if _use_log1p:
        sample_df["fwd_plot"] = np.log1p(sample_df["fwd_alt_count"])
        sample_df["rev_plot"] = np.log1p(sample_df["rev_alt_count"])
        _x_title = "Forward alt reads"
        _y_title = "Reverse alt reads"
    else:
        sample_df["fwd_plot"] = sample_df["fwd_alt_count"]
        sample_df["rev_plot"] = sample_df["rev_alt_count"]
        _x_title = "Forward alt reads"
        _y_title = "Reverse alt reads"

    max_val = max(
        float(sample_df["fwd_plot"].max()) if len(sample_df) > 0 else 50,
        float(sample_df["rev_plot"].max()) if len(sample_df) > 0 else 50,
        1.0,
    )

    _x_lin = np.arange(0, int(sample_df["fwd_alt_count"].max()) + 1, dtype=float)
    _z = 1.96
    _s_lo = (-_z + np.sqrt(_z**2 + 8 * _x_lin)) / 2
    _s_hi = ( _z + np.sqrt(_z**2 + 8 * _x_lin)) / 2
    _rev_min_lin = np.maximum(_s_lo**2 - _x_lin, 0.0)
    _rev_max_lin = _s_hi**2 - _x_lin

    if _use_log1p:
        _ci_band = pd.DataFrame({
            "fwd":     np.log1p(_x_lin),
            "rev_min": np.log1p(_rev_min_lin),
            "rev_max": np.log1p(_rev_max_lin),
        })
        _diag = pd.DataFrame({"fwd": [0.0, max_val], "rev": [0.0, max_val]})
    else:
        _ci_band = pd.DataFrame({
            "fwd":     _x_lin,
            "rev_min": _rev_min_lin,
            "rev_max": _rev_max_lin,
        })
        _diag = pd.DataFrame({"fwd": [0.0, max_val], "rev": [0.0, max_val]})

    ci_lower = (
        alt.Chart(_ci_band)
        .mark_line(color="steelblue", opacity=0.6, strokeDash=[3, 3], tooltip=None)
        .encode(alt.X("fwd:Q"), alt.Y("rev_min:Q"))
    )
    ci_upper = (
        alt.Chart(_ci_band)
        .mark_line(color="steelblue", opacity=0.6, strokeDash=[3, 3], tooltip=None)
        .encode(alt.X("fwd:Q"), alt.Y("rev_max:Q"))
    )
    diag_line = (
        alt.Chart(_diag)
        .mark_line(strokeDash=[6, 4], color="gray", opacity=0.7)
        .encode(
            alt.X("fwd:Q"),
            alt.Y("rev:Q"),
        )
    )
    _sb_n_pts = len(sample_df)
    _sb_title = (
        f"Strand Bias (log1p scale) — solid: perfect balance; dashed: 95% CI under Binomial(n, 0.5) — showing {_sb_n_pts:,} alt loci"
        if _use_log1p else
        f"Strand Bias — solid: perfect balance; dashed: 95% CI under Binomial(n, 0.5) — showing {_sb_n_pts:,} alt loci"
    )

    if _use_log1p:
        _max_linear = max(
            int(sample_df["fwd_alt_count"].max()),
            int(sample_df["rev_alt_count"].max()),
        )
        _tick_vals = [np.log1p(v) for v in _log1p_ticks_linear if v <= _max_linear * 1.1]
        _log1p_axis = alt.Axis(
            values=_tick_vals,
            labelExpr="format(exp(datum.value) - 1, 'd')",
        )
        _enc_x = alt.X("fwd_plot:Q", title=_x_title, axis=_log1p_axis,
                        scale=alt.Scale(domain=[0, max_val]))
        _enc_y = alt.Y("rev_plot:Q", title=_y_title, axis=_log1p_axis,
                        scale=alt.Scale(domain=[0, max_val]))
    else:
        _enc_x = alt.X("fwd_plot:Q", title=_x_title)
        _enc_y = alt.Y("rev_plot:Q", title=_y_title)

    _color_field, _color_title, _color_scale = {
        "Variant type":   ("variant_type:N",  "Variant type",   alt.Scale()),
        "Sample":         ("sample_id:N",      "Sample",         alt.Scale()),
        "Batch":          ("batch:N",          "Batch",          alt.Scale()),
        "Label 1":        ("label1:N",         "Label 1",        alt.Scale()),
        "Label 2":        ("label2:N",         "Label 2",        alt.Scale()),
        "Label 3":        ("label3:N",         "Label 3",        alt.Scale()),
        "Timepoint":      ("timepoint:N",      "Timepoint",      alt.Scale()),
        "On target":      ("on_target:N",      "On target",
                           alt.Scale(domain=[True, False], range=["#2ca02c", "#d62728"])),
        "Called variant": ("variant_called:N", "Called variant",
                           alt.Scale(domain=[True, False], range=["#2ca02c", "#d62728"])),
    }[_color_by]

    _sb_point_sel = alt.selection_point(
        name="sb_select",
        fields=["sample_id", "chrom", "pos", "ref_allele", "alt_allele"],
        on="click",
        toggle="event.shiftKey",
    )

    scatter = (
        alt.Chart(sample_df)
        .mark_point(size=40)
        .encode(
            _enc_x,
            _enc_y,
            alt.Color(_color_field, title=_color_title, scale=_color_scale),
            opacity=alt.condition(_sb_point_sel, alt.value(1.0), alt.value(0.35)),
            size=alt.condition(_sb_point_sel, alt.value(80), alt.value(30)),
            tooltip=(
                ["sample_id", "chrom", "pos", "ref_allele", "alt_allele",
                 "variant_type", "fwd_alt_count", "rev_alt_count", "vaf"]
                + (["on_target"] if ctx.has_data("on_target") else [])
                + (["variant_called"] if ctx.has_data("variant_called") else [])
                + (["gene"] if _genes_available else [])
            ),
        )
        .add_params(_sb_point_sel)
        .properties(title=_sb_title, height=350)
    )
    _sb_data_key = hashlib.md5(ctx.where.encode()).hexdigest()[:8]
    sb_event = st.altair_chart(
        (ci_lower + ci_upper + diag_line + scatter).resolve_scale(color="independent"),
        width="stretch",
        on_select="rerun",
        key=f"strand_bias_scatter_{_sb_data_key}",
    )

    sb_pts = (sb_event.selection or {}).get("sb_select", [])
    if sb_pts:
        _sb_or_clauses = " OR ".join(
            f"(sample_id = '{p['sample_id']}' AND chrom = '{p['chrom']}' "
            f"AND pos = {int(p['pos'])} AND ref_allele = '{p['ref_allele']}' "
            f"AND alt_allele = '{p['alt_allele']}')"
            for p in sb_pts
            if all(k in p for k in ["sample_id", "chrom", "pos", "ref_allele", "alt_allele"])
        )
        if _sb_or_clauses:
            _sb_sel_df = ctx.con.execute(f"""
                SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                       pos + 1 AS pos_display
                FROM {ctx.table_expr}
                WHERE {_sb_or_clauses}
                ORDER BY sample_id, chrom, pos
            """).df()

            n_pts = len(_sb_sel_df)
            n_smp = _sb_sel_df["sample_id"].nunique()
            st.caption(
                f"{n_pts} selected loci across {n_smp} sample(s) — "
                "shift-click to select multiple points"
            )
            st.dataframe(_sb_sel_df[ctx.table_cols], width="stretch", hide_index=True)
            ctx.igv_buttons(
                [f"({_sb_or_clauses})"],
                _sb_sel_df,
                key=f"sb_sel_{'_'.join(str(int(p['pos'])) for p in sb_pts[:5] if 'pos' in p)}",
            )
    else:
        st.caption("Click a point to select it; shift-click to select multiple.")
