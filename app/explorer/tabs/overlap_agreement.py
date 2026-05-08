"""Overlap agreement tab — alt-agree vs alt-disagree scatter with concordance histogram and drill-down."""
from __future__ import annotations

import altair as alt
import hashlib
import numpy as np
import streamlit as st

from explorer.tab_context import TabContext

LABEL = "🤝 Overlap agreement"

_SAMPLE_THRESHOLD = 5_000


def render(ctx: TabContext) -> None:
    if not ctx.has_data("overlap_alt_agree"):
        st.info(
            "Overlap agreement columns (`overlap_alt_agree`, `overlap_alt_disagree`, "
            "`overlap_ref_agree`) are not present in this dataset. Re-run "
            "`geac collect` without `--no-overlap` to populate them.",
            icon="ℹ️",
        )
        return

    _opt_cols = ", ".join(filter(None, [
        "batch"          if ctx.has_data("batch")          else None,
        "on_target"      if ctx.has_data("on_target")      else None,
        "variant_called" if ctx.has_data("variant_called") else None,
        "gene"           if ctx.has_data("gene")           else None,
    ]))

    _sql_base = f"""
        SELECT
            sample_id, chrom, pos, ref_allele, alt_allele, variant_type,
            overlap_alt_agree, overlap_alt_disagree,
            ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
            ROUND(
                overlap_alt_agree * 1.0
                    / NULLIF(overlap_alt_agree + overlap_alt_disagree, 0),
                4
            ) AS concordance
            {(', ' + _opt_cols) if _opt_cols else ''}
        FROM {ctx.table_expr}
        WHERE {ctx.where}
          AND (overlap_alt_agree + overlap_alt_disagree) > 0
    """

    _needs_sample = ctx.total_count > _SAMPLE_THRESHOLD
    if _needs_sample:
        st.warning(
            f"**{ctx.total_count:,} loci** exceed the scatter plot threshold of "
            f"{_SAMPLE_THRESHOLD:,}. Showing a random sample of {_SAMPLE_THRESHOLD:,} points.",
            icon="⚠️",
        )
        _show_all = st.checkbox("Show all loci (may be slow)", value=False, key="oa_show_all")
        if _show_all:
            st.warning("Loading all loci — this may take a moment.", icon="🐢")
            df = ctx.con.execute(_sql_base).df()
        else:
            df = ctx.con.execute(
                f"{_sql_base} USING SAMPLE reservoir({_SAMPLE_THRESHOLD}) REPEATABLE(42)"
            ).df()
    else:
        df = ctx.con.execute(_sql_base).df()

    if df.empty:
        st.warning("No overlapping pairs found under current filters.", icon="🔎")
        return

    _color_options = ["Variant type", "Sample"]
    if ctx.has_data("batch"):
        _color_options.append("Batch")
    if ctx.has_data("on_target"):
        _color_options.append("On target")
    if ctx.has_data("variant_called"):
        _color_options.append("Called variant")

    _col1, _col2 = st.columns(2)
    _color_by = _col1.radio("Color by", _color_options, horizontal=True, key="oa_color_by")
    _scale = _col2.radio("Axis scale", ["Linear", "log1p"], horizontal=True, key="oa_scale")
    _use_log = _scale == "log1p"

    _color_field, _color_title, _color_scale = {
        "Variant type":   ("variant_type:N",  "Variant type",   alt.Scale()),
        "Sample":         ("sample_id:N",      "Sample",         alt.Scale()),
        "Batch":          ("batch:N",          "Batch",          alt.Scale()),
        "On target":      ("on_target:N",      "On target",
                           alt.Scale(domain=[True, False], range=["#2ca02c", "#d62728"])),
        "Called variant": ("variant_called:N", "Called variant",
                           alt.Scale(domain=[True, False], range=["#2ca02c", "#d62728"])),
    }[_color_by]

    _log1p_ticks = [1, 2, 5, 10, 25, 50, 100, 250, 500, 1000]

    if _use_log:
        df["agree_plot"] = np.log1p(df["overlap_alt_agree"])
        df["disagree_plot"] = np.log1p(df["overlap_alt_disagree"])
        _max_val = max(float(df["agree_plot"].max()), float(df["disagree_plot"].max()), 1.0)
        _max_linear = max(int(df["overlap_alt_agree"].max()), int(df["overlap_alt_disagree"].max()))
        _tick_vals = [np.log1p(v) for v in _log1p_ticks if v <= _max_linear * 1.1]
        _log_axis = alt.Axis(values=_tick_vals, labelExpr="format(exp(datum.value) - 1, 'd')")
        _enc_x = alt.X("agree_plot:Q",    title="Pairs agreeing on alt",       axis=_log_axis, scale=alt.Scale(domain=[0, _max_val]))
        _enc_y = alt.Y("disagree_plot:Q", title="Pairs disagreeing (alt vs ref)", axis=_log_axis, scale=alt.Scale(domain=[0, _max_val]))
    else:
        df["agree_plot"] = df["overlap_alt_agree"]
        df["disagree_plot"] = df["overlap_alt_disagree"]
        _max_val = max(float(df["agree_plot"].max()), float(df["disagree_plot"].max()), 1.0)
        _enc_x = alt.X("agree_plot:Q",    title="Pairs agreeing on alt")
        _enc_y = alt.Y("disagree_plot:Q", title="Pairs disagreeing (alt vs ref)")

    _diag_df = __import__("pandas").DataFrame({"x": [0.0, _max_val], "y": [0.0, _max_val]})
    _diag = (
        alt.Chart(_diag_df)
        .mark_line(strokeDash=[6, 4], color="gray", opacity=0.5)
        .encode(alt.X("x:Q"), alt.Y("y:Q"))
    )

    _sel = alt.selection_point(
        name="oa_select",
        fields=["sample_id", "chrom", "pos", "ref_allele", "alt_allele"],
        on="click",
        toggle="event.shiftKey",
    )

    _tooltip = (
        ["sample_id", "chrom", "pos", "ref_allele", "alt_allele", "variant_type",
         "overlap_alt_agree", "overlap_alt_disagree", "concordance", "vaf"]
        + (["on_target"] if ctx.has_data("on_target") else [])
        + (["variant_called"] if ctx.has_data("variant_called") else [])
        + (["gene"] if ctx.has_data("gene") else [])
    )

    _n_pts = len(df)
    _scatter = (
        alt.Chart(df)
        .mark_point(size=40)
        .encode(
            _enc_x,
            _enc_y,
            alt.Color(_color_field, title=_color_title, scale=_color_scale),
            opacity=alt.condition(_sel, alt.value(1.0), alt.value(0.35)),
            size=alt.condition(_sel, alt.value(80), alt.value(30)),
            tooltip=_tooltip,
        )
        .add_params(_sel)
        .properties(
            title=f"Overlap agreement — agree vs disagree ({_n_pts:,} loci with ≥1 overlapping pair)",
            height=350,
        )
    )

    _oa_data_key = hashlib.md5(ctx.where.encode()).hexdigest()[:8]
    _oa_event = st.altair_chart(
        (_diag + _scatter).resolve_scale(color="independent"),
        width="stretch",
        on_select="rerun",
        key=f"oa_scatter_{_oa_data_key}",
    )

    _sel_pts = (_oa_event.selection or {}).get("oa_select", [])
    if _sel_pts:
        _or_clauses = " OR ".join(
            f"(sample_id = '{p['sample_id']}' AND chrom = '{p['chrom']}' "
            f"AND pos = {int(p['pos'])} AND ref_allele = '{p['ref_allele']}' "
            f"AND alt_allele = '{p['alt_allele']}')"
            for p in _sel_pts
            if all(k in p for k in ["sample_id", "chrom", "pos", "ref_allele", "alt_allele"])
        )
        if _or_clauses:
            _sel_df = ctx.con.execute(f"""
                SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                       pos + 1 AS pos_display
                FROM {ctx.table_expr}
                WHERE {_or_clauses}
                ORDER BY sample_id, chrom, pos
            """).df()
            n = len(_sel_df)
            n_smp = _sel_df["sample_id"].nunique()
            st.caption(
                f"{n} selected loci across {n_smp} sample(s) — "
                "shift-click to select multiple points"
            )
            st.dataframe(_sel_df[ctx.table_cols], width="stretch", hide_index=True)
            ctx.igv_buttons(
                [f"({_or_clauses})"],
                _sel_df,
                key=f"oa_sel_{'_'.join(str(int(p['pos'])) for p in _sel_pts[:5] if 'pos' in p)}",
            )
    else:
        st.caption("Click a point to select it; shift-click to select multiple.")

    st.divider()
    _render_concordance_histogram(df)


def _render_concordance_histogram(df) -> None:
    st.subheader("Concordance distribution")
    st.caption(
        "Concordance = agree / (agree + disagree). "
        "Values near 1.0 mean both reads consistently support the alt; "
        "values near 0.0 suggest the alt is seen on only one strand or one read of a pair."
    )

    _hist = (
        alt.Chart(df.dropna(subset=["concordance"]))
        .mark_bar()
        .encode(
            alt.X("concordance:Q", bin=alt.Bin(maxbins=40), title="Concordance"),
            alt.Y("count():Q", title="Loci"),
            tooltip=[alt.Tooltip("concordance:Q", bin=True), alt.Tooltip("count():Q", title="Loci")],
        )
        .properties(height=200)
    )
    st.altair_chart(_hist, width="stretch")
