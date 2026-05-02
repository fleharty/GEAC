"""Pipeline comparison tab — A/B compares the same sample under two `pipeline` values.

Builds a single FULL OUTER JOIN dataset on (sample_id, chrom, pos, alt_allele)
and shares it across six views: concordance summary, VAF correlation,
unique-loci table, unique-loci characterization, side-by-side SBS96 spectra,
and depth comparison. The join result is cached in `st.session_state` keyed
on the two WHERE fragments to avoid recomputing on every Streamlit rerun.
"""
from __future__ import annotations

import altair as alt
import pandas as pd
import streamlit as st

from explorer.sbs96 import strat_sbs96_chart, to_spec96_strat
from explorer.tab_context import TabContext
from pipeline_compare_helpers import (
    build_unique_pipeline_characterization_df,
    summarize_unique_pipeline_groups,
)

LABEL = "🔀 Pipeline comparison"


def render(ctx: TabContext) -> None:
    if not ctx.data_source.is_duckdb:
        st.info(
            "Pipeline comparison requires a merged DuckDB file. "
            "Run `geac collect` twice with different `--pipeline` values and the same "
            "`--sample-id`, then `geac merge` both Parquets into one DuckDB."
        )
        return

    _pc_pipelines = ctx.data_source.distinct_values("pipeline")
    if len(_pc_pipelines) < 2:
        st.info(
            "Pipeline comparison requires at least 2 distinct `pipeline` values in the database. "
            f"This database contains only: {', '.join(str(p) for p in _pc_pipelines)}. "
            "Re-run `geac collect` with a different `--pipeline` value for the same sample "
            "and merge both Parquets into one DuckDB."
        )
        return

    st.caption(
        "Compare the same sample processed through two different pipelines. "
        "All active sidebar filters (chromosome, sample, VAF, etc.) are applied to both pipelines."
    )

    _pc_col1, _pc_col2 = st.columns(2)
    _pc_pipe_a = _pc_col1.selectbox(
        "Pipeline A", _pc_pipelines, index=0, key="pipe_cmp_a"
    )
    _pc_pipe_b_opts = [p for p in _pc_pipelines if p != _pc_pipe_a]
    _pc_pipe_b = _pc_col2.selectbox(
        "Pipeline B", _pc_pipe_b_opts, index=0, key="pipe_cmp_b"
    )

    # Strip any pipeline_sel filter from conditions — this tab manages
    # pipeline selection internally via its own A/B selectors.
    _pc_base_conds = [c for c in ctx.conditions if not c.startswith("pipeline IN")]
    _pc_wa = " AND ".join(
        _pc_base_conds + [f"pipeline = '{ctx.sql_str(str(_pc_pipe_a))}'"]
    ) if _pc_base_conds else f"pipeline = '{ctx.sql_str(str(_pc_pipe_a))}'"
    _pc_wb = " AND ".join(
        _pc_base_conds + [f"pipeline = '{ctx.sql_str(str(_pc_pipe_b))}'"]
    ) if _pc_base_conds else f"pipeline = '{ctx.sql_str(str(_pc_pipe_b))}'"

    _pc_df = _build_or_get_join_df(ctx, _pc_wa, _pc_wb)

    if _pc_df.empty:
        st.info("No records for either pipeline under the current filters.")
        return

    _pc_label_map = {
        "shared": "Shared",
        "only_a": f"Only {_pc_pipe_a}",
        "only_b": f"Only {_pc_pipe_b}",
    }
    _pc_df["concordance_label"] = _pc_df["concordance"].map(_pc_label_map)
    _pc_shared = _pc_df[_pc_df["concordance"] == "shared"]

    _render_concordance_summary(_pc_df, _pc_pipe_a, _pc_pipe_b)
    st.divider()
    _render_vaf_correlation(ctx, _pc_shared, _pc_pipe_a, _pc_pipe_b)
    st.divider()
    _render_unique_loci_table(_pc_df, _pc_pipe_a, _pc_pipe_b)
    st.divider()
    _render_unique_characterization(ctx, _pc_df, _pc_pipe_a, _pc_pipe_b)
    st.divider()
    _render_sbs96_side_by_side(ctx, _pc_wa, _pc_wb, _pc_pipe_a, _pc_pipe_b)
    st.divider()
    _render_depth_comparison(ctx, _pc_shared, _pc_pipe_a, _pc_pipe_b)


def _build_or_get_join_df(ctx: TabContext, wa: str, wb: str) -> pd.DataFrame:
    # Cache on the two WHERE fragments — this full-outer-join scans
    # table_expr twice and runs on every rerun (Streamlit evaluates all
    # tab bodies eagerly, not just the active tab).
    _cache_key = ("pc_df", wa, wb)
    if st.session_state.get("_pc_cache_key") != _cache_key:
        st.session_state["_pc_cache_key"] = _cache_key
        with ctx.timed("pc_df FULL OUTER JOIN [cache miss]"):
            _result = ctx.con.execute(f"""
                WITH a AS (
                    SELECT sample_id, chrom, pos, alt_allele, variant_type,
                           ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                           total_depth, alt_count,
                           trinuc_context, ref_allele
                    FROM {ctx.table_expr}
                    WHERE {wa}
                ),
                b AS (
                    SELECT sample_id, chrom, pos, alt_allele, variant_type,
                           ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                           total_depth, alt_count,
                           trinuc_context, ref_allele
                    FROM {ctx.table_expr}
                    WHERE {wb}
                )
                SELECT
                    COALESCE(a.sample_id,    b.sample_id)    AS sample_id,
                    COALESCE(a.chrom,        b.chrom)        AS chrom,
                    COALESCE(a.pos,          b.pos)          AS pos,
                    COALESCE(a.alt_allele,   b.alt_allele)   AS alt_allele,
                    COALESCE(a.variant_type, b.variant_type) AS variant_type,
                    a.vaf         AS vaf_a,
                    b.vaf         AS vaf_b,
                    a.total_depth AS depth_a,
                    b.total_depth AS depth_b,
                    a.alt_count   AS alt_count_a,
                    b.alt_count   AS alt_count_b,
                    COALESCE(a.trinuc_context, b.trinuc_context) AS trinuc_context,
                    COALESCE(a.ref_allele, b.ref_allele) AS ref_allele,
                    CASE
                        WHEN a.chrom IS NOT NULL AND b.chrom IS NOT NULL THEN 'shared'
                        WHEN a.chrom IS NOT NULL                         THEN 'only_a'
                        ELSE                                                  'only_b'
                    END AS concordance
                FROM a
                FULL OUTER JOIN b
                    ON  a.sample_id  = b.sample_id
                    AND a.chrom      = b.chrom
                    AND a.pos        = b.pos
                    AND a.alt_allele = b.alt_allele
            """).df()
        st.session_state["_pc_df_cache"] = _result
    return st.session_state["_pc_df_cache"]


def _render_concordance_summary(_pc_df: pd.DataFrame, pipe_a, pipe_b) -> None:
    st.subheader("Locus concordance")

    _n_shared = int((_pc_df["concordance"] == "shared").sum())
    _n_only_a = int((_pc_df["concordance"] == "only_a").sum())
    _n_only_b = int((_pc_df["concordance"] == "only_b").sum())
    _m1, _m2, _m3 = st.columns(3)
    _m1.metric("Shared loci", f"{_n_shared:,}")
    _m2.metric(f"Only {pipe_a}", f"{_n_only_a:,}")
    _m3.metric(f"Only {pipe_b}", f"{_n_only_b:,}")

    _conc_counts = (
        _pc_df.groupby(["concordance_label", "variant_type"])
        .size()
        .reset_index(name="n_loci")
    )
    _conc_order = ["Shared", f"Only {pipe_a}", f"Only {pipe_b}"]
    _conc_chart = (
        alt.Chart(_conc_counts)
        .mark_bar()
        .encode(
            x=alt.X("concordance_label:N", title=None, sort=_conc_order),
            y=alt.Y("n_loci:Q", title="Loci"),
            color=alt.Color("variant_type:N", title="Variant type"),
            tooltip=["concordance_label", "variant_type", "n_loci"],
        )
        .properties(height=280)
    )
    st.altair_chart(_conc_chart, width="stretch")


def _render_vaf_correlation(ctx: TabContext, _pc_shared: pd.DataFrame, pipe_a, pipe_b) -> None:
    st.subheader("VAF correlation (shared loci)")

    if _pc_shared.empty:
        st.info("No shared loci for current filters.")
        return

    _diag = (
        alt.Chart(pd.DataFrame({"x": [0.0, 1.0], "y": [0.0, 1.0]}))
        .mark_line(color="gray", strokeDash=[4, 4], opacity=0.6)
        .encode(x="x:Q", y="y:Q")
    )
    _vaf_sel = alt.selection_point(
        name="pc_vaf_select",
        fields=["sample_id", "chrom", "pos", "alt_allele"],
        on="click",
        toggle="event.shiftKey",
    )
    _vaf_scatter = (
        alt.Chart(_pc_shared)
        .mark_circle(size=40)
        .encode(
            x=alt.X("vaf_a:Q", title=f"VAF — {pipe_a}",
                    scale=alt.Scale(domain=[0.0, 1.0])),
            y=alt.Y("vaf_b:Q", title=f"VAF — {pipe_b}",
                    scale=alt.Scale(domain=[0.0, 1.0])),
            color=alt.Color("variant_type:N", title="Variant type"),
            opacity=alt.condition(_vaf_sel, alt.value(1.0), alt.value(0.4)),
            size=alt.condition(_vaf_sel, alt.value(120), alt.value(40)),
            tooltip=[
                "sample_id", "chrom", "pos", "alt_allele", "variant_type",
                alt.Tooltip("vaf_a:Q", title=f"VAF ({pipe_a})", format=".4f"),
                alt.Tooltip("vaf_b:Q", title=f"VAF ({pipe_b})", format=".4f"),
            ],
        )
        .add_params(_vaf_sel)
        .properties(width=450, height=400)
        .interactive()
    )
    _vaf_event = st.altair_chart(
        _diag + _vaf_scatter,
        width="stretch",
        on_select="rerun",
        key="pc_vaf_scatter",
    )

    _vaf_pts = (_vaf_event.selection or {}).get("pc_vaf_select", [])
    if _vaf_pts:
        _vaf_or = " OR ".join(
            f"(sample_id = '{ctx.sql_str(p['sample_id'])}' AND chrom = '{ctx.sql_str(p['chrom'])}' "
            f"AND pos = {int(p['pos'])} AND alt_allele = '{ctx.sql_str(p['alt_allele'])}')"
            for p in _vaf_pts
            if all(k in p for k in ["sample_id", "chrom", "pos", "alt_allele"])
        )
        if _vaf_or:
            _vaf_sel_df = ctx.con.execute(f"""
                SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
                FROM {ctx.table_expr}
                WHERE ({_vaf_or})
                ORDER BY sample_id, chrom, pos, pipeline
            """).df()
            _show_cols = [c for c in ctx.table_cols + ["pipeline"] if c in _vaf_sel_df.columns]
            st.caption(
                f"{len(_vaf_sel_df)} records for {len(_vaf_pts)} selected "
                f"loci (both pipelines shown) — shift-click to select multiple"
            )
            st.dataframe(_vaf_sel_df[_show_cols], width="stretch", hide_index=True)
            ctx.igv_buttons(
                [f"({_vaf_or})"],
                _vaf_sel_df,
                key=f"pc_vaf_{'_'.join(str(int(p['pos'])) for p in _vaf_pts[:5] if 'pos' in p)}",
                use_global_filters=False,
            )
    else:
        st.caption("Click a point to select it; shift-click to select multiple.")

    if len(_pc_shared) >= 2:
        _r = _pc_shared["vaf_a"].corr(_pc_shared["vaf_b"])
        st.caption(
            f"Pearson r = {_r:.4f} across {len(_pc_shared):,} shared loci. "
            "Points off the diagonal are discordant calls."
        )


def _render_unique_loci_table(_pc_df: pd.DataFrame, pipe_a, pipe_b) -> None:
    st.subheader("Loci unique to one pipeline")
    st.caption(
        "One pipeline calls these loci; the other does not. "
        "Highest-priority candidates for manual review."
    )

    _uniq = _pc_df[_pc_df["concordance"] != "shared"].copy()
    _uniq["pipeline"] = _uniq["concordance"].map({"only_a": pipe_a, "only_b": pipe_b})
    _filter = st.radio(
        "Show",
        ["Both", f"Only {pipe_a}", f"Only {pipe_b}"],
        horizontal=True,
        key="pipe_cmp_uniq_filter",
    )
    if _filter == f"Only {pipe_a}":
        _show = _uniq[_uniq["concordance"] == "only_a"]
    elif _filter == f"Only {pipe_b}":
        _show = _uniq[_uniq["concordance"] == "only_b"]
    else:
        _show = _uniq

    _cols = [c for c in [
        "pipeline", "sample_id", "chrom", "pos", "alt_allele",
        "variant_type", "vaf_a", "vaf_b", "depth_a", "depth_b",
    ] if c in _show.columns]

    if _show.empty:
        st.success("No unique loci for current filter selection.")
    else:
        st.dataframe(
            _show[_cols]
            .rename(columns={
                "vaf_a":   f"vaf_{pipe_a}",
                "vaf_b":   f"vaf_{pipe_b}",
                "depth_a": f"depth_{pipe_a}",
                "depth_b": f"depth_{pipe_b}",
            })
            .sort_values(["pipeline", "chrom", "pos"]),
            width="stretch",
            hide_index=True,
            key="pipe_cmp_uniq_table",
        )


def _render_unique_characterization(ctx: TabContext, _pc_df: pd.DataFrame, pipe_a, pipe_b) -> None:
    st.subheader("Unique loci characterization")
    st.caption(
        "These plots compare loci called by only one pipeline under the current filters. "
        "They are intended to help explain whether unique calls differ in VAF, support, "
        "or mutational context between the two pipelines."
    )

    _uniq = _pc_df[_pc_df["concordance"] != "shared"].copy()
    _uniq["pipeline"] = _uniq["concordance"].map({"only_a": pipe_a, "only_b": pipe_b})

    _uniq_cmp = build_unique_pipeline_characterization_df(_uniq, str(pipe_a), str(pipe_b))
    _uniq_summary = summarize_unique_pipeline_groups(_uniq_cmp, str(pipe_a), str(pipe_b))
    _group_a = f"Only {pipe_a}"
    _group_b = f"Only {pipe_b}"
    _uniq_a = _uniq_cmp[_uniq_cmp["group"] == _group_a]
    _uniq_b = _uniq_cmp[_uniq_cmp["group"] == _group_b]

    if _uniq_a.empty or _uniq_b.empty:
        st.info(
            "Unique-loci characterization needs at least one unique locus in each pipeline "
            "under the current filters."
        )
        return

    def _fmt_metric(v, fmt):
        return "NA" if v is None or pd.isna(v) else format(v, fmt)

    _um = _uniq_summary["metrics"]
    _mcols = st.columns(6)
    _mcols[0].metric(f"{pipe_a} unique loci", f"{_um[_group_a]['count']:,}")
    _mcols[1].metric(f"{pipe_a} median VAF", _fmt_metric(_um[_group_a]["median_vaf"], ".4f"))
    _mcols[2].metric(f"{pipe_a} median depth", _fmt_metric(_um[_group_a]["median_depth"], ".1f"))
    _mcols[3].metric(f"{pipe_b} unique loci", f"{_um[_group_b]['count']:,}")
    _mcols[4].metric(f"{pipe_b} median VAF", _fmt_metric(_um[_group_b]["median_vaf"], ".4f"))
    _mcols[5].metric(f"{pipe_b} median depth", _fmt_metric(_um[_group_b]["median_depth"], ".1f"))
    st.caption(_uniq_summary["summary"])

    _dist_cols = st.columns(2)
    with _dist_cols[0]:
        st.altair_chart(
            _overlay_dist_chart(
                _uniq_cmp.dropna(subset=["vaf"]),
                "vaf", "VAF distribution", "VAF", domain=[0, 1],
            ),
            width="stretch",
        )
    with _dist_cols[1]:
        st.altair_chart(
            _overlay_dist_chart(
                _uniq_cmp.dropna(subset=["depth"]),
                "depth", "Total depth distribution", "Total depth",
            ),
            width="stretch",
        )

    st.altair_chart(
        _overlay_dist_chart(
            _uniq_cmp.dropna(subset=["alt_count"]),
            "alt_count", "Alt-count distribution", "Alt count",
        ),
        width="stretch",
    )

    st.markdown("**Unique-only SBS96 spectrum**")
    if not ctx.has_data("trinuc_context"):
        st.info(
            "SBS96 spectrum requires the `trinuc_context` column. "
            "Re-run `geac collect` with a reference FASTA."
        )
        return

    def _uniq_sbs96(group_name: str):
        _raw = (
            _uniq_cmp[
                (_uniq_cmp["group"] == group_name)
                & (_uniq_cmp["variant_type"] == "SNV")
                & (_uniq_cmp["trinuc_context"].notna())
            ][["trinuc_context", "ref_allele", "alt_allele"]]
            .groupby(["trinuc_context", "ref_allele", "alt_allele"])
            .size()
            .reset_index(name="count")
        )
        return to_spec96_strat(_raw)

    _u96_a, _u_total_a = _uniq_sbs96(_group_a)
    _u96_b, _u_total_b = _uniq_sbs96(_group_b)
    if _u96_a is None or _u96_b is None:
        st.info("Insufficient unique SNV data for unique-only SBS96 comparison.")
        return

    _sc1, _sc2 = st.columns(2)
    with _sc1:
        st.markdown(f"**{pipe_a} unique loci** ({_u_total_a:,} SNVs)")
        st.altair_chart(strat_sbs96_chart(_u96_a, f"Only {pipe_a}"), width="stretch")
    with _sc2:
        st.markdown(f"**{pipe_b} unique loci** ({_u_total_b:,} SNVs)")
        st.altair_chart(strat_sbs96_chart(_u96_b, f"Only {pipe_b}"), width="stretch")


def _overlay_dist_chart(
    df: pd.DataFrame,
    value_col: str,
    title: str,
    x_title: str,
    key_col: str = "group",
    domain=None,
):
    _counts = df.groupby(key_col).size()
    if not _counts.empty and int(_counts.min()) >= 5:
        _density_kwargs = {
            "as_": [value_col, "density"],
            "groupby": [key_col],
            "steps": 100,
        }
        if domain is not None:
            _density_kwargs["extent"] = domain
        return (
            alt.Chart(df)
            .transform_density(value_col, **_density_kwargs)
            .mark_line(strokeWidth=2.5)
            .encode(
                x=alt.X(
                    f"{value_col}:Q",
                    title=x_title,
                    scale=alt.Scale(domain=domain) if domain else alt.Scale(),
                ),
                y=alt.Y("density:Q", title="Density", stack=None),
                color=alt.Color(f"{key_col}:N", title="Unique set"),
                tooltip=[
                    alt.Tooltip(f"{value_col}:Q", format=".4f" if value_col == "vaf" else ".2f"),
                    f"{key_col}:N",
                    alt.Tooltip("density:Q", format=".3f"),
                ],
            )
            .properties(title=title, height=260)
        )
    return (
        alt.Chart(df)
        .mark_bar(opacity=0.45)
        .encode(
            x=alt.X(
                f"{value_col}:Q",
                title=x_title,
                bin=alt.Bin(maxbins=40),
                scale=alt.Scale(domain=domain) if domain else alt.Scale(),
            ),
            y=alt.Y("count():Q", title="Count", stack=None),
            color=alt.Color(f"{key_col}:N", title="Unique set"),
            tooltip=[
                f"{key_col}:N",
                alt.Tooltip("count():Q", title="Count"),
            ],
        )
        .properties(title=title, height=260)
    )


def _render_sbs96_side_by_side(ctx: TabContext, wa: str, wb: str, pipe_a, pipe_b) -> None:
    st.subheader("SBS96 error spectrum")

    if not ctx.has_data("trinuc_context"):
        st.info(
            "SBS96 spectrum requires the `trinuc_context` column. "
            "Re-run `geac collect` with a reference FASTA."
        )
        return

    def _sbs96(pipe_where: str):
        _raw = ctx.con.execute(f"""
            SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
            FROM {ctx.table_expr}
            WHERE {pipe_where}
              AND variant_type = 'SNV'
              AND trinuc_context IS NOT NULL
              AND length(trinuc_context) = 3
            GROUP BY trinuc_context, ref_allele, alt_allele
        """).df()
        return to_spec96_strat(_raw)

    _s96_a, _total_a = _sbs96(wa)
    _s96_b, _total_b = _sbs96(wb)

    if _s96_a is None or _s96_b is None:
        st.info("Insufficient SNV data for spectrum comparison.")
        return

    _sc1, _sc2 = st.columns(2)
    with _sc1:
        st.markdown(f"**{pipe_a}** ({_total_a:,} SNVs)")
        st.altair_chart(strat_sbs96_chart(_s96_a, str(pipe_a)), width="stretch")
    with _sc2:
        st.markdown(f"**{pipe_b}** ({_total_b:,} SNVs)")
        st.altair_chart(strat_sbs96_chart(_s96_b, str(pipe_b)), width="stretch")


def _render_depth_comparison(ctx: TabContext, _pc_shared: pd.DataFrame, pipe_a, pipe_b) -> None:
    st.subheader("Depth comparison (shared loci)")

    if _pc_shared.empty:
        st.info("No shared loci for current filters.")
        return

    _depth_max = int(max(_pc_shared["depth_a"].max(), _pc_shared["depth_b"].max()))
    _diag = (
        alt.Chart(pd.DataFrame({"x": [0, _depth_max], "y": [0, _depth_max]}))
        .mark_line(color="gray", strokeDash=[4, 4], opacity=0.6)
        .encode(x="x:Q", y="y:Q")
    )
    _depth_sel = alt.selection_point(
        name="pc_depth_select",
        fields=["sample_id", "chrom", "pos", "alt_allele"],
        on="click",
        toggle="event.shiftKey",
    )
    _depth_scatter = (
        alt.Chart(_pc_shared)
        .mark_circle(size=40)
        .encode(
            x=alt.X("depth_a:Q", title=f"Total depth — {pipe_a}"),
            y=alt.Y("depth_b:Q", title=f"Total depth — {pipe_b}"),
            color=alt.Color("variant_type:N", title="Variant type"),
            opacity=alt.condition(_depth_sel, alt.value(1.0), alt.value(0.4)),
            size=alt.condition(_depth_sel, alt.value(120), alt.value(40)),
            tooltip=[
                "sample_id", "chrom", "pos", "alt_allele", "variant_type",
                alt.Tooltip("depth_a:Q", title=f"Depth ({pipe_a})"),
                alt.Tooltip("depth_b:Q", title=f"Depth ({pipe_b})"),
            ],
        )
        .add_params(_depth_sel)
        .properties(width=450, height=350)
        .interactive()
    )
    _depth_event = st.altair_chart(
        _diag + _depth_scatter,
        width="stretch",
        on_select="rerun",
        key="pc_depth_scatter",
    )

    _depth_pts = (_depth_event.selection or {}).get("pc_depth_select", [])
    if _depth_pts:
        _depth_or = " OR ".join(
            f"(sample_id = '{ctx.sql_str(p['sample_id'])}' AND chrom = '{ctx.sql_str(p['chrom'])}' "
            f"AND pos = {int(p['pos'])} AND alt_allele = '{ctx.sql_str(p['alt_allele'])}')"
            for p in _depth_pts
            if all(k in p for k in ["sample_id", "chrom", "pos", "alt_allele"])
        )
        if _depth_or:
            _depth_sel_df = ctx.con.execute(f"""
                SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
                FROM {ctx.table_expr}
                WHERE ({_depth_or})
                ORDER BY sample_id, chrom, pos, pipeline
            """).df()
            _show_cols = [c for c in ctx.table_cols + ["pipeline"] if c in _depth_sel_df.columns]
            st.caption(
                f"{len(_depth_sel_df)} records for {len(_depth_pts)} selected "
                f"loci (both pipelines shown) — shift-click to select multiple"
            )
            st.dataframe(_depth_sel_df[_show_cols], width="stretch", hide_index=True)
            ctx.igv_buttons(
                [f"({_depth_or})"],
                _depth_sel_df,
                key=f"pc_depth_{'_'.join(str(int(p['pos'])) for p in _depth_pts[:5] if 'pos' in p)}",
                use_global_filters=False,
            )
    else:
        st.caption("Click a point to select it; shift-click to select multiple.")

    st.caption(
        f"Points above the diagonal: higher depth in {pipe_b}. "
        f"Points below: higher depth in {pipe_a}. "
        "Systematic offset suggests different duplicate-collapsing or overlap behaviour."
    )
