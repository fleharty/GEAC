"""Read-type comparison tab — A/B compares the same sample under two `read_type` values.

Designed for evaluating duplex consensus processing: which loci are removed,
which are retained, and how the VAF / strand-balance / SBS96 distributions
shift between raw and consensus read types. Six views share one cached
FULL OUTER JOIN dataset, just like the Pipeline comparison tab.
"""
from __future__ import annotations

import altair as alt
import pandas as pd
import streamlit as st

from explorer.sbs96 import strat_sbs96_chart, to_spec96_strat
from explorer.tab_context import TabContext

LABEL = "🔬 Read-type comparison"


def render(ctx: TabContext) -> None:
    if not ctx.data_source.is_duckdb:
        st.info(
            "Read-type comparison requires a merged DuckDB file. "
            "Run `geac collect` separately for each BAM (raw, simplex, duplex) with the same "
            "`--sample-id`, then `geac merge` all Parquets into one DuckDB."
        )
        return

    _rt_read_types = ctx.data_source.distinct_values("read_type")
    if len(_rt_read_types) < 2:
        st.info(
            "Read-type comparison requires at least 2 distinct `read_type` values in the database. "
            f"This database contains only: {', '.join(str(r) for r in _rt_read_types)}. "
            "Re-run `geac collect` with a different `--read-type` value for the same sample "
            "and merge both Parquets into one DuckDB."
        )
        return

    st.caption(
        "Compare the same sample processed at two different read types (e.g. raw vs duplex). "
        "All active sidebar filters (chromosome, sample, VAF, etc.) are applied to both. "
        "The goal is to quantify what duplex consensus processing removes vs retains."
    )

    _rt_col1, _rt_col2 = st.columns(2)
    _rt_a = _rt_col1.selectbox(
        "Read type A", _rt_read_types, index=0, key="rt_cmp_a"
    )
    _rt_b_opts = [r for r in _rt_read_types if r != _rt_a]
    _rt_b = _rt_col2.selectbox(
        "Read type B", _rt_b_opts, index=0, key="rt_cmp_b"
    )

    _rt_wa = _with_read_type(ctx, _rt_a)
    _rt_wb = _with_read_type(ctx, _rt_b)

    _rt_df = _build_or_get_join_df(ctx, _rt_wa, _rt_wb)

    if _rt_df.empty:
        st.info("No records for either read type under the current filters.")
        return

    _rt_label_map = {
        "shared": "Shared",
        "only_a": f"Only {_rt_a}",
        "only_b": f"Only {_rt_b}",
    }
    _rt_df["concordance_label"] = _rt_df["concordance"].map(_rt_label_map)
    _rt_shared = _rt_df[_rt_df["concordance"] == "shared"]

    _render_concordance_summary(_rt_df, _rt_a, _rt_b)
    st.divider()
    _render_vaf_distribution(_rt_df, _rt_a, _rt_b)
    st.divider()
    _render_vaf_correlation(_rt_shared, _rt_a, _rt_b)
    st.divider()
    _render_strand_balance(_rt_df, _rt_a, _rt_b)
    st.divider()
    _render_sbs96_side_by_side(ctx, _rt_wa, _rt_wb, _rt_a, _rt_b)
    st.divider()
    _render_unique_loci_table(_rt_df, _rt_a, _rt_b)


def _build_or_get_join_df(ctx: TabContext, wa: str, wb: str) -> pd.DataFrame:
    # Cache on the two WHERE fragments — Streamlit evaluates all tab
    # bodies on every rerun, so without caching this runs on every
    # widget interaction anywhere in the app.
    _cache_key = ("rt_df", wa, wb)
    if _rt_cache_get("key") != _cache_key:
        _rt_cache_set("key", _cache_key)
        with ctx.timed("rt_df FULL OUTER JOIN [cache miss]"):
            _result = ctx.con.execute(f"""
                WITH a AS (
                    SELECT sample_id, chrom, pos, alt_allele, variant_type,
                           ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                           total_depth,
                           fwd_alt_count, rev_alt_count,
                           trinuc_context, ref_allele
                    FROM {ctx.table_expr}
                    WHERE {wa}
                ),
                b AS (
                    SELECT sample_id, chrom, pos, alt_allele, variant_type,
                           ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                           total_depth,
                           fwd_alt_count, rev_alt_count
                    FROM {ctx.table_expr}
                    WHERE {wb}
                )
                SELECT
                    COALESCE(a.sample_id,    b.sample_id)    AS sample_id,
                    COALESCE(a.chrom,        b.chrom)        AS chrom,
                    COALESCE(a.pos,          b.pos)          AS pos,
                    COALESCE(a.alt_allele,   b.alt_allele)   AS alt_allele,
                    COALESCE(a.variant_type, b.variant_type) AS variant_type,
                    a.vaf            AS vaf_a,
                    b.vaf            AS vaf_b,
                    a.total_depth    AS depth_a,
                    b.total_depth    AS depth_b,
                    a.fwd_alt_count  AS fwd_alt_a,
                    a.rev_alt_count  AS rev_alt_a,
                    b.fwd_alt_count  AS fwd_alt_b,
                    b.rev_alt_count  AS rev_alt_b,
                    a.trinuc_context,
                    a.ref_allele,
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
        _rt_cache_set("df", _result)
    return _rt_cache_get("df")


def _with_read_type(ctx: TabContext, read_type: str) -> str:
    read_filter = f"read_type = '{ctx.sql_str(str(read_type))}'"
    if not ctx.conditions:
        return read_filter
    return " AND ".join(ctx.conditions + [read_filter])


def _rt_cache_get(name: str):
    return st.session_state.get(f"_rt_{name}_cache")


def _rt_cache_set(name: str, value) -> None:
    st.session_state[f"_rt_{name}_cache"] = value


def _render_concordance_summary(_rt_df: pd.DataFrame, rt_a, rt_b) -> None:
    st.subheader("Locus concordance")
    st.caption(
        "Loci unique to the less-filtered read type were likely removed by consensus "
        "processing. Shared loci are retained across both read types."
    )

    _n_shared = int((_rt_df["concordance"] == "shared").sum())
    _n_only_a = int((_rt_df["concordance"] == "only_a").sum())
    _n_only_b = int((_rt_df["concordance"] == "only_b").sum())
    _m1, _m2, _m3 = st.columns(3)
    _m1.metric("Shared loci", f"{_n_shared:,}")
    _m2.metric(f"Only {rt_a}", f"{_n_only_a:,}")
    _m3.metric(f"Only {rt_b}", f"{_n_only_b:,}")

    _conc_counts = (
        _rt_df.groupby(["concordance_label", "variant_type"])
        .size()
        .reset_index(name="n_loci")
    )
    _conc_order = ["Shared", f"Only {rt_a}", f"Only {rt_b}"]
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


def _render_vaf_distribution(_rt_df: pd.DataFrame, rt_a, rt_b) -> None:
    st.subheader("VAF distribution")
    st.caption(
        "Duplex processing typically shifts the VAF distribution by removing "
        "low-VAF artefacts that arise from PCR or sequencing errors."
    )

    _vaf_a = _rt_df[_rt_df["concordance"].isin(["shared", "only_a"])][
        ["vaf_a", "variant_type"]
    ].dropna(subset=["vaf_a"]).rename(columns={"vaf_a": "vaf"})
    _vaf_a["read_type"] = str(rt_a)
    _vaf_b = _rt_df[_rt_df["concordance"].isin(["shared", "only_b"])][
        ["vaf_b", "variant_type"]
    ].dropna(subset=["vaf_b"]).rename(columns={"vaf_b": "vaf"})
    _vaf_b["read_type"] = str(rt_b)
    _vaf_both = pd.concat([_vaf_a, _vaf_b], ignore_index=True)

    if _vaf_both.empty:
        return

    _vaf_chart = (
        alt.Chart(_vaf_both)
        .transform_density(
            "vaf",
            as_=["vaf", "density"],
            groupby=["read_type"],
            extent=[0, 1],
            steps=100,
        )
        .mark_area(opacity=0.4)
        .encode(
            x=alt.X("vaf:Q", title="VAF", scale=alt.Scale(domain=[0, 1])),
            y=alt.Y("density:Q", title="Density", stack=None),
            color=alt.Color("read_type:N", title="Read type"),
            tooltip=[
                alt.Tooltip("vaf:Q", format=".3f"),
                "read_type:N",
                alt.Tooltip("density:Q", format=".3f"),
            ],
        )
        .properties(height=280)
    )
    st.altair_chart(_vaf_chart, width="stretch")


def _render_vaf_correlation(_rt_shared: pd.DataFrame, rt_a, rt_b) -> None:
    st.subheader("VAF correlation (shared loci)")

    if _rt_shared.empty:
        st.info("No shared loci for current filters.")
        return

    _vaf_diag = (
        alt.Chart(pd.DataFrame({"x": [0.0, 1.0], "y": [0.0, 1.0]}))
        .mark_line(color="gray", strokeDash=[4, 4], opacity=0.6)
        .encode(x="x:Q", y="y:Q")
    )
    _vaf_scatter = (
        alt.Chart(_rt_shared)
        .mark_circle(opacity=0.5, size=40)
        .encode(
            x=alt.X("vaf_a:Q", title=f"VAF — {rt_a}",
                    scale=alt.Scale(domain=[0.0, 1.0])),
            y=alt.Y("vaf_b:Q", title=f"VAF — {rt_b}",
                    scale=alt.Scale(domain=[0.0, 1.0])),
            color=alt.Color("variant_type:N", title="Variant type"),
            tooltip=[
                "sample_id", "chrom", "pos", "alt_allele", "variant_type",
                alt.Tooltip("vaf_a:Q", title=f"VAF ({rt_a})", format=".4f"),
                alt.Tooltip("vaf_b:Q", title=f"VAF ({rt_b})", format=".4f"),
            ],
        )
        .properties(width=450, height=380)
        .interactive()
    )
    st.altair_chart(_vaf_diag + _vaf_scatter, width="stretch")

    if len(_rt_shared) >= 2:
        _r = _rt_shared["vaf_a"].corr(_rt_shared["vaf_b"])
        st.caption(
            f"Pearson r = {_r:.4f} across {len(_rt_shared):,} shared loci. "
            "Points below the diagonal indicate higher VAF in read type B."
        )


def _render_strand_balance(_rt_df: pd.DataFrame, rt_a, rt_b) -> None:
    st.subheader("Strand balance")
    st.caption(
        "Fraction of alt-supporting reads on the forward strand. "
        "Values near 0.5 indicate balanced strand support; "
        "values near 0 or 1 suggest strand artefacts."
    )

    _rows = []
    for _lbl, _fw_col, _rv_col in [
        (str(rt_a), "fwd_alt_a", "rev_alt_a"),
        (str(rt_b), "fwd_alt_b", "rev_alt_b"),
    ]:
        _sub = _rt_df[[_fw_col, _rv_col, "variant_type"]].dropna()
        _total_alt = _sub[_fw_col] + _sub[_rv_col]
        _sub = _sub[_total_alt > 0].copy()
        _total_alt = _total_alt[_total_alt > 0]
        if not _sub.empty:
            _sub["frac_fwd"] = _sub[_fw_col] / _total_alt
            _sub["read_type"] = _lbl
            _rows.append(_sub[["frac_fwd", "variant_type", "read_type"]])

    if not _rows:
        st.info("No strand balance data available under current filters.")
        return

    _strand_df = pd.concat(_rows, ignore_index=True)
    _strand_chart = (
        alt.Chart(_strand_df)
        .transform_density(
            "frac_fwd",
            as_=["frac_fwd", "density"],
            groupby=["read_type"],
            extent=[0, 1],
            steps=80,
        )
        .mark_area(opacity=0.4)
        .encode(
            x=alt.X("frac_fwd:Q", title="Fraction of alt reads on forward strand",
                    scale=alt.Scale(domain=[0, 1])),
            y=alt.Y("density:Q", title="Density", stack=None),
            color=alt.Color("read_type:N", title="Read type"),
            tooltip=[
                alt.Tooltip("frac_fwd:Q", format=".3f", title="Frac forward"),
                "read_type:N",
            ],
        )
        .properties(height=250)
    )
    st.altair_chart(_strand_chart, width="stretch")


def _render_sbs96_side_by_side(ctx: TabContext, wa: str, wb: str, rt_a, rt_b) -> None:
    st.subheader("SBS96 error spectrum")

    if not ctx.has_data("trinuc_context"):
        st.info(
            "SBS96 spectrum requires the `trinuc_context` column. "
            "Re-run `geac collect` with a reference FASTA."
        )
        return

    def _sbs96(rt_where: str):
        _raw = ctx.con.execute(f"""
            SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
            FROM {ctx.table_expr}
            WHERE {rt_where}
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
        st.markdown(f"**{rt_a}** ({_total_a:,} SNVs)")
        st.altair_chart(strat_sbs96_chart(_s96_a, str(rt_a)), width="stretch")
    with _sc2:
        st.markdown(f"**{rt_b}** ({_total_b:,} SNVs)")
        st.altair_chart(strat_sbs96_chart(_s96_b, str(rt_b)), width="stretch")


def _render_unique_loci_table(_rt_df: pd.DataFrame, rt_a, rt_b) -> None:
    st.subheader("Loci unique to one read type")
    st.caption(
        "These loci are called in one read type but absent in the other. "
        "Loci only in the less-filtered read type are candidates for "
        "consensus-removed artefacts."
    )

    _uniq = _rt_df[_rt_df["concordance"] != "shared"].copy()
    _uniq["read_type"] = _uniq["concordance"].map(
        {"only_a": str(rt_a), "only_b": str(rt_b)}
    )
    _filter = st.radio(
        "Show",
        ["Both", f"Only {rt_a}", f"Only {rt_b}"],
        horizontal=True,
        key="rt_cmp_uniq_filter",
    )
    if _filter == f"Only {rt_a}":
        _show = _uniq[_uniq["concordance"] == "only_a"]
    elif _filter == f"Only {rt_b}":
        _show = _uniq[_uniq["concordance"] == "only_b"]
    else:
        _show = _uniq

    _cols = [c for c in [
        "read_type", "sample_id", "chrom", "pos", "alt_allele",
        "variant_type", "vaf_a", "vaf_b", "depth_a", "depth_b",
    ] if c in _show.columns]

    if _show.empty:
        st.success("No unique loci for current filter selection.")
    else:
        st.dataframe(
            _show[_cols]
            .rename(columns={
                "vaf_a":   f"vaf_{rt_a}",
                "vaf_b":   f"vaf_{rt_b}",
                "depth_a": f"depth_{rt_a}",
                "depth_b": f"depth_{rt_b}",
            })
            .sort_values(["read_type", "chrom", "pos"]),
            width="stretch",
            hide_index=True,
            key="rt_cmp_uniq_table",
        )
