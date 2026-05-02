"""VAF distribution tab — per-variant-type VAF histograms plus gnomAD het bait-bias check."""
from __future__ import annotations

import altair as alt
import pandas as pd
import streamlit as st

from explorer.tab_context import TabContext

LABEL = "📊 VAF distribution"


def render(ctx: TabContext) -> None:
    for vtype, color in [
        ("SNV",       "#4c78a8"),
        ("insertion", "#f58518"),
        ("deletion",  "#e45756"),
    ]:
        with ctx.timed(f"vaf_histogram_{vtype}"):
            counts = ctx.con.execute(f"""
                SELECT
                    FLOOR(ROUND(alt_count * 1.0 / total_depth, 4) * 50) / 50.0 AS vaf_bin,
                    FLOOR(ROUND(alt_count * 1.0 / total_depth, 4) * 50) / 50.0 + 0.02 AS vaf_bin_end,
                    COUNT(*) AS count
                FROM {ctx.table_expr}
                WHERE {ctx.where} AND variant_type = '{vtype}'
                  AND total_depth > 0
                  AND alt_count <= total_depth
                GROUP BY vaf_bin, vaf_bin_end
                HAVING vaf_bin IS NOT NULL AND vaf_bin >= 0.0
                ORDER BY vaf_bin
            """).df()
            _vaf_stats = ctx.con.execute(f"""
                SELECT
                    AVG(alt_count * 1.0 / total_depth)    AS mean_vaf,
                    MEDIAN(alt_count * 1.0 / total_depth) AS median_vaf
                FROM {ctx.table_expr}
                WHERE {ctx.where} AND variant_type = '{vtype}'
                  AND total_depth > 0
                  AND alt_count <= total_depth
            """).fetchone()

        if counts.empty:
            st.info(f"No {vtype}s in current selection.")
        else:
            _sel_name = f"bar_click_{vtype}"
            sel_param = alt.selection_point(
                name=_sel_name,
                fields=["vaf_bin", "vaf_bin_end"],
                on="click",
            )
            chart = (
                alt.Chart(counts)
                .mark_bar(color=color)
                .encode(
                    alt.X("vaf_bin:Q",     title="VAF", scale=alt.Scale(domain=[0, 1])),
                    alt.X2("vaf_bin_end:Q"),
                    alt.Y("count:Q",       title="Count"),
                    opacity=alt.condition(sel_param, alt.value(1.0), alt.value(0.4)),
                    tooltip=[
                        alt.Tooltip("vaf_bin:Q",     title="Bin start", format=".3f"),
                        alt.Tooltip("vaf_bin_end:Q", title="Bin end",   format=".3f"),
                        alt.Tooltip("count:Q",       title="Count"),
                    ],
                )
                .add_params(sel_param)
                .properties(title=f"{vtype} VAF Distribution", height=300)
            )
            _vtype_count = int(counts["count"].sum())
            event = st.altair_chart(chart, width="stretch", on_select="rerun", key=f"vaf_chart_{vtype}")
            _mean_vaf = _vaf_stats[0] if _vaf_stats and _vaf_stats[0] is not None else None
            _median_vaf = _vaf_stats[1] if _vaf_stats and _vaf_stats[1] is not None else None
            _stats_txt = (
                f"mean VAF {_mean_vaf:.4f}, median VAF {_median_vaf:.4f}. "
                if _mean_vaf is not None and _median_vaf is not None
                else ""
            )
            st.caption(
                f"{_vtype_count:,} alt-allele records (one per unique alt allele observed "
                f"at a locus in a sample). {_stats_txt}Click a bar to drill down."
            )

            pts = (event.selection or {}).get(_sel_name, [])
            if pts:
                bin_start = pts[0].get("vaf_bin")
                bin_end   = pts[0].get("vaf_bin_end")
                if bin_start is not None and bin_end is not None:
                    sel = ctx.query_records([
                        f"variant_type = '{vtype}'",
                        f"ROUND(alt_count * 1.0 / total_depth, 4) >= {bin_start}",
                        f"ROUND(alt_count * 1.0 / total_depth, 4) < {bin_end}",
                    ])
                    st.caption(
                        f"{len(sel):,} {vtype} records with VAF in "
                        f"[{bin_start:.3f}, {bin_end:.3f})"
                    )
                    st.dataframe(sel[ctx.table_cols], width="stretch")
                    ctx.igv_buttons([
                        f"variant_type = '{vtype}'",
                        f"ROUND(alt_count * 1.0 / total_depth, 4) >= {bin_start}",
                        f"ROUND(alt_count * 1.0 / total_depth, 4) < {bin_end}",
                    ], sel, key=f"vaf_{vtype}_{bin_start}")

    # ── Allele fraction at gnomAD het sites (bait-bias proxy) ────────────────
    st.divider()
    st.subheader("Allele fraction at gnomAD het sites (bait-bias proxy)")
    st.caption(
        "At gnomAD-annotated sites where the observed allele fraction is 30–70%, "
        "the sample is likely a true germline heterozygote. Bait-bias would cause "
        "indel-containing fragments to be under-captured relative to reference-only "
        "fragments, pulling the allele fraction below 0.5 for insertions and deletions "
        "while SNVs remain near 0.5. "
        "The dashed line marks the expected allele fraction of 0.5."
    )
    if "gnomad_af" not in ctx.schema_cols:
        st.info(
            "This section requires `gnomad_af` annotation in the alt-bases table. "
            "Run `geac collect` with `--gnomad` to produce this, then re-open the database."
        )
        return

    with ctx.timed("bait_bias_gnomad_het"):
        _bb_df = ctx.con.execute(f"""
            SELECT
                variant_type,
                MIN(alt_count * 1.0 / total_depth)                                              AS min_vaf,
                PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY alt_count * 1.0 / total_depth)     AS q1,
                PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY alt_count * 1.0 / total_depth)     AS median_vaf,
                PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY alt_count * 1.0 / total_depth)     AS q3,
                MAX(alt_count * 1.0 / total_depth)                                              AS max_vaf,
                COUNT(*)                                                                        AS n
            FROM {ctx.table_expr}
            WHERE {ctx.where}
              AND gnomad_af IS NOT NULL
              AND total_depth > 0
              AND alt_count * 1.0 / total_depth BETWEEN 0.30 AND 0.70
            GROUP BY variant_type
        """).df()

    if _bb_df.empty:
        st.info(
            "No gnomAD-annotated het-confirmed records under current filters. "
            "Ensure `--gnomad` was provided during collection."
        )
        return

    _bb_df["iqr"]   = _bb_df["q3"] - _bb_df["q1"]
    _bb_df["lower"] = (_bb_df["q1"] - 1.5 * _bb_df["iqr"]).clip(lower=_bb_df["min_vaf"])
    _bb_df["upper"] = (_bb_df["q3"] + 1.5 * _bb_df["iqr"]).clip(upper=_bb_df["max_vaf"])

    _bb_base    = alt.Chart(_bb_df)
    _bb_whisker = _bb_base.mark_rule().encode(
        x=alt.X("variant_type:N", title="Variant type", axis=alt.Axis(labelAngle=0)),
        y=alt.Y("lower:Q",     title="Allele fraction", scale=alt.Scale(domain=[0.25, 0.75])),
        y2=alt.Y2("upper:Q"),
        color=alt.Color("variant_type:N", legend=None),
    )
    _bb_box = _bb_base.mark_bar(size=40).encode(
        x=alt.X("variant_type:N"),
        y=alt.Y("q1:Q"),
        y2=alt.Y2("q3:Q"),
        color=alt.Color("variant_type:N", legend=None),
        tooltip=[
            "variant_type",
            alt.Tooltip("n:Q",          title="N",           format=","),
            alt.Tooltip("median_vaf:Q", title="Median VAF",  format=".3f"),
            alt.Tooltip("q1:Q",         title="Q1",          format=".3f"),
            alt.Tooltip("q3:Q",         title="Q3",          format=".3f"),
            alt.Tooltip("lower:Q",      title="Whisker low",  format=".3f"),
            alt.Tooltip("upper:Q",      title="Whisker high", format=".3f"),
        ],
    )
    _bb_median = _bb_base.mark_tick(color="white", size=40, thickness=2).encode(
        x=alt.X("variant_type:N"),
        y=alt.Y("median_vaf:Q"),
    )
    _bb_ref_line = (
        alt.Chart(pd.DataFrame({"y": [0.5]}))
        .mark_rule(color="gray", strokeDash=[4, 4])
        .encode(y="y:Q")
    )
    st.altair_chart(
        (_bb_whisker + _bb_box + _bb_median + _bb_ref_line)
        .properties(
            height=300,
            title="Allele fraction at gnomAD het sites by variant type",
        )
        .interactive(),
        use_container_width=True,
    )
    _bb_counts = {r["variant_type"]: int(r["n"]) for _, r in _bb_df.iterrows()}
    st.caption(
        "gnomAD-annotated sites, observed allele fraction 30–70%. "
        + "  ".join(f"{vt}: n={_bb_counts.get(vt, 0):,}" for vt in ["SNV", "insertion", "deletion"] if vt in _bb_counts)
    )
