"""Cohort tab — per-sample summary table plus four QC charts.

Drops any `sample_id IN (...)` clause from the WHERE so every sample in the
merged DuckDB always appears, even when the sidebar focuses on one sample.
The summary dataframe supports row-click to focus all other tabs on the
selected sample.
"""
from __future__ import annotations

import altair as alt
import pandas as pd
import streamlit as st

from explorer.tab_context import TabContext

LABEL = "👥 Cohort"


def render(ctx: TabContext) -> None:
    if not ctx.path.endswith(".duckdb"):
        st.info("Cohort view is available when loading a merged DuckDB file (`geac merge` output).")
        return

    st.subheader("Per-sample summary")
    st.caption(
        "Applies all active sidebar filters except the sample selection. "
        "Click a row to focus all other tabs on that sample."
    )

    # Build a where clause without the sample_sel condition so all samples
    # always appear in the cohort summary table.
    _cohort_conditions = [c for c in ctx.conditions if not c.startswith("sample_id IN")]
    _cohort_where = " AND ".join(_cohort_conditions) if _cohort_conditions else "true"

    # When batch is present each row is (sample_id, batch); combine into a
    # single display label so every chart has one series per unique unit.
    _has_batch = ctx.has_data("batch")
    _id_sql = "sample_id || ' / ' || batch" if _has_batch else "sample_id"
    _group_by = "sample_id, batch" if _has_batch else "sample_id"

    _has_overlap = ctx.has_data("overlap_alt_agree")
    _overlap_col = """
        ROUND(
            SUM(overlap_alt_agree) * 1.0
                / NULLIF(SUM(overlap_alt_agree + overlap_alt_disagree), 0),
            4
        ) AS overlap_concordance,
    """ if _has_overlap else "NULL AS overlap_concordance,"

    _cohort_stats = ctx.con.execute(f"""
        SELECT
            {_id_sql}                                            AS sample_label,
            sample_id,
            {'batch,' if _has_batch else ''}
            COUNT(*) FILTER (WHERE variant_type = 'SNV')        AS n_snv,
            COUNT(*) FILTER (WHERE variant_type = 'insertion')   AS n_insertion,
            COUNT(*) FILTER (WHERE variant_type = 'deletion')    AS n_deletion,
            ROUND(AVG(total_depth), 1)                           AS mean_depth,
            ROUND(AVG(alt_count * 1.0 / total_depth), 6)        AS mean_vaf,
            ROUND(
                AVG(fwd_alt_count * 1.0
                    / NULLIF(fwd_alt_count + rev_alt_count, 0)),
                4
            ) AS strand_balance,
            {_overlap_col}
            COUNT(*) AS n_total
        FROM {ctx.table_expr}
        WHERE {_cohort_where}
        GROUP BY {_group_by}
        ORDER BY {_group_by}
    """).df()

    if _cohort_stats.empty:
        st.warning("No records match the current filters.", icon="🔎")
        return

    _render_summary_table(_cohort_stats)

    st.divider()
    _render_strand_balance_scatter(ctx, _id_sql, _group_by, _cohort_where)

    st.divider()
    _render_alt_loci_vs_basequal(ctx, _has_batch, _cohort_where)

    _render_sbs6_stacked_bar(ctx, _id_sql, _group_by, _cohort_where)


def _render_summary_table(_cohort_stats: pd.DataFrame) -> None:
    _event = st.dataframe(
        _cohort_stats,
        width="stretch",
        on_select="rerun",
        selection_mode="single-row",
        hide_index=True,
        key="cohort_data_table",
    )

    _sel = (_event.selection or {}).get("rows", [])
    if _sel:
        _row    = _cohort_stats.iloc[_sel[0]]
        _sample = _row["sample_id"]
        _label  = _row["sample_label"]
        st.caption(
            f"Focused on **{_label}** — "
            "click button below to filter all tabs to this sample."
        )
        if st.button(f"Filter all tabs to {_label}"):
            st.session_state["sample_sel"] = [_sample]
            st.rerun()


def _render_strand_balance_scatter(
    ctx: TabContext, id_sql: str, group_by: str, cohort_where: str
) -> None:
    st.subheader("Strand balance by sample")

    _strand_stats = ctx.con.execute(f"""
        SELECT
            {id_sql} AS sample_label,
            ROUND(AVG(alt_count * 1.0 / total_depth), 6) AS mean_vaf,
            ROUND(
                AVG(fwd_alt_count * 1.0
                    / NULLIF(fwd_alt_count + rev_alt_count, 0)),
                4
            ) AS mean_strand_balance,
            COUNT(*) AS n_loci
        FROM {ctx.table_expr}
        WHERE {cohort_where} AND variant_type = 'SNV'
        GROUP BY {group_by}
    """).df()

    if _strand_stats.empty:
        st.info("No SNVs in current selection.")
        return

    _chart = (
        alt.Chart(_strand_stats)
        .mark_circle(size=80, opacity=0.85)
        .encode(
            alt.X("mean_strand_balance:Q",
                  title="Mean strand balance (0.5 = perfect)",
                  scale=alt.Scale(domain=[0, 1])),
            alt.Y("mean_vaf:Q", title="Mean VAF",
                  scale=alt.Scale(zero=False)),
            alt.Color("sample_label:N", title="Sample"),
            alt.Size("n_loci:Q", title="SNV loci",
                     scale=alt.Scale(range=[40, 300])),
            tooltip=[
                "sample_label:N",
                alt.Tooltip("mean_strand_balance:Q", format=".4f",
                            title="Mean strand balance"),
                alt.Tooltip("mean_vaf:Q", format=".6f", title="Mean VAF"),
                alt.Tooltip("n_loci:Q", format=",", title="SNV loci"),
            ],
        )
        .properties(title="Strand Balance vs Mean VAF (per sample)", height=350)
    )
    _ref = (
        alt.Chart(pd.DataFrame({"x": [0.5]}))
        .mark_rule(strokeDash=[4, 4], color="gray", opacity=0.6)
        .encode(alt.X("x:Q"))
    )
    st.altair_chart((_chart + _ref).properties(height=350), width="stretch")
    st.caption(
        "Each dot is one sample. x = mean strand balance (0.5 = perfect), "
        "y = mean VAF. Outliers in either axis may indicate a problematic sample."
    )


def _render_alt_loci_vs_basequal(
    ctx: TabContext, has_batch: bool, cohort_where: str
) -> None:
    st.subheader("Alt loci count vs mean base quality (per sample)")

    _label_sql  = "ab.sample_id || ' / ' || ab.batch" if has_batch else "ab.sample_id"
    _group_sql  = "ab.sample_id, ab.batch"            if has_batch else "ab.sample_id"
    _batch_sel  = "batch,"                             if has_batch else ""

    _df = ctx.con.execute(f"""
        SELECT
            {_label_sql} AS sample_label,
            COUNT(DISTINCT CONCAT(ab.chrom, ':', ab.pos, ':', ab.alt_allele)) AS n_alt_loci,
            ROUND(AVG(ar.base_qual), 2) AS mean_base_qual,
            COUNT(ar.rowid) AS n_reads
        FROM (
            SELECT sample_id, {_batch_sel} chrom, pos, alt_allele
            FROM {ctx.table_expr}
            WHERE {cohort_where}
        ) ab
        INNER JOIN alt_reads ar
            ON  ab.sample_id  = ar.sample_id
            AND ab.chrom      = ar.chrom
            AND ab.pos        = ar.pos
            AND ab.alt_allele = ar.alt_allele
        WHERE ar.base_qual IS NOT NULL
        GROUP BY {_group_sql}
    """).df()

    if _df.empty:
        st.info("No base quality data available (alt_reads table may be absent).")
        return

    _chart = (
        alt.Chart(_df)
        .mark_circle(size=80, opacity=0.85)
        .encode(
            alt.X("n_alt_loci:Q", title="Number of alt loci",
                  scale=alt.Scale(zero=True)),
            alt.Y("mean_base_qual:Q", title="Mean base quality (Phred)",
                  scale=alt.Scale(zero=False)),
            alt.Color("sample_label:N", title="Sample"),
            alt.Size("n_reads:Q", title="Alt-supporting reads",
                     scale=alt.Scale(range=[40, 300])),
            tooltip=[
                "sample_label:N",
                alt.Tooltip("n_alt_loci:Q", format=",", title="Alt loci"),
                alt.Tooltip("mean_base_qual:Q", format=".1f", title="Mean base qual"),
                alt.Tooltip("n_reads:Q", format=",", title="Alt-supporting reads"),
            ],
        )
        .properties(height=350, title="Alt loci count vs mean base quality (per sample)")
    )
    st.altair_chart(_chart, width="stretch")
    st.caption(
        "Each dot is one sample. x = number of distinct alt loci, "
        "y = mean base quality across all alt-supporting reads. "
        "Samples with many loci but low base quality may be artefact-driven."
    )


def _render_sbs6_stacked_bar(
    ctx: TabContext, id_sql: str, group_by: str, cohort_where: str
) -> None:
    st.subheader("SNV Count by Sample (SBS6 breakdown)")

    _df = ctx.con.execute(f"""
        SELECT
            {id_sql} AS sample_label,
            CASE
                WHEN ref_allele IN ('C','G') AND alt_allele IN ('A','T') THEN 'C>A'
                WHEN ref_allele IN ('C','G') AND alt_allele IN ('G','C') THEN 'C>G'
                WHEN ref_allele IN ('C','G') AND alt_allele IN ('T','A') THEN 'C>T'
                WHEN ref_allele IN ('T','A') AND alt_allele IN ('A','T') THEN 'T>A'
                WHEN ref_allele IN ('T','A') AND alt_allele IN ('C','G') THEN 'T>C'
                WHEN ref_allele IN ('T','A') AND alt_allele IN ('G','C') THEN 'T>G'
                ELSE 'other'
            END AS substitution,
            COUNT(*) AS n_snv
        FROM {ctx.table_expr}
        WHERE {cohort_where} AND variant_type = 'SNV'
        GROUP BY {group_by}, substitution
        ORDER BY {group_by}, substitution
    """).df()

    if _df.empty:
        st.info("No SNVs in current selection.")
        return

    _color_scale = alt.Scale(
        domain=["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"],
        range=["#1BBDEB", "#808080", "#E22926", "#CBCACB", "#A1CE63", "#EDB5C0"],
    )
    _chart = (
        alt.Chart(_df)
        .mark_bar()
        .encode(
            alt.X("sample_label:N", title="Sample", sort="-y",
                  axis=alt.Axis(labelAngle=-45, labelLimit=200)),
            alt.Y("n_snv:Q", title="SNV count", stack="zero"),
            alt.Color("substitution:N", title="Substitution",
                      scale=_color_scale),
            alt.Tooltip(["sample_label:N", "substitution:N", "n_snv:Q"]),
        )
        .properties(height=350, title="SNV count per sample colored by SBS6 substitution type")
    )
    _total = int(_df["n_snv"].sum())
    st.altair_chart(_chart, width="stretch")
    st.caption(
        f"{_total:,} SNVs across {_df['sample_label'].nunique():,} samples. "
        "Each bar shows the total SNV count for one sample, stacked by the six SBS "
        "substitution classes (C>A, C>G, C>T, T>A, T>C, T>G). Samples are sorted by "
        "total SNV count. A dominant C>T signal can indicate UV or FFPE damage; elevated "
        "C>A is associated with oxidative damage or smoking; a relatively flat distribution "
        "is typical of background noise."
    )
