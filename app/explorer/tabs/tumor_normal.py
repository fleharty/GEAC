"""Tumor/Normal tab — joins alt_bases to normal_evidence to classify each tumor alt locus."""
from __future__ import annotations

import altair as alt
import pandas as pd
import streamlit as st

from explorer.tab_context import TabContext

LABEL = "🔬 Tumor/Normal"

_CLS_COLOR = {
    "Somatic candidate":   "#2ca02c",
    "Germline-like":       "#d62728",
    "Artifact-like":       "#ff7f0e",
    "No normal coverage":  "#aec7e8",
    "No normal data":      "#c7c7c7",
}
_CLS_ORDER = [
    "Somatic candidate", "Artifact-like", "Germline-like",
    "No normal coverage", "No normal data",
]


def render(ctx: TabContext) -> None:
    if not ctx.has_normal_evidence:
        st.info(
            "No `normal_evidence` table found in this database. "
            "Run `geac annotate-normal` on each tumor/normal pair and include the "
            "`.normal_evidence.parquet` files when running `geac merge`."
        )
        return

    st.subheader("Normal evidence at tumor loci")
    st.caption(
        "Joins `alt_bases` to `normal_evidence` to show, for every tumor alt locus, "
        "how much support for the same allele was observed in the paired normal. "
        "All active sidebar filters are applied to the tumor loci."
    )

    # For each tumor locus, look up two pieces of info via LEFT JOINs:
    #   normal_depth     — from the anchor row (normal_alt_allele IS NULL)
    #   normal_alt_count — count of reads in normal supporting the same allele
    _tn_df = ctx.con.execute(f"""
        WITH ne_anchor AS (
            SELECT tumor_sample_id, chrom, pos, tumor_alt_allele,
                   normal_sample_id, normal_depth
            FROM normal_evidence
            WHERE normal_alt_allele IS NULL
        ),
        ne_match AS (
            SELECT tumor_sample_id, chrom, pos, tumor_alt_allele,
                   SUM(normal_alt_count) AS normal_alt_count
            FROM normal_evidence
            WHERE normal_alt_allele = tumor_alt_allele
            GROUP BY tumor_sample_id, chrom, pos, tumor_alt_allele
        )
        SELECT
            ab.sample_id                                       AS tumor_sample_id,
            ab.chrom,
            ab.pos,
            ab.alt_allele                                      AS tumor_alt_allele,
            ab.variant_type,
            ROUND(ab.alt_count * 1.0 / ab.total_depth, 4)    AS tumor_vaf,
            ab.alt_count                                       AS tumor_alt_count,
            ab.total_depth                                     AS tumor_depth,
            a.normal_sample_id,
            COALESCE(a.normal_depth, 0)                       AS normal_depth,
            COALESCE(m.normal_alt_count, 0)                   AS normal_alt_count,
            CASE
                WHEN a.normal_depth IS NULL OR a.normal_depth = 0 THEN 0.0
                ELSE ROUND(COALESCE(m.normal_alt_count, 0) * 1.0 / a.normal_depth, 4)
            END AS normal_vaf,
            CASE
                WHEN a.normal_depth IS NULL            THEN 'No normal data'
                WHEN a.normal_depth = 0                THEN 'No normal coverage'
                WHEN COALESCE(m.normal_alt_count, 0) = 0 THEN 'Somatic candidate'
                WHEN COALESCE(m.normal_alt_count, 0) * 1.0 / a.normal_depth >= 0.2
                                                       THEN 'Germline-like'
                ELSE                                        'Artifact-like'
            END AS classification
        FROM (SELECT * FROM {ctx.table_expr} WHERE {ctx.where}) ab
        LEFT JOIN ne_anchor a ON ab.sample_id = a.tumor_sample_id
            AND ab.chrom = a.chrom AND ab.pos = a.pos
            AND ab.alt_allele = a.tumor_alt_allele
        LEFT JOIN ne_match m ON ab.sample_id = m.tumor_sample_id
            AND ab.chrom = m.chrom AND ab.pos = m.pos
            AND ab.alt_allele = m.tumor_alt_allele
    """).df()

    if _tn_df.empty:
        st.info("No records match the current filters.")
        return

    _render_classification_summary(_tn_df)
    st.divider()
    _render_vaf_scatter(_tn_df)
    st.divider()
    _render_normal_depth_distribution(_tn_df)
    st.divider()
    _render_data_table(_tn_df)


def _render_classification_summary(_tn_df: pd.DataFrame) -> None:
    _cls_counts = (
        _tn_df["classification"]
        .value_counts()
        .rename_axis("classification")
        .reset_index(name="n_loci")
    )
    _cls_total = len(_tn_df)
    _cls_counts["pct"] = (_cls_counts["n_loci"] / _cls_total * 100).round(1)

    _cls_chart = (
        alt.Chart(_cls_counts)
        .mark_bar()
        .encode(
            alt.X("n_loci:Q", title="Loci"),
            alt.Y("classification:N", sort=_CLS_ORDER, title=None),
            alt.Color("classification:N",
                      scale=alt.Scale(
                          domain=list(_CLS_COLOR.keys()),
                          range=list(_CLS_COLOR.values()),
                      ),
                      legend=None),
            tooltip=[
                alt.Tooltip("classification:N"),
                alt.Tooltip("n_loci:Q", title="Loci"),
                alt.Tooltip("pct:Q", title="%", format=".1f"),
            ],
        )
        .properties(height=200, title="Loci by normal-evidence classification")
    )
    st.altair_chart(_cls_chart, width="stretch")
    st.caption(
        "**Somatic candidate**: no alt allele detected in the normal (normal_depth > 0).  "
        "**Germline-like**: same allele seen in normal at ≥ 20 % frequency.  "
        "**Artifact-like**: same allele seen in normal at < 20 % frequency (possible sequencing noise in normal).  "
        "**No normal coverage**: locus has no reads in the normal BAM.  "
        "**No normal data**: locus was not found in the normal_evidence table."
    )


def _render_vaf_scatter(_tn_df: pd.DataFrame) -> None:
    st.subheader("Tumor VAF vs Normal VAF")

    _scatter_df = _tn_df[_tn_df["classification"] != "No normal data"].copy()

    if _scatter_df.empty:
        st.info("No records with normal evidence data.")
        return

    _tn_scatter = (
        alt.Chart(_scatter_df)
        .mark_circle(opacity=0.6, size=40)
        .encode(
            alt.X("tumor_vaf:Q", title="Tumor VAF",
                  scale=alt.Scale(domain=[0, 1])),
            alt.Y("normal_vaf:Q", title="Normal VAF",
                  scale=alt.Scale(domain=[0, 1])),
            alt.Color("classification:N",
                      scale=alt.Scale(
                          domain=list(_CLS_COLOR.keys()),
                          range=list(_CLS_COLOR.values()),
                      )),
            tooltip=[
                alt.Tooltip("tumor_sample_id:N", title="Tumor sample"),
                alt.Tooltip("normal_sample_id:N", title="Normal sample"),
                alt.Tooltip("chrom:N"),
                alt.Tooltip("pos:Q"),
                alt.Tooltip("tumor_alt_allele:N", title="Alt allele"),
                alt.Tooltip("tumor_vaf:Q", title="Tumor VAF", format=".4f"),
                alt.Tooltip("normal_vaf:Q", title="Normal VAF", format=".4f"),
                alt.Tooltip("tumor_alt_count:Q", title="Tumor alt count"),
                alt.Tooltip("normal_alt_count:Q", title="Normal alt count"),
                alt.Tooltip("normal_depth:Q", title="Normal depth"),
                alt.Tooltip("classification:N"),
            ],
        )
        .properties(height=500)
    )
    _diag_df = pd.DataFrame({"x": [0, 1], "y": [0, 1]})
    _diag_line = (
        alt.Chart(_diag_df)
        .mark_line(strokeDash=[4, 4], color="grey", opacity=0.5)
        .encode(x="x:Q", y="y:Q")
    )
    st.altair_chart((_diag_line + _tn_scatter).properties(height=500), width="stretch")
    st.caption(
        "Each point is one tumor alt locus. Points near the diagonal suggest germline variants. "
        "Points on the left edge (low tumor VAF, no normal support) are somatic candidates."
    )


def _render_normal_depth_distribution(_tn_df: pd.DataFrame) -> None:
    st.subheader("Normal depth at tumor loci")

    _depth_df = _tn_df[_tn_df["classification"] != "No normal data"][["normal_depth"]].copy()
    if _depth_df.empty:
        return

    _nd_chart = (
        alt.Chart(_depth_df)
        .mark_bar(color="#4c78a8")
        .encode(
            alt.X("normal_depth:Q", bin=alt.Bin(maxbins=50), title="Normal depth"),
            alt.Y("count():Q", title="Loci"),
            tooltip=[
                alt.Tooltip("normal_depth:Q", title="Normal depth", bin=True),
                alt.Tooltip("count():Q", title="Loci"),
            ],
        )
        .properties(height=300)
    )
    st.altair_chart(_nd_chart, width="stretch")
    st.caption(
        "Distribution of normal read depth at tumor alt loci. "
        "Low-depth positions may have unreliable normal evidence."
    )


def _render_data_table(_tn_df: pd.DataFrame) -> None:
    with st.expander("Normal evidence data table"):
        _tbl_cols = [c for c in [
            "tumor_sample_id", "normal_sample_id", "chrom", "pos",
            "tumor_alt_allele", "variant_type",
            "tumor_vaf", "tumor_alt_count", "tumor_depth",
            "normal_vaf", "normal_alt_count", "normal_depth",
            "classification",
        ] if c in _tn_df.columns]
        st.dataframe(
            _tn_df[_tbl_cols].sort_values(
                ["classification", "tumor_vaf"], ascending=[True, False]
            ),
            width="stretch",
            hide_index=True,
        )
