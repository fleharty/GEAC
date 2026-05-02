"""Panel of Normals tab — joins alt_bases to pon_evidence to flag artefacts/germline.

Each tumor alt locus is classified by how frequently the same allele appears
across the PoN cohort. Loci absent from the PoN are somatic candidates;
loci common in the PoN are likely recurrent artefacts or germline variants.
"""
from __future__ import annotations

import altair as alt
import streamlit as st

from explorer.tab_context import TabContext

LABEL = "🛡️ Panel of Normals"

_PON_CLS_COLOR = {
    "PoN clean":     "#2ca02c",
    "Rare in PoN":   "#ff7f0e",
    "Common in PoN": "#d62728",
    "No PoN data":   "#c7c7c7",
}
_PON_CLS_ORDER = ["PoN clean", "Rare in PoN", "Common in PoN", "No PoN data"]


def render(ctx: TabContext) -> None:
    if not ctx.has_pon_evidence:
        st.info(
            "No `pon_evidence` table found in this database. "
            "To build one: run `geac collect` on each normal sample, "
            "`geac merge` the results into a PoN DuckDB, then run "
            "`geac annotate-pon --tumor-parquet <tumor.parquet> --pon-db <pon.duckdb>` "
            "and include the `.pon_evidence.parquet` output when running `geac merge` "
            "for the cohort."
        )
        return

    st.subheader("Panel of Normals evidence at tumor loci")
    st.caption(
        "Joins `alt_bases` to `pon_evidence` to show, for every tumor alt locus, "
        "how frequently the same allele was observed across the Panel of Normals. "
        "Loci common in the PoN are likely sequencing artefacts or germline variants. "
        "All active sidebar filters are applied to the tumor loci."
    )

    _pon_df = ctx.con.execute(f"""
        SELECT
            ab.sample_id                                      AS tumor_sample_id,
            ab.chrom,
            ab.pos,
            ab.alt_allele                                     AS tumor_alt_allele,
            ab.variant_type,
            ROUND(ab.alt_count * 1.0 / ab.total_depth, 4)   AS tumor_vaf,
            COALESCE(pe.n_pon_samples, 0)                    AS n_pon_samples,
            pe.pon_total_samples,
            CASE
                WHEN pe.pon_total_samples IS NULL OR pe.pon_total_samples = 0 THEN NULL
                ELSE ROUND(COALESCE(pe.n_pon_samples, 0) * 1.0 / pe.pon_total_samples, 4)
            END AS pon_sample_fraction,
            pe.max_pon_vaf,
            pe.mean_pon_vaf,
            CASE
                WHEN pe.pon_total_samples IS NULL THEN 'No PoN data'
                WHEN COALESCE(pe.n_pon_samples, 0) = 0 THEN 'PoN clean'
                WHEN COALESCE(pe.n_pon_samples, 0) * 1.0 / pe.pon_total_samples >= 0.1
                    THEN 'Common in PoN'
                ELSE 'Rare in PoN'
            END AS pon_classification
        FROM (SELECT * FROM {ctx.table_expr} WHERE {ctx.where}) ab
        LEFT JOIN pon_evidence pe
               ON pe.tumor_sample_id = ab.sample_id
              AND pe.chrom           = ab.chrom
              AND pe.pos             = ab.pos
              AND pe.tumor_alt_allele = ab.alt_allele
    """).df()

    if _pon_df.empty:
        st.info("No records match the current filters.")
        return

    _pon_total = (
        int(_pon_df["pon_total_samples"].dropna().iloc[0])
        if not _pon_df["pon_total_samples"].dropna().empty
        else 0
    )
    _render_summary_metrics(_pon_df, _pon_total)
    st.divider()
    _render_classification_chart(_pon_df)
    st.divider()
    _render_vaf_scatter(_pon_df, _pon_total)
    st.divider()
    _render_max_pon_vaf_distribution(_pon_df)
    st.divider()
    _render_data_table(_pon_df)


def _render_summary_metrics(_pon_df, _pon_total: int) -> None:
    _n_clean   = int((_pon_df["pon_classification"] == "PoN clean").sum())
    _n_rare    = int((_pon_df["pon_classification"] == "Rare in PoN").sum())
    _n_common  = int((_pon_df["pon_classification"] == "Common in PoN").sum())
    _n_nodata  = int((_pon_df["pon_classification"] == "No PoN data").sum())

    _pc1, _pc2, _pc3, _pc4, _pc5 = st.columns(5)
    _pc1.metric("PoN samples",   f"{_pon_total:,}")
    _pc2.metric("PoN clean",     f"{_n_clean:,}",  help="No alt allele seen in any PoN sample")
    _pc3.metric("Rare in PoN",   f"{_n_rare:,}",   help="Alt seen in < 10 % of PoN samples")
    _pc4.metric("Common in PoN", f"{_n_common:,}", help="Alt seen in ≥ 10 % of PoN samples")
    _pc5.metric("No PoN data",   f"{_n_nodata:,}", help="Locus not in pon_evidence table")


def _render_classification_chart(_pon_df) -> None:
    _counts = (
        _pon_df["pon_classification"]
        .value_counts()
        .rename_axis("classification")
        .reset_index(name="n_loci")
    )
    _total = len(_pon_df)
    _counts["pct"] = (_counts["n_loci"] / _total * 100).round(1)

    _chart = (
        alt.Chart(_counts)
        .mark_bar()
        .encode(
            alt.X("n_loci:Q", title="Loci"),
            alt.Y("classification:N", sort=_PON_CLS_ORDER, title=None),
            alt.Color("classification:N",
                      scale=alt.Scale(
                          domain=list(_PON_CLS_COLOR.keys()),
                          range=list(_PON_CLS_COLOR.values()),
                      ),
                      legend=None),
            tooltip=[
                alt.Tooltip("classification:N"),
                alt.Tooltip("n_loci:Q", title="Loci"),
                alt.Tooltip("pct:Q", title="%", format=".1f"),
            ],
        )
        .properties(height=180, title="Loci by PoN classification")
    )
    st.altair_chart(_chart, width="stretch")


def _render_vaf_scatter(_pon_df, _pon_total: int) -> None:
    st.subheader("Tumor VAF vs PoN sample fraction")

    _scatter_df = _pon_df[_pon_df["pon_classification"] != "No PoN data"].copy()
    if _scatter_df.empty:
        return

    _scatter = (
        alt.Chart(_scatter_df)
        .mark_circle(opacity=0.6, size=40)
        .encode(
            alt.X("tumor_vaf:Q", title="Tumor VAF",
                  scale=alt.Scale(domain=[0, 1])),
            alt.Y("pon_sample_fraction:Q",
                  title=f"PoN sample fraction (N={_pon_total})",
                  scale=alt.Scale(domain=[0, 1])),
            alt.Color("pon_classification:N",
                      scale=alt.Scale(
                          domain=list(_PON_CLS_COLOR.keys()),
                          range=list(_PON_CLS_COLOR.values()),
                      )),
            tooltip=[
                alt.Tooltip("tumor_sample_id:N", title="Tumor sample"),
                alt.Tooltip("chrom:N"),
                alt.Tooltip("pos:Q"),
                alt.Tooltip("tumor_alt_allele:N", title="Alt allele"),
                alt.Tooltip("variant_type:N", title="Variant type"),
                alt.Tooltip("tumor_vaf:Q", title="Tumor VAF", format=".4f"),
                alt.Tooltip("n_pon_samples:Q", title="PoN samples with alt"),
                alt.Tooltip("pon_sample_fraction:Q", title="PoN fraction", format=".3f"),
                alt.Tooltip("max_pon_vaf:Q", title="Max PoN VAF", format=".4f"),
                alt.Tooltip("mean_pon_vaf:Q", title="Mean PoN VAF", format=".4f"),
            ],
        )
        .properties(height=500)
    )
    st.altair_chart(_scatter, width="stretch")
    st.caption(
        "Each point is one tumor alt locus. "
        "Somatic candidates cluster near Y = 0 (absent from PoN). "
        "Recurrent artefacts and germline variants appear at higher Y values. "
        f"PoN fraction threshold for 'Common in PoN' is ≥ 10 % "
        f"({max(1, round(_pon_total * 0.1))} / {_pon_total} samples)."
    )


def _render_max_pon_vaf_distribution(_pon_df) -> None:
    st.subheader("Max PoN VAF distribution")

    _vaf_df = _pon_df[_pon_df["max_pon_vaf"].notna()][["max_pon_vaf"]].copy()
    if _vaf_df.empty:
        st.info("No loci with PoN alt evidence in current selection.")
        return

    _chart = (
        alt.Chart(_vaf_df)
        .mark_bar(color="#d62728")
        .encode(
            alt.X("max_pon_vaf:Q", bin=alt.Bin(maxbins=50),
                  title="Max PoN VAF (highest alt frequency seen in any PoN sample)"),
            alt.Y("count():Q", title="Loci"),
            tooltip=[
                alt.Tooltip("max_pon_vaf:Q", title="Max PoN VAF", bin=True),
                alt.Tooltip("count():Q", title="Loci"),
            ],
        )
        .properties(height=280)
    )
    st.altair_chart(_chart, width="stretch")
    st.caption(
        "Distribution of the highest alt VAF seen across all PoN samples at each locus "
        "(restricted to loci with ≥ 1 PoN sample showing the alt). "
        "Low max VAF suggests sequencing noise; high max VAF suggests a germline polymorphism."
    )


def _render_data_table(_pon_df) -> None:
    with st.expander("PoN evidence data table"):
        _cols = [c for c in [
            "tumor_sample_id", "chrom", "pos", "tumor_alt_allele",
            "variant_type", "tumor_vaf",
            "n_pon_samples", "pon_total_samples", "pon_sample_fraction",
            "max_pon_vaf", "mean_pon_vaf", "pon_classification",
        ] if c in _pon_df.columns]
        st.dataframe(
            _pon_df[_cols].sort_values("pon_sample_fraction", ascending=False),
            width="stretch",
            hide_index=True,
        )
