"""Duplex/Simplex tab — family-size, AB/BA strand, and family-stratified spectrum plots.

All plots in this tab require the per-read `alt_reads` table, which is produced
when `geac collect` is run with `--reads-output`. Without it, the tab shows an
explanatory banner and no plots.
"""
from __future__ import annotations

import altair as alt
import streamlit as st

from explorer.sbs96 import strat_sbs96_chart, to_spec96_strat
from explorer.tab_context import TabContext

LABEL = "🔁 Duplex/Simplex"


def render(ctx: TabContext) -> None:
    if not ctx.has_alt_reads:
        st.info(
            "Per-read detail table not available. "
            "Re-run `geac collect` with `--reads-output` and merge the resulting "
            "`.reads.parquet` files to enable this tab."
        )
        return

    st.caption(
        "All plots reflect alt-supporting reads linked to loci that pass the current sidebar filters."
    )

    _render_family_size_distribution(ctx)
    _render_family_size_vs_vaf(ctx)
    _render_family_size_stratified_spectrum(ctx)
    _render_cohort_artefact_family_size(ctx)
    _render_ab_ba_strand_counts(ctx)


def _render_family_size_distribution(ctx: TabContext) -> None:
    st.subheader("Family size distribution")
    _fs_ctrl_col1, _fs_ctrl_col2 = st.columns(2)
    _fs_color_options = ["All samples (aggregate)", "Sample"]
    if ctx.alt_reads_has("batch"):
        _fs_color_options.append("Batch")
    if ctx.alt_reads_has("label1"):
        _fs_color_options.append("Label 1")
    if ctx.alt_reads_has("label2"):
        _fs_color_options.append("Label 2")
    if ctx.alt_reads_has("label3"):
        _fs_color_options.append("Label 3")
    if ctx.alt_reads_has("timepoint"):
        _fs_color_options.append("Timepoint")
    _fs_color_by = _fs_ctrl_col1.radio(
        "Color by", _fs_color_options,
        horizontal=True, key="fs_color_by",
    )
    _fs_y_mode = _fs_ctrl_col2.radio(
        "Y axis", ["Fraction", "Count"],
        horizontal=True, key="fs_y_mode",
    )
    _fs_by_sample = _fs_color_by == "Sample"
    _fs_by_batch  = _fs_color_by == "Batch"
    _fs_by_label  = _fs_color_by in ("Label 1", "Label 2", "Label 3", "Timepoint")
    _fs_lbl_col   = {"Label 1": "label1", "Label 2": "label2", "Label 3": "label3", "Timepoint": "timepoint"}.get(_fs_color_by)
    _fs_normalize = _fs_y_mode == "Fraction"

    _fs_group_col = (
        "ar.sample_id"        if _fs_by_sample else
        "ar.batch"            if _fs_by_batch  else
        f"ar.{_fs_lbl_col}"   if _fs_by_label  else
        None
    )
    _fs_select = f"{_fs_group_col}, " if _fs_group_col else ""
    _fs_group  = f"{_fs_group_col}, " if _fs_group_col else ""

    _fs_df = ctx.con.execute(f"""
        SELECT {_fs_select}ar.family_size, COUNT(*) AS n_reads
        FROM {ctx.r_join}
        WHERE ar.family_size IS NOT NULL
        GROUP BY {_fs_group}ar.family_size
        ORDER BY {_fs_group}ar.family_size
    """).df()

    if _fs_df.empty:
        st.info("No family size data (fgbio cD tag absent in this dataset).")
        return

    _fs_x_min = int(_fs_df["family_size"].min())
    _fs_x_max = int(_fs_df["family_size"].max())
    if _fs_x_min < _fs_x_max:
        _fs_x_range = st.slider(
            "X-axis range (family size)",
            min_value=_fs_x_min,
            max_value=_fs_x_max,
            value=(_fs_x_min, _fs_x_max),
            key="fs_x_range",
        )
        _fs_df = _fs_df[
            (_fs_df["family_size"] >= _fs_x_range[0])
            & (_fs_df["family_size"] <= _fs_x_range[1])
        ]

    _fs_label_col = (
        "sample_id"  if _fs_by_sample else
        "batch"      if _fs_by_batch  else
        _fs_lbl_col  if _fs_by_label  else
        None
    )
    if _fs_normalize:
        if _fs_label_col:
            _fs_df["y_val"] = _fs_df.groupby(_fs_label_col)["n_reads"].transform(
                lambda x: x / x.sum()
            )
        else:
            _fs_df["y_val"] = _fs_df["n_reads"] / _fs_df["n_reads"].sum()
        _y_field = "y_val:Q"
        _y_title = "Fraction of alt-supporting reads"
        _y_fmt = ".3f"
    else:
        _fs_df["y_val"] = _fs_df["n_reads"]
        _y_field = "y_val:Q"
        _y_title = "Alt-supporting reads"
        _y_fmt = "d"

    _fs_enc = dict(
        x=alt.X("family_size:Q", title="Family size (cD tag)", bin=False),
        y=alt.Y(_y_field, title=_y_title),
        tooltip=[
            *([f"{_fs_label_col}:N"] if _fs_label_col else []),
            "family_size:Q",
            alt.Tooltip(_y_field, format=_y_fmt, title=_y_title),
        ],
    )
    if _fs_label_col:
        _fs_enc["color"] = alt.Color(f"{_fs_label_col}:N", scale=alt.Scale(scheme="tableau10"))
        _fs_chart = (
            alt.Chart(_fs_df)
            .mark_line(point=True, opacity=0.8)
            .encode(**_fs_enc)
            .properties(height=300)
        )
    else:
        _fs_chart = (
            alt.Chart(_fs_df)
            .mark_line(point=True, color="#4c78a8")
            .encode(**_fs_enc)
            .properties(height=300)
        )
    st.altair_chart(_fs_chart, width="stretch")
    _fs_norm_note = (
        "Fraction mode normalizes each batch independently."
        if _fs_by_batch else
        "Fraction mode normalizes each sample independently, allowing shape comparison across samples with different read counts."
    )
    st.caption(f"Artefacts are enriched in singletons (family_size = 1). {_fs_norm_note}")


def _render_family_size_vs_vaf(ctx: TabContext) -> None:
    st.subheader("Family size vs VAF (per locus)")
    _fsvaf_df = ctx.con.execute(f"""
        SELECT
            ab.sample_id,
            ab.chrom,
            ab.pos,
            ab.alt_allele,
            ROUND(ab.alt_count * 1.0 / ab.total_depth, 4) AS vaf,
            AVG(ar.family_size) AS mean_family_size,
            COUNT(*)            AS n_reads
        FROM (SELECT * FROM {ctx.table_expr} WHERE {ctx.where}) ab
        INNER JOIN alt_reads ar
            ON  ab.sample_id  = ar.sample_id
            AND ab.chrom      = ar.chrom
            AND ab.pos        = ar.pos
            AND ab.alt_allele = ar.alt_allele
        WHERE ar.family_size IS NOT NULL
        GROUP BY ab.sample_id, ab.chrom, ab.pos, ab.alt_allele, ab.alt_count, ab.total_depth
    """).df()

    if _fsvaf_df.empty:
        st.info("No data with family size available for the current selection.")
        return

    _fsvaf_plot_df = (
        _fsvaf_df
        .sort_values(["sample_id", "chrom", "pos", "alt_allele"])
        .head(3000)
        .reset_index(drop=True)
    )
    _fsvaf_sel = alt.selection_point(
        name="fsvaf_select",
        fields=["sample_id", "chrom", "pos", "alt_allele"],
        on="click",
        toggle="event.shiftKey",
    )
    _fsvaf_chart = (
        alt.Chart(_fsvaf_plot_df)
        .mark_point(filled=True, size=60)
        .encode(
            alt.X("mean_family_size:Q", title="Mean family size of alt reads"),
            alt.Y("vaf:Q", title="VAF", scale=alt.Scale(domain=[0, 1])),
            alt.Color("sample_id:N", title="Sample",
                      scale=alt.Scale(scheme="tableau10")),
            opacity=alt.condition(_fsvaf_sel, alt.value(1.0), alt.value(0.35)),
            size=alt.condition(_fsvaf_sel, alt.value(110), alt.value(55)),
            tooltip=[
                alt.Tooltip("sample_id:N"),
                alt.Tooltip("chrom:N"),
                alt.Tooltip("pos:Q"),
                alt.Tooltip("alt_allele:N"),
                alt.Tooltip("vaf:Q", format=".4f"),
                alt.Tooltip("mean_family_size:Q", format=".1f", title="Mean family size"),
                alt.Tooltip("n_reads:Q", title="Alt reads"),
            ],
        )
        .add_params(_fsvaf_sel)
        .properties(height=350)
    )
    _fsvaf_event = st.altair_chart(
        _fsvaf_chart,
        width="stretch",
        on_select="rerun",
        key="fsvaf_scatter",
    )
    st.caption(
        "True low-VAF variants should have reasonable mean family sizes. "
        "Artefacts at low VAF tend to cluster at low family size."
    )
    _fsvaf_pts = (_fsvaf_event.selection or {}).get("fsvaf_select", [])
    if _fsvaf_pts:
        _fsvaf_or = " OR ".join(
            f"(sample_id = '{ctx.sql_str(p['sample_id'])}' AND chrom = '{ctx.sql_str(p['chrom'])}' "
            f"AND pos = {int(p['pos'])} AND alt_allele = '{ctx.sql_str(p['alt_allele'])}')"
            for p in _fsvaf_pts
            if all(k in p for k in ["sample_id", "chrom", "pos", "alt_allele"])
        )
        if _fsvaf_or:
            _fsvaf_sel_df = ctx.con.execute(f"""
                SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                       pos + 1 AS pos_display
                FROM {ctx.table_expr}
                WHERE ({_fsvaf_or})
                ORDER BY sample_id, chrom, pos, alt_allele
            """).df()
            st.caption(
                f"{len(_fsvaf_sel_df)} selected loci across "
                f"{_fsvaf_sel_df['sample_id'].nunique()} sample(s) — "
                "shift-click to select multiple"
            )
            st.dataframe(_fsvaf_sel_df[ctx.table_cols], width="stretch", hide_index=True)
            ctx.igv_buttons(
                [f"({_fsvaf_or})"],
                _fsvaf_sel_df,
                key=f"fsvaf_{'_'.join(str(int(p['pos'])) for p in _fsvaf_pts[:5] if 'pos' in p)}",
            )
    else:
        st.caption("Click a point to select it; shift-click to select multiple.")


def _render_family_size_stratified_spectrum(ctx: TabContext) -> None:
    if not ctx.reads.fs_has_data:
        return

    st.divider()
    st.subheader("Family-size stratified Spectrum")
    st.caption(
        "Singleton reads (family_size = 1) vs multi-member families (family_size > 1). "
        "Singletons are enriched for sequencing errors; differences in profile shape "
        "reveal the true variant signal from the error process."
    )

    _locus_fs_filter = f"AND {ctx.reads_where}" if ctx.reads.active else ""
    _fs_strat_raw = ctx.con.execute(f"""
        WITH locus_fs AS (
            SELECT sample_id, chrom, pos, alt_allele,
                   MEDIAN(family_size) AS median_fs
            FROM alt_reads
            WHERE family_size IS NOT NULL {_locus_fs_filter}
            GROUP BY sample_id, chrom, pos, alt_allele
        )
        SELECT
            _t.trinuc_context, _t.ref_allele, _t.alt_allele,
            CASE WHEN COALESCE(lfs.median_fs, 1) <= 1
                 THEN 'singleton' ELSE 'multi' END AS fs_group,
            COUNT(*) AS count
        FROM (SELECT * FROM {ctx.table_expr} WHERE {ctx.where}) _t
        LEFT JOIN locus_fs lfs
            ON  lfs.sample_id  = _t.sample_id
            AND lfs.chrom      = _t.chrom
            AND lfs.pos        = _t.pos
            AND lfs.alt_allele = _t.alt_allele
        WHERE _t.variant_type = 'SNV'
          AND _t.trinuc_context IS NOT NULL
          AND length(_t.trinuc_context) = 3
        GROUP BY _t.trinuc_context, _t.ref_allele, _t.alt_allele, fs_group
    """).df()

    _fs_sing_raw  = _fs_strat_raw[_fs_strat_raw["fs_group"] == "singleton"].drop(columns="fs_group")
    _fs_multi_raw = _fs_strat_raw[_fs_strat_raw["fs_group"] == "multi"].drop(columns="fs_group")

    _fs_sing_s96,  _n_sing  = to_spec96_strat(_fs_sing_raw)
    _fs_multi_s96, _n_multi = to_spec96_strat(_fs_multi_raw)

    _fc1, _fc2 = st.columns(2)
    with _fc1:
        if _fs_sing_s96 is not None:
            st.altair_chart(
                strat_sbs96_chart(_fs_sing_s96, f"Singleton (family_size = 1, n={_n_sing:,})"),
                width="stretch",
            )
        else:
            st.info("No singleton loci in current selection.")
    with _fc2:
        if _fs_multi_s96 is not None:
            st.altair_chart(
                strat_sbs96_chart(_fs_multi_s96, f"Multi-member (family_size > 1, n={_n_multi:,})"),
                width="stretch",
            )
        else:
            st.info("No multi-member loci in current selection.")


def _render_cohort_artefact_family_size(ctx: TabContext) -> None:
    if not ctx.path.endswith(".duckdb"):
        return

    st.divider()
    st.subheader("Cohort artefact vs rare variant: family size comparison")
    _n_samples_total = ctx.con.execute(
        f"SELECT COUNT(DISTINCT sample_id) FROM {ctx.table_expr} WHERE {ctx.where}"
    ).fetchone()[0]

    if _n_samples_total < 2:
        st.info("Need at least 2 samples for cohort artefact comparison.")
        return

    _cohort_fs_df = ctx.con.execute(f"""
        WITH _base AS (
            SELECT sample_id,
                   chrom,
                   CAST(pos AS BIGINT) AS pos,
                   alt_allele
            FROM {ctx.table_expr}
            WHERE {ctx.where}
        ),
        locus_counts AS (
            SELECT chrom, pos, alt_allele,
                   COUNT(DISTINCT sample_id) AS n_samples_with_alt
            FROM _base
            GROUP BY chrom, pos, alt_allele
        ),
        labeled AS (
            SELECT
                CASE
                    WHEN lc.n_samples_with_alt = 1 THEN 'Seen in 1 sample'
                    WHEN lc.n_samples_with_alt <= 3 THEN 'Seen in 2–3 samples'
                    ELSE 'Seen in 4+ samples'
                END AS cohort_freq,
                ar.family_size
            FROM alt_reads ar
            INNER JOIN _base _filt
                ON  ar.sample_id  = _filt.sample_id
                AND ar.chrom      = _filt.chrom
                AND ar.pos        = _filt.pos
                AND ar.alt_allele = _filt.alt_allele
            INNER JOIN locus_counts lc
                ON  ar.chrom      = lc.chrom
                AND ar.pos        = lc.pos
                AND ar.alt_allele = lc.alt_allele
            WHERE ar.family_size IS NOT NULL
        )
        SELECT
            cohort_freq,
            PERCENTILE_CONT(0.0)  WITHIN GROUP (ORDER BY family_size) AS min_val,
            PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY family_size) AS q1,
            PERCENTILE_CONT(0.5)  WITHIN GROUP (ORDER BY family_size) AS median,
            PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY family_size) AS q3,
            PERCENTILE_CONT(1.0)  WITHIN GROUP (ORDER BY family_size) AS max_val,
            COUNT(*) AS n_reads
        FROM labeled
        GROUP BY cohort_freq
    """).df()

    if _cohort_fs_df.empty:
        st.info("No family size data for cohort comparison.")
        return

    _cf_order = ["Seen in 1 sample", "Seen in 2–3 samples", "Seen in 4+ samples"]
    _cf_box = (
        alt.Chart(_cohort_fs_df)
        .mark_bar(size=40)
        .encode(
            alt.X("cohort_freq:N", title="Cohort frequency", sort=_cf_order),
            alt.Y("q1:Q", title="Family size", scale=alt.Scale(zero=True)),
            alt.Y2("q3:Q"),
            alt.Color("cohort_freq:N", legend=None,
                      scale=alt.Scale(scheme="tableau10")),
            tooltip=[
                alt.Tooltip("cohort_freq:N", title="Group"),
                alt.Tooltip("min_val:Q", title="Min"),
                alt.Tooltip("q1:Q", title="Q1"),
                alt.Tooltip("median:Q", title="Median"),
                alt.Tooltip("q3:Q", title="Q3"),
                alt.Tooltip("max_val:Q", title="Max"),
                alt.Tooltip("n_reads:Q", title="Reads"),
            ],
        )
    )
    _cf_whisker = (
        alt.Chart(_cohort_fs_df)
        .mark_rule()
        .encode(
            alt.X("cohort_freq:N", sort=_cf_order),
            alt.Y("min_val:Q"),
            alt.Y2("max_val:Q"),
            alt.Color("cohort_freq:N", legend=None,
                      scale=alt.Scale(scheme="tableau10")),
        )
    )
    _cf_median = (
        alt.Chart(_cohort_fs_df)
        .mark_tick(color="white", thickness=2, size=40)
        .encode(
            alt.X("cohort_freq:N", sort=_cf_order),
            alt.Y("median:Q"),
        )
    )
    _cf_chart = (_cf_whisker + _cf_box + _cf_median).properties(height=300)
    st.altair_chart(_cf_chart, width="stretch")
    st.caption(
        "Cohort artefacts (seen in many samples) tend to have lower family sizes "
        "than rare variants, confirming they are sequencing noise rather than "
        "recurrent true variants."
    )


def _render_ab_ba_strand_counts(ctx: TabContext) -> None:
    _ab_has_data = ctx.con.execute(
        "SELECT COUNT(*) FROM alt_reads WHERE ab_count IS NOT NULL LIMIT 1"
    ).fetchone()[0] > 0
    if not _ab_has_data:
        return

    st.divider()
    st.subheader("AB vs BA strand counts")
    st.caption(
        "`ab_count` (aD tag) and `ba_count` (bD tag) are the number of raw reads from "
        "the AB (top) and BA (bottom) strands that contributed to each consensus read. "
        "Points on the diagonal indicate balanced duplex support from both strands. "
        "Points on the axes (ba_count = 0 or ab_count = 0) came from only one strand."
    )

    _ab_heat_df = ctx.con.execute(f"""
        SELECT
            ar.ab_count,
            ar.ba_count,
            COUNT(*) AS n_reads
        FROM {ctx.r_join}
        WHERE ar.ab_count IS NOT NULL AND ar.ba_count IS NOT NULL
        GROUP BY ar.ab_count, ar.ba_count
    """).df()

    if _ab_heat_df.empty:
        st.info("No AB/BA data in current selection.")
        return

    _ab_sel = alt.selection_point(
        name="ab_ba_click", fields=["ab_count", "ba_count"], on="click"
    )
    _ab_heat_chart = (
        alt.Chart(_ab_heat_df)
        .mark_rect()
        .encode(
            alt.X("ab_count:O", title="AB strand count (aD tag)"),
            alt.Y("ba_count:O", title="BA strand count (bD tag)"),
            alt.Color("n_reads:Q", title="Alt-supporting reads",
                      scale=alt.Scale(scheme="blues")),
            opacity=alt.condition(_ab_sel, alt.value(1.0), alt.value(0.4)),
            tooltip=[
                alt.Tooltip("ab_count:O", title="AB count"),
                alt.Tooltip("ba_count:O", title="BA count"),
                alt.Tooltip("n_reads:Q", title="Reads"),
            ],
        )
        .add_params(_ab_sel)
        .properties(height=400)
    )
    _ev_ab = st.altair_chart(
        _ab_heat_chart, width="stretch", on_select="rerun",
        key="ab_ba_heatmap",
    )
    st.caption(
        "Colour intensity = number of alt-supporting reads with each (ab_count, ba_count) combination. "
        "Reads on the diagonal have balanced strand support; reads on the axes came from one strand only. "
        "Click a cell to drill down."
    )

    _pts_ab = (_ev_ab.selection or {}).get("ab_ba_click", [])
    if _pts_ab:
        _ab_or = " OR ".join(
            f"(ab_count = {int(p['ab_count'])} AND ba_count = {int(p['ba_count'])})"
            for p in _pts_ab
        )
        _ab_cond = (
            f"(sample_id, chrom, pos, alt_allele) IN ("
            f"SELECT sample_id, chrom, pos, alt_allele "
            f"FROM alt_reads WHERE {_ab_or})"
        )
        _ab_sel_df = ctx.query_records([_ab_cond])
        _pairs_str = ", ".join(
            f"({int(p['ab_count'])},{int(p['ba_count'])})"
            for p in _pts_ab
        )
        st.caption(f"{len(_ab_sel_df):,} loci · selected (AB, BA): {_pairs_str}")
        st.dataframe(_ab_sel_df[ctx.table_cols], width="stretch")
        _ab_key = "ab_ba_" + "_".join(
            f"{int(p['ab_count'])}x{int(p['ba_count'])}"
            for p in _pts_ab
        )
        ctx.igv_buttons([_ab_cond], _ab_sel_df, key=_ab_key)
