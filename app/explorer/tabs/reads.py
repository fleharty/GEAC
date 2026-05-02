"""Reads tab extracted from the explorer monolith."""
from __future__ import annotations

import altair as alt
import pandas as pd
import streamlit as st

from explorer.tab_context import TabContext
from read_context_helpers import compute_locus_n_asymmetry_precomputed

LABEL = "📖 Reads"


def render(ctx: TabContext) -> None:
    alt_reads_cols = ctx.alt_reads_cols
    has_alt_reads = ctx.has_alt_reads
    _alt_reads_has_pipeline = "pipeline" in alt_reads_cols
    _alt_reads_has_read_type = "read_type" in alt_reads_cols
    _alt_reads_has_batch = "batch" in alt_reads_cols
    _alt_reads_has_label1 = "label1" in alt_reads_cols
    _alt_reads_has_label2 = "label2" in alt_reads_cols
    _alt_reads_has_label3 = "label3" in alt_reads_cols
    _alt_reads_has_timepoint = "timepoint" in alt_reads_cols
    _alt_reads_has_frag_gc = "frag_gc" in alt_reads_cols

    def _has_alt_reads_cols(*cols: str) -> bool:
        return all(col in alt_reads_cols for col in cols)

    def _read_series_config(color_by: str, show_read_pair: bool, read_expr_name: str):
        """Return SQL fragments for the common sample/batch/label/read grouping modes."""
        by_sample = color_by == "Sample"
        by_batch = color_by == "Batch"
        by_label = color_by in ("Label 1", "Label 2", "Label 3", "Timepoint")
        lbl_col = {"Label 1": "label1", "Label 2": "label2", "Label 3": "label3", "Timepoint": "timepoint"}.get(color_by)
        if by_batch and show_read_pair:
            return _r_join, f"ar.batch || ' ' || {read_expr_name} AS label, ", "label, ", "label"
        if by_label and show_read_pair:
            return _r_join, f"ar.{lbl_col} || ' ' || {read_expr_name} AS label, ", "label, ", "label"
        if by_sample and show_read_pair:
            return _r_join, f"ar.sample_id || ' ' || {read_expr_name} AS label, ", "label, ", "label"
        if show_read_pair:
            return _r_join, f"{read_expr_name} AS read, ", "read, ", "read"
        if by_batch:
            return _r_join, "ar.batch, ", "ar.batch, ", "batch"
        if by_label:
            return _r_join, f"ar.{lbl_col}, ", f"ar.{lbl_col}, ", lbl_col
        if by_sample:
            return _r_join, "ar.sample_id, ", "ar.sample_id, ", "sample_id"
        return _r_join, "", "", None

    def _n_table_or(df: pd.DataFrame) -> str:
        return " OR ".join(
            f"(sample_id = '{_sql_str(r['sample_id'])}' AND chrom = '{_sql_str(r['chrom'])}' "
            f"AND pos = {int(r['pos'])} AND alt_allele = '{_sql_str(r['alt_allele'])}')"
            for _, r in df.iterrows()
        )

    def _nctx_cache_key() -> tuple[str, str, str]:
        return ("nctx", where, _r_reads_filter)

    def _nasym_cache_key(use_fast: bool) -> tuple[str, str, str, bool]:
        return ("nasym_df", where, _r_reads_filter, use_fast)

    def _nrich_cache_key(min_frac: float) -> tuple[str, str, str, float]:
        return ("nrich_df", where, _r_reads_filter, min_frac)

    def _cache_get(key: str, default=None):
        return st.session_state.get(key, default)

    def _cache_set(key: str, value) -> None:
        st.session_state[key] = value

    def _n_flags() -> dict[str, object]:
        return {
            "min_delta_n_frac": st.session_state.get("min_delta_n_frac", 0.0),
            "nctx_enabled": st.session_state.get("nctx_enabled", False),
            "nasym_enabled": st.session_state.get("nasym_enabled", False),
            "nrich_enabled": st.session_state.get("nrich_enabled", False),
            "_use_fast_nasym": st.session_state.get("_use_fast_nasym", False),
        }

    def _insert_size_state() -> dict[str, object]:
        return {
            "_cached_has_insert_size": st.session_state.get("_cached_has_insert_size"),
            "_read_len_median": st.session_state.get("_cached_read_len_median"),
            "_cached_has_frag_gc": st.session_state.get("_cached_has_frag_gc"),
        }

    _IS_MIN, _IS_MAX = 20, 500
    _r_join = ctx.r_join
    _r_reads_filter = ctx.reads_where
    _reads_active = ctx.reads.active
    _fs_has_data = ctx.reads.fs_has_data
    _fs_lo = ctx.reads.fs_lo
    _fs_hi = ctx.reads.fs_hi
    _fs_max = ctx.reads.fs_max
    _mq_lo = ctx.reads.mq_lo
    _mq_hi = ctx.reads.mq_hi
    _mq_max = ctx.reads.mq_max
    _is_lo = ctx.reads.is_lo
    _is_hi = ctx.reads.is_hi
    _gc_has_data = ctx.reads.gc_has_data
    _gc_lo = ctx.reads.gc_lo
    _gc_hi = ctx.reads.gc_hi
    _n_state = _n_flags()
    min_delta_n_frac = _n_state["min_delta_n_frac"]
    nctx_enabled = _n_state["nctx_enabled"]
    nasym_enabled = _n_state["nasym_enabled"]
    nrich_enabled = _n_state["nrich_enabled"]
    _use_fast_nasym = _n_state["_use_fast_nasym"]
    _density_color = st.session_state.get("density_color", "All samples (aggregate)")
    table_expr = ctx.table_expr
    where = ctx.where
    con = ctx.con
    _sql_str = ctx.sql_str
    igv_buttons = ctx.igv_buttons
    _schema_cols = ctx.schema_cols
    _TIMINGS = []
    _timed = ctx.timed
    _ins_state = _insert_size_state()
    _cached_has_insert_size = _ins_state["_cached_has_insert_size"]
    _has_alt_reads = has_alt_reads
    _read_len_median = _ins_state["_read_len_median"]

    if not _has_alt_reads:
        st.info(
            "Per-read detail table not available. "
            "Re-run `geac collect` with `--reads-output` and merge the resulting "
            "`.reads.parquet` files to enable this tab."
        )
    else:
        st.caption(
            "All plots reflect alt-supporting reads linked to loci that pass the current sidebar filters."
        )

        # ── Row 1: Read position bias ──────────────────────────────────────────
        with st.expander("Read position bias", expanded=True):
            _dfe_ctrl1, _dfe_ctrl2, _dfe_ctrl3 = st.columns([3, 2, 1])
            _dfe_color_options = ["All samples (aggregate)", "Sample"]
            if _alt_reads_has_batch:
                _dfe_color_options.append("Batch")
            if _alt_reads_has_pipeline:
                _dfe_color_options.append("Pipeline")
            if _alt_reads_has_label1:
                _dfe_color_options.append("Label 1")
            if _alt_reads_has_label2:
                _dfe_color_options.append("Label 2")
            if _alt_reads_has_label3:
                _dfe_color_options.append("Label 3")
            if _alt_reads_has_timepoint:
                _dfe_color_options.append("Timepoint")
            _dfe_color_by = _dfe_ctrl1.radio(
                "Color by", _dfe_color_options,
                horizontal=True, key="dfe_color_by",
            )
            _dfe_y_mode = _dfe_ctrl2.radio(
                "Y axis", ["Fraction", "Count"],
                horizontal=True, key="dfe_y_mode",
            )
            _dfe_by_read   = _dfe_ctrl3.checkbox("Show R1/R2", value=False, key="dfe_show_r1r2")
            _dfe_normalize = _dfe_y_mode == "Fraction"
            _DFE_READ_EXPR = "CASE WHEN ar.is_read1 THEN 'R1' ELSE 'R2' END"
            _dfe_source, _dfe_select_expr, _dfe_group_expr, _dfe_label_col = _read_series_config(
                _dfe_color_by, _dfe_by_read, _DFE_READ_EXPR
            )

            _dfe_df = con.execute(f"""
                SELECT {_dfe_select_expr}ar.cycle, COUNT(*) AS n_reads
                FROM {_dfe_source}
                GROUP BY {_dfe_group_expr}ar.cycle
                ORDER BY {_dfe_group_expr}ar.cycle
            """).df()

            if _dfe_df.empty:
                st.info("No data.")
            else:
                if _dfe_normalize:
                    if _dfe_label_col:
                        _dfe_df["y_val"] = _dfe_df.groupby(_dfe_label_col)["n_reads"].transform(
                            lambda x: x / x.sum()
                        )
                    else:
                        _dfe_df["y_val"] = _dfe_df["n_reads"] / _dfe_df["n_reads"].sum()
                    _dfe_y_field = "y_val:Q"
                    _dfe_y_title = "Fraction of alt-supporting reads"
                    _dfe_y_fmt = ".3f"
                else:
                    _dfe_df["y_val"] = _dfe_df["n_reads"]
                    _dfe_y_field = "y_val:Q"
                    _dfe_y_title = "Alt-supporting reads"
                    _dfe_y_fmt = "d"

                _dfe_enc = dict(
                    x=alt.X("cycle:Q", title="Cycle (1-based)", bin=False),
                    y=alt.Y(_dfe_y_field, title=_dfe_y_title),
                    tooltip=[
                        *([f"{_dfe_label_col}:N"] if _dfe_label_col else []),
                        alt.Tooltip("cycle:Q", title="Cycle"),
                        alt.Tooltip(_dfe_y_field, format=_dfe_y_fmt, title=_dfe_y_title),
                    ],
                )
                if _dfe_label_col:
                    _dfe_enc["color"] = alt.Color(f"{_dfe_label_col}:N", scale=alt.Scale(scheme="tableau10"))
                    _dfe_chart = (
                        alt.Chart(_dfe_df)
                        .mark_line(point=True, opacity=0.8)
                        .encode(**_dfe_enc)
                        .properties(height=300)
                    )
                else:
                    _dfe_chart = (
                        alt.Chart(_dfe_df)
                        .mark_line(point=True, color="#f58518")
                        .encode(**_dfe_enc)
                        .properties(height=300)
                    )
                st.altair_chart(_dfe_chart, width="stretch")
                st.caption(
                    "A spike at high cycle numbers indicates alt-supporting reads clustered at read ends — "
                    "a red flag for alignment artefacts or damaged bases."
                )

        # ── Row 2: Base qual vs dist from read end scatter ────────────────────
        with st.expander("Mean base quality by cycle", expanded=True):
            _bq_ctrl1, _bq_ctrl2 = st.columns([4, 1])
            _bq_color_options = ["All samples (aggregate)", "Sample"]
            if _alt_reads_has_batch:
                _bq_color_options.append("Batch")
            if _alt_reads_has_label1:
                _bq_color_options.append("Label 1")
            if _alt_reads_has_label2:
                _bq_color_options.append("Label 2")
            if _alt_reads_has_label3:
                _bq_color_options.append("Label 3")
            if _alt_reads_has_timepoint:
                _bq_color_options.append("Timepoint")
            _bq_color_by = _bq_ctrl1.radio(
                "Color by", _bq_color_options,
                horizontal=True, key="bq_color_by",
            )
            _bq_by_read   = _bq_ctrl2.checkbox("Show R1/R2", value=False, key="bq_show_r1r2")
            _BQ_READ_EXPR = "CASE WHEN ar.is_read1 THEN 'R1' ELSE 'R2' END"
            _bq_source, _bq_select_expr, _bq_group_expr, _bq_label_col = _read_series_config(
                _bq_color_by, _bq_by_read, _BQ_READ_EXPR
            )

            _bq_df = con.execute(f"""
                SELECT
                    {_bq_select_expr}ar.cycle,
                    ROUND(AVG(ar.base_qual), 2) AS mean_base_qual,
                    COUNT(*) AS n_reads
                FROM {_bq_source}
                GROUP BY {_bq_group_expr}ar.cycle
                ORDER BY {_bq_group_expr}ar.cycle
            """).df()

            if _bq_df.empty:
                st.info("No data.")
            else:
                _bq_enc = dict(
                    x=alt.X("cycle:Q", title="Cycle (1-based)"),
                    y=alt.Y("mean_base_qual:Q", title="Mean base quality (Phred)",
                            scale=alt.Scale(zero=False)),
                    tooltip=[
                        *([f"{_bq_label_col}:N"] if _bq_label_col else []),
                        alt.Tooltip("cycle:Q", title="Cycle"),
                        alt.Tooltip("mean_base_qual:Q", format=".1f", title="Mean base qual"),
                        alt.Tooltip("n_reads:Q", title="Reads"),
                    ],
                )
                if _bq_label_col:
                    _bq_enc["color"] = alt.Color(f"{_bq_label_col}:N", scale=alt.Scale(scheme="tableau10"))
                    _bq_chart = (
                        alt.Chart(_bq_df)
                        .mark_line(point=True, opacity=0.8)
                        .encode(**_bq_enc)
                        .properties(height=350)
                    )
                else:
                    _bq_chart = (
                        alt.Chart(_bq_df)
                        .mark_line(point=True, color="#f58518")
                        .encode(**_bq_enc)
                        .properties(height=350)
                    )
                st.altair_chart(_bq_chart, width="stretch")
                st.caption(
                    "A drop in mean base quality at high cycle numbers (late in the read) "
                    "indicates that alt-supporting reads at those positions may be artefacts. "
                    "Quality is averaged arithmetically over Phred scores (not over error probabilities)."
                )

        # ── Row 3: N context around alt-supporting reads ──────────────────────
        st.subheader("N context around alt-supporting reads")
        if not _has_alt_reads_cols(
            "n_before_alt",
            "n_after_alt",
            "n_n_before_alt",
            "n_n_after_alt",
            "leading_n_run_len",
            "trailing_n_run_len",
        ):
            st.info(
                "N-context analysis requires a newer `alt_reads` schema. "
                "Re-run `geac collect --reads-output` and `geac merge` to enable these plots."
            )
        else:
            _nctx_enabled = st.checkbox(
                "Enable N-context read distribution",
                value=False,
                key="nctx_enabled",
                help="Scans all alt-supporting reads to build N-context distributions. "
                     "Can be slow on large cohorts — enable when you need it.",
            )
            if not _nctx_enabled:
                st.info(
                    "N-context read distribution is disabled by default. "
                    "Tick the checkbox above to run it."
                )
            else:
                _nctx_color_mode = st.radio(
                    "Group by",
                    ["All reads", "R1/R2"],
                    horizontal=True,
                    key="nctx_group_by",
                )
                # Cache on the filter strings — this query scans alt_reads and runs
                # on every rerun because Streamlit evaluates tab bodies eagerly.
                # Pre-aggregate in DuckDB to avoid sending millions of read-level
                # rows to the browser. Two small queries replace the former full
                # alt_reads fetch: one for scalar metrics, one for histogram bins
                # used to draw the density curve (~200 rows max).
                _nctx_cache_key = _nctx_cache_key()
                if st.session_state.get("_nctx_cache_key") != _nctx_cache_key:
                    st.session_state["_nctx_cache_key"] = _nctx_cache_key
                    with _timed("nctx scalars [cache miss]"):
                        st.session_state["_nctx_scalars_cache"] = con.execute(f"""
                            SELECT
                                AVG(CASE WHEN ar.n_before_alt > 0
                                         THEN ar.n_n_before_alt * 1.0 / ar.n_before_alt
                                         ELSE NULL END)                              AS mean_frac_n_before,
                                AVG(CASE WHEN ar.n_after_alt  > 0
                                         THEN ar.n_n_after_alt  * 1.0 / ar.n_after_alt
                                         ELSE NULL END)                              AS mean_frac_n_after,
                                AVG(CASE WHEN ar.n_n_before_alt > 0 THEN 1.0 ELSE 0.0 END) AS frac_any_n_before,
                                AVG(CASE WHEN ar.n_n_after_alt  > 0 THEN 1.0 ELSE 0.0 END) AS frac_any_n_after
                            FROM {_r_join}
                        """).fetchone()
                    with _timed("nctx density bins [cache miss]"):
                        st.session_state["_nctx_density_cache"] = con.execute(f"""
                            WITH sides AS (
                                SELECT
                                    CASE WHEN ar.is_read1 THEN 'R1' ELSE 'R2' END AS read_group,
                                    'Before alt'                                   AS side,
                                    CASE WHEN ar.n_before_alt > 0
                                         THEN ar.n_n_before_alt * 1.0 / ar.n_before_alt
                                         ELSE NULL END                             AS frac_n
                                FROM {_r_join}
                                UNION ALL
                                SELECT
                                    CASE WHEN ar.is_read1 THEN 'R1' ELSE 'R2' END AS read_group,
                                    'After alt'                                    AS side,
                                    CASE WHEN ar.n_after_alt > 0
                                         THEN ar.n_n_after_alt * 1.0 / ar.n_after_alt
                                         ELSE NULL END                             AS frac_n
                                FROM {_r_join}
                            ),
                            binned AS (
                                SELECT
                                    side, read_group,
                                    ROUND(FLOOR(frac_n / 0.001) * 0.001, 4) AS frac_n_bin,
                                    COUNT(*)                                  AS n
                                FROM sides
                                WHERE frac_n IS NOT NULL AND frac_n BETWEEN 0 AND 0.05
                                GROUP BY side, read_group, frac_n_bin
                            )
                            SELECT side, read_group, frac_n_bin AS frac_n, n
                            FROM binned
                            ORDER BY side, read_group, frac_n
                        """).df()
                else:
                    _TIMINGS.append(("nctx [cache hit]", 0.0))

                _nctx_scalars    = st.session_state["_nctx_scalars_cache"]
                _nctx_density_df = st.session_state["_nctx_density_cache"]

                if _nctx_density_df.empty and (
                    _nctx_scalars is None or all(v is None for v in _nctx_scalars)
                ):
                    st.info("No read-context data available under current filters.")
                else:
                    _nctx_mean_before, _nctx_mean_after, _nctx_any_before, _nctx_any_after = (
                        float(v) if v is not None else float("nan")
                        for v in _nctx_scalars
                    )
                    _nctx_mean_delta = _nctx_mean_after - _nctx_mean_before

                    _nctx_metrics = st.columns(4)
                    _nctx_metrics[0].metric("Mean frac N before", f"{_nctx_mean_before:.3f}")
                    _nctx_metrics[1].metric(
                        "Mean frac N after", f"{_nctx_mean_after:.3f}",
                        delta=f"{_nctx_mean_delta:+.3f} vs before",
                    )
                    _nctx_metrics[2].metric("Reads with any N before", f"{_nctx_any_before:.1%}")
                    _nctx_metrics[3].metric("Reads with any N after",  f"{_nctx_any_after:.1%}")

                    _nctx_any_df = pd.DataFrame([
                        {"side": "Before alt", "fraction": _nctx_any_before},
                        {"side": "After alt",  "fraction": _nctx_any_after},
                    ])
                    st.altair_chart(
                        alt.Chart(_nctx_any_df)
                        .mark_bar(opacity=0.85)
                        .encode(
                            x=alt.X("side:N", title=None),
                            y=alt.Y(
                                "fraction:Q",
                                title="Fraction of alt-supporting reads with any N",
                                scale=alt.Scale(domain=[0, 1]),
                            ),
                            color=alt.Color("side:N", title=None),
                            tooltip=[
                                "side:N",
                                alt.Tooltip("fraction:Q", title="Fraction", format=".3%"),
                            ],
                        )
                        .properties(title="Fraction of alt-supporting reads with any N", height=220),
                        width="stretch",
                    )

                    if not _nctx_density_df.empty:
                        # Build density curves from pre-binned counts (bin width = 0.001).
                        if _nctx_color_mode == "R1/R2":
                            _density_plot_df = _nctx_density_df.copy()
                            _totals = _density_plot_df.groupby(["side", "read_group"])["n"].transform("sum")
                            _density_plot_df["density"] = _density_plot_df["n"] / (_totals * 0.001)
                            _density_plot_df["group"] = (
                                _density_plot_df["side"] + " / " + _density_plot_df["read_group"]
                            )
                            _density_color = alt.Color("group:N", title=None)
                        else:
                            _density_plot_df = (
                                _nctx_density_df
                                .groupby(["side", "frac_n"])["n"]
                                .sum()
                                .reset_index()
                            )
                            _totals = _density_plot_df.groupby("side")["n"].transform("sum")
                            _density_plot_df["density"] = _density_plot_df["n"] / (_totals * 0.001)
                            _density_color = alt.Color("side:N", title=None)

                        st.altair_chart(
                            alt.Chart(_density_plot_df)
                            .mark_line(strokeWidth=3)
                            .encode(
                                x=alt.X(
                                    "frac_n:Q",
                                    title="Fraction of available bases that are N",
                                    scale=alt.Scale(domain=[0, 0.05]),
                                ),
                                y=alt.Y("density:Q", title="Density"),
                                color=_density_color,
                                tooltip=[
                                    alt.Tooltip("side:N"),
                                    alt.Tooltip("frac_n:Q", title="Fraction N", format=".3f"),
                                    alt.Tooltip("density:Q", title="Density", format=".3f"),
                                ],
                            )
                            .properties(title="Leading vs trailing N burden", height=320),
                            width="stretch",
                        )

                    st.caption(
                        "Each read is normalized by the number of available bases on that side of the alt. "
                        "For example, 1 `N` out of 1 base after the alt contributes 1.00, while 5 `N`s out of 100 "
                        "bases after the alt contributes 0.05. The curves show the density of those per-read fractions "
                        "before versus after the alt."
                    )

        # ── Row 3b: N-asymmetry locus discovery ───────────────────────────────
        st.subheader("N-asymmetry locus discovery")
        if not _has_alt_reads_cols(
            "n_before_alt", "n_after_alt", "n_n_before_alt", "n_n_after_alt"
        ):
            st.info(
                "N-asymmetry analysis requires a newer `alt_reads` schema. "
                "Re-run `geac collect --reads-output` and `geac merge` to enable this section."
            )
        else:
            _nasym_enabled = st.checkbox(
                "Enable N-asymmetry locus discovery",
                value=False,
                key="nasym_enabled",
                help="Computes per-locus N-asymmetry scores. Fast when precomputed columns are "
                     "present (v0.4.7+ cohorts); slower on older cohorts that require a "
                     "full alt_reads GROUP BY.",
            )
            if not _nasym_enabled:
                st.info(
                    "N-asymmetry locus discovery is disabled by default. "
                    "Tick the checkbox above to run it."
                )
            else:
                st.caption(
                    "Each point is one locus. The x-axis shows the mean fraction of upstream bases "
                    "that are N; the y-axis shows the same for downstream bases. Loci in the upper-left "
                    "have many Ns **after** the alt base and few before — a hallmark of sequencing "
                    "error enrichment at late cycles. Use the sidebar **Min N-asymmetry score** slider "
                    "to highlight the most extreme sites."
                )
                # Fast path: read precomputed per-locus N-asymmetry columns directly
                # from alt_bases (added in v0.4.7) when no read-level filters are
                # active. This avoids the multi-million-row alt_reads GROUP BY
                # entirely. Falls back to the slow path when read filters change
                # the per-read denominator.
                _use_fast_nasym = (
                    not _reads_active
                    and "mean_delta_n_frac" in _schema_cols
                    and "n_alt_reads_with_n_ctx" in _schema_cols
                )
                # Cache keyed on the actual filter strings — stable across rerenders,
                # unlike hash() which uses a randomised seed per process.
                _nasym_cache_key = _nasym_cache_key(_use_fast_nasym)
                if st.session_state.get("_nasym_cache_key") != _nasym_cache_key:
                    st.session_state["_nasym_cache_key"] = _nasym_cache_key
                    _nasym_label = "nasym_df fast [cache miss]" if _use_fast_nasym else "nasym_df slow GROUP BY [cache miss]"
                    with _timed(_nasym_label):
                        if _use_fast_nasym:
                            st.session_state["_nasym_df_cache"] = (
                                compute_locus_n_asymmetry_precomputed(con, table_expr, where)
                            )
                        else:
                            st.session_state["_nasym_df_cache"] = compute_locus_n_asymmetry(
                                con, _r_join
                            )
                    st.session_state["_nasym_locus_cache_key"] = None  # invalidate dependent caches
                    st.session_state["_nasym_tbl_cache_key"] = None
                else:
                    _TIMINGS.append(("nasym_df [cache hit]", 0.0))
                _nasym_df = st.session_state["_nasym_df_cache"]

                if _nasym_df.empty:
                    st.info("No loci available under current filters.")
                else:
                    # When the fast path is in use, alt_bases gives one row per
                    # (sample, locus, alt, pipeline), so we merge per-pipeline.
                    # The slow path collapses across pipelines via the alt_reads
                    # GROUP BY, so we merge on the locus key only.
                    _nasym_merge_keys = ["sample_id", "chrom", "pos", "alt_allele"]
                    if _use_fast_nasym:
                        _nasym_merge_keys.append("pipeline")
                    _nasym_locus_select_extra = ", pipeline" if _use_fast_nasym else ""

                    # Merge in VAF and variant_type — cached separately on the locus WHERE clause.
                    _nasym_locus_cache_key = ("nasym_locus", where, _use_fast_nasym)
                    if st.session_state.get("_nasym_locus_cache_key") != _nasym_locus_cache_key:
                        st.session_state["_nasym_locus_cache_key"] = _nasym_locus_cache_key
                        st.session_state["_nasym_locus_cache"] = con.execute(f"""
                            SELECT sample_id, chrom, pos, alt_allele, variant_type{_nasym_locus_select_extra},
                                   ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                                   alt_count, total_depth
                            FROM {table_expr}
                            WHERE {where}
                        """).df()
                    _nasym_locus = st.session_state["_nasym_locus_cache"]
                    _nasym_df = _nasym_df.merge(
                        _nasym_locus,
                        on=_nasym_merge_keys,
                        how="left",
                    )
                    _nasym_df = _nasym_df.dropna(subset=["mean_delta_n_frac"])

                    # Apply the sidebar threshold filter.
                    if min_delta_n_frac > 0.0:
                        _nasym_df = _nasym_df[_nasym_df["mean_delta_n_frac"] >= min_delta_n_frac]

                    if _nasym_df.empty:
                        st.info(
                            f"No loci with N-asymmetry score ≥ {min_delta_n_frac:.2f} under current filters. "
                            "Try lowering the sidebar threshold."
                        )
                    else:
                        # ── Scatter: frac_n_before vs frac_n_after ─────────────
                        _nasym_sel = alt.selection_point(
                            name="nasym_select",
                            fields=["sample_id", "chrom", "pos", "alt_allele"],
                            on="click",
                            toggle="event.shiftKey",
                        )
                        _nasym_diag_max = float(_nasym_df[["mean_frac_n_before", "mean_frac_n_after"]].max().max())
                        _nasym_diag = (
                            alt.Chart(
                                pd.DataFrame({"x": [0.0, _nasym_diag_max], "y": [0.0, _nasym_diag_max]})
                            )
                            .mark_line(color="gray", strokeDash=[4, 4], opacity=0.5)
                            .encode(x="x:Q", y="y:Q")
                        )
                        _nasym_scatter = (
                            alt.Chart(_nasym_df)
                            .mark_circle(size=50)
                            .encode(
                                x=alt.X(
                                    "mean_frac_n_before:Q",
                                    title="Mean frac N before alt",
                                    scale=alt.Scale(zero=True),
                                ),
                                y=alt.Y(
                                    "mean_frac_n_after:Q",
                                    title="Mean frac N after alt",
                                    scale=alt.Scale(zero=True),
                                ),
                                color=alt.Color("variant_type:N", title="Variant type"),
                                opacity=alt.condition(_nasym_sel, alt.value(1.0), alt.value(0.4)),
                                size=alt.condition(_nasym_sel, alt.value(150), alt.value(50)),
                                tooltip=[
                                    "sample_id", "chrom", "pos", "alt_allele", "variant_type",
                                    alt.Tooltip("vaf:Q", format=".4f"),
                                    alt.Tooltip("mean_frac_n_before:Q", title="Mean frac N before", format=".4f"),
                                    alt.Tooltip("mean_frac_n_after:Q", title="Mean frac N after", format=".4f"),
                                    alt.Tooltip("mean_delta_n_frac:Q", title="Delta (after−before)", format=".4f"),
                                    alt.Tooltip("frac_reads_asymmetric:Q", title="Frac reads asymmetric", format=".3f"),
                                    alt.Tooltip("n_alt_reads:Q", title="Alt reads"),
                                ],
                            )
                            .add_params(_nasym_sel)
                            .properties(height=380, title="N-asymmetry scatter (upper-left = high downstream N)")
                            .interactive()
                        )

                        # ── Histogram of mean_delta_n_frac ─────────────────────
                        _nasym_hist = (
                            alt.Chart(_nasym_df)
                            .mark_bar(opacity=0.8)
                            .encode(
                                x=alt.X(
                                    "mean_delta_n_frac:Q",
                                    bin=alt.Bin(maxbins=40),
                                    title="N-asymmetry score (mean frac N after − before)",
                                ),
                                y=alt.Y("count():Q", title="Loci"),
                                tooltip=[
                                    alt.Tooltip("mean_delta_n_frac:Q", title="Score", format=".3f"),
                                    alt.Tooltip("count():Q", title="Loci"),
                                ],
                            )
                            .properties(height=380, title="Distribution of N-asymmetry score")
                        )

                        _nasym_col_scatter, _nasym_col_hist = st.columns(2)
                        with _nasym_col_scatter:
                            _nasym_event = st.altair_chart(
                                _nasym_diag + _nasym_scatter,
                                on_select="rerun",
                                key="nasym_scatter",
                            )
                        with _nasym_col_hist:
                            st.altair_chart(_nasym_hist, use_container_width=True)

                        # ── Ranked table + IGV ─────────────────────────────────
                        st.subheader("Most N-asymmetric loci")
                        _nasym_display = (
                            _nasym_df
                            .sort_values("mean_delta_n_frac", ascending=False)
                            [["sample_id", "chrom", "pos", "alt_allele", "variant_type",
                              "vaf", "alt_count", "total_depth", "n_alt_reads",
                              "mean_frac_n_before", "mean_frac_n_after",
                              "mean_delta_n_frac", "frac_reads_asymmetric"]]
                            .reset_index(drop=True)
                        )
                        st.dataframe(_nasym_display, width="stretch", hide_index=True)

                        # Build a WHERE fragment covering all loci in the ranked table.
                        # Cache on the filter + threshold — this query rebuilds a huge
                        # OR clause and re-scans table_expr on every rerun otherwise.
                        _nasym_tbl_cache_key = (
                            "nasym_table_df", where, _r_reads_filter, float(min_delta_n_frac)
                        )
                        if (
                            st.session_state.get("_nasym_tbl_cache_key") != _nasym_tbl_cache_key
                            or "_nasym_tbl_df_cache" not in st.session_state
                        ):
                            _nasym_table_or = _n_table_or(_nasym_display)
                            st.session_state["_nasym_tbl_cache_key"] = _nasym_tbl_cache_key
                            st.session_state["_nasym_tbl_or_cache"] = _nasym_table_or
                            st.session_state["_nasym_tbl_df_cache"] = con.execute(f"""
                                SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
                                FROM {table_expr}
                                WHERE ({_nasym_table_or})
                                ORDER BY chrom, pos
                            """).df()
                        _nasym_table_or = st.session_state["_nasym_tbl_or_cache"]
                        _nasym_table_df = st.session_state["_nasym_tbl_df_cache"]
                        igv_buttons(
                            [f"({_nasym_table_or})"],
                            _nasym_table_df,
                            key="nasym_table_igv",
                            use_global_filters=False,
                        )

                        # ── IGV for scatter selection ──────────────────────────
                        _nasym_pts = (_nasym_event.selection or {}).get("nasym_select", [])
                        if _nasym_pts:
                            _nasym_or = " OR ".join(
                                f"(sample_id = '{_sql_str(p['sample_id'])}' AND chrom = '{_sql_str(p['chrom'])}' "
                                f"AND pos = {int(p['pos'])} AND alt_allele = '{_sql_str(p['alt_allele'])}')"
                                for p in _nasym_pts
                                if all(k in p for k in ["sample_id", "chrom", "pos", "alt_allele"])
                            )
                            if _nasym_or:
                                _nasym_sel_df = con.execute(f"""
                                    SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
                                    FROM {table_expr}
                                    WHERE ({_nasym_or})
                                """).df()
                                st.caption(
                                    f"{len(_nasym_pts)} loci selected in scatter — shift-click to add more."
                                )
                                igv_buttons(
                                    [f"({_nasym_or})"],
                                    _nasym_sel_df,
                                    key=f"nasym_scatter_{'_'.join(str(int(p['pos'])) for p in _nasym_pts[:5] if 'pos' in p)}",
                                    use_global_filters=False,
                                )
                        else:
                            st.caption(
                                "Click a point in the scatter to load just those loci in IGV. "
                                "Shift-click to select multiple."
                            )

        # ── Row 3c: N-rich supporting loci ────────────────────────────────────
        st.subheader("N-rich supporting loci")
        if not _has_alt_reads_cols("n_n_before_alt", "n_n_after_alt", "trailing_n_run_len"):
            st.info(
                "N-rich loci analysis requires a newer `alt_reads` schema. "
                "Re-run `geac collect --reads-output` and `geac merge` to enable this section."
            )
        else:
            _nrich_enabled = st.checkbox(
                "Enable N-rich loci analysis",
                value=False,
                key="nrich_enabled",
                help=(
                    "Identifies loci where alt-supporting reads are systematically enriched "
                    "for N bases — a signal that the apparent alt support may be dominated by "
                    "low-confidence or damaged read context. Runs a GROUP BY on alt_reads; "
                    "may be slow on large cohorts."
                ),
            )
            if not _nrich_enabled:
                st.info(
                    "N-rich loci analysis is disabled by default. "
                    "Tick the checkbox above to run it."
                )
            else:
                _nrich_min_frac = st.slider(
                    "Min fraction of alt reads with any N",
                    min_value=0.0,
                    max_value=1.0,
                    value=0.1,
                    step=0.05,
                    help=(
                        "Show only loci where at least this fraction of alt-supporting reads "
                        "carry one or more N bases in the tracked read context "
                        "(n_n_before_alt + n_n_after_alt > 0)."
                    ),
                )

                _nrich_cache_key = _nrich_cache_key(_nrich_min_frac)
                if st.session_state.get("_nrich_cache_key") != _nrich_cache_key:
                    st.session_state["_nrich_cache_key"] = _nrich_cache_key
                    st.session_state["_nrich_locus_cache_key"] = None  # invalidate dependent cache
                    with _timed("nrich_df GROUP BY [cache miss]"):
                        st.session_state["_nrich_df_cache"] = con.execute(f"""
                            SELECT
                                ar.sample_id,
                                ar.chrom,
                                ar.pos,
                                ar.alt_allele,
                                COUNT(*)                                                     AS n_alt_reads,
                                AVG(CASE WHEN ar.n_n_before_alt + ar.n_n_after_alt > 0
                                         THEN 1.0 ELSE 0.0 END)                             AS frac_reads_with_any_n,
                                AVG(CASE WHEN ar.trailing_n_run_len > 0
                                         THEN 1.0 ELSE 0.0 END)                             AS frac_reads_with_trailing_n,
                                AVG(ar.trailing_n_run_len)                                   AS mean_trailing_n_run_len,
                                AVG(CASE WHEN (ar.n_before_alt + ar.n_after_alt) > 0
                                         THEN (ar.n_n_before_alt + ar.n_n_after_alt) * 1.0
                                              / (ar.n_before_alt + ar.n_after_alt)
                                         ELSE NULL END)                                      AS mean_total_n_frac
                            FROM {_r_join}
                            GROUP BY ar.sample_id, ar.chrom, ar.pos, ar.alt_allele
                        """).df()
                else:
                    _TIMINGS.append(("nrich_df [cache hit]", 0.0))
                _nrich_df = st.session_state["_nrich_df_cache"]

                if _nrich_df.empty:
                    st.info("No loci available under current filters.")
                else:
                    # Merge in VAF, variant_type, alt_count, total_depth from the locus table.
                    _nrich_locus_cache_key = ("nrich_locus", where)
                    if st.session_state.get("_nrich_locus_cache_key") != _nrich_locus_cache_key:
                        st.session_state["_nrich_locus_cache_key"] = _nrich_locus_cache_key
                        st.session_state["_nrich_locus_cache"] = con.execute(f"""
                            SELECT sample_id, chrom, pos, alt_allele, variant_type,
                                   ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                                   alt_count, total_depth
                            FROM {table_expr}
                            WHERE {where}
                        """).df()
                    _nrich_df = _nrich_df.merge(
                        st.session_state["_nrich_locus_cache"],
                        on=["sample_id", "chrom", "pos", "alt_allele"],
                        how="left",
                    )

                    # Apply the threshold filter.
                    _nrich_df = _nrich_df[_nrich_df["frac_reads_with_any_n"] >= _nrich_min_frac]

                    if _nrich_df.empty:
                        st.info(
                            f"No loci where ≥{_nrich_min_frac:.0%} of alt reads carry N. "
                            "Try lowering the threshold."
                        )
                    else:
                        _nrich_display = (
                            _nrich_df
                            .sort_values("frac_reads_with_any_n", ascending=False)
                            [["sample_id", "chrom", "pos", "alt_allele", "variant_type",
                              "vaf", "alt_count", "total_depth", "n_alt_reads",
                              "frac_reads_with_any_n", "frac_reads_with_trailing_n",
                              "mean_trailing_n_run_len", "mean_total_n_frac"]]
                            .reset_index(drop=True)
                        )
                        st.caption(
                            f"{len(_nrich_display):,} loci where ≥{_nrich_min_frac:.0%} of alt reads "
                            "carry N — sorted by fraction of reads with any N, descending."
                        )
                        st.dataframe(_nrich_display, width="stretch", hide_index=True)

                        # IGV buttons for the full N-rich loci set.
                        # Cache the OR clause and backing df on filters + threshold so it
                        # doesn't rebuild the huge OR string on every widget interaction.
                        _nrich_tbl_cache_key = (
                            "nrich_tbl", where, _r_reads_filter, _nrich_min_frac
                        )
                        if (
                            st.session_state.get("_nrich_tbl_cache_key") != _nrich_tbl_cache_key
                            or "_nrich_tbl_df_cache" not in st.session_state
                        ):
                            _nrich_table_or = _n_table_or(_nrich_display)
                            st.session_state["_nrich_tbl_cache_key"] = _nrich_tbl_cache_key
                            st.session_state["_nrich_tbl_or_cache"] = _nrich_table_or
                            st.session_state["_nrich_tbl_df_cache"] = con.execute(f"""
                                SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                                       pos + 1 AS pos_display
                                FROM {table_expr}
                                WHERE ({_nrich_table_or})
                                ORDER BY chrom, pos
                            """).df()
                        igv_buttons(
                            [f"({st.session_state['_nrich_tbl_or_cache']})"],
                            st.session_state["_nrich_tbl_df_cache"],
                            key="nrich_table_igv",
                            use_global_filters=False,
                        )

        # ── Row 4: Insert size distribution ───────────────────────────────────
        if "_cached_has_insert_size" not in st.session_state:
            st.session_state["_cached_has_insert_size"] = (
                con.execute("SELECT COUNT(*) FROM alt_reads WHERE insert_size IS NOT NULL LIMIT 1").fetchone()[0] > 0
            )
        if st.session_state["_cached_has_insert_size"]:
            # Median read length — used for gap correction in both insert size plots.
            # Fragments longer than 2×read_len have a coverage gap; the correction
            # weights each count by 1 / min(1, 2R/L) to recover the unbiased distribution.
            _read_len_median = con.execute(
                "SELECT MEDIAN(read_length) FROM alt_reads WHERE read_length IS NOT NULL"
            ).fetchone()[0] or 0
            _gap_threshold = 2 * _read_len_median  # bp where gap effect begins

            if _gap_threshold > 0:
                st.info(
                    f"**Gap correction**: for paired-end reads of length {int(_read_len_median)} bp, "
                    f"fragments longer than {int(_gap_threshold)} bp have a gap between R1 and R2 where neither "
                    "read covers the position. Longer fragments are therefore underrepresented in the raw count, "
                    "producing a kink at that threshold. **Frequency** mode divides each bin by its capture "
                    "probability \u2014 min(1, 2R\u2009/\u2009L) \u2014 to recover the unbiased fragment size distribution. "
                    "Switch to **Count** to see the raw read counts."
                )

            st.subheader("Insert size distribution")
            _ins_color_options = ["All samples (aggregate)", "Sample"]
            if _alt_reads_has_batch:
                _ins_color_options.append("Batch")
            if _alt_reads_has_label1:
                _ins_color_options.append("Label 1")
            if _alt_reads_has_label2:
                _ins_color_options.append("Label 2")
            if _alt_reads_has_label3:
                _ins_color_options.append("Label 3")
            if _alt_reads_has_timepoint:
                _ins_color_options.append("Timepoint")
            _ins_color_by = st.radio(
                "Color by", _ins_color_options,
                horizontal=True, key="ins_color_by",
            )
            _ins_source, _ins_select_expr, _ins_group_expr, _ins_label_col = _read_series_config(
                _ins_color_by, False, "CASE WHEN ar.is_read1 THEN 'R1' ELSE 'R2' END"
            )

            _ins_df = con.execute(f"""
                SELECT
                    {_ins_select_expr}ar.insert_size,
                    COUNT(*) AS n_reads
                FROM {_ins_source}
                WHERE ar.insert_size BETWEEN {_is_lo} AND {_is_hi}
                GROUP BY {_ins_group_expr}ar.insert_size
                ORDER BY {_ins_group_expr}ar.insert_size
            """).df()

            if not _ins_df.empty:
                _ins_y_mode = st.radio(
                    "Y axis",
                    ["Frequency", "Count"],
                    horizontal=True,
                    key="ins_y_mode",
                )
                _ins_use_corrected = _ins_y_mode == "Frequency"

                # Apply gap correction: weight = min(1, 2R/L); corrected = raw / weight
                if _ins_use_corrected and _read_len_median > 0:
                    _ins_df["n_eff"] = _ins_df["n_reads"] / _ins_df["insert_size"].apply(
                        lambda x: min(1.0, 2.0 * _read_len_median / x) if x > 0 else 1.0
                    )
                else:
                    _ins_df["n_eff"] = _ins_df["n_reads"].astype(float)

                if _ins_use_corrected:
                    if _ins_label_col:
                        _ins_totals = _ins_df.groupby(_ins_label_col)["n_eff"].transform("sum")
                    else:
                        _ins_totals = _ins_df["n_eff"].sum()
                    _ins_df["frequency"] = _ins_df["n_eff"] / _ins_totals

                _ins_y_field = "frequency" if _ins_use_corrected else "n_reads"
                _ins_y_title = "Frequency" if _ins_use_corrected else "Alt-supporting reads"
                _ins_enc = dict(
                    x=alt.X("insert_size:Q", title="Insert size (bp)"),
                    y=alt.Y(f"{_ins_y_field}:Q", title=_ins_y_title,
                            **({"axis": alt.Axis(format=".3f")} if _ins_use_corrected else {})),
                    tooltip=[
                        *([f"{_ins_label_col}:N"] if _ins_label_col else []),
                        alt.Tooltip("insert_size:Q", title="Insert size (bp)"),
                        alt.Tooltip(f"{_ins_y_field}:Q",
                                    title=_ins_y_title,
                                    **({"format": ".4f"} if _ins_use_corrected else {})),
                    ],
                )
                if _ins_label_col:
                    _ins_enc["color"] = alt.Color(f"{_ins_label_col}:N", scale=alt.Scale(scheme="tableau10"))
                    _ins_chart = (
                        alt.Chart(_ins_df)
                        .mark_line(opacity=0.8)
                        .encode(**_ins_enc)
                        .properties(height=300)
                    )
                else:
                    _ins_chart = (
                        alt.Chart(_ins_df)
                        .mark_line(color="#4c78a8")
                        .encode(**_ins_enc)
                        .properties(height=300)
                    )
                st.altair_chart(_ins_chart, width="stretch")
                _ins_caption = (
                    "Insert size distribution of alt-supporting reads. "
                    "A shift toward shorter inserts can indicate adapter contamination or artefacts."
                )
                if _gap_threshold > 0:
                    _ins_caption += (
                        f" A kink near {int(_gap_threshold)} bp (2× read length) is expected: "
                        "fragments longer than this have a coverage gap between R1 and R2, "
                        "so not every fragment is captured at every locus. "
                        "Frequency mode divides each bin by its capture probability "
                        "min(1, 2R/L), recovering the unbiased fragment size distribution."
                    )
                st.caption(_ins_caption)

        # ── Row 3b: Insert size by AF class (germline vs somatic) ─────────────
        if st.session_state.get("_cached_has_insert_size", False):
            st.subheader("Insert size by allele frequency class")
            _af_ins_color_options = ["All samples (aggregate)", "Sample"]
            if _alt_reads_has_batch:
                _af_ins_color_options.append("Batch")
            if _alt_reads_has_label1:
                _af_ins_color_options.append("Label 1")
            if _alt_reads_has_label2:
                _af_ins_color_options.append("Label 2")
            if _alt_reads_has_label3:
                _af_ins_color_options.append("Label 3")
            if _alt_reads_has_timepoint:
                _af_ins_color_options.append("Timepoint")
            _af_ins_ctrl1, _af_ins_ctrl2 = st.columns(2)
            _af_ins_color_by = _af_ins_ctrl1.radio(
                "Color by", _af_ins_color_options,
                horizontal=True, key="af_ins_color_by",
            )
            _af_ins_y_mode = _af_ins_ctrl2.radio(
                "Y axis", ["Frequency", "Count"],
                horizontal=True, key="af_ins_y_mode",
            )
            _af_ins_source, _af_ins_extra_select, _af_ins_extra_group, _af_ins_group_col = _read_series_config(
                _af_ins_color_by, False, "CASE WHEN ar.is_read1 THEN 'R1' ELSE 'R2' END"
            )

            _af_ins_df = con.execute(f"""
                SELECT
                    {_af_ins_extra_select}ar.insert_size,
                    CASE
                        WHEN (_locus.alt_count * 1.0 / _locus.total_depth) > 0.30
                            THEN 'Likely germline (VAF > 30%)'
                        ELSE 'Likely somatic (VAF ≤ 30%)'
                    END AS af_class,
                    COUNT(*) AS n_reads
                FROM {_r_join}
                INNER JOIN (
                    SELECT DISTINCT sample_id, chrom, pos, alt_allele, alt_count, total_depth
                    FROM {table_expr}
                    WHERE {where}
                ) _locus ON  ar.sample_id  = _locus.sample_id
                         AND ar.chrom      = _locus.chrom
                         AND ar.pos        = _locus.pos
                         AND ar.alt_allele = _locus.alt_allele
                WHERE ar.insert_size BETWEEN {_is_lo} AND {_is_hi}
                GROUP BY {_af_ins_extra_group}ar.insert_size, af_class
                ORDER BY {_af_ins_extra_group}ar.insert_size
            """).df()

            if _af_ins_df.empty:
                st.info("No data.")
            else:
                _af_ins_use_corrected = _af_ins_y_mode == "Frequency"

                # When grouping by sample/batch, combine group + af_class into one label
                # so color encodes the group and line style implicitly encodes AF class.
                if _af_ins_group_col:
                    _af_ins_df["series"] = (
                        _af_ins_df[_af_ins_group_col].astype(str)
                        + " — "
                        + _af_ins_df["af_class"]
                    )
                    _af_ins_color_field = "series"
                    _af_ins_color_enc = alt.Color("series:N", title=_af_ins_color_by,
                                                   scale=alt.Scale(scheme="tableau10"))
                else:
                    _af_ins_color_field = "af_class"
                    _af_ins_color_enc = alt.Color("af_class:N", title="AF class",
                                                   scale=alt.Scale(
                                                       domain=["Likely germline (VAF > 30%)", "Likely somatic (VAF ≤ 30%)"],
                                                       range=["#e45756", "#4c78a8"],
                                                   ))

                # Apply gap correction
                if _af_ins_use_corrected and _read_len_median > 0:
                    _af_ins_df["n_eff"] = _af_ins_df["n_reads"] / _af_ins_df["insert_size"].apply(
                        lambda x: min(1.0, 2.0 * _read_len_median / x) if x > 0 else 1.0
                    )
                else:
                    _af_ins_df["n_eff"] = _af_ins_df["n_reads"].astype(float)

                if _af_ins_use_corrected:
                    _norm_key = "series" if _af_ins_group_col else "af_class"
                    _af_totals = _af_ins_df.groupby(_norm_key)["n_eff"].transform("sum")
                    _af_ins_df["frequency"] = _af_ins_df["n_eff"] / _af_totals

                _af_ins_y_title = "Frequency" if _af_ins_use_corrected else "Alt-supporting reads"
                _af_ins_y_field = "frequency" if _af_ins_use_corrected else "n_reads"
                _af_ins_y = (
                    alt.Y(f"{_af_ins_y_field}:Q", title=_af_ins_y_title,
                          axis=alt.Axis(format=".3f"))
                    if _af_ins_use_corrected else
                    alt.Y("n_reads:Q", title="Alt-supporting reads")
                )
                _af_ins_tooltip = [
                    *([f"{_af_ins_group_col}:N"] if _af_ins_group_col else []),
                    alt.Tooltip("insert_size:Q", title="Insert size (bp)"),
                    alt.Tooltip("af_class:N", title="AF class"),
                    alt.Tooltip(_af_ins_y_field + ":Q",
                                title=_af_ins_y_title,
                                **({"format": ".4f"} if _af_ins_use_corrected else {})),
                ]
                _af_ins_chart = (
                    alt.Chart(_af_ins_df)
                    .mark_line(opacity=0.85)
                    .encode(
                        alt.X("insert_size:Q", title="Insert size (bp)"),
                        _af_ins_y,
                        _af_ins_color_enc,
                        tooltip=_af_ins_tooltip,
                    )
                    .properties(height=300)
                )
                st.altair_chart(_af_ins_chart, width="stretch")
                _af_ins_caption = (
                    "Insert size distributions split by allele frequency. "
                    "Each series is normalised independently so lines are directly comparable "
                    "regardless of how many reads fall in each group. "
                    "A shift in one class toward very short or very long inserts suggests artefacts in that group."
                )
                if _gap_threshold > 0:
                    _af_ins_caption += (
                        f" Select 'Frequency' to remove the coverage-gap bias near {int(_gap_threshold)} bp."
                    )
                st.caption(_af_ins_caption)

        # ── Fragment GC distribution ──────────────────────────────────────────
        if _alt_reads_has_frag_gc:
            if "_cached_has_frag_gc" not in st.session_state:
                st.session_state["_cached_has_frag_gc"] = (
                    con.execute("SELECT COUNT(*) FROM alt_reads WHERE frag_gc IS NOT NULL LIMIT 1").fetchone()[0] > 0
                )
            if st.session_state["_cached_has_frag_gc"]:
                st.subheader("Fragment GC distribution")
                _gc_df = con.execute(f"""
                    SELECT
                        ROUND(ar.frag_gc, 2) AS frag_gc_bin,
                        COUNT(*) AS n_reads
                    FROM {_r_join}
                    WHERE ar.frag_gc IS NOT NULL
                    GROUP BY frag_gc_bin
                    ORDER BY frag_gc_bin
                """).df()

                if _gc_df.empty:
                    st.info("No fragment GC data in current selection.")
                else:
                    _gc_chart = (
                        alt.Chart(_gc_df)
                        .mark_bar()
                        .encode(
                            alt.X("frag_gc_bin:Q", title="Fragment GC fraction",
                                  scale=alt.Scale(domain=[0.0, 1.0])),
                            alt.Y("n_reads:Q", title="Alt-supporting reads"),
                            tooltip=[
                                alt.Tooltip("frag_gc_bin:Q", title="GC fraction", format=".2f"),
                                alt.Tooltip("n_reads:Q", title="Reads"),
                            ],
                        )
                        .properties(height=250)
                    )
                    st.altair_chart(_gc_chart, width="stretch")
                    st.caption(
                        "GC fraction of the inferred fragment (reference bases from read start to read start + |TLEN|). "
                        "A bimodal or shifted distribution relative to the target panel GC may indicate GC bias."
                    )

                    # GC by family size
                    if _fs_has_data:
                        _gc_fs_sizes = [2] + list(range(5, 81, 5))
                        _gc_fs_in    = ", ".join(str(s) for s in _gc_fs_sizes)
                        _gc_fs_df = con.execute(f"""
                            SELECT
                                ROUND(ar.frag_gc, 2)  AS frag_gc_bin,
                                ar.family_size,
                                COUNT(*)              AS n_reads
                            FROM {_r_join}
                            WHERE ar.frag_gc IS NOT NULL
                              AND ar.family_size IN ({_gc_fs_in})
                            GROUP BY frag_gc_bin, ar.family_size
                            ORDER BY frag_gc_bin, ar.family_size
                        """).df()

                        if not _gc_fs_df.empty:
                            # Normalise each family-size series independently
                            _gc_fs_df["frac"] = _gc_fs_df.groupby("family_size")["n_reads"].transform(
                                lambda x: x / x.sum()
                            )
                            _gc_fs_df["family_size"] = _gc_fs_df["family_size"].astype(str)

                            _gc_fs_chart = (
                                alt.Chart(_gc_fs_df)
                                .mark_line(opacity=0.85)
                                .encode(
                                    alt.X("frag_gc_bin:Q", title="Fragment GC fraction",
                                          scale=alt.Scale(domain=[0.0, 1.0])),
                                    alt.Y("frac:Q", title="Fraction of reads"),
                                    alt.Color("family_size:O",
                                              title="Family size",
                                              scale=alt.Scale(scheme="viridis"),
                                              sort=[str(s) for s in _gc_fs_sizes]),
                                    tooltip=[
                                        alt.Tooltip("family_size:O", title="Family size"),
                                        alt.Tooltip("frag_gc_bin:Q", title="GC fraction", format=".2f"),
                                        alt.Tooltip("frac:Q", title="Fraction", format=".4f"),
                                    ],
                                )
                                .properties(height=250)
                            )
                            st.altair_chart(_gc_fs_chart, width="stretch")
                            st.caption(
                                "Fragment GC distribution for family sizes 2, 5, 10, 15, … 80, "
                                "each series normalised independently. A consistent GC shape "
                                "across family sizes indicates GC bias is not driven by consensus depth."
                            )

                            # GC bin vs family size (transpose view)
                            _gc_bins = [round(v * 0.1, 1) for v in range(1, 9)]  # 0.1 … 0.8
                            _gc_bin_in = ", ".join(str(b) for b in _gc_bins)
                            _gc_bin_df = con.execute(f"""
                                SELECT
                                    ar.family_size,
                                    ROUND(FLOOR(ar.frag_gc * 10) / 10, 1) AS gc_bin,
                                    COUNT(*) AS n_reads
                                FROM {_r_join}
                                WHERE ar.frag_gc IS NOT NULL
                                  AND ar.family_size IN ({_gc_fs_in})
                                  AND ROUND(FLOOR(ar.frag_gc * 10) / 10, 1) IN ({_gc_bin_in})
                                GROUP BY ar.family_size, gc_bin
                                ORDER BY ar.family_size, gc_bin
                            """).df()

                            if not _gc_bin_df.empty:
                                _gc_bin_df["frac"] = _gc_bin_df.groupby("gc_bin")["n_reads"].transform(
                                    lambda x: x / x.sum()
                                )
                                _gc_bin_df["gc_bin"] = _gc_bin_df["gc_bin"].apply(lambda v: f"{v:.1f}")

                                _gc_bin_chart = (
                                    alt.Chart(_gc_bin_df)
                                    .mark_line(point=True, opacity=0.85)
                                    .encode(
                                        alt.X("family_size:O", title="Family size",
                                              sort=[str(s) for s in _gc_fs_sizes]),
                                        alt.Y("frac:Q", title="Fraction of reads"),
                                        alt.Color("gc_bin:O",
                                                  title="GC bin",
                                                  scale=alt.Scale(scheme="spectral"),
                                                  sort=[f"{b:.1f}" for b in _gc_bins]),
                                        tooltip=[
                                            alt.Tooltip("gc_bin:O",       title="GC bin"),
                                            alt.Tooltip("family_size:O",  title="Family size"),
                                            alt.Tooltip("frac:Q",         title="Fraction", format=".4f"),
                                        ],
                                    )
                                    .properties(height=250)
                                )
                                st.altair_chart(_gc_bin_chart, width="stretch")
                                st.caption(
                                    "Each line is a GC bin (0.1–0.9 in steps of 0.1), normalised "
                                    "independently across family sizes. A line that rises or falls "
                                    "with family size indicates that GC content is associated with "
                                    "consensus depth."
                                )

                            # Mean family size per GC bin
                            _gc_mean_df = con.execute(f"""
                                SELECT
                                    ROUND(FLOOR(ar.frag_gc * 10) / 10, 1) AS gc_bin,
                                    AVG(ar.family_size)                    AS mean_family_size,
                                    COUNT(*)                               AS n_reads
                                FROM {_r_join}
                                WHERE ar.frag_gc IS NOT NULL
                                  AND ar.family_size IS NOT NULL
                                  AND ROUND(FLOOR(ar.frag_gc * 10) / 10, 1) IN ({_gc_bin_in})
                                GROUP BY gc_bin
                                ORDER BY gc_bin
                            """).df()

                            if not _gc_mean_df.empty:
                                _gc_mean_df["gc_bin"] = _gc_mean_df["gc_bin"].apply(lambda v: f"{v:.1f}")
                                _gc_mean_chart = (
                                    alt.Chart(_gc_mean_df)
                                    .mark_bar()
                                    .encode(
                                        alt.X("gc_bin:O", title="GC bin",
                                              sort=[f"{b:.1f}" for b in _gc_bins]),
                                        alt.Y("mean_family_size:Q", title="Mean family size"),
                                        alt.Color("gc_bin:O",
                                                  title="GC bin",
                                                  scale=alt.Scale(scheme="spectral"),
                                                  sort=[f"{b:.1f}" for b in _gc_bins],
                                                  legend=None),
                                        tooltip=[
                                            alt.Tooltip("gc_bin:O",            title="GC bin"),
                                            alt.Tooltip("mean_family_size:Q",  title="Mean family size", format=".2f"),
                                            alt.Tooltip("n_reads:Q",           title="Reads"),
                                        ],
                                    )
                                    .properties(height=250)
                                )
                                st.altair_chart(_gc_mean_chart, width="stretch")
                                st.caption(
                                    "Mean family size for each GC bin. A systematic difference "
                                    "across bins suggests GC-dependent capture efficiency affecting "
                                    "consensus depth."
                                )

        # ── Row 4: Mapping quality distribution ───────────────────────────────
        st.subheader("Mapping quality distribution")
        _mq_df = con.execute(f"""
            SELECT
                map_qual,
                CASE
                    WHEN homopolymer_len >= 5 OR str_len >= 6 THEN 'Repetitive'
                    ELSE 'Non-repetitive'
                END AS locus_type,
                COUNT(*) AS n_reads
            FROM {_r_join}
            INNER JOIN (
                SELECT DISTINCT sample_id, chrom, pos, alt_allele,
                    COALESCE(homopolymer_len, 0) AS homopolymer_len,
                    COALESCE(str_len, 0) AS str_len
                FROM {table_expr}
                WHERE {where}
            ) _locus ON  ar.sample_id  = _locus.sample_id
                     AND ar.chrom      = _locus.chrom
                     AND ar.pos        = _locus.pos
                     AND ar.alt_allele = _locus.alt_allele
            GROUP BY map_qual, locus_type
            ORDER BY map_qual
        """).df()

        if _mq_df.empty:
            st.info("No data.")
        else:
            _mq_chart = (
                alt.Chart(_mq_df)
                .mark_bar(opacity=0.7)
                .encode(
                    alt.X("map_qual:Q", title="Mapping quality (MAPQ)", bin=False),
                    alt.Y("n_reads:Q", title="Alt-supporting reads", stack="zero"),
                    alt.Color("locus_type:N", title="Locus type",
                              scale=alt.Scale(
                                  domain=["Repetitive", "Non-repetitive"],
                                  range=["#e45756", "#4c78a8"],
                              )),
                    alt.Tooltip(["map_qual:Q", "locus_type:N", "n_reads:Q"]),
                )
                .properties(height=300)
            )
            st.altair_chart(_mq_chart, width="stretch")
            st.caption(
                "Stacked by locus type (repetitive = homopolymer ≥ 5 or STR length ≥ 6). "
                "Low MAPQ at repetitive loci indicates multi-mapping artefacts."
            )
