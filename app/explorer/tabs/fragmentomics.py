"""Fragmentomics tab — 4-mer end-motif and GC-content analysis of cfDNA fragments.

Three views: (1) end-motif frequency vs insert size with optional smoothing
and sample/batch grouping; (2) FFT power spectrum to surface the ~10.5 bp
nucleosome-positioning periodicity; (3) end-motif frequency by GC content
plus a fragment-level GC distribution. Sample, motif end (5′/3′), insert-size
range, and grouping dimension are picked in-tab; chrom/region/batch/label/
timepoint filters come from the sidebar via `ctx.fragments_where_parts`.
"""
from __future__ import annotations

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from explorer.tab_context import TabContext

LABEL = "🧬 Fragmentomics"


def render(ctx: TabContext) -> None:
    if not ctx.has_fragments:
        st.info(
            "No `fragments` view found in this database. "
            "Run `geac fragments` on each sample and include the "
            "`.fragments.parquet` files when running `geac merge`, or use "
            "`run_fragments=true` in `geac_cohort.wdl`."
        )
        return

    st.subheader("End-motif frequency by insert size")
    st.caption(
        "4-mer end-motif frequency at each insert size. "
        "Motifs are centered on the fragment cut site [cut−2, cut+2) using the reference sequence. "
        "Color encodes the fraction of fragments at each insert size carrying that motif."
    )

    cols = st.columns([1, 1, 1, 1, 1, 1, 1])
    with cols[0]:
        end = st.selectbox("End", ["5′", "3′"], key="fm_end")
    with cols[1]:
        top_n = st.slider("Top N motifs", 5, 64, 16, key="fm_top_n")
    with cols[2]:
        is_min = st.number_input("Min insert size", value=50, step=10, key="fm_is_min")
    with cols[3]:
        is_max = st.number_input("Max insert size", value=400, step=10, key="fm_is_max")
    with cols[4]:
        smooth = st.slider("Smoothing window (bp)", 1, 51, 11, step=2, key="fm_smooth")
    with cols[5]:
        n_columns = st.slider("Columns", 1, 8, 4, key="fm_columns")
    with cols[6]:
        group_by = st.radio("Group by", ["Totals", "Sample", "Batch"], key="fm_by_sample")
        group_dim, group_label = _group_dim_for(group_by)

    motif_col = "end_motif_5p" if end == "5′" else "end_motif_3p"

    with ctx.timed("fm_sample_list"):
        samples = ctx.con.execute(
            "SELECT DISTINCT sample_id FROM fragments ORDER BY sample_id"
        ).df()["sample_id"].tolist()

    if not st.session_state.get("fm_samples"):
        st.session_state["fm_samples"] = samples[:min(len(samples), 6)]
    sel_samples = st.multiselect("Samples", samples, key="fm_samples")
    if not sel_samples:
        st.info("Select at least one sample.")
        return

    sample_list = ", ".join(f"'{s}'" for s in sel_samples)
    where_parts = [
        f"{motif_col} IS NOT NULL",
        f"insert_size BETWEEN {int(is_min)} AND {int(is_max)}",
        f"sample_id IN ({sample_list})",
    ] + list(ctx.fragments_where_parts)
    where_sql = " AND ".join(where_parts)

    with ctx.timed("fm_top_motifs"):
        top_df = ctx.con.execute(f"""
            SELECT {motif_col} AS motif, COUNT(*) AS n
            FROM fragments
            WHERE {where_sql}
            GROUP BY motif
            ORDER BY n DESC
            LIMIT {top_n}
        """).df()

    if top_df.empty:
        st.info("No fragments match the current filters.")
        return

    top_motifs = top_df["motif"].tolist()
    motif_list = ", ".join(f"'{m}'" for m in top_motifs)

    df_raw, df_smoothed = _query_motif_by_insert_size(
        ctx, motif_col, where_sql, motif_list, group_dim, smooth
    )
    _render_motif_facets(df_smoothed, top_motifs, group_dim, group_label, end, n_columns)

    _render_power_spectrum(df_raw, top_motifs, group_dim, group_label, end, n_columns)

    _render_motif_by_gc(
        ctx, motif_col, where_sql, motif_list, top_motifs,
        group_dim, group_label, end,
    )

    _render_fragment_gc_distribution(ctx, where_sql, group_dim, group_label)


def _group_dim_for(group_by: str) -> tuple[str | None, str]:
    if group_by == "Totals":
        return None, group_by
    if group_by == "Sample":
        return "sample_id", group_by
    if group_by == "Batch":
        return "batch", group_by
    raise ValueError(f"Unknown grouping mode: {group_by}")


def _query_motif_by_insert_size(
    ctx: TabContext, motif_col: str, where_sql: str, motif_list: str,
    group_dim, smooth: int,
):
    with ctx.timed("fm_query"):
        if group_dim:
            df = ctx.con.execute(f"""
                WITH counts AS (
                    SELECT
                        {group_dim},
                        {motif_col}      AS motif,
                        insert_size,
                        COUNT(*)         AS n
                    FROM fragments
                    WHERE {where_sql}
                      AND {motif_col} IN ({motif_list})
                    GROUP BY {group_dim}, motif, insert_size
                )
                SELECT
                    {group_dim},
                    motif,
                    insert_size,
                    n,
                    n * 1.0 / SUM(n) OVER (PARTITION BY {group_dim}, insert_size) AS freq
                FROM counts
                ORDER BY {group_dim}, insert_size, motif
            """).df()
        else:
            df = ctx.con.execute(f"""
                WITH counts AS (
                    SELECT
                        {motif_col}      AS motif,
                        insert_size,
                        COUNT(*)         AS n
                    FROM fragments
                    WHERE {where_sql}
                      AND {motif_col} IN ({motif_list})
                    GROUP BY motif, insert_size
                )
                SELECT
                    motif,
                    insert_size,
                    n,
                    n * 1.0 / SUM(n) OVER (PARTITION BY insert_size) AS freq
                FROM counts
                ORDER BY insert_size, motif
            """).df()

    sort_cols  = [group_dim, "motif", "insert_size"] if group_dim else ["motif", "insert_size"]
    group_cols = [group_dim, "motif"] if group_dim else ["motif"]
    raw = df.copy()
    if smooth > 1:
        df = (
            df.sort_values(sort_cols)
            .assign(freq=lambda d: d.groupby(group_cols)["freq"]
                    .transform(lambda s: s.rolling(smooth, center=True, min_periods=1).mean()))
        )
    return raw, df


def _render_motif_facets(
    df: pd.DataFrame, top_motifs, group_dim, group_label, end: str, n_columns: int,
) -> None:
    encode = dict(
        x=alt.X("insert_size:Q", title="Insert size (bp)"),
        y=alt.Y("freq:Q", title="Fraction"),
        tooltip=[
            alt.Tooltip("motif:N",       title="Motif"),
            alt.Tooltip("insert_size:Q", title="Insert size"),
            alt.Tooltip("freq:Q",        title="Fraction", format=".4f"),
            alt.Tooltip("n:Q",           title="Count"),
        ],
    )
    if group_dim:
        encode["color"] = alt.Color(f"{group_dim}:N", title=group_label)
        encode["tooltip"].insert(0, alt.Tooltip(f"{group_dim}:N", title=group_label))

    facet_width = max(80, (700 // n_columns) - 20)
    chart = (
        alt.Chart(df)
        .mark_line()
        .encode(**encode)
        .properties(width=facet_width, height=120)
        .facet(
            facet=alt.Facet("motif:N", title=f"{end} end motif", sort=top_motifs),
            columns=n_columns,
        )
        .resolve_scale(y="independent")
    )
    st.altair_chart(chart, use_container_width=True)


def _render_power_spectrum(
    df_raw: pd.DataFrame, top_motifs, group_dim, group_label, end: str, n_columns: int,
) -> None:
    st.subheader("Power spectrum of end-motif frequency (insert size axis)")
    st.caption(
        "FFT power spectrum of end-motif frequency along the insert size axis. "
        "A peak near **10.5 bp** (red line) indicates rotational nucleosome "
        "positioning signal. Computed on unsmoothed data; the smoothing window "
        "above does not affect this plot. Restrict the insert size range to the "
        "mononucleosome window (100–300 bp) for best signal-to-noise."
    )

    fc1, fc2 = st.columns(2)
    with fc1:
        fft_is_min = st.number_input(
            "Min insert size for FFT", value=100, step=10, key="fm_fft_is_min"
        )
    with fc2:
        fft_is_max = st.number_input(
            "Max insert size for FFT", value=300, step=10, key="fm_fft_is_max"
        )

    records = []
    groups = (
        [((gval, m), grp)
         for (gval, m), grp in df_raw.groupby([group_dim, "motif"])]
        if group_dim
        else [((None, m), grp) for m, grp in df_raw.groupby("motif")]
    )
    for (gval, motif), grp in groups:
        sub = (
            grp[(grp["insert_size"] >= fft_is_min) & (grp["insert_size"] <= fft_is_max)]
            .sort_values("insert_size")
        )
        if len(sub) < 20:
            continue
        is_vals   = sub["insert_size"].values
        freq_vals = sub["freq"].values
        grid   = np.arange(is_vals[0], is_vals[-1] + 1)
        signal = np.interp(grid, is_vals, freq_vals)
        signal -= signal.mean()
        signal *= np.hanning(len(signal))
        power  = np.abs(np.fft.rfft(signal)) ** 2
        freqs  = np.fft.rfftfreq(len(grid), d=1.0)
        for f, p in zip(freqs[1:], power[1:]):
            period = 1.0 / f
            if 5.0 <= period <= 50.0:
                rec = {"motif": motif, "period_bp": round(period, 3), "power": p}
                if gval is not None:
                    rec[group_dim] = gval
                records.append(rec)

    if not records:
        st.info("Not enough data in the selected insert size range for FFT.")
        return

    fft_df = pd.DataFrame(records)
    encode = dict(
        x=alt.X("period_bp:Q", title="Period (bp)"),
        y=alt.Y("power:Q",     title="Power"),
        tooltip=[
            alt.Tooltip("motif:N",     title="Motif"),
            alt.Tooltip("period_bp:Q", title="Period (bp)", format=".2f"),
            alt.Tooltip("power:Q",     title="Power",       format=".3e"),
        ],
    )
    if group_dim:
        encode["color"] = alt.Color(f"{group_dim}:N", title=group_label)
        encode["tooltip"].insert(0, alt.Tooltip(f"{group_dim}:N", title=group_label))

    facet_width = max(80, (700 // n_columns) - 20)
    line = alt.Chart(fft_df).mark_line().encode(**encode)
    rule = (
        alt.Chart(fft_df)
        .mark_rule(color="red", strokeDash=[4, 2])
        .encode(x=alt.datum(10.5))
    )
    chart = (
        alt.layer(line, rule)
        .properties(width=facet_width, height=120)
        .facet(
            facet=alt.Facet("motif:N", title=f"{end} end motif", sort=top_motifs),
            columns=n_columns,
        )
        .resolve_scale(y="independent")
    )
    st.altair_chart(chart, use_container_width=True)


def _render_motif_by_gc(
    ctx: TabContext, motif_col: str, where_sql: str, motif_list: str, top_motifs,
    group_dim, group_label, end: str,
) -> None:
    st.subheader("End-motif frequency by GC content")
    st.caption(
        "4-mer end-motif frequency across fragment GC content bins. "
        "Fraction is normalized within each GC bin."
    )

    gc1, gc2, gc3 = st.columns([1, 1, 1])
    with gc1:
        gc_bin = st.select_slider(
            "GC bin size", options=[0.01, 0.02, 0.05, 0.10], value=0.05,
            key="fm_gc_bin",
        )
    with gc2:
        gc_smooth = st.slider(
            "Smoothing window (bins)", 1, 11, 3, step=2, key="fm_gc_smooth",
        )
    with gc3:
        gc_columns = st.slider("Columns", 1, 8, 4, key="fm_gc_columns")

    gc_where = where_sql + " AND gc_content IS NOT NULL"

    with ctx.timed("fm_gc_query"):
        if group_dim:
            df = ctx.con.execute(f"""
                WITH counts AS (
                    SELECT
                        {group_dim},
                        {motif_col}                                  AS motif,
                        ROUND(gc_content / {gc_bin}) * {gc_bin}     AS gc_bin,
                        COUNT(*)                                     AS n
                    FROM fragments
                    WHERE {gc_where}
                      AND {motif_col} IN ({motif_list})
                    GROUP BY {group_dim}, motif, gc_bin
                )
                SELECT
                    {group_dim}, motif, gc_bin, n,
                    n * 1.0 / SUM(n) OVER (PARTITION BY {group_dim}, gc_bin) AS freq
                FROM counts
                ORDER BY {group_dim}, gc_bin, motif
            """).df()
        else:
            df = ctx.con.execute(f"""
                WITH counts AS (
                    SELECT
                        {motif_col}                                  AS motif,
                        ROUND(gc_content / {gc_bin}) * {gc_bin}     AS gc_bin,
                        COUNT(*)                                     AS n
                    FROM fragments
                    WHERE {gc_where}
                      AND {motif_col} IN ({motif_list})
                    GROUP BY motif, gc_bin
                )
                SELECT
                    motif, gc_bin, n,
                    n * 1.0 / SUM(n) OVER (PARTITION BY gc_bin) AS freq
                FROM counts
                ORDER BY gc_bin, motif
            """).df()

    if df.empty:
        st.info("No fragments with GC content data match the current filters.")
        return

    sort_cols  = [group_dim, "motif", "gc_bin"] if group_dim else ["motif", "gc_bin"]
    group_cols = [group_dim, "motif"] if group_dim else ["motif"]
    if gc_smooth > 1:
        df = (
            df.sort_values(sort_cols)
            .assign(freq=lambda d: d.groupby(group_cols)["freq"]
                    .transform(lambda s: s.rolling(gc_smooth, center=True, min_periods=1).mean()))
        )

    encode = dict(
        x=alt.X("gc_bin:Q", title="GC content",
                scale=alt.Scale(domain=[0, 1]),
                axis=alt.Axis(format=".0%")),
        y=alt.Y("freq:Q", title="Fraction"),
        tooltip=[
            alt.Tooltip("motif:N",  title="Motif"),
            alt.Tooltip("gc_bin:Q", title="GC bin", format=".2f"),
            alt.Tooltip("freq:Q",   title="Fraction", format=".4f"),
            alt.Tooltip("n:Q",      title="Count"),
        ],
    )
    if group_dim:
        encode["color"] = alt.Color(f"{group_dim}:N", title=group_label)
        encode["tooltip"].insert(0, alt.Tooltip(f"{group_dim}:N", title=group_label))

    facet_width = max(80, (700 // gc_columns) - 20)
    chart = (
        alt.Chart(df)
        .mark_line()
        .encode(**encode)
        .properties(width=facet_width, height=120)
        .facet(
            facet=alt.Facet("motif:N", title=f"{end} end motif", sort=top_motifs),
            columns=gc_columns,
        )
        .resolve_scale(y="independent")
    )
    st.altair_chart(chart, use_container_width=True)


def _render_fragment_gc_distribution(
    ctx: TabContext, where_sql: str, group_dim, group_label,
) -> None:
    st.subheader("Fragment GC content distribution")
    st.caption(
        "Distribution of fragment-level GC content (computed from the reference "
        "over the full fragment span). Useful for assessing library prep GC bias "
        "and cfDNA GC enrichment. Note: the end-motif by GC plot above is "
        "reference-confounded; this plot shows real sequencing data."
    )

    fgc_bin = st.select_slider(
        "Bin size", options=[0.01, 0.02, 0.05], value=0.02,
        key="fm_frag_gc_bin",
    )

    with ctx.timed("fm_frag_gc_query"):
        if group_dim:
            df = ctx.con.execute(f"""
                SELECT
                    {group_dim},
                    ROUND(gc_content / {fgc_bin}) * {fgc_bin} AS gc_bin,
                    COUNT(*) AS n
                FROM fragments
                WHERE {where_sql} AND gc_content IS NOT NULL
                GROUP BY {group_dim}, gc_bin
                ORDER BY {group_dim}, gc_bin
            """).df()
        else:
            df = ctx.con.execute(f"""
                SELECT
                    ROUND(gc_content / {fgc_bin}) * {fgc_bin} AS gc_bin,
                    COUNT(*) AS n
                FROM fragments
                WHERE {where_sql} AND gc_content IS NOT NULL
                GROUP BY gc_bin
                ORDER BY gc_bin
            """).df()

    if df.empty:
        return

    encode = dict(
        x=alt.X("gc_bin:Q", title="GC content",
                scale=alt.Scale(domain=[0, 1]),
                axis=alt.Axis(format=".0%")),
        y=alt.Y("n:Q", title="Fragment count"),
        tooltip=[
            alt.Tooltip("gc_bin:Q", title="GC bin",   format=".2f"),
            alt.Tooltip("n:Q",      title="Fragments", format=","),
        ],
    )
    if group_dim:
        encode["color"]   = alt.Color(f"{group_dim}:N", title=group_label)
        encode["opacity"] = alt.value(0.5)
        encode["tooltip"].insert(0, alt.Tooltip(f"{group_dim}:N", title=group_label))

    chart = (
        alt.Chart(df)
        .mark_area(opacity=0.5 if not group_dim else 0.35, line=True)
        .encode(**encode)
        .properties(height=250)
    )
    st.altair_chart(chart, use_container_width=True)
