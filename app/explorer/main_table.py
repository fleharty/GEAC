from __future__ import annotations

import duckdb
import pandas as pd
import streamlit as st


def render_summary_metrics(
    stats: pd.DataFrame,
    filtered_stats: pd.DataFrame,
    pct_called: str,
    filtered_pct_called: str,
) -> None:
    st.caption("Overall")
    c1, c2, c3, c4, c5, c6 = st.columns(6)
    c1.metric("Alt records", f"{int(stats['n_records'][0]):,}")
    c2.metric("Samples", f"{int(stats['n_samples'][0]):,}")
    c3.metric("Total alt bases", f"{int(stats['total_alt_bases'][0]):,}")
    c4.metric("Mean VAF", str(stats["mean_vaf"][0]))
    c5.metric("Mean depth", str(stats["mean_depth"][0]))
    c6.metric("% variant called", pct_called)

    st.caption("Filtered")
    c1, c2, c3, c4, c5, c6 = st.columns(6)
    c1.metric("Alt records", f"{int(filtered_stats['n_records'][0]):,}")
    c2.metric("Samples", f"{int(filtered_stats['n_samples'][0]):,}")
    c3.metric("Total alt bases", f"{int(filtered_stats['total_alt_bases'][0]):,}")
    c4.metric("Mean VAF", str(filtered_stats["mean_vaf"][0]))
    c5.metric("Mean depth", str(filtered_stats["mean_depth"][0]))
    c6.metric("% variant called", filtered_pct_called)


def render_records_table(
    df: pd.DataFrame,
    total_count: int,
    table_cols: list[str],
    *,
    igv_buttons,
    key: str,
):
    with st.expander("Data table", expanded=True):
        caption = (
            f"Showing {len(df):,} of {total_count:,} rows. Increase the table row limit above to see more."
            if len(df) < total_count
            else f"Showing all {total_count:,} rows."
        )
        st.caption(caption)
        event = st.dataframe(
            df[table_cols],
            width="stretch",
            on_select="rerun",
            selection_mode="single-row",
            key=f"{key}_data_table",
        )
        igv_buttons([], df, key=key)
        return event


def render_position_drilldown(
    *,
    con: duckdb.DuckDBPyConnection,
    table_expr: str,
    schema_cols: set[str],
    chrom: str,
    pos: int,
    selected_alt: str,
    has_alt_reads: bool,
    sql_str,
    igv_buttons,
) -> None:
    drill_df = con.execute(
        f"""
        SELECT
            sample_id,
            ref_allele,
            alt_allele,
            variant_type,
            alt_count,
            ref_count,
            total_depth,
            ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
            fwd_alt_count,
            rev_alt_count,
            overlap_alt_agree,
            overlap_alt_disagree,
            variant_called,
            variant_filter
        FROM {table_expr}
        WHERE chrom = '{chrom}' AND pos = {pos}
        ORDER BY sample_id, alt_allele
        """
    ).df()

    locus_cols = ["ref_allele"] + [
        c
        for c in [
            "gene",
            "on_target",
            "homopolymer_len",
            "str_period",
            "str_len",
            "trinuc_context",
        ]
        if c in schema_cols
    ]
    locus_row = con.execute(
        f"""
        SELECT {", ".join(locus_cols)}
        FROM {table_expr}
        WHERE chrom = '{chrom}' AND pos = {pos}
        LIMIT 1
        """
    ).df()

    st.subheader(f"🔍 Position drill-down: {chrom}:{pos}")

    match_alt = st.checkbox(
        f"Same alt allele only ({selected_alt})",
        value=False,
        key=f"drill_match_alt_{chrom}_{pos}",
    )
    if match_alt:
        drill_df = drill_df[drill_df["alt_allele"] == selected_alt].reset_index(drop=True)

    info_cols = st.columns(6)
    info_cols[0].metric("Ref allele", str(locus_row["ref_allele"].iloc[0]))
    info_cols[1].metric("Samples with alt", str(drill_df["sample_id"].nunique()))
    if "gene" in locus_row.columns:
        info_cols[2].metric("Gene", str(locus_row["gene"].iloc[0] or "intergenic"))
    if "on_target" in locus_row.columns:
        info_cols[3].metric("On target", str(locus_row["on_target"].iloc[0]))
    if "homopolymer_len" in locus_row.columns:
        info_cols[4].metric("Homopolymer len", str(locus_row["homopolymer_len"].iloc[0]))
    if "trinuc_context" in locus_row.columns:
        info_cols[5].metric("Trinuc context", str(locus_row["trinuc_context"].iloc[0] or ""))

    st.dataframe(drill_df, width="stretch", hide_index=True)

    drill_igv_conditions = [f"chrom = '{chrom}'", f"pos = {pos}"]
    if match_alt:
        drill_igv_conditions.append(f"alt_allele = '{sql_str(selected_alt)}'")
    igv_buttons(
        drill_igv_conditions,
        drill_df,
        key=f"drill_{chrom}_{pos}",
        use_global_filters=False,
    )

    if not has_alt_reads:
        return

    reads_alt_clause = f" AND alt_allele = '{sql_str(selected_alt)}'" if match_alt else ""
    reads_df = con.execute(
        f"""
        SELECT
            sample_id,
            alt_allele,
            cycle,
            read_length,
            is_read1,
            ab_count,
            ba_count,
            family_size,
            base_qual,
            map_qual
        FROM alt_reads
        WHERE chrom = '{chrom}' AND pos = {pos}{reads_alt_clause}
        ORDER BY sample_id, alt_allele, family_size DESC NULLS LAST
        """
    ).df()

    if reads_df.empty:
        st.caption("No per-read detail available for this locus.")
        return

    st.markdown(
        f"**Per-read detail** — {len(reads_df):,} alt-supporting reads across {reads_df['sample_id'].nunique()} sample(s)"
    )

    reads_summary = (
        reads_df.groupby(["sample_id", "alt_allele"])
        .agg(
            n_reads=("cycle", "count"),
            median_cycle=("cycle", "median"),
            median_family_size=("family_size", "median"),
            min_family_size=("family_size", "min"),
            max_family_size=("family_size", "max"),
            mean_base_qual=("base_qual", "mean"),
        )
        .reset_index()
        .round(1)
    )
    st.caption("Summary by sample / allele")
    st.dataframe(reads_summary, width="stretch", hide_index=True)

    with st.expander("📖 Individual reads"):
        st.dataframe(reads_df, width="stretch", hide_index=True)
