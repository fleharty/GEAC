"""MNV candidates tab — cross-locus fragment_id join to find co-occurring substitutions."""
from __future__ import annotations

import altair as alt
import streamlit as st

from explorer.tab_context import TabContext

LABEL = "🔗 MNV candidates"


def render(ctx: TabContext) -> None:
    if not ctx.has_alt_reads:
        st.info(
            "No `alt_reads` table found. Re-run `geac collect --reads-output` to enable MNV detection.",
            icon="ℹ️",
        )
        return

    if not ctx.alt_reads_has("fragment_id"):
        st.info(
            "`fragment_id` column is absent from `alt_reads`. Re-run `geac collect --reads-output` "
            "with a build that includes fragment ID hashing (v0.4.25+) to enable MNV detection.",
            icon="ℹ️",
        )
        return

    # ── Controls ─────────────────────────────────────────────────────────────
    _c1, _c2, _c3 = st.columns(3)
    _min_dist = _c1.slider(
        "Min distance between substitutions (bp)",
        min_value=1, max_value=20, value=1, step=1,
        key="mnv_min_dist",
        help="Only pairs of alt loci separated by at least this many base pairs are considered.",
    )
    _max_dist = _c2.slider(
        "Max distance between substitutions (bp)",
        min_value=1, max_value=20, value=3, step=1,
        key="mnv_max_dist",
        help="Only pairs of alt loci within this many base pairs are considered MNV candidates.",
    )
    _min_co = _c3.slider(
        "Min co-occurring reads",
        min_value=1, max_value=20, value=2, step=1,
        key="mnv_min_co",
        help="Only MNV pairs seen together on at least this many fragments are shown.",
    )

    # ── Query ─────────────────────────────────────────────────────────────────
    # filtered_reads: alt reads restricted to loci passing the main sidebar filters,
    # with the active per-read filter applied (via ctx.r_join).
    _sql = f"""
        WITH filtered_reads AS (
            SELECT ar.fragment_id, ar.sample_id, ar.chrom, ar.pos, ar.alt_allele
            FROM {ctx.r_join}
        ),
        co_occurring AS (
            SELECT
                a.sample_id,
                a.chrom,
                a.pos        AS pos1,
                a.alt_allele AS alt_allele1,
                b.pos        AS pos2,
                b.alt_allele AS alt_allele2,
                COUNT(*)     AS co_count
            FROM filtered_reads a
            JOIN filtered_reads b
                ON  a.fragment_id = b.fragment_id
                AND a.sample_id   = b.sample_id
                AND a.chrom       = b.chrom
                AND a.pos         < b.pos
                AND (b.pos - a.pos) BETWEEN {_min_dist} AND {_max_dist}
            GROUP BY a.sample_id, a.chrom, a.pos, a.alt_allele, b.pos, b.alt_allele
            HAVING COUNT(*) >= {_min_co}
        ),
        locus_info AS (
            SELECT DISTINCT sample_id, chrom, pos, ref_allele, alt_allele, alt_count
            FROM {ctx.table_expr}
            WHERE {ctx.where}
        )
        SELECT
            c.sample_id,
            c.chrom,
            c.pos1,
            c.pos1 + 1                                              AS pos1_display,
            l1.ref_allele                                           AS ref_allele1,
            c.alt_allele1,
            c.pos2,
            c.pos2 + 1                                              AS pos2_display,
            l2.ref_allele                                           AS ref_allele2,
            c.alt_allele2,
            c.co_count,
            ROUND(c.co_count * 1.0 / NULLIF(l1.alt_count, 0), 4)  AS frac_cooccurring,
            l1.ref_allele || l2.ref_allele                          AS dinuc_ref,
            c.alt_allele1 || c.alt_allele2                         AS dinuc_alt,
            (c.pos2 - c.pos1)                                       AS distance_bp
        FROM co_occurring c
        JOIN locus_info l1
            ON  c.sample_id  = l1.sample_id
            AND c.chrom      = l1.chrom
            AND c.pos1       = l1.pos
            AND c.alt_allele1 = l1.alt_allele
        JOIN locus_info l2
            ON  c.sample_id  = l2.sample_id
            AND c.chrom      = l2.chrom
            AND c.pos2       = l2.pos
            AND c.alt_allele2 = l2.alt_allele
        ORDER BY c.co_count DESC, c.sample_id, c.chrom, c.pos1
    """

    with ctx.timed("mnv_candidates"):
        df = ctx.con.execute(_sql).df()

    if df.empty:
        st.info(
            f"No MNV candidates found with ≥{_min_co} co-occurring reads "
            f"separated by {_min_dist}–{_max_dist} bp under the current filters.",
            icon="🔎",
        )
        return

    st.caption(
        f"**{len(df):,} MNV candidate pair(s)** — substitutions sharing ≥{_min_co} fragment(s) "
        f"separated by {_min_dist}–{_max_dist} bp. Rows are sorted by co-occurrence count."
    )

    # ── Display table ─────────────────────────────────────────────────────────
    _display_cols = [
        "sample_id", "chrom", "pos1_display", "ref_allele1", "alt_allele1",
        "pos2_display", "ref_allele2", "alt_allele2",
        "co_count", "frac_cooccurring", "dinuc_ref", "dinuc_alt", "distance_bp",
    ]
    _col_cfg = {
        "pos1_display":       st.column_config.NumberColumn("pos1 (1-based)", format="%d"),
        "pos2_display":       st.column_config.NumberColumn("pos2 (1-based)", format="%d"),
        "frac_cooccurring":   st.column_config.NumberColumn("frac cooccurring", format="%.4f"),
        "co_count":           st.column_config.NumberColumn("co-occurring reads"),
        "distance_bp":        st.column_config.NumberColumn("distance (bp)"),
        "dinuc_ref":          st.column_config.TextColumn("dinuc ref"),
        "dinuc_alt":          st.column_config.TextColumn("dinuc alt"),
    }

    _sel = st.dataframe(
        df[_display_cols],
        column_config=_col_cfg,
        hide_index=True,
        use_container_width=True,
        on_select="rerun",
        selection_mode="multi-row",
        key="mnv_table",
    )

    # ── IGV drill-down ────────────────────────────────────────────────────────
    _sel_rows = (_sel.selection or {}).get("rows", [])
    if _sel_rows:
        _sel_df = df.iloc[_sel_rows]

        # Build OR clauses covering both positions of each selected pair so IGV
        # can show the reads at either locus.
        _or_clauses = []
        for _, row in _sel_df.iterrows():
            sid = ctx.sql_str(str(row["sample_id"]))
            chrom = ctx.sql_str(str(row["chrom"]))
            _or_clauses.append(
                f"(sample_id = {sid} AND chrom = {chrom} AND pos = {int(row['pos1'])})"
            )
            _or_clauses.append(
                f"(sample_id = {sid} AND chrom = {chrom} AND pos = {int(row['pos2'])})"
            )

        _locus_df = ctx.con.execute(f"""
            SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                   pos + 1 AS pos_display
            FROM {ctx.table_expr}
            WHERE ({" OR ".join(_or_clauses)})
            ORDER BY sample_id, chrom, pos
        """).df()

        n_pairs = len(_sel_rows)
        st.caption(f"{n_pairs} pair(s) selected — showing both loci for IGV navigation")
        st.dataframe(_locus_df[ctx.table_cols], use_container_width=True, hide_index=True)
        ctx.igv_buttons(
            [f"({' OR '.join(_or_clauses)})"],
            _locus_df,
            key=f"mnv_igv_{'_'.join(str(df.iloc[i]['pos1']) for i in _sel_rows[:3])}",
        )
    else:
        st.caption("Select rows to navigate to those loci in IGV.")

    # ── Dinucleotide context bar chart ────────────────────────────────────────
    if len(df) >= 2:
        st.divider()
        _render_dinuc_chart(df)


def _render_dinuc_chart(df) -> None:
    st.subheader("Dinucleotide context")
    st.caption(
        "Substitution pairs grouped by their reference dinucleotide context "
        "(ref1+ref2 → alt1+alt2). APOBEC signature: CC→TT or similar clustered transitions."
    )

    _dinuc = (
        df.groupby(["dinuc_ref", "dinuc_alt"])["co_count"]
        .sum()
        .reset_index()
        .rename(columns={"co_count": "total_co_count"})
        .sort_values("total_co_count", ascending=False)
    )
    _dinuc["context"] = _dinuc["dinuc_ref"] + "→" + _dinuc["dinuc_alt"]

    _chart = (
        alt.Chart(_dinuc)
        .mark_bar()
        .encode(
            alt.X("context:N", sort="-y", title="Dinucleotide context (ref→alt)"),
            alt.Y("total_co_count:Q", title="Total co-occurring reads"),
            tooltip=["context:N", "dinuc_ref:N", "dinuc_alt:N",
                     alt.Tooltip("total_co_count:Q", title="Total co-occurring reads")],
        )
        .properties(height=220)
    )
    st.altair_chart(_chart, use_container_width=True)
