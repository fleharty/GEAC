import streamlit as st
import duckdb
import altair as alt
import pandas as pd
import geac_config

st.set_page_config(page_title="GEAC Coverage Explorer", layout="wide")
st.title("GEAC Coverage Explorer")
st.markdown(
    "Per-position depth quality metrics from `geac coverage`. "
    "Accepts a merged cohort DuckDB (with a `coverage` table) or a single `.coverage.parquet` file."
)

# ── Project config ─────────────────────────────────────────────────────────────
_cfg = geac_config.load()

# ── File input ─────────────────────────────────────────────────────────────────
path = st.text_input(
    "Data file path",
    value=_cfg.get("data", ""),
    placeholder="/path/to/cohort.duckdb  or  sample.coverage.parquet",
)

if not path or not path.strip():
    st.info("Enter a coverage DuckDB or Parquet file path above to begin.")
    st.stop()

path = path.strip()

try:
    if path.endswith(".duckdb"):
        con = duckdb.connect(path, read_only=True)
        table_expr = "coverage"
        try:
            con.execute("SELECT 1 FROM coverage LIMIT 1")
        except Exception:
            st.error(
                "No `coverage` table found in this DuckDB. "
                "Run `geac merge` with `.coverage.parquet` files to create it."
            )
            st.stop()
    else:
        con = duckdb.connect()
        table_expr = f"read_parquet('{path}', union_by_name=true)"
except Exception as e:
    st.error(f"Could not open file: {e}")
    st.stop()

# ── Detect optional columns ────────────────────────────────────────────────────
_cols = {
    row[0]
    for row in con.execute(f"DESCRIBE SELECT * FROM {table_expr} LIMIT 0").fetchall()
}
_has_gene      = "gene" in _cols
_has_on_target = "on_target" in _cols

# ── Sidebar filters ────────────────────────────────────────────────────────────
_hdr_col, _btn_col = st.sidebar.columns([2, 1])
_hdr_col.header("Filters")
_FILTER_KEYS = ["sample_sel", "chrom_sel", "gene_text", "on_target_sel"]
if _btn_col.button("Clear", help="Reset all filters"):
    for k in _FILTER_KEYS:
        if k in st.session_state:
            del st.session_state[k]
    st.rerun()

all_samples = (
    con.execute(f"SELECT DISTINCT sample_id FROM {table_expr} ORDER BY sample_id")
    .df()["sample_id"]
    .tolist()
)
sample_sel = st.sidebar.multiselect(
    "Samples", all_samples, default=all_samples, key="sample_sel"
)

all_chroms = (
    con.execute(f"SELECT DISTINCT chrom FROM {table_expr} ORDER BY chrom")
    .df()["chrom"]
    .tolist()
)
chrom_sel = st.sidebar.selectbox("Chromosome", ["All"] + all_chroms, key="chrom_sel")

if _has_gene:
    gene_text = st.sidebar.text_input("Gene (partial match)", "", key="gene_text")
else:
    gene_text = ""

if _has_on_target:
    on_target_sel = st.sidebar.selectbox(
        "Target bases", ["All", "On target", "Off target"], key="on_target_sel"
    )
else:
    on_target_sel = "All"


# ── Filter helpers ─────────────────────────────────────────────────────────────
def _filter_clauses(extra: list[str] | None = None) -> list[str]:
    clauses: list[str] = []
    if sample_sel:
        ids = ", ".join(f"'{s.replace(chr(39), chr(39)*2)}'" for s in sample_sel)
        clauses.append(f"sample_id IN ({ids})")
    if chrom_sel != "All":
        c = chrom_sel.replace("'", "''")
        clauses.append(f"chrom = '{c}'")
    if gene_text.strip():
        g = gene_text.strip().replace("'", "''")
        clauses.append(f"LOWER(gene) LIKE '%{g.lower()}%'")
    if on_target_sel == "On target":
        clauses.append("on_target = true")
    elif on_target_sel == "Off target":
        clauses.append("on_target = false")
    if extra:
        clauses.extend(extra)
    return clauses


def where(extra: list[str] | None = None) -> str:
    clauses = _filter_clauses(extra)
    return ("WHERE " + " AND ".join(clauses)) if clauses else ""


# ── Tabs ───────────────────────────────────────────────────────────────────────
tab_summary, tab_depth, tab_gc, tab_low = st.tabs(
    ["Summary", "Depth Distribution", "GC Bias", "Low Coverage"]
)

# ──────────────────────────────────────────────────────────────────────────────
# Tab 1 — Per-sample summary
# ──────────────────────────────────────────────────────────────────────────────
with tab_summary:
    st.subheader("Per-sample coverage summary")

    insert_col = (
        "ROUND(AVG(mean_insert_size) FILTER (WHERE n_insert_size_obs > 0), 1)"
        if "n_insert_size_obs" in _cols
        else "NULL"
    )

    summary_df = con.execute(f"""
        SELECT
            sample_id,
            COUNT(*)                                AS n_positions,
            ROUND(AVG(total_depth),   1)            AS mean_depth,
            ROUND(MEDIAN(total_depth), 1)           AS median_depth,
            ROUND(MIN(total_depth * 1.0),  0)       AS min_depth,
            ROUND(MAX(total_depth * 1.0),  0)       AS max_depth,
            ROUND(AVG(frac_dup),      3)            AS mean_frac_dup,
            ROUND(AVG(frac_overlap),  3)            AS mean_frac_overlap,
            ROUND(AVG(frac_mapq0),    3)            AS mean_frac_mapq0,
            {insert_col}                            AS mean_insert_size
        FROM {table_expr}
        {where()}
        GROUP BY sample_id
        ORDER BY sample_id
    """).df()

    st.dataframe(summary_df, use_container_width=True)

    if len(summary_df) > 0:
        st.markdown("**Mean depth per sample**")
        bar = (
            alt.Chart(summary_df)
            .mark_bar()
            .encode(
                x=alt.X("mean_depth:Q", title="Mean depth"),
                y=alt.Y("sample_id:N", sort="-x", title="Sample"),
                tooltip=["sample_id", "mean_depth", "median_depth", "n_positions"],
            )
            .properties(height=max(120, 28 * len(summary_df)))
        )
        st.altair_chart(bar, use_container_width=True)

        qc_metrics = summary_df.melt(
            id_vars="sample_id",
            value_vars=["mean_frac_dup", "mean_frac_overlap", "mean_frac_mapq0"],
            var_name="metric",
            value_name="value",
        )
        _metric_labels = {
            "mean_frac_dup":     "Frac dup",
            "mean_frac_overlap": "Frac overlap",
            "mean_frac_mapq0":   "Frac MAPQ=0",
        }
        qc_metrics["metric"] = qc_metrics["metric"].map(_metric_labels)

        st.markdown("**QC fractions per sample** (lower is better for dup and MAPQ=0)")
        qc_chart = (
            alt.Chart(qc_metrics)
            .mark_bar()
            .encode(
                x=alt.X("value:Q", title="Fraction", scale=alt.Scale(domain=[0, 1])),
                y=alt.Y("sample_id:N", title="Sample"),
                color=alt.Color("metric:N", title="Metric"),
                row=alt.Row("metric:N", title=""),
                tooltip=["sample_id", "metric", alt.Tooltip("value:Q", format=".3f")],
            )
            .properties(height=max(60, 22 * len(summary_df)))
            .resolve_scale(x="independent")
        )
        st.altair_chart(qc_chart, use_container_width=True)

# ──────────────────────────────────────────────────────────────────────────────
# Tab 2 — Depth distribution
# ──────────────────────────────────────────────────────────────────────────────
with tab_depth:
    st.subheader("Depth distribution")
    depth_cap = st.slider("Cap depth display at", 10, 1000, 200, step=10)

    depth_df = con.execute(f"""
        SELECT
            sample_id,
            LEAST(total_depth, {depth_cap}) AS depth_bin,
            COUNT(*) AS n_positions
        FROM {table_expr}
        {where()}
        GROUP BY sample_id, depth_bin
        ORDER BY sample_id, depth_bin
    """).df()

    if depth_df.empty:
        st.info("No data for current filters.")
    else:
        chart = (
            alt.Chart(depth_df)
            .mark_line(opacity=0.8)
            .encode(
                x=alt.X("depth_bin:Q", title=f"Total depth (capped at {depth_cap})"),
                y=alt.Y("n_positions:Q", title="Number of positions"),
                color=alt.Color("sample_id:N", title="Sample"),
                tooltip=["sample_id", "depth_bin", "n_positions"],
            )
            .properties(width=700, height=350)
            .interactive()
        )
        st.altair_chart(chart, use_container_width=True)

        # Fraction at depth thresholds
        st.markdown("**Fraction of positions at or above depth threshold**")
        thresholds = [1, 10, 20, 30, 50, 100]
        rows = []
        for sid, grp in depth_df.groupby("sample_id"):
            total = grp["n_positions"].sum()
            row = {"sample_id": sid}
            cumsum = grp.sort_values("depth_bin")["n_positions"].cumsum()
            for t in thresholds:
                below = grp.loc[grp["depth_bin"] < t, "n_positions"].sum()
                row[f">={t}x"] = round((total - below) / total, 3) if total > 0 else 0.0
            rows.append(row)
        st.dataframe(pd.DataFrame(rows), use_container_width=True)

# ──────────────────────────────────────────────────────────────────────────────
# Tab 3 — GC bias
# ──────────────────────────────────────────────────────────────────────────────
with tab_gc:
    st.subheader("GC bias")
    st.caption(
        "Mean depth per GC content bin (5% bins). "
        "A U-shaped or monotone curve indicates GC-dependent depth variation."
    )

    gc_df = con.execute(f"""
        SELECT
            sample_id,
            ROUND(gc_content * 20) / 20     AS gc_bin,
            ROUND(AVG(total_depth), 2)       AS mean_depth,
            COUNT(*)                         AS n_positions
        FROM {table_expr}
        {where()}
        GROUP BY sample_id, gc_bin
        HAVING COUNT(*) >= 5
        ORDER BY sample_id, gc_bin
    """).df()

    if gc_df.empty:
        st.info("No data for current filters.")
    else:
        gc_chart = (
            alt.Chart(gc_df)
            .mark_line(point=True, opacity=0.8)
            .encode(
                x=alt.X(
                    "gc_bin:Q",
                    title="GC content (5% bins)",
                    scale=alt.Scale(domain=[0, 1]),
                    axis=alt.Axis(format=".0%"),
                ),
                y=alt.Y("mean_depth:Q", title="Mean depth"),
                color=alt.Color("sample_id:N", title="Sample"),
                tooltip=[
                    "sample_id",
                    alt.Tooltip("gc_bin:Q", title="GC bin", format=".0%"),
                    alt.Tooltip("mean_depth:Q", title="Mean depth", format=".1f"),
                    alt.Tooltip("n_positions:Q", title="Positions"),
                ],
            )
            .properties(width=700, height=350)
            .interactive()
        )
        st.altair_chart(gc_chart, use_container_width=True)

        # Optionally overlay frac_mapq0 on a second y-axis (shown as separate chart)
        mapq_df = con.execute(f"""
            SELECT
                sample_id,
                ROUND(gc_content * 20) / 20     AS gc_bin,
                ROUND(AVG(frac_mapq0), 4)        AS mean_frac_mapq0,
                COUNT(*)                         AS n_positions
            FROM {table_expr}
            {where()}
            GROUP BY sample_id, gc_bin
            HAVING COUNT(*) >= 5
            ORDER BY sample_id, gc_bin
        """).df()

        with st.expander("Frac MAPQ=0 by GC bin (mappability signal)"):
            mapq_chart = (
                alt.Chart(mapq_df)
                .mark_line(point=True, opacity=0.8)
                .encode(
                    x=alt.X(
                        "gc_bin:Q",
                        title="GC content (5% bins)",
                        scale=alt.Scale(domain=[0, 1]),
                        axis=alt.Axis(format=".0%"),
                    ),
                    y=alt.Y(
                        "mean_frac_mapq0:Q",
                        title="Mean frac MAPQ=0",
                        scale=alt.Scale(domain=[0, 1]),
                    ),
                    color=alt.Color("sample_id:N", title="Sample"),
                    tooltip=[
                        "sample_id",
                        alt.Tooltip("gc_bin:Q", title="GC bin", format=".0%"),
                        alt.Tooltip("mean_frac_mapq0:Q", title="Frac MAPQ=0", format=".3f"),
                    ],
                )
                .properties(width=700, height=250)
                .interactive()
            )
            st.altair_chart(mapq_chart, use_container_width=True)

# ──────────────────────────────────────────────────────────────────────────────
# Tab 4 — Low coverage
# ──────────────────────────────────────────────────────────────────────────────
with tab_low:
    st.subheader("Undercovered positions")

    col1, col2 = st.columns(2)
    depth_threshold = col1.slider("Depth threshold", 1, 200, 20)
    min_samples_frac = col2.slider(
        "Min fraction of samples below threshold", 0.0, 1.0, 0.5, step=0.05,
        help="Only show positions below threshold in this fraction of selected samples.",
    )

    n_selected = len(sample_sel) if sample_sel else 1
    min_samples_count = max(1, int(min_samples_frac * n_selected))

    gene_col  = "gene,"       if _has_gene else ""
    gene_sel  = "gene,"       if _has_gene else ""

    low_df = con.execute(f"""
        SELECT
            chrom,
            pos,
            {gene_col}
            COUNT(DISTINCT sample_id)          AS n_samples_below,
            ROUND(AVG(total_depth), 1)         AS mean_depth,
            ROUND(MIN(total_depth * 1.0), 0)   AS min_depth,
            ROUND(AVG(frac_mapq0), 3)          AS mean_frac_mapq0,
            ROUND(AVG(gc_content), 3)          AS mean_gc_content
        FROM {table_expr}
        {where([f"total_depth < {depth_threshold}"])}
        GROUP BY chrom, pos {(',' + gene_sel.rstrip(',')) if _has_gene else ''}
        HAVING COUNT(DISTINCT sample_id) >= {min_samples_count}
        ORDER BY n_samples_below DESC, mean_depth
        LIMIT 1000
    """).df()

    n_low = len(low_df)
    st.caption(
        f"Positions with `total_depth < {depth_threshold}` in ≥ {min_samples_frac:.0%} of samples "
        f"({min_samples_count}/{n_selected}). Showing up to 1000 rows."
    )

    if low_df.empty:
        st.success("No positions below threshold for the current filters.")
    else:
        st.dataframe(low_df, use_container_width=True)

        if _has_gene and "gene" in low_df.columns:
            st.markdown("**Affected positions per gene**")
            by_gene = (
                low_df.groupby("gene", dropna=False)
                .agg(n_positions=("pos", "count"), mean_depth=("mean_depth", "mean"))
                .reset_index()
                .sort_values("n_positions", ascending=False)
                .head(30)
            )
            gene_bar = (
                alt.Chart(by_gene)
                .mark_bar()
                .encode(
                    x=alt.X("n_positions:Q", title="Positions below threshold"),
                    y=alt.Y("gene:N", sort="-x", title="Gene"),
                    color=alt.Color(
                        "mean_depth:Q",
                        scale=alt.Scale(scheme="reds", reverse=True),
                        title="Mean depth",
                    ),
                    tooltip=["gene", "n_positions", alt.Tooltip("mean_depth:Q", format=".1f")],
                )
                .properties(height=max(120, 22 * len(by_gene)))
            )
            st.altair_chart(gene_bar, use_container_width=True)
