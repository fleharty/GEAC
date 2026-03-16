import streamlit as st
import duckdb
import altair as alt
import pandas as pd

st.set_page_config(page_title="GEAC Explorer", layout="wide")
st.title("GEAC Explorer")
st.markdown(
    "**Genomic Evidence Atlas of Cohorts** — inspect alt base metrics from "
    "per-sample Parquet files or a merged cohort DuckDB database."
)

# ── File input ────────────────────────────────────────────────────────────────
path = st.text_input(
    "Data file path",
    placeholder="/path/to/sample.parquet  or  cohort.duckdb",
)

if not path or not path.strip():
    st.info("Enter a Parquet or DuckDB file path above to begin.")
    st.stop()

path = path.strip()

@st.cache_resource
def open_connection(p: str):
    if p.endswith(".duckdb"):
        return duckdb.connect(p, read_only=True), "alt_bases"
    else:
        con = duckdb.connect()
        return con, f"read_parquet('{p}')"

try:
    con, table_expr = open_connection(path)
except Exception as e:
    st.error(f"Could not open file: {e}")
    st.stop()

# ── Summary stats ─────────────────────────────────────────────────────────────
stats = con.execute(f"""
    SELECT
        COUNT(*)                              AS n_records,
        COUNT(DISTINCT sample_id)             AS n_samples,
        COUNT(DISTINCT chrom)                 AS n_chroms,
        COUNT(DISTINCT chrom || ':' || pos)   AS n_positions,
        SUM(alt_count)                        AS total_alt_reads,
        ROUND(AVG(alt_count * 1.0 / total_depth), 4) AS mean_vaf
    FROM {table_expr}
""").df()

c1, c2, c3, c4, c5, c6 = st.columns(6)
c1.metric("Alt records",      f"{int(stats['n_records'][0]):,}")
c2.metric("Samples",          f"{int(stats['n_samples'][0]):,}")
c3.metric("Chromosomes",      f"{int(stats['n_chroms'][0]):,}")
c4.metric("Positions",        f"{int(stats['n_positions'][0]):,}")
c5.metric("Total alt reads",  f"{int(stats['total_alt_reads'][0]):,}")
c6.metric("Mean VAF",         str(stats["mean_vaf"][0]))

# ── Filters (sidebar) ─────────────────────────────────────────────────────────
st.sidebar.header("Filters")

chroms = con.execute(f"SELECT DISTINCT chrom FROM {table_expr} ORDER BY chrom").df()["chrom"].tolist()
samples = con.execute(f"SELECT DISTINCT sample_id FROM {table_expr} ORDER BY sample_id").df()["sample_id"].tolist()

chrom_sel = st.sidebar.selectbox("Chromosome", ["All"] + chroms)
sample_sel = st.sidebar.multiselect("Samples (blank = all)", samples)
variant_sel = st.sidebar.multiselect(
    "Variant type",
    ["SNV", "insertion", "deletion", "MNV"],
    default=["SNV", "insertion", "deletion", "MNV"],
)
vaf_range = st.sidebar.slider("VAF range", 0.0, 1.0, (0.0, 1.0), step=0.01)
min_alt = st.sidebar.number_input("Min alt count", min_value=1, max_value=10000, value=1, step=1)

# ── Filtered query ────────────────────────────────────────────────────────────
conditions = [
    f"alt_count >= {min_alt}",
    f"alt_count * 1.0 / total_depth BETWEEN {vaf_range[0]} AND {vaf_range[1]}",
]
if chrom_sel != "All":
    conditions.append(f"chrom = '{chrom_sel}'")
if sample_sel:
    s_list = ", ".join(f"'{s}'" for s in sample_sel)
    conditions.append(f"sample_id IN ({s_list})")
if variant_sel:
    t_list = ", ".join(f"'{t}'" for t in variant_sel)
    conditions.append(f"variant_type IN ({t_list})")

where = " AND ".join(conditions)
df = con.execute(f"""
    SELECT *,
           ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
    FROM {table_expr}
    WHERE {where}
    LIMIT 50000
""").df()

if len(df) == 0:
    st.warning("No records match the current filters.")
    st.stop()

cap_msg = " (capped at 50,000 — refine filters to see more)" if len(df) == 50000 else ""
st.info(f"**{len(df):,}** records{cap_msg}")

# ── Data table ────────────────────────────────────────────────────────────────
with st.expander("Data table", expanded=True):
    st.dataframe(
        df[[
            "sample_id", "chrom", "pos", "ref_allele", "alt_allele",
            "variant_type", "vaf", "alt_count", "ref_count", "total_depth",
            "fwd_alt_count", "rev_alt_count", "overlap_alt_agree",
            "overlap_alt_disagree", "variant_called", "variant_filter",
        ]],
        use_container_width=True,
    )

# ── Plots ─────────────────────────────────────────────────────────────────────
tab1, tab2, tab3, tab4 = st.tabs(["VAF distribution", "Error spectrum", "Strand bias", "Overlap agreement"])

with tab1:
    vaf_col1, vaf_col2, vaf_col3 = st.columns(3)
    for col, vtype, color in [
        (vaf_col1, "SNV",       "#4c78a8"),
        (vaf_col2, "insertion", "#f58518"),
        (vaf_col3, "deletion",  "#e45756"),
    ]:
        subset = df[df["variant_type"] == vtype]
        with col:
            if len(subset) == 0:
                st.info(f"No {vtype}s in current selection.")
            else:
                chart = (
                    alt.Chart(subset)
                    .mark_bar(opacity=0.8, color=color)
                    .encode(
                        alt.X("vaf:Q", bin=alt.Bin(maxbins=50), title="VAF",
                              scale=alt.Scale(domain=[0, 1])),
                        alt.Y("count():Q", title="Count"),
                        tooltip=["count():Q"],
                    )
                    .properties(title=f"{vtype} VAF Distribution", height=300)
                )
                st.altair_chart(chart, use_container_width=True)

with tab2:
    snvs = df[df["variant_type"] == "SNV"].copy()
    if len(snvs) > 0:
        spec = snvs.groupby(["ref_allele", "alt_allele"]).size().reset_index(name="count")
        spec["substitution"] = spec["ref_allele"] + ">" + spec["alt_allele"]
        chart = (
            alt.Chart(spec)
            .mark_bar()
            .encode(
                alt.X("substitution:N", sort="-y", title="Substitution"),
                alt.Y("count:Q", title="Count"),
                alt.Color("substitution:N", legend=None),
                tooltip=["substitution:N", "count:Q"],
            )
            .properties(title="SNV Error Spectrum", height=350)
        )
        st.altair_chart(chart, use_container_width=True)
    else:
        st.info("No SNVs in current selection.")

with tab3:
    sample_df = df.sample(min(2000, len(df)))
    chart = (
        alt.Chart(sample_df)
        .mark_point(opacity=0.5, size=30)
        .encode(
            alt.X("fwd_alt_count:Q", title="Forward alt reads"),
            alt.Y("rev_alt_count:Q", title="Reverse alt reads"),
            alt.Color("variant_type:N", title="Variant type"),
            tooltip=["sample_id", "chrom", "pos", "ref_allele", "alt_allele",
                     "fwd_alt_count", "rev_alt_count", "vaf"],
        )
        .properties(title="Strand Bias (up to 2,000 points)", height=350)
    )
    st.altair_chart(chart, use_container_width=True)

with tab4:
    ov = df[df["overlap_depth"] > 0].copy()
    if len(ov) > 0:
        ov["agree_frac"] = (
            ov["overlap_alt_agree"]
            / (ov["overlap_alt_agree"] + ov["overlap_alt_disagree"]).clip(lower=1)
        )
        chart = (
            alt.Chart(ov)
            .mark_bar(opacity=0.8)
            .encode(
                alt.X("agree_frac:Q", bin=alt.Bin(maxbins=20), title="Overlap agreement fraction"),
                alt.Y("count():Q", title="Count"),
                tooltip=["count():Q"],
            )
            .properties(title="Overlap Agreement Fraction", height=350)
        )
        st.altair_chart(chart, use_container_width=True)
    else:
        st.info("No overlapping fragments in current selection.")
