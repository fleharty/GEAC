import io
import zipfile
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
        COUNT(*)                                        AS n_records,
        COUNT(DISTINCT sample_id)                       AS n_samples,
        SUM(alt_count)                                  AS total_alt_bases,
        ROUND(AVG(alt_count * 1.0 / total_depth), 4)   AS mean_vaf,
        ROUND(AVG(total_depth), 1)                      AS mean_depth,
        COUNT(*) FILTER (WHERE variant_called IS NOT NULL) AS n_annotated,
        COUNT(*) FILTER (WHERE variant_called = true)   AS n_called
    FROM {table_expr}
""").df()

n_annotated = int(stats["n_annotated"][0])
n_called    = int(stats["n_called"][0])
pct_called  = f"{100 * n_called / n_annotated:.1f}%" if n_annotated > 0 else "N/A"

# ── Filters (sidebar) ─────────────────────────────────────────────────────────
st.sidebar.header("Filters")

chroms = con.execute(f"SELECT DISTINCT chrom FROM {table_expr} ORDER BY chrom").df()["chrom"].tolist()
samples = con.execute(f"SELECT DISTINCT sample_id FROM {table_expr} ORDER BY sample_id").df()["sample_id"].tolist()

chrom_sel = st.sidebar.selectbox("Chromosome", ["All"] + chroms)
sample_sel = st.sidebar.multiselect("Samples (blank = all)", samples)

_genes_available = con.execute(f"SELECT COUNT(*) FROM {table_expr} WHERE gene IS NOT NULL").fetchone()[0] > 0
genes = con.execute(f"SELECT DISTINCT gene FROM {table_expr} WHERE gene IS NOT NULL ORDER BY gene").df()["gene"].tolist() if _genes_available else []
gene_sel = st.sidebar.multiselect("Gene (blank = all)", genes) if _genes_available else []
if not _genes_available:
    st.sidebar.caption("Gene filter unavailable — run geac collect with --gene-annotations to enable.")
variant_sel = st.sidebar.multiselect(
    "Variant type",
    ["SNV", "insertion", "deletion", "MNV"],
    default=["SNV", "insertion", "deletion", "MNV"],
)
vaf_range = st.sidebar.slider("VAF range", 0.0, 1.0, (0.0, 1.0), step=0.01)
min_alt = st.sidebar.number_input("Min alt count", min_value=1, max_value=10000, value=1, step=1)
variant_called_sel = st.sidebar.selectbox("Variant called", ["All", "Yes", "No", "Unknown (no VCF/TSV)"])
on_target_sel = st.sidebar.selectbox("Target bases", ["All", "On target", "Off target"])
min_depth = st.sidebar.number_input("Min depth (0 = no minimum)", min_value=0, value=0, step=1)
max_depth = st.sidebar.number_input("Max depth (0 = no maximum)", min_value=0, value=0, step=1)

_limit_options = [100, 500, 1000, 5000, 10000, 50000, "All"]
_limit_sel = st.sidebar.selectbox("Display limit (rows)", _limit_options, index=0)
display_limit = None if _limit_sel == "All" else int(_limit_sel)

# ── IGV integration (sidebar) ─────────────────────────────────────────────────
st.sidebar.divider()
st.sidebar.header("IGV Integration")

import os as _os
_default_manifest = _os.path.join(_os.path.dirname(_os.path.abspath(path)), "manifest.tsv")

manifest_path = st.sidebar.text_input(
    "Manifest file (optional)",
    value=_default_manifest,
    help="Tab-separated file with columns: sample_id, bam_path, bai_path",
)
genome = st.sidebar.selectbox("Genome", ["hg19", "hg38", "mm10", "mm39", "other"])
if genome == "other":
    genome = st.sidebar.text_input("Genome ID", value="hg38")

@st.cache_data
def load_manifest(p: str) -> dict:
    mdf = pd.read_csv(p.strip(), sep="\t")
    result = {}
    for row in mdf.itertuples(index=False):
        bai = str(row.bai_path) if hasattr(row, "bai_path") and pd.notna(row.bai_path) else None
        result[str(row.sample_id)] = {"bam": str(row.bam_path), "bai": bai}
    return result

manifest = {}
if manifest_path and manifest_path.strip():
    try:
        manifest = load_manifest(manifest_path.strip())
        st.sidebar.success(f"{len(manifest):,} samples loaded from manifest")
        with st.sidebar.expander("Manifest sample IDs"):
            st.write(sorted(manifest.keys()))
    except Exception as e:
        st.sidebar.error(f"Could not load manifest: {e}")

# ── IGV helper functions ───────────────────────────────────────────────────────
def make_bed(df: pd.DataFrame) -> str:
    # For deletions the alt_allele is e.g. "-ACGT", so the deleted bases span
    # pos+1 .. pos+del_len.  Extend the BED end to cover the full deleted region
    # so IGV highlights the right coordinates.  For all other variant types a
    # single-base interval (pos, pos+1) is correct.
    def _end(row) -> int:
        alt = str(row["alt_allele"])
        if row["variant_type"] == "deletion" and alt.startswith("-"):
            return int(row["pos"]) + len(alt)   # anchor + deleted bases
        return int(row["pos"]) + 1

    tmp = df.copy()
    tmp["bed_end"] = tmp.apply(_end, axis=1)
    # Where multiple records share the same locus, take the largest end coord.
    positions = (
        tmp.groupby(["chrom", "pos"])["bed_end"]
        .max()
        .reset_index()
        .sort_values(["chrom", "pos"])
    )
    lines = [
        f"{row.chrom}\t{int(row.pos)}\t{int(row.bed_end)}"
        for row in positions.itertuples(index=False)
    ]
    return "\n".join(lines) + "\n"


def make_igv_session(df: pd.DataFrame, manifest: dict, genome: str) -> str:
    sample_ids = df["sample_id"].unique().tolist()
    first = df.sort_values(["chrom", "pos"]).iloc[0]
    locus = f"{first['chrom']}:{max(0, int(first['pos']) - 100)}-{int(first['pos']) + 100}"

    resources, tracks = [], []
    for sid in sample_ids:
        entry = manifest.get(str(sid))
        if entry:
            bam, bai = entry["bam"], entry["bai"]
            index_attr = f' index="{bai}"' if bai else ""
            resources.append(f'        <Resource path="{bam}" name="{sid}"{index_attr}/>')
            tracks.append(f'        <Track id="{bam}" name="{sid}"/>')

    resources.append('        <Resource path="positions.bed" name="Selected positions"/>')
    tracks.append('        <Track id="positions.bed" name="Selected positions" color="255,0,0" height="40"/>')

    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        f'<Session genome="{genome}" locus="{locus}" version="8">\n'
        '    <Resources>\n'
        + "\n".join(resources) + "\n"
        '    </Resources>\n'
        '    <Tracks>\n'
        + "\n".join(tracks) + "\n"
        '    </Tracks>\n'
        '</Session>\n'
    )


IGV_CAP = 5

_IGV_CHUNK = 10_000

def igv_buttons(extra_conditions: list[str], display_df: pd.DataFrame, key: str):
    """Render IGV Prepare + Download buttons with chunked progress.

    extra_conditions — additional SQL WHERE fragments (on top of global filters)
    display_df       — already-fetched display DataFrame (used only to enumerate sample_ids)
    key              — unique widget key prefix
    """
    if not manifest:
        st.caption("Add a manifest in the sidebar to enable IGV session download.")
        return

    sample_ids = display_df["sample_id"].unique().tolist()
    n = len(sample_ids)
    cap_samples = sample_ids[:IGV_CAP]

    missing = [sid for sid in sample_ids if str(sid) not in manifest]
    if missing:
        st.warning(
            f"Sample(s) not found in manifest — no BAM track will be added for: "
            f"{', '.join(str(s) for s in missing)}"
        )

    if n > IGV_CAP:
        st.warning(
            f"{n} samples in this selection. IGV session capped at {IGV_CAP}. "
            "Check the box below to override (you're on your own for 50,000 BAMs!)."
        )
        if st.checkbox(f"Load all {n} samples", key=f"{key}_override"):
            cap_samples = sample_ids

    if st.button("Prepare IGV session", key=f"{key}_prepare"):
        w = " AND ".join(conditions + extra_conditions)
        estimated = con.execute(
            f"SELECT COUNT(*) FROM {table_expr} WHERE {w}"
        ).fetchone()[0]

        pbar = st.progress(0, text=f"Querying 0 / ~{estimated:,} records…")

        cursor = con.execute(f"""
            SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
            FROM {table_expr}
            WHERE {w}
        """)
        col_names = [d[0] for d in cursor.description]
        chunks, fetched = [], 0
        while True:
            rows = cursor.fetchmany(_IGV_CHUNK)
            if not rows:
                break
            chunks.append(pd.DataFrame(rows, columns=col_names))
            fetched += len(rows)
            pct = min(fetched / estimated, 1.0) if estimated > 0 else 1.0
            pbar.progress(pct, text=f"Querying {fetched:,} / ~{estimated:,} records…")

        pbar.progress(1.0, text=f"Done — {fetched:,} records fetched.")

        full_df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
        igv_df = full_df[full_df["sample_id"].isin(cap_samples)]
        bed = make_bed(igv_df)
        session = make_igv_session(igv_df, manifest, genome)
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("session.xml", session)
            zf.writestr("positions.bed", bed)
        st.session_state[f"{key}_igv"] = buf.getvalue()

    if f"{key}_igv" in st.session_state:
        st.download_button(
            label="Download IGV session (.zip)",
            data=st.session_state[f"{key}_igv"],
            file_name="igv_session.zip",
            mime="application/zip",
            key=f"{key}_dl",
            help="Extract both files to the same folder, then open session.xml in IGV.",
        )


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
if min_depth > 0:
    conditions.append(f"total_depth >= {min_depth}")
if max_depth > 0:
    conditions.append(f"total_depth <= {max_depth}")
if variant_called_sel == "Yes":
    conditions.append("variant_called = true")
elif variant_called_sel == "No":
    conditions.append("variant_called = false")
elif variant_called_sel == "Unknown (no VCF/TSV)":
    conditions.append("variant_called IS NULL")
if on_target_sel == "On target":
    conditions.append("on_target = true")
elif on_target_sel == "Off target":
    conditions.append("on_target = false")
if gene_sel:
    g_list = ", ".join(f"'{g}'" for g in gene_sel)
    conditions.append(f"gene IN ({g_list})")

where = " AND ".join(conditions)

# ── Summary stats display ──────────────────────────────────────────────────────
fstats = con.execute(f"""
    SELECT
        COUNT(*)                                        AS n_records,
        COUNT(DISTINCT sample_id)                       AS n_samples,
        SUM(alt_count)                                  AS total_alt_bases,
        ROUND(AVG(alt_count * 1.0 / total_depth), 4)   AS mean_vaf,
        ROUND(AVG(total_depth), 1)                      AS mean_depth,
        COUNT(*) FILTER (WHERE variant_called IS NOT NULL) AS n_annotated,
        COUNT(*) FILTER (WHERE variant_called = true)   AS n_called
    FROM {table_expr}
    WHERE {where}
""").df()

fn_annotated = int(fstats["n_annotated"][0])
fn_called    = int(fstats["n_called"][0])
fpct_called  = f"{100 * fn_called / fn_annotated:.1f}%" if fn_annotated > 0 else "N/A"

st.caption("Overall")
c1, c2, c3, c4, c5, c6 = st.columns(6)
c1.metric("Alt records",      f"{int(stats['n_records'][0]):,}")
c2.metric("Samples",          f"{int(stats['n_samples'][0]):,}")
c3.metric("Total alt bases",  f"{int(stats['total_alt_bases'][0]):,}")
c4.metric("Mean VAF",         str(stats["mean_vaf"][0]))
c5.metric("Mean depth",       str(stats["mean_depth"][0]))
c6.metric("% variant called", pct_called)

st.caption("Filtered")
c1, c2, c3, c4, c5, c6 = st.columns(6)
c1.metric("Alt records",      f"{int(fstats['n_records'][0]):,}")
c2.metric("Samples",          f"{int(fstats['n_samples'][0]):,}")
c3.metric("Total alt bases",  f"{int(fstats['total_alt_bases'][0]):,}")
c4.metric("Mean VAF",         str(fstats["mean_vaf"][0]))
c5.metric("Mean depth",       str(fstats["mean_depth"][0]))
c6.metric("% variant called", fpct_called)

def query_records(extra: list[str] = [], limit: int | None = display_limit) -> pd.DataFrame:
    """Query records with current filters plus any extra conditions.
    limit=None fetches all rows (used for IGV); otherwise applies LIMIT."""
    w = " AND ".join(conditions + extra)
    limit_clause = f"LIMIT {limit}" if limit is not None else ""
    return con.execute(f"""
        SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
        FROM {table_expr}
        WHERE {w}
        {limit_clause}
    """).df()

total_count = con.execute(f"SELECT COUNT(*) FROM {table_expr} WHERE {where}").fetchone()[0]

if total_count == 0:
    st.warning("No records match the current filters.")
    st.stop()

df = query_records()

cap_msg = (
    f" (showing {len(df):,} of {total_count:,} — plots and IGV drill-downs use full dataset)"
    if len(df) < total_count else ""
)
st.info(f"**{total_count:,}** records{cap_msg}")

# ── Data table ────────────────────────────────────────────────────────────────
_table_cols = [
    "sample_id", "chrom", "pos", "ref_allele", "alt_allele",
    "variant_type", "vaf", "alt_count", "ref_count", "total_depth",
    "fwd_alt_count", "rev_alt_count", "overlap_alt_agree",
    "overlap_alt_disagree", "variant_called", "variant_filter", "on_target", "gene",
]

with st.expander("Data table", expanded=True):
    st.dataframe(df[_table_cols], use_container_width=True)
    igv_buttons([], df, key="main")

# ── Plots ─────────────────────────────────────────────────────────────────────
tab1, tab2, tab3, tab4 = st.tabs(["VAF distribution", "Error spectrum", "Strand bias", "Overlap agreement"])

with tab1:
    for vtype, color in [
        ("SNV",       "#4c78a8"),
        ("insertion", "#f58518"),
        ("deletion",  "#e45756"),
    ]:
        counts = con.execute(f"""
            SELECT
                FLOOR(ROUND(alt_count * 1.0 / total_depth, 4) * 50) / 50.0 AS vaf_bin,
                FLOOR(ROUND(alt_count * 1.0 / total_depth, 4) * 50) / 50.0 + 0.02 AS vaf_bin_end,
                COUNT(*) AS count
            FROM {table_expr}
            WHERE {where} AND variant_type = '{vtype}'
            GROUP BY vaf_bin, vaf_bin_end
            ORDER BY vaf_bin
        """).df()

        if counts.empty:
            st.info(f"No {vtype}s in current selection.")
        else:
            sel_param = alt.selection_point(
                name="bar_click",
                fields=["vaf_bin", "vaf_bin_end"],
                on="click",
            )
            chart = (
                alt.Chart(counts)
                .mark_bar(color=color)
                .encode(
                    alt.X("vaf_bin:Q",     title="VAF", scale=alt.Scale(domain=[0, 1])),
                    alt.X2("vaf_bin_end:Q"),
                    alt.Y("count:Q",       title="Count"),
                    opacity=alt.condition(sel_param, alt.value(1.0), alt.value(0.4)),
                    tooltip=[
                        alt.Tooltip("vaf_bin:Q",     title="Bin start", format=".3f"),
                        alt.Tooltip("vaf_bin_end:Q", title="Bin end",   format=".3f"),
                        alt.Tooltip("count:Q",       title="Count"),
                    ],
                )
                .add_params(sel_param)
                .properties(title=f"{vtype} VAF Distribution", height=300)
            )
            event = st.altair_chart(chart, use_container_width=True, on_select="rerun")

            pts = (event.selection or {}).get("bar_click", [])
            if pts:
                bin_start = pts[0].get("vaf_bin")
                bin_end   = pts[0].get("vaf_bin_end")
                if bin_start is not None and bin_end is not None:
                    sel = query_records([
                        f"variant_type = '{vtype}'",
                        f"ROUND(alt_count * 1.0 / total_depth, 4) >= {bin_start}",
                        f"ROUND(alt_count * 1.0 / total_depth, 4) < {bin_end}",
                    ])
                    st.caption(
                        f"{len(sel):,} {vtype} records with VAF in "
                        f"[{bin_start:.3f}, {bin_end:.3f})"
                    )
                    st.dataframe(sel[_table_cols], use_container_width=True)
                    igv_buttons([
                        f"variant_type = '{vtype}'",
                        f"ROUND(alt_count * 1.0 / total_depth, 4) >= {bin_start}",
                        f"ROUND(alt_count * 1.0 / total_depth, 4) < {bin_end}",
                    ], sel, key=f"vaf_{vtype}_{bin_start}")

with tab2:
    spec = con.execute(f"""
        SELECT ref_allele || '>' || alt_allele AS substitution, COUNT(*) AS count
        FROM {table_expr}
        WHERE {where} AND variant_type = 'SNV'
        GROUP BY substitution
        ORDER BY count DESC
    """).df()

    if spec.empty:
        st.info("No SNVs in current selection.")
    else:
        sel_param = alt.selection_point(name="bar_click", fields=["substitution"], on="click")
        chart = (
            alt.Chart(spec)
            .mark_bar()
            .encode(
                alt.X("substitution:N", sort="-y", title="Substitution"),
                alt.Y("count:Q", title="Count"),
                alt.Color("substitution:N", legend=None),
                opacity=alt.condition(sel_param, alt.value(1.0), alt.value(0.4)),
                tooltip=["substitution:N", "count:Q"],
            )
            .add_params(sel_param)
            .properties(title="SNV Error Spectrum", height=350)
        )
        event = st.altair_chart(chart, use_container_width=True, on_select="rerun")

        pts = (event.selection or {}).get("bar_click", [])
        if pts:
            sub = pts[0].get("substitution")
            if sub:
                ref, alt_allele = sub.split(">")
                sel = query_records([
                    "variant_type = 'SNV'",
                    f"ref_allele = '{ref}'",
                    f"alt_allele = '{alt_allele}'",
                ])
                st.caption(f"{len(sel):,} records with substitution {sub}")
                st.dataframe(sel[_table_cols], use_container_width=True)
                igv_buttons([
                    "variant_type = 'SNV'",
                    f"ref_allele = '{ref}'",
                    f"alt_allele = '{alt_allele}'",
                ], sel, key=f"spectrum_{sub}")

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
