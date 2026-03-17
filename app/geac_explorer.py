import io
import zipfile
import numpy as np
import streamlit as st
import duckdb
import altair as alt
import pandas as pd
from scipy.optimize import nnls

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
        return con, f"read_parquet('{p}', union_by_name=true)"

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

_schema_cols = set(con.execute(f"DESCRIBE SELECT * FROM {table_expr} LIMIT 0").df()["column_name"].tolist())

def _has_data(col: str) -> bool:
    """True iff col exists in the schema AND has at least one non-null value."""
    if col not in _schema_cols:
        return False
    return con.execute(f"SELECT COUNT(*) FROM {table_expr} WHERE {col} IS NOT NULL").fetchone()[0] > 0

_genes_available = _has_data("gene")
if _genes_available:
    gene_text = st.sidebar.text_input("Gene (partial match, blank = all)", "")
else:
    gene_text = ""
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

_repeat_cols_present = _has_data("homopolymer_len")
if _repeat_cols_present:
    homopolymer_range = st.sidebar.slider("Homopolymer length range", 0, 20, (0, 20), step=1)
    str_len_range     = st.sidebar.slider("STR length range",         0, 50, (0, 50), step=1)
else:
    homopolymer_range = (0, 20)
    str_len_range     = (0, 50)
    st.sidebar.caption("Repeat filters unavailable — run geac collect with a newer build to enable.")
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
if "on_target" in _schema_cols:
    if on_target_sel == "On target":
        conditions.append("on_target = true")
    elif on_target_sel == "Off target":
        conditions.append("on_target = false")
if gene_text.strip() and "gene" in _schema_cols:
    _gene_escaped = gene_text.strip().replace("'", "''")
    conditions.append(f"gene ILIKE '%{_gene_escaped}%'")

if _repeat_cols_present:
    conditions.append(f"homopolymer_len BETWEEN {homopolymer_range[0]} AND {homopolymer_range[1]}")
    conditions.append(f"str_len BETWEEN {str_len_range[0]} AND {str_len_range[1]}")

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
    _tbl_event = st.dataframe(
        df[_table_cols],
        use_container_width=True,
        on_select="rerun",
        selection_mode="single-row",
    )
    igv_buttons([], df, key="main")

# ── Position-level drill-down ──────────────────────────────────────────────────
_selected_rows = (_tbl_event.selection or {}).get("rows", [])
if _selected_rows:
    _row = df.iloc[_selected_rows[0]]
    _chrom, _pos = _row["chrom"], int(_row["pos"])

    # Query ALL samples/alleles at this locus, ignoring current filters.
    _drill_df = con.execute(f"""
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
        WHERE chrom = '{_chrom}' AND pos = {_pos}
        ORDER BY sample_id, alt_allele
    """).df()

    # Extract locus-level annotations from the first row (same for all samples).
    _locus_cols = ["ref_allele"] + [
        c for c in ["gene", "on_target", "homopolymer_len", "str_period", "str_len", "trinuc_context"]
        if c in _schema_cols
    ]
    _locus_row = con.execute(f"""
        SELECT {", ".join(_locus_cols)}
        FROM {table_expr}
        WHERE chrom = '{_chrom}' AND pos = {_pos}
        LIMIT 1
    """).df()

    st.subheader(f"Position drill-down: {_chrom}:{_pos}")

    # Locus-level info as metrics
    _info_cols = st.columns(6)
    _info_cols[0].metric("Ref allele", str(_locus_row["ref_allele"].iloc[0]))
    _info_cols[1].metric("Samples with alt", str(_drill_df["sample_id"].nunique()))
    if "gene" in _locus_row.columns:
        _info_cols[2].metric("Gene", str(_locus_row["gene"].iloc[0] or "intergenic"))
    if "on_target" in _locus_row.columns:
        _info_cols[3].metric("On target", str(_locus_row["on_target"].iloc[0]))
    if "homopolymer_len" in _locus_row.columns:
        _info_cols[4].metric("Homopolymer len", str(_locus_row["homopolymer_len"].iloc[0]))
    if "trinuc_context" in _locus_row.columns:
        _info_cols[5].metric("Trinuc context", str(_locus_row["trinuc_context"].iloc[0] or ""))

    st.dataframe(_drill_df, use_container_width=True, hide_index=True)
    igv_buttons(
        [f"chrom = '{_chrom}'", f"pos = {_pos}"],
        _drill_df,
        key=f"drill_{_chrom}_{_pos}",
    )

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
    _COMP = str.maketrans('ACGT', 'TGCA')
    _SBS_MUT_TYPES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    _SBS_COLORS    = {
        "C>A": "#1BBDEB", "C>G": "#231F20", "C>T": "#E22926",
        "T>A": "#CBCACB", "T>C": "#97D54C", "T>G": "#ECC6C5",
    }
    _SBS_ORDER = [
        f"{b5}[{mt}]{b3}"
        for mt in _SBS_MUT_TYPES
        for b5 in "ACGT"
        for b3 in "ACGT"
    ]

    _trinuc_available = _has_data("trinuc_context")

    if _trinuc_available:
        raw = con.execute(f"""
            SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
            FROM {table_expr}
            WHERE {where} AND variant_type = 'SNV' AND trinuc_context IS NOT NULL
              AND length(trinuc_context) = 3
            GROUP BY trinuc_context, ref_allele, alt_allele
        """).df()

        if raw.empty:
            st.info("No SNVs with trinucleotide context in current selection.")
        else:
            def _sbs(row):
                ctx, r, a = row["trinuc_context"], row["ref_allele"], row["alt_allele"]
                if not all(b in 'ACGT' for b in (ctx + r + a)):
                    return None
                if r in ('A', 'G'):
                    ctx = ctx[::-1].translate(_COMP)
                    r = r.translate(_COMP)
                    a = a.translate(_COMP)
                return f"{ctx[0]}[{r}>{a}]{ctx[2]}"

            raw["sbs_label"] = raw.apply(_sbs, axis=1)
            raw = raw.dropna(subset=["sbs_label"])
            raw["mut_type"] = raw["sbs_label"].str.extract(r'\[([A-Z]>[A-Z])\]')[0]

            agg = raw.groupby(["sbs_label", "mut_type"], as_index=False)["count"].sum()
            # Ensure all 96 contexts are present (fill missing with 0)
            full = pd.DataFrame({
                "sbs_label": _SBS_ORDER,
                "mut_type":  [lbl[2:5] for lbl in _SBS_ORDER],
            })
            spec96 = full.merge(agg, on=["sbs_label", "mut_type"], how="left")
            spec96["count"] = spec96["count"].fillna(0).astype(int)

            sel_param = alt.selection_point(name="bar_click", fields=["sbs_label"], on="click")

            _sub_charts = []
            for _mt in _SBS_MUT_TYPES:
                _sub = spec96[spec96["mut_type"] == _mt].copy()
                _order = [lbl for lbl in _SBS_ORDER if f"[{_mt}]" in lbl]
                _c = (
                    alt.Chart(_sub)
                    .mark_bar(color=_SBS_COLORS[_mt])
                    .encode(
                        alt.X("sbs_label:N", sort=_order, title=None,
                              axis=alt.Axis(labelAngle=-90, labelFontSize=8)),
                        alt.Y("count:Q", title="Count"),
                        opacity=alt.condition(sel_param, alt.value(1.0), alt.value(0.4)),
                        tooltip=["sbs_label:N", "count:Q"],
                    )
                    .add_params(sel_param)
                    .properties(
                        title=alt.TitleParams(_mt, color=_SBS_COLORS[_mt], fontSize=13, fontWeight="bold"),
                        width=150, height=140,
                    )
                )
                _sub_charts.append(_c)

            chart = (
                alt.concat(*_sub_charts, columns=3)
                .resolve_scale(y="shared")
                .properties(title=alt.TitleParams("SNV Trinucleotide Spectrum (SBS96)", fontSize=15))
            )
            event = st.altair_chart(chart, use_container_width=True, on_select="rerun")

            pts = (event.selection or {}).get("bar_click", [])
            if pts:
                clicked = pts[0].get("sbs_label")
                if clicked:
                    matching = raw[raw["sbs_label"] == clicked][["trinuc_context", "ref_allele", "alt_allele"]]
                    if not matching.empty:
                        or_clauses = " OR ".join(
                            f"(trinuc_context = '{r.trinuc_context}' AND ref_allele = '{r.ref_allele}' AND alt_allele = '{r.alt_allele}')"
                            for r in matching.itertuples(index=False)
                        )
                        extra_cond = f"variant_type = 'SNV' AND ({or_clauses})"
                        sel = query_records([extra_cond])
                        st.caption(f"{len(sel):,} records with context {clicked}")
                        st.dataframe(sel[_table_cols], use_container_width=True)
                        igv_buttons([extra_cond], sel, key=f"sbs_{clicked}")

            # ── COSMIC signature known-cause annotations ──────────────────────
            _SBS_ETIOLOGY = {
                # Age / endogenous
                "SBS1":  "Age-related CpG deamination (5-methylcytosine → T)",
                "SBS5":  "Age-related, unknown mechanism",
                # APOBEC
                "SBS2":  "APOBEC (C>T at TC context)",
                "SBS13": "APOBEC (C>G at TC context)",
                # Mismatch repair deficiency / MSI
                "SBS6":  "Mismatch repair deficiency (MMRd/MSI)",
                "SBS14": "MMRd + POLE mutation",
                "SBS15": "Mismatch repair deficiency (MMRd/MSI)",
                "SBS20": "Mismatch repair deficiency (MMRd/MSI)",
                "SBS21": "Mismatch repair deficiency (MMRd/MSI)",
                "SBS26": "Mismatch repair deficiency (MMRd/MSI)",
                "SBS44": "Mismatch repair deficiency (MMRd/MSI)",
                # POLE / proofreading
                "SBS10a": "POLE proofreading exonuclease mutation",
                "SBS10b": "POLE proofreading exonuclease mutation",
                "SBS10c": "POLD1 proofreading mutation",
                "SBS10d": "POLD1 proofreading mutation",
                "SBS28": "POLE mutation (mechanism unclear)",
                # UV light
                "SBS7a": "UV light (C>T at dipyrimidines)",
                "SBS7b": "UV light (C>T at dipyrimidines)",
                "SBS7c": "UV light",
                "SBS7d": "UV light",
                # Tobacco / smoking
                "SBS4":  "Tobacco smoking (BPDE adducts, C>A)",
                "SBS29": "Tobacco chewing",
                # Chemotherapy
                "SBS25": "Chemotherapy (mechanism unclear)",
                "SBS31": "Platinum chemotherapy",
                "SBS35": "Platinum chemotherapy",
                "SBS86": "Chemotherapy (unknown agent)",
                "SBS87": "Thiopurine chemotherapy",
                # Environmental carcinogens
                "SBS22": "Aristolochic acid exposure",
                "SBS24": "Aflatoxin exposure",
                # HR deficiency
                "SBS3":  "Homologous recombination deficiency (BRCA1/2)",
                # Oxidative damage / sequencing artifacts
                "SBS18": "Oxidative damage (8-oxoG, C>A)",
                "SBS36": "Base excision repair deficiency (MUTYH)",
                "SBS58": "Oxidative damage artifact (8-oxoG) — common sequencing artifact",
                # Other / tissue-specific
                "SBS8":  "Unknown — late replication timing",
                "SBS9":  "Polymerase η somatic hypermutation",
                "SBS11": "Temozolomide chemotherapy",
                "SBS16": "Unknown — liver-specific, associated with alcohol",
                "SBS17a": "Unknown — esophageal/gastric enrichment (T>G)",
                "SBS17b": "Unknown — esophageal/gastric enrichment (T>G)",
                "SBS19": "Unknown",
                "SBS23": "Unknown",
                "SBS30": "Base excision repair deficiency (NTHL1)",
                "SBS33": "Unknown",
                "SBS34": "Unknown",
                "SBS37": "Unknown",
                "SBS38": "Indirect UV damage",
                "SBS39": "Unknown",
                "SBS40": "Unknown — age-related",
                "SBS41": "Unknown",
                "SBS84": "AID activity (immune/lymphoid)",
                "SBS85": "AID activity (immune/lymphoid)",
            }

            # ── COSMIC signature decomposition (NNLS) ─────────────────────────
            st.divider()
            st.subheader("COSMIC Signature Decomposition")
            st.caption(
                "Provide a COSMIC SBS matrix file (tab-separated, 96 contexts × N signatures). "
                "Download from the COSMIC website under Mutational Signatures → Downloads "
                "(e.g. COSMIC_v3.4_SBS_GRCh37.txt)."
            )

            cosmic_path = st.text_input(
                "COSMIC SBS matrix path",
                placeholder="/path/to/COSMIC_v3.4_SBS_GRCh37.txt",
                key="cosmic_path",
            )

            @st.cache_data
            def _load_cosmic(p: str) -> pd.DataFrame:
                df_c = pd.read_csv(p, sep="\t", index_col=0)
                return df_c

            if cosmic_path and cosmic_path.strip():
                try:
                    cosmic_df = _load_cosmic(cosmic_path.strip())
                    # Align rows to our canonical SBS96 order
                    cosmic_aligned = cosmic_df.reindex(_SBS_ORDER)
                    missing = cosmic_aligned.isna().any(axis=1).sum()
                    if missing > 0:
                        st.warning(
                            f"{missing} context(s) not found in COSMIC matrix — "
                            "check that the file uses the standard A[C>A]A format."
                        )
                    else:
                        W = cosmic_aligned.values.astype(float)  # 96 × N
                        obs = (
                            spec96.set_index("sbs_label")["count"]
                            .reindex(_SBS_ORDER)
                            .fillna(0)
                            .values.astype(float)
                        )

                        h, _residual = nnls(W, obs)
                        total = h.sum()
                        h_norm = h / total if total > 0 else h

                        sig_df = pd.DataFrame({
                            "signature": cosmic_aligned.columns.tolist(),
                            "exposure":  h_norm,
                        })
                        sig_df["etiology"] = sig_df["signature"].map(
                            lambda s: _SBS_ETIOLOGY.get(s, "")
                        )
                        sig_df = sig_df[sig_df["exposure"] > 0].sort_values(
                            "exposure", ascending=False
                        ).reset_index(drop=True)

                        # ── Goodness-of-fit ───────────────────────────────
                        reconstructed = W @ h
                        cos_sim = (
                            float(np.dot(obs, reconstructed))
                            / (np.linalg.norm(obs) * np.linalg.norm(reconstructed) + 1e-12)
                        )
                        residual_pct = float(np.linalg.norm(obs - reconstructed)) / (float(obs.sum()) + 1e-12) * 100

                        fit_col1, fit_col2 = st.columns(2)
                        fit_col1.metric(
                            "Cosine similarity",
                            f"{cos_sim:.4f}",
                            help="1.0 = perfect reconstruction. Values above 0.95 indicate a good fit.",
                        )
                        fit_col2.metric(
                            "Residual (% of counts)",
                            f"{residual_pct:.1f}%",
                            help="L2 norm of unexplained counts as a percentage of total SNV count. Lower is better.",
                        )

                        top_n_sig = st.slider(
                            "Top signatures to display", 3, min(20, len(sig_df)), 4,
                            key="top_n_sig",
                        )
                        top_df = sig_df.head(top_n_sig)

                        # Fix color domain to full ranked list so colors
                        # don't shift when the slider changes.
                        _all_sigs = list(sig_df["signature"])
                        sig_chart = (
                            alt.Chart(top_df)
                            .mark_bar()
                            .encode(
                                alt.X("signature:N",
                                      sort=list(top_df["signature"]),
                                      title="Signature"),
                                alt.Y("exposure:Q",
                                      title="Exposure (proportion)",
                                      axis=alt.Axis(format=".0%")),
                                alt.Color("signature:N",
                                          scale=alt.Scale(domain=_all_sigs),
                                          legend=None),
                                tooltip=[
                                    "signature:N",
                                    alt.Tooltip("exposure:Q", format=".2%", title="Exposure"),
                                    alt.Tooltip("etiology:N", title="Etiology"),
                                ],
                            )
                            .properties(
                                title=f"Top {top_n_sig} COSMIC SBS Signatures (NNLS fit)",
                                height=300,
                            )
                        )
                        st.altair_chart(sig_chart, use_container_width=True)

                        display = top_df.copy()
                        display["exposure"] = display["exposure"].map("{:.2%}".format)
                        st.dataframe(display, use_container_width=True, hide_index=True)

                except Exception as exc:
                    st.error(f"Failed to load COSMIC matrix: {exc}")

    else:
        # Fallback: simple ref>alt spectrum for older Parquet files
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
    max_val = max(
        int(sample_df["fwd_alt_count"].max()) if len(sample_df) > 0 else 50,
        int(sample_df["rev_alt_count"].max()) if len(sample_df) > 0 else 50,
        1,
    )

    # 95% CI band using normal approximation to Binomial(n, 0.5).
    # For a point (fwd=x, rev=y), total n = x+y. Under H0 (equal strand probability)
    # fwd ~ N(n/2, n/4). The point is in the 95% CI iff lo(n) <= x <= hi(n).
    #
    # For each x, we solve for the rev range [rev_min, rev_max]:
    #   rev_min: x = hi(x+y)  →  s² + z·s − 2x = 0  →  s = (−z + √(z²+8x))/2
    #   rev_max: x = lo(x+y)  →  s² − z·s − 2x = 0  →  s = ( z + √(z²+8x))/2
    # where s = √(x+y) and n = s².
    _x = np.arange(0, max_val + 1, dtype=float)
    _z = 1.96
    _s_lo = (-_z + np.sqrt(_z**2 + 8 * _x)) / 2
    _s_hi = ( _z + np.sqrt(_z**2 + 8 * _x)) / 2
    _ci_band = pd.DataFrame({
        "fwd":     _x,
        "rev_min": np.maximum(_s_lo**2 - _x, 0.0),
        "rev_max": _s_hi**2 - _x,
    })

    _diag = pd.DataFrame({"fwd": [0.0, float(max_val)], "rev": [0.0, float(max_val)]})

    ci_area = (
        alt.Chart(_ci_band)
        .mark_area(opacity=0.25, color="steelblue")
        .encode(
            alt.X("fwd:Q"),
            alt.Y("rev_min:Q"),
            alt.Y2("rev_max:Q"),
        )
    )
    diag_line = (
        alt.Chart(_diag)
        .mark_line(strokeDash=[6, 4], color="gray", opacity=0.7)
        .encode(
            alt.X("fwd:Q"),
            alt.Y("rev:Q"),
        )
    )
    scatter = (
        alt.Chart(sample_df)
        .mark_point(opacity=0.5, size=30)
        .encode(
            alt.X("fwd_alt_count:Q", title="Forward alt reads"),
            alt.Y("rev_alt_count:Q", title="Reverse alt reads"),
            alt.Color("variant_type:N", title="Variant type"),
            tooltip=(
                ["sample_id", "chrom", "pos", "ref_allele", "alt_allele",
                 "fwd_alt_count", "rev_alt_count", "vaf"]
                + (["gene"] if _genes_available else [])
            ),
        )
        .properties(title="Strand Bias — dashed line: perfect balance; shaded band: 95% CI under Binomial(n, 0.5)", height=350)
    )
    st.altair_chart((ci_area + diag_line + scatter).resolve_scale(color="independent"), use_container_width=True)

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
