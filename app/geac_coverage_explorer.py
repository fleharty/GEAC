import streamlit as st
import duckdb
import altair as alt
import pandas as pd
import geac_config
import json
import streamlit.components.v1 as components
from igv_helpers import load_manifest, resolve_index_uri
from explorer import COVERAGE_FILTER_STATE, GEAC_VERSION, DataSource

st.set_page_config(page_title="GEAC Coverage Explorer", layout="wide")
st.title("GEAC Coverage Explorer")
st.markdown(
    "Per-position depth quality metrics from `geac coverage`. "
    "Accepts a merged cohort DuckDB (with a `coverage` table) or a single `.coverage.parquet` file."
)

st.sidebar.caption(f"geac v{GEAC_VERSION}")

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
    data_source = DataSource.open_coverage(path)
except Exception as e:
    st.error(f"Could not open file: {e}")
    st.stop()

con = data_source.con
table_expr = data_source.table_expr

# ── Version metadata (DuckDB only) ────────────────────────────────────────────
if data_source.is_duckdb:
    _db_version = data_source.db_version
    _db_created = data_source.db_created
    if _db_version is not None and _db_version != GEAC_VERSION:
        st.warning(
            f"Version mismatch: database was created with geac v{_db_version}, "
            f"but this Explorer expects v{GEAC_VERSION}. "
            "Results may be incomplete or columns may differ.",
            icon="⚠️",
        )
    if _db_created is not None:
        st.sidebar.caption(str(_db_created)[:10])

_missing_required_cols = data_source.required_columns_missing()
if _missing_required_cols:
    st.warning(
        "This dataset is missing required `coverage` columns expected by the current Explorer: "
        + ", ".join(sorted(_missing_required_cols)),
        icon="⚠️",
    )

# ── Detect optional columns ────────────────────────────────────────────────────
_cols = set(data_source.schema_cols)
_has_gene      = "gene" in _cols
_has_on_target = "on_target" in _cols
_has_bin_n     = "bin_n" in _cols
_has_intervals = data_source.is_duckdb and "coverage_intervals" in data_source.available_tables

# ── Sidebar filters ────────────────────────────────────────────────────────────
_hdr_col, _btn_col = st.sidebar.columns([2, 1])
_hdr_col.header("Filters")
if _btn_col.button("Clear", help="Reset all filters"):
    COVERAGE_FILTER_STATE.clear(st.session_state)
    st.rerun()

all_samples = data_source.distinct_values("sample_id")
sample_sel = st.sidebar.multiselect(
    "Samples", all_samples, default=all_samples, key="sample_sel"
)

all_chroms = data_source.distinct_values("chrom")
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

# ── Manifest (IGV) ─────────────────────────────────────────────────────────────
st.sidebar.divider()
st.sidebar.header("IGV Integration")
_default_manifest = _cfg.get("manifest", "")
manifest_path = st.sidebar.text_input(
    "Manifest file",
    value=_default_manifest,
    placeholder="/path/to/manifest.tsv",
    help="TSV with columns: collaborator_sample_id, duplex_output_bam, duplex_output_bam_index",
)
_manifest: dict = {}
if manifest_path and manifest_path.strip():
    try:
        _manifest = load_manifest(manifest_path.strip())
        st.sidebar.success(f"{len(_manifest):,} samples loaded")
    except Exception as _e:
        st.sidebar.error(f"Could not load manifest: {_e}")

if data_source.is_duckdb:
    st.sidebar.divider()
    with st.sidebar.expander("Advanced", expanded=False):
        st.caption("Database metadata")
        _meta_header = data_source.metadata_header()
        if _meta_header.empty:
            st.caption("No `geac_metadata` table found.")
        else:
            _meta_row = _meta_header.iloc[0]
            for _col in _meta_header.columns:
                _val = _meta_row[_col]
                if pd.isna(_val):
                    _display = "NULL"
                elif isinstance(_val, float):
                    _display = f"{_val:g}"
                else:
                    _display = str(_val)
                st.caption(f"{_col}: {_display}")

        _meta_inputs = data_source.metadata_inputs()
        if not _meta_inputs.empty:
            st.caption("Merged inputs")
            st.dataframe(
                _meta_inputs,
                width="stretch",
                hide_index=True,
                key="coverage_advanced_metadata_inputs",
            )

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


# ── IGV.js component ───────────────────────────────────────────────────────────
@st.cache_data(ttl=2700, show_spinner=False)  # refresh well before the 60-min token expiry
def _gcs_access_token() -> tuple[str | None, str | None]:
    """Return (token, error_message) from ADC.

    error_message is None on success.  On failure it is a human-readable string
    suitable for display in the UI.  Two distinct failure modes are reported:
    - 'not_installed': google-auth package is absent
    - other: credentials exist but refresh failed
    """
    try:
        import google.auth
        import google.auth.transport.requests
    except ModuleNotFoundError:
        return None, "not_installed"
    try:
        credentials, project = google.auth.default(
            scopes=["https://www.googleapis.com/auth/devstorage.read_only"]
        )
        credentials.refresh(google.auth.transport.requests.Request())
        token = credentials.token
        print(f"[IGV] credential type : {type(credentials).__name__}", flush=True)
        print(f"[IGV] project         : {project}", flush=True)
        print(f"[IGV] token obtained  : {bool(token)}", flush=True)
        print(f"[IGV] token prefix    : {token[:20] if token else 'None'}...", flush=True)
        print(f"[IGV] token expiry    : {getattr(credentials, 'expiry', 'unknown')}", flush=True)
        return token, None
    except Exception as e:
        print(f"[IGV] token error: {e}", flush=True)
        return None, str(e)
def _igv_tracks(
    sample_ids: tuple[str, ...],
    manifest_json: str,
    oauth_token: str | None = None,
    track_height: int = 250,
) -> list[dict]:
    """Build IGV track dicts for *sample_ids*, passing URIs (including gs://) directly."""
    manifest = json.loads(manifest_json)
    tracks = []
    for sid in sample_ids:
        entry = manifest.get(str(sid))
        if not entry:
            continue
        bam_uri = entry["bam"]
        index_uri = resolve_index_uri(bam_uri, entry.get("bai"))
        track: dict = {
            "name": sid,
            "url": bam_uri,
            "height": track_height,
            "colorBy": "strand",
        }
        if index_uri:
            track["indexURL"] = index_uri
        if oauth_token:
            track["oauthToken"] = oauth_token
        tracks.append(track)
    return tracks


def _igv_html(locus: str, tracks: list[dict], genome: str = "hg38", oauth_token: str | None = None) -> str:
    tracks_json = json.dumps(tracks, indent=2)
    token_js = f'igv.setOauthToken("{oauth_token}", "storage.googleapis.com");' if oauth_token else ""
    oauth_option = f'"oauthToken": "{oauth_token}",' if oauth_token else ""
    # Explicit hg38 reference ensures CRAM decoding works regardless of what
    # reference path is embedded in the CRAM header.
    reference = {
        "id": "hg38",
        "fastaURL": "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa",
        "indexURL": "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa.fai",
        "cytobandURL": "https://s3.amazonaws.com/igv.org.genomes/hg38/cytoBandIdeo.txt.gz",
    }
    reference_json = json.dumps(reference)
    return f"""<!DOCTYPE html>
<html>
<head>
<style>
  body {{ margin: 0; padding: 4px; font-family: sans-serif; }}
  #igv-div {{ width: 100%; }}
</style>
<script src="https://cdn.jsdelivr.net/npm/igv@3.0.2/dist/igv.min.js"></script>
</head>
<body>
<div id="igv-div"></div>
<script>
  {token_js}
  igv.createBrowser(document.getElementById("igv-div"), {{
    {oauth_option}
    reference: {reference_json},
    locus: "{locus}",
    tracks: {tracks_json}
  }});
</script>
</body>
</html>"""


# ── Tabs ───────────────────────────────────────────────────────────────────────
_tab_names = ["Summary", "Depth Distribution", "GC Bias", "Low Coverage", "Depth Profile", "IGV"]
if _has_intervals:
    _tab_names.append("Intervals")
_tabs = st.tabs(_tab_names)
tab_summary, tab_depth, tab_gc, tab_low, tab_profile, tab_igv = _tabs[:6]
tab_intervals = _tabs[6] if _has_intervals else None

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

    st.dataframe(summary_df, width="stretch")

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
        st.altair_chart(bar, width="stretch")

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
        st.altair_chart(qc_chart, width="stretch")

# ──────────────────────────────────────────────────────────────────────────────
# Tab 2 — Depth distribution
# ──────────────────────────────────────────────────────────────────────────────
with tab_depth:
    st.subheader("Depth distribution")
    depth_cap = st.slider("Cap depth display at", 10, 1000, 200, step=10)

    _loci_expr = "SUM(bin_n)" if _has_bin_n else "COUNT(*)"

    depth_df = con.execute(f"""
        SELECT
            sample_id,
            LEAST(total_depth, {depth_cap}) AS depth_bin,
            {_loci_expr} AS n_loci
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
                y=alt.Y("n_loci:Q", title="Number of loci"),
                color=alt.Color("sample_id:N", title="Sample"),
                tooltip=["sample_id", "depth_bin", "n_loci"],
            )
            .properties(width=700, height=350)
            .interactive()
        )
        st.altair_chart(chart, width="stretch")

        # Fraction at depth thresholds
        st.markdown("**Fraction of loci at or above depth threshold**")
        thresholds = [1, 10, 20, 30, 50, 100]
        rows = []
        for sid, grp in depth_df.groupby("sample_id"):
            total = grp["n_loci"].sum()
            row = {"sample_id": sid}
            for t in thresholds:
                below = grp.loc[grp["depth_bin"] < t, "n_loci"].sum()
                row[f">={t}x"] = round((total - below) / total, 3) if total > 0 else 0.0
            rows.append(row)
        st.dataframe(pd.DataFrame(rows), width="stretch")

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
        st.altair_chart(gc_chart, width="stretch")

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
            st.altair_chart(mapq_chart, width="stretch")

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

    _bin_n_col = "MAX(bin_n) AS bin_n," if _has_bin_n else ""

    low_df = con.execute(f"""
        SELECT
            chrom,
            pos,
            {gene_col}
            {_bin_n_col}
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
        _low_event = st.dataframe(
            low_df,
            width="stretch",
            on_select="rerun",
            selection_mode="single-row",
            key="low_coverage_table",
        )
        _low_rows = (_low_event.selection or {}).get("rows", [])
        if _low_rows:
            _row = low_df.iloc[_low_rows[0]]
            _end = int(_row.get("end", int(_row["pos"]) + 1)) if "end" in _row else int(_row["pos"]) + 1
            st.session_state["_igv_locus"] = f"{_row['chrom']}:{int(_row['pos'])}-{_end}"
            st.caption(f"Locus set to {st.session_state['_igv_locus']} — switch to IGV tab to view.")

        if _has_gene and "gene" in low_df.columns:
            st.markdown("**Affected loci per gene**")
            _loci_agg = ("bin_n", "sum") if _has_bin_n else ("pos", "count")
            by_gene = (
                low_df.groupby("gene", dropna=False)
                .agg(n_loci=_loci_agg, mean_depth=("mean_depth", "mean"))
                .reset_index()
                .sort_values("n_loci", ascending=False)
                .head(30)
            )
            _gene_sel = alt.selection_point(fields=["gene"])
            gene_bar = (
                alt.Chart(by_gene)
                .mark_bar()
                .encode(
                    x=alt.X("n_loci:Q", title="Loci below threshold"),
                    y=alt.Y("gene:N", sort="-x", title="Gene"),
                    color=alt.Color(
                        "mean_depth:Q",
                        scale=alt.Scale(scheme="reds", reverse=True),
                        title="Mean depth",
                    ),
                    opacity=alt.condition(_gene_sel, alt.value(1.0), alt.value(0.4)),
                    tooltip=["gene", "n_loci", alt.Tooltip("mean_depth:Q", format=".1f")],
                )
                .add_params(_gene_sel)
                .properties(height=max(120, 22 * len(by_gene)))
            )
            event = st.altair_chart(
                gene_bar, width="stretch", on_select="rerun", key="low_coverage_gene_bar"
            )
            # event.selection is an AttributeDictionary keyed by the Altair param name;
            # find the first non-empty list value regardless of key.
            # Only update persisted state when the event actually carries a selection;
            # leave it alone on reruns where selection is absent (e.g. filter change).
            pts = []
            if event.selection:
                for v in event.selection.values():
                    if isinstance(v, list) and v:
                        pts = v
                        break
                # A non-empty selection dict with no populated list means the user
                # explicitly deselected (clicked the same bar again).
                if not pts:
                    st.session_state.pop("_low_selected_gene", None)
            if pts:
                st.session_state["_low_selected_gene"] = pts[0].get("gene")

            selected_gene = st.session_state.get("_low_selected_gene")
            if selected_gene and selected_gene in low_df["gene"].values:
                st.markdown(f"**Records for {selected_gene}**")
                st.dataframe(
                    low_df[low_df["gene"] == selected_gene].reset_index(drop=True),
                    width="stretch",
                )

            # ── Gene coverage summary table ────────────────────────────────
            st.markdown("**Gene coverage summary**")
            st.caption("All genes ranked by mean depth. Click a row then use the button to open it in the Depth Profile.")
            _gene_summary = con.execute(f"""
                SELECT
                    gene,
                    COUNT(DISTINCT pos)                                          AS n_positions,
                    ROUND(AVG(total_depth), 1)                                  AS mean_depth,
                    ROUND(MIN(total_depth * 1.0), 0)                            AS min_depth,
                    ROUND(AVG(mean_mapq), 1)                                    AS mean_mapq,
                    ROUND(AVG(gc_content), 3)                                   AS mean_gc,
                    SUM(CASE WHEN total_depth < {depth_threshold} THEN 1 ELSE 0 END) AS n_loci_below,
                    ROUND(
                        100.0 * SUM(CASE WHEN total_depth < {depth_threshold} THEN 1 ELSE 0 END)
                        / NULLIF(COUNT(*), 0), 1
                    )                                                            AS pct_below
                FROM {table_expr}
                {where(["gene IS NOT NULL", "gene <> ''"])}
                GROUP BY gene
                ORDER BY mean_depth
            """).df()

            _gs_event = st.dataframe(
                _gene_summary,
                use_container_width=True,
                hide_index=True,
                on_select="rerun",
                selection_mode="single-row",
                key="low_gene_summary_table",
            )
            _gs_rows = (_gs_event.selection or {}).get("rows", [])
            if _gs_rows:
                _gs_gene = _gene_summary.iloc[_gs_rows[0]]["gene"]
                if st.button(f"Open **{_gs_gene}** in Depth Profile →", key="open_in_profile"):
                    st.session_state["prof_gene"] = _gs_gene
                    st.rerun()

# ──────────────────────────────────────────────────────────────────────────────
# Tab 5 — Depth Profile
# ──────────────────────────────────────────────────────────────────────────────
with tab_profile:
    st.subheader("Depth Profile")
    st.caption(
        "Aggregate depth across a gene or genomic region across all selected samples. "
        "The band shows the min–max range; the shaded area shows the IQR (p25–p75); "
        "the line shows the mean."
    )

    _prof_mode = "gene" if _has_gene else "region"

    # ACMG SF v3.2 gene list (81 genes)
    _ACMG_GENES = frozenset({
        "ACTA2", "ACTC1", "APC", "APOB", "ATP7B", "BMPR1A", "BRCA1", "BRCA2",
        "CACNA1S", "COL3A1", "DSC2", "DSG2", "DSP", "FBN1", "FBN2", "FLNC",
        "GLA", "HFE", "HNF1A", "KCNQ1", "KCNH2", "LDLR", "LMNA", "MAX",
        "MEN1", "MLH1", "MSH2", "MSH6", "MUTYH", "MYH11", "MYH7", "MYL2",
        "MYL3", "MYBPC3", "NF2", "OTC", "PALB2", "PCSK9", "PKP2", "PMS2",
        "POLD1", "POLE", "PTEN", "RB1", "RET", "RPE65", "RYR1", "RYR2",
        "SCN5A", "SCN1A", "SDHAF2", "SDHB", "SDHC", "SDHD", "SMAD3", "SMAD4",
        "STK11", "TGFBR1", "TGFBR2", "TMEM43", "TNNI3", "TNNT2", "TP53",
        "TPM1", "TSC1", "TSC2", "TTN", "TTR", "VHL", "WT1",
        # Additional SF v3.2 genes
        "EZH2", "ANKRD26", "CEBPA", "DDX41", "ELANE", "GATA2", "RUNX1",
        "SRP72", "TERC", "TERT", "TINF2",
    })

    _prof_col1, _prof_col2 = st.columns([2, 1])
    with _prof_col1:
        if _has_gene:
            _all_genes = sorted(
                con.execute(
                    f"SELECT DISTINCT gene FROM {table_expr} "
                    f"WHERE gene IS NOT NULL AND gene <> '' ORDER BY gene"
                ).df()["gene"].tolist()
            )
            _acmg_only = st.checkbox(
                "ACMG genes only", value=False, key="prof_acmg_only",
                help="Filter gene list to ACMG Secondary Findings v3.2 genes present in this dataset",
            )
            _gene_choices = (
                [g for g in _all_genes if g in _ACMG_GENES] if _acmg_only else _all_genes
            )
            if not _gene_choices:
                st.info("No ACMG genes found in this dataset.")
                st.stop()
            else:
                # If a gene was sent here from Low Coverage but isn't in the
                # current filtered list, clear the ACMG filter automatically.
                _incoming = st.session_state.get("prof_gene")
                if _incoming and _incoming not in _gene_choices:
                    st.session_state["prof_acmg_only"] = False
                    _gene_choices = _all_genes
                _prof_gene = st.selectbox("Gene", _gene_choices, key="prof_gene")
            _prof_chrom_val = None
            _prof_start_val = None
            _prof_end_val = None
        else:
            _prof_gene = None
            _prof_chrom_val = st.selectbox(
                "Chromosome", all_chroms, key="prof_chrom"
            )
    with _prof_col2:
        if not _has_gene:
            _prof_start_val = st.number_input(
                "Start", min_value=0, value=0, step=1000, key="prof_start"
            )
            _prof_end_val = st.number_input(
                "End", min_value=1, value=1_000_000, step=1000, key="prof_end"
            )

    # Build WHERE clause for profile region
    _prof_clauses: list[str] = []
    _ids = ", ".join(f"'{s.replace(chr(39), chr(39)*2)}'" for s in sample_sel) if sample_sel else ""
    if _ids:
        _prof_clauses.append(f"sample_id IN ({_ids})")
    if _has_gene and _prof_gene:
        _pg = _prof_gene.replace("'", "''")
        _prof_clauses.append(f"gene = '{_pg}'")
    elif _prof_chrom_val:
        _pc = _prof_chrom_val.replace("'", "''")
        _prof_clauses.append(f"chrom = '{_pc}'")
        _prof_clauses.append(f"pos >= {int(_prof_start_val)}")
        _prof_clauses.append(f"pos < {int(_prof_end_val)}")

    _prof_where = ("WHERE " + " AND ".join(_prof_clauses)) if _prof_clauses else ""

    # Count positions before fetching
    _prof_n_pos = con.execute(
        f"SELECT COUNT(DISTINCT pos) FROM {table_expr} {_prof_where}"
    ).fetchone()[0]

    if _prof_n_pos == 0:
        st.info("No data for the selected region/gene.")
    elif _prof_n_pos > 50_000:
        st.warning(
            f"This region spans {_prof_n_pos:,} distinct positions. "
            "Consider using a smaller region or re-running `geac coverage` with a larger `--bin-size`."
        )
    else:
        _prof_df = con.execute(
            f"""
            SELECT
                pos,
                AVG(total_depth)   AS mean_depth,
                MIN(total_depth)   AS min_depth,
                MAX(total_depth)   AS max_depth,
                PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY total_depth) AS p25_depth,
                PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY total_depth) AS p75_depth,
                COUNT(DISTINCT sample_id) AS n_samples,
                AVG(mean_mapq)     AS mean_mapq,
                AVG(frac_mapq0)    AS mean_frac_mapq0,
                AVG(frac_low_mapq) AS mean_frac_low_mapq,
                AVG(gc_content)    AS mean_gc_content
            FROM {table_expr}
            {_prof_where}
            GROUP BY pos
            ORDER BY pos
            """
        ).df()

        _ctl_col1, _ctl_col2 = st.columns([1, 2])
        with _ctl_col1:
            _color_by = st.selectbox(
                "Color line by",
                ["Mean MAPQ", "GC content"],
                key="prof_color_by",
            )
        with _ctl_col2:
            _show_samples = st.checkbox(
                "Show individual sample lines", value=False, key="prof_show_samples"
            )

        if _color_by == "Mean MAPQ":
            _color_field = "mean_mapq"
            _color_title = "Mean MAPQ"
            _color_scale = alt.Scale(scheme="redyellowgreen", domainMid=30)
            _color_fmt   = ".1f"
        else:
            _color_field = "mean_gc_content"
            _color_title = "GC content"
            _color_scale = alt.Scale(scheme="blueorange", domainMid=0.5, domain=[0, 1])
            _color_fmt   = ".1%"

        # ── Aggregate chart ────────────────────────────────────────────────
        _base = alt.Chart(_prof_df)

        _band_minmax = (
            _base.mark_area(opacity=0.12, color="#4c78a8")
            .encode(
                x=alt.X("pos:Q", title="Position"),
                y=alt.Y("min_depth:Q", title="Depth"),
                y2=alt.Y2("max_depth:Q"),
                tooltip=["pos:Q", alt.Tooltip("min_depth:Q", format=".1f"),
                         alt.Tooltip("max_depth:Q", format=".1f")],
            )
        )

        _band_iqr = (
            _base.mark_area(opacity=0.30, color="#4c78a8")
            .encode(
                x="pos:Q",
                y=alt.Y("p25_depth:Q"),
                y2=alt.Y2("p75_depth:Q"),
            )
        )

        _line_mean = (
            _base.mark_line(strokeWidth=2.5)
            .encode(
                x="pos:Q",
                y=alt.Y("mean_depth:Q"),
                color=alt.Color(
                    f"{_color_field}:Q",
                    scale=_color_scale,
                    legend=alt.Legend(title=_color_title),
                ),
                tooltip=["pos:Q",
                         alt.Tooltip("mean_depth:Q", format=".1f", title="Mean depth"),
                         alt.Tooltip("mean_mapq:Q", format=".1f", title="Mean MAPQ"),
                         alt.Tooltip("mean_frac_low_mapq:Q", format=".1%", title="Frac low MAPQ"),
                         alt.Tooltip("mean_gc_content:Q", format=".1%", title="GC content"),
                         "n_samples:Q"],
            )
        )

        _prof_chart = _band_minmax + _band_iqr + _line_mean

        # ── Optional per-sample lines ──────────────────────────────────────
        if _show_samples and sample_sel:
            _samp_df = con.execute(
                f"""
                SELECT pos, sample_id, total_depth AS depth
                FROM {table_expr}
                {_prof_where}
                ORDER BY pos
                """
            ).df()
            _samp_chart = (
                alt.Chart(_samp_df)
                .mark_line(opacity=0.5, strokeWidth=1)
                .encode(
                    x="pos:Q",
                    y=alt.Y("depth:Q"),
                    color=alt.Color("sample_id:N", legend=alt.Legend(title="Sample")),
                    tooltip=["pos:Q", "sample_id:N",
                             alt.Tooltip("depth:Q", format=".1f")],
                )
            )
            _prof_chart = _prof_chart + _samp_chart

        # ── Exon / interval shaded bands ──────────────────────────────────
        if _has_intervals and _has_gene and _prof_gene:
            _pg2 = _prof_gene.replace("'", "''")
            _ivl_cols = con.execute("PRAGMA table_info(coverage_intervals)").df()["name"].tolist()
            _has_feature_type = "feature_type" in _ivl_cols
            _has_exon_number  = "exon_number"  in _ivl_cols

            _ivl_select = ", ".join(filter(None, [
                "start", "end", "interval_name",
                "feature_type" if _has_feature_type else None,
                "exon_number"  if _has_exon_number  else None,
            ]))
            _ivl_bounds = con.execute(
                f"""
                SELECT DISTINCT {_ivl_select}
                FROM coverage_intervals
                WHERE gene = '{_pg2}'
                ORDER BY start
                """
            ).df()

            if not _ivl_bounds.empty:
                # Fill missing feature_type so color encoding always has a value
                if "feature_type" not in _ivl_bounds.columns:
                    _ivl_bounds["feature_type"] = "Exon"
                else:
                    _ivl_bounds["feature_type"] = (
                        _ivl_bounds["feature_type"].fillna("Exon")
                    )

                _ft_color = alt.Color(
                    "feature_type:N",
                    scale=alt.Scale(
                        domain=["CDS", "Exon", "UTR", "Gene"],
                        range=["#4c78a8", "#72b7b2", "#f58518", "#b279a2"],
                    ),
                    legend=alt.Legend(title="Feature type"),
                )

                _tooltip_fields = ["interval_name:N", "start:Q", "end:Q", "feature_type:N"]
                if "exon_number" in _ivl_bounds.columns:
                    _tooltip_fields.append("exon_number:Q")

                # Shaded rect bands spanning the full y extent
                _bands = (
                    alt.Chart(_ivl_bounds)
                    .mark_rect(opacity=0.15)
                    .encode(
                        x=alt.X("start:Q"),
                        x2=alt.X2("end:Q"),
                        color=_ft_color,
                        tooltip=_tooltip_fields,
                    )
                )

                # Thin border lines at start/end of each interval for crispness
                _borders = (
                    alt.Chart(_ivl_bounds)
                    .mark_rule(opacity=0.4, strokeWidth=1)
                    .encode(
                        x=alt.X("start:Q"),
                        color=_ft_color,
                    )
                )

                # Exon number labels along the top of the chart
                if "exon_number" in _ivl_bounds.columns:
                    _mid = _ivl_bounds.copy()
                    _mid["mid"] = (_mid["start"] + _mid["end"]) / 2
                    _labels = (
                        alt.Chart(_mid)
                        .mark_text(dy=-8, fontSize=9, align="center", color="#555")
                        .encode(
                            x=alt.X("mid:Q"),
                            text=alt.Text("exon_number:Q"),
                            tooltip=_tooltip_fields,
                        )
                    )
                    _prof_chart = _bands + _borders + _labels + _prof_chart
                else:
                    _prof_chart = _bands + _borders + _prof_chart

        st.altair_chart(
            _prof_chart.resolve_scale(color="independent").properties(height=320).interactive(),
            use_container_width=True,
        )

        # ── Interval-level summary ─────────────────────────────────────────
        if _has_intervals and _has_gene and _prof_gene:
            _pg3 = _prof_gene.replace("'", "''")
            _ivl_df = con.execute(
                f"""
                SELECT
                    ci.interval_name,
                    ci.start,
                    ci.end,
                    AVG(c.total_depth)  AS mean_depth,
                    MIN(c.total_depth)  AS min_depth,
                    MAX(c.total_depth)  AS max_depth,
                    PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY c.total_depth) AS p25_depth,
                    PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY c.total_depth) AS p75_depth,
                    COUNT(DISTINCT c.sample_id) AS n_samples
                FROM coverage_intervals ci
                JOIN {table_expr} c
                    ON c.pos >= ci.start AND c.pos < ci.end
                    AND c.gene = ci.gene
                WHERE ci.gene = '{_pg3}'
                {"AND c.sample_id IN (" + _ids + ")" if sample_sel else ""}
                GROUP BY ci.interval_name, ci.start, ci.end
                ORDER BY ci.start
                """
            ).df()
            if not _ivl_df.empty:
                with st.expander("Interval summary", expanded=False):
                    _ivl_bar = (
                        alt.Chart(_ivl_df)
                        .mark_bar()
                        .encode(
                            x=alt.X("interval_name:N", sort=alt.SortField("start"),
                                    title="Interval"),
                            y=alt.Y("mean_depth:Q", title="Mean depth"),
                            color=alt.Color("mean_depth:Q",
                                            scale=alt.Scale(scheme="blues"),
                                            legend=None),
                            tooltip=["interval_name:N", "start:Q", "end:Q",
                                     alt.Tooltip("mean_depth:Q", format=".1f"),
                                     alt.Tooltip("min_depth:Q", format=".1f"),
                                     alt.Tooltip("max_depth:Q", format=".1f"),
                                     "n_samples:Q"],
                        )
                        .properties(height=200)
                    )
                    st.altair_chart(_ivl_bar, use_container_width=True)
                    st.dataframe(_ivl_df.round(1), use_container_width=True, hide_index=True)

        # ── Per-sample stats ───────────────────────────────────────────────
        with st.expander("Per-sample statistics", expanded=False):
            _samp_stats = con.execute(
                f"""
                SELECT
                    sample_id,
                    COUNT(DISTINCT pos)   AS n_positions,
                    AVG(total_depth)            AS mean_depth,
                    MIN(total_depth)            AS min_depth,
                    MAX(total_depth)            AS max_depth,
                    PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY total_depth) AS p25_depth,
                    PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY total_depth) AS p75_depth
                FROM {table_expr}
                {_prof_where}
                GROUP BY sample_id
                ORDER BY sample_id
                """
            ).df()
            st.dataframe(_samp_stats.round(1), use_container_width=True, hide_index=True)

# ──────────────────────────────────────────────────────────────────────────────
# Tab 6 — IGV
# ──────────────────────────────────────────────────────────────────────────────
with tab_igv:
    st.subheader("IGV.js Genome Browser")

    # Locus input — pre-populated from Low Coverage row click if available
    _default_locus = st.session_state.get("_igv_locus", "chr20")
    igv_locus = st.text_input(
        "Locus",
        value=_default_locus,
        placeholder="e.g. chr20:1000000-1010000",
        key="igv_locus_input",
    )
    if igv_locus:
        st.session_state["_igv_locus"] = igv_locus

    if not _manifest:
        st.info("Add a manifest file in the sidebar to load BAM/CRAM tracks.")
    else:
        # Use only samples currently selected in the sidebar
        _igv_samples = tuple(s for s in sample_sel if str(s) in _manifest)
        _missing = [s for s in sample_sel if str(s) not in _manifest]
        if _missing:
            st.warning(f"Not in manifest (no track): {', '.join(str(s) for s in _missing)}")

        if not _igv_samples:
            st.info("No selected samples found in manifest.")
        else:
            _track_height = st.slider(
                "Track height (px)", min_value=100, max_value=600,
                value=250, step=50, key="igv_track_height",
            )

            if st.button("Load IGV", type="primary"):
                st.session_state["_igv_loaded"] = True

            if st.session_state.get("_igv_loaded"):
                _needs_gcs = any(
                    str(_manifest.get(str(s), {}).get("bam", "")).startswith("gs://")
                    for s in _igv_samples
                )
                _token: str | None = None
                if _needs_gcs:
                    _token, _token_err = _gcs_access_token()
                    if _token_err == "not_installed":
                        st.warning(
                            "GCS BAMs detected but the `google-auth` package is not installed. "
                            "Install it with `pip install google-auth` and restart the app.",
                            icon="⚠️",
                        )
                    elif _token_err:
                        st.warning(f"GCS auth failed: {_token_err}", icon="⚠️")
                    else:
                        st.caption("Auth token: ✓ obtained")
                tracks = _igv_tracks(
                    _igv_samples, json.dumps(_manifest),
                    oauth_token=_token, track_height=_track_height,
                )
                # IGV browser chrome (toolbar + cytoband + reference track) is ~180px;
                # each BAM/CRAM track adds track_height + ~30px for its header.
                component_height = 180 + (_track_height + 30) * len(tracks)
                components.html(
                    _igv_html(igv_locus or "chr20", tracks, oauth_token=_token),
                    height=component_height,
                    scrolling=False,
                )

# ──────────────────────────────────────────────────────────────────────────────
# Tab 6 — Intervals  (DuckDB only, gated on coverage_intervals table)
# ──────────────────────────────────────────────────────────────────────────────
if _has_intervals and tab_intervals is not None:
    with tab_intervals:
        st.subheader("Per-interval coverage")

        # ── Detect columns available in coverage_intervals ────────────────────
        _iv_cols = set(
            con.execute("DESCRIBE SELECT * FROM coverage_intervals LIMIT 0")
            .df()["column_name"]
            .tolist()
        )
        _iv_has_gene = "gene" in _iv_cols

        # Sample filter clause for coverage_intervals (no chrom/gene/on_target filter)
        def _iv_where(extra: list[str] | None = None) -> str:
            clauses: list[str] = []
            if sample_sel:
                ids = ", ".join(f"'{s.replace(chr(39), chr(39)*2)}'" for s in sample_sel)
                clauses.append(f"sample_id IN ({ids})")
            if extra:
                clauses.extend(extra)
            return ("WHERE " + " AND ".join(clauses)) if clauses else ""

        # ── View 1: Undercovered intervals ────────────────────────────────────
        st.markdown("#### Undercovered intervals")
        iv_col1, iv_col2 = st.columns(2)
        iv_depth_thresh = iv_col1.slider(
            "Mean depth threshold", 1, 500, 30, key="iv_depth_thresh",
            help="Flag intervals whose mean depth falls below this value.",
        )
        iv_min_frac = iv_col2.slider(
            "Min fraction of samples below threshold", 0.0, 1.0, 0.5, step=0.05,
            key="iv_min_frac",
            help="Only show intervals undercovered in this fraction of selected samples.",
        )
        iv_n_selected = max(1, len(sample_sel)) if sample_sel else 1
        iv_min_count = max(1, int(iv_min_frac * iv_n_selected))

        _iv_gene_col = "gene," if _iv_has_gene else ""
        _iv_gene_grp = ", gene" if _iv_has_gene else ""

        undercov_df = con.execute(f"""
            SELECT
                chrom,
                start,
                "end",
                interval_name,
                {_iv_gene_col}
                COUNT(DISTINCT sample_id)              AS n_samples_below,
                ROUND(AVG(mean_depth), 1)              AS mean_depth,
                ROUND(MIN(mean_depth), 1)              AS min_mean_depth,
                ROUND(AVG(frac_at_20x), 3)             AS mean_frac_at_20x,
                ROUND(AVG(frac_at_30x), 3)             AS mean_frac_at_30x,
                ROUND(AVG(mean_gc_content), 3)         AS mean_gc_content,
                ROUND(AVG(mean_frac_mapq0), 3)         AS mean_frac_mapq0
            FROM coverage_intervals
            {_iv_where([f"mean_depth < {iv_depth_thresh}"])}
            GROUP BY chrom, start, "end", interval_name{_iv_gene_grp}
            HAVING COUNT(DISTINCT sample_id) >= {iv_min_count}
            ORDER BY n_samples_below DESC, mean_depth
            LIMIT 500
        """).df()

        st.caption(
            f"Intervals with `mean_depth < {iv_depth_thresh}` in ≥ {iv_min_frac:.0%} of samples "
            f"({iv_min_count}/{iv_n_selected}). Showing up to 500 rows."
        )
        if undercov_df.empty:
            st.success("No undercovered intervals for the current filters.")
        else:
            st.dataframe(undercov_df, width="stretch", hide_index=True)

        st.divider()

        # ── View 2: GC bias scatter per interval ──────────────────────────────
        st.markdown("#### GC bias: mean depth vs GC content per interval")
        st.caption(
            "Each point is one interval (averaged across selected samples). "
            "Colour = mean fraction of reads with MAPQ=0 (mappability proxy)."
        )

        gc_scatter_df = con.execute(f"""
            SELECT
                interval_name,
                {_iv_gene_col}
                ROUND(AVG(mean_gc_content), 3)  AS mean_gc_content,
                ROUND(AVG(mean_depth), 2)        AS mean_depth,
                ROUND(AVG(mean_frac_mapq0), 4)   AS mean_frac_mapq0,
                COUNT(DISTINCT sample_id)        AS n_samples
            FROM coverage_intervals
            {_iv_where()}
            GROUP BY interval_name{_iv_gene_grp}
            ORDER BY mean_gc_content
        """).df()

        if gc_scatter_df.empty:
            st.info("No data for current sample selection.")
        else:
            _tooltip = [
                "interval_name",
                alt.Tooltip("mean_gc_content:Q", title="GC content", format=".2f"),
                alt.Tooltip("mean_depth:Q", title="Mean depth", format=".1f"),
                alt.Tooltip("mean_frac_mapq0:Q", title="Frac MAPQ=0", format=".3f"),
                "n_samples",
            ]
            if _iv_has_gene:
                _tooltip.insert(1, "gene")

            gc_scatter = (
                alt.Chart(gc_scatter_df)
                .mark_circle(opacity=0.7, size=60)
                .encode(
                    x=alt.X(
                        "mean_gc_content:Q",
                        title="Mean GC content",
                        scale=alt.Scale(domain=[0, 1]),
                        axis=alt.Axis(format=".0%"),
                    ),
                    y=alt.Y("mean_depth:Q", title="Mean depth"),
                    color=alt.Color(
                        "mean_frac_mapq0:Q",
                        title="Frac MAPQ=0",
                        scale=alt.Scale(scheme="orangered"),
                    ),
                    tooltip=_tooltip,
                )
                .properties(width=700, height=400)
                .interactive()
            )
            st.altair_chart(gc_scatter, width="stretch")

        st.divider()

        # ── View 3: Per-exon heatmap ──────────────────────────────────────────
        if _iv_has_gene:
            st.markdown("#### Per-exon coverage heatmap (frac at 30×)")
            st.caption(
                "Rows = genes, columns = interval names. "
                "Colour = fraction of bases covered at ≥ 30×, averaged across selected samples. "
                "Only genes with ≥ 2 intervals are shown."
            )

            heatmap_col1, heatmap_col2 = st.columns(2)
            hm_top_n = heatmap_col1.slider(
                "Max genes to display", 5, 100, 30, key="hm_top_n",
            )
            hm_sort_by = heatmap_col2.selectbox(
                "Sort genes by", ["Min frac_at_30x (worst first)", "Gene name"],
                key="hm_sort_by",
            )

            heatmap_df = con.execute(f"""
                SELECT
                    gene,
                    interval_name,
                    ROUND(AVG(frac_at_30x), 3) AS frac_at_30x
                FROM coverage_intervals
                {_iv_where(["gene IS NOT NULL"])}
                GROUP BY gene, interval_name
            """).df()

            if heatmap_df.empty:
                st.info("No gene/interval data for current selection.")
            else:
                # Keep only genes with >= 2 intervals
                gene_counts = heatmap_df.groupby("gene")["interval_name"].nunique()
                multi_genes = gene_counts[gene_counts >= 2].index
                heatmap_df = heatmap_df[heatmap_df["gene"].isin(multi_genes)]

                if heatmap_df.empty:
                    st.info("No genes with multiple intervals found.")
                else:
                    # Sort genes
                    if hm_sort_by.startswith("Min"):
                        gene_min = (
                            heatmap_df.groupby("gene")["frac_at_30x"].min()
                            .sort_values()
                            .head(hm_top_n)
                        )
                    else:
                        gene_min = (
                            heatmap_df.groupby("gene")["frac_at_30x"].min()
                            .loc[sorted(heatmap_df["gene"].unique())]
                            .head(hm_top_n)
                        )
                    heatmap_df = heatmap_df[heatmap_df["gene"].isin(gene_min.index)]
                    gene_order = list(gene_min.index)

                    heatmap_chart = (
                        alt.Chart(heatmap_df)
                        .mark_rect()
                        .encode(
                            x=alt.X("interval_name:N", title="Interval", sort=None),
                            y=alt.Y(
                                "gene:N",
                                title="Gene",
                                sort=gene_order,
                            ),
                            color=alt.Color(
                                "frac_at_30x:Q",
                                title="Frac at 30×",
                                scale=alt.Scale(
                                    domain=[0, 1],
                                    scheme="redyellowgreen",
                                ),
                            ),
                            tooltip=[
                                "gene",
                                "interval_name",
                                alt.Tooltip("frac_at_30x:Q", title="Frac at 30×", format=".3f"),
                            ],
                        )
                        .properties(
                            height=max(150, 22 * len(gene_order)),
                            width=700,
                        )
                    )
                    st.altair_chart(heatmap_chart, width="stretch")
        else:
            st.info(
                "Per-exon heatmap requires a `gene` column in `coverage_intervals`. "
                "Re-run `geac coverage` with `--gene-annotations` to populate gene names."
            )
