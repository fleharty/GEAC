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
def _gcs_access_token() -> str | None:
    """Return a GCS read-only access token from ADC, or None if unavailable."""
    try:
        import google.auth
        import google.auth.transport.requests
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
        return token
    except Exception as e:
        print(f"[IGV] token error: {e}", flush=True)
        return None
def _igv_tracks(sample_ids: tuple[str, ...], manifest_json: str, oauth_token: str | None = None) -> list[dict]:
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
            "height": 100,
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
tab_summary, tab_depth, tab_gc, tab_low, tab_igv = st.tabs(
    ["Summary", "Depth Distribution", "GC Bias", "Low Coverage", "IGV"]
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

# ──────────────────────────────────────────────────────────────────────────────
# Tab 5 — IGV
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
            if st.button("Load IGV", type="primary"):
                st.session_state["_igv_loaded"] = True

            if st.session_state.get("_igv_loaded"):
                _token = _gcs_access_token()
                st.caption(f"Auth token: {'✓ obtained' if _token else '✗ not available'}")
                tracks = _igv_tracks(_igv_samples, json.dumps(_manifest), oauth_token=_token)
                component_height = 300 + 110 * len(tracks)
                components.html(
                    _igv_html(igv_locus or "chr20", tracks, oauth_token=_token),
                    height=component_height,
                    scrolling=False,
                )
