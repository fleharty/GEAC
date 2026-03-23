import io
import zipfile
import numpy as np
import streamlit as st
import duckdb
import altair as alt
import pandas as pd
from scipy.optimize import nnls
from igv_helpers import query_distinct_samples
import geac_config

st.set_page_config(page_title="GEAC Explorer", layout="wide")
st.title("GEAC Explorer")
st.markdown(
    "**Genomic Evidence Atlas of Cohorts** — inspect alt base metrics from "
    "per-sample Parquet files or a merged cohort DuckDB database."
)

# ── Project config (geac.toml or --config flag) ───────────────────────────────
_cfg = geac_config.load()

# ── File input ────────────────────────────────────────────────────────────────
path = st.text_input(
    "Data file path",
    value=_cfg.get("data", ""),
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

# Detect whether an alt_reads table is available (only possible in DuckDB mode).
_has_alt_reads = False
if path.endswith(".duckdb"):
    try:
        con.execute("SELECT 1 FROM alt_reads LIMIT 1")
        _has_alt_reads = True
    except Exception:
        _has_alt_reads = False

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
_FILTER_KEYS = [
    "chrom_sel", "sample_sel", "batch_sel", "gene_text", "variant_sel", "vaf_range",
    "min_alt", "min_fwd_alt", "min_rev_alt",
    "min_overlap_agree", "min_overlap_disagree",
    "variant_called_sel", "variant_filter_sel", "on_target_sel",
    "homopolymer_range", "str_len_range", "min_depth", "max_depth",
    "table_limit_sel",
]

_hdr_col, _btn_col = st.sidebar.columns([2, 1])
_hdr_col.header("Filters")
if _btn_col.button("Clear all", help="Reset all filters to defaults"):
    st.session_state["chrom_sel"]          = "All"
    st.session_state["sample_sel"]         = []
    st.session_state["batch_sel"]          = []
    st.session_state["gene_text"]          = ""
    st.session_state["variant_sel"]        = ["SNV", "insertion", "deletion"]
    st.session_state["vaf_range"]          = (0.0, 1.0)
    st.session_state["min_alt"]            = 1
    st.session_state["min_fwd_alt"]        = 0
    st.session_state["min_rev_alt"]        = 0
    st.session_state["min_overlap_agree"]   = 0
    st.session_state["min_overlap_disagree"] = 0
    st.session_state["variant_called_sel"]  = "All"
    st.session_state["variant_filter_sel"] = []
    st.session_state["on_target_sel"]      = "All"
    st.session_state["homopolymer_range"]  = (0, 20)
    st.session_state["str_len_range"]      = (0, 50)
    st.session_state["min_depth"]          = 0
    st.session_state["max_depth"]          = 0
    st.session_state["table_limit_sel"]    = 500
    st.session_state.pop("family_size_range", None)
    st.session_state.pop("dist_from_end_range", None)
    st.session_state.pop("map_qual_range", None)
    st.session_state["fs_exclude_mode"]  = False
    st.session_state["dfe_exclude_mode"] = False
    st.session_state["mq_exclude_mode"]  = False
    st.rerun()

chroms = con.execute(f"SELECT DISTINCT chrom FROM {table_expr} ORDER BY chrom").df()["chrom"].tolist()
samples = con.execute(f"SELECT DISTINCT sample_id FROM {table_expr} ORDER BY sample_id").df()["sample_id"].tolist()

chrom_sel = st.sidebar.selectbox("Chromosome", ["All"] + chroms, key="chrom_sel")
if "sample_sel" not in st.session_state:
    st.session_state["sample_sel"] = []
sample_sel = st.sidebar.multiselect("Samples (blank = all)", samples, key="sample_sel")

_schema_cols = set(con.execute(f"DESCRIBE SELECT * FROM {table_expr} LIMIT 0").df()["column_name"].tolist())

def _has_data(col: str) -> bool:
    """True iff col exists in the schema AND has at least one non-null value."""
    if col not in _schema_cols:
        return False
    return con.execute(f"SELECT COUNT(*) FROM {table_expr} WHERE {col} IS NOT NULL").fetchone()[0] > 0

if _has_data("batch"):
    _batches = con.execute(f"SELECT DISTINCT batch FROM {table_expr} WHERE batch IS NOT NULL ORDER BY batch").df()["batch"].tolist()
    if "batch_sel" not in st.session_state:
        st.session_state["batch_sel"] = []
    batch_sel = st.sidebar.multiselect("Batch (blank = all)", _batches, key="batch_sel")
else:
    batch_sel = []

_genes_available = _has_data("gene")
if _genes_available:
    gene_text = st.sidebar.text_input("Gene (exact match, blank = all)", "", key="gene_text")
else:
    gene_text = ""
    st.sidebar.caption("Gene filter unavailable — run geac collect with --gene-annotations to enable.")
if "variant_sel" not in st.session_state:
    st.session_state["variant_sel"] = ["SNV", "insertion", "deletion"]
variant_sel = st.sidebar.multiselect(
    "Variant type",
    ["SNV", "insertion", "deletion"],
    key="variant_sel",
)
vaf_range = st.sidebar.slider("VAF range", 0.0, 1.0, (0.0, 1.0), step=0.01, key="vaf_range")
min_alt = st.sidebar.number_input("Min alt count", min_value=1, max_value=10000, value=1, step=1, key="min_alt")
min_fwd_alt = st.sidebar.number_input("Min fwd alt count (0 = no minimum)", min_value=0, max_value=10000, value=0, step=1, key="min_fwd_alt")
min_rev_alt = st.sidebar.number_input("Min rev alt count (0 = no minimum)", min_value=0, max_value=10000, value=0, step=1, key="min_rev_alt")
min_overlap_agree    = st.sidebar.number_input("Min overlap alt agree (0 = no minimum)",    min_value=0, max_value=10000, value=0, step=1, key="min_overlap_agree")
min_overlap_disagree = st.sidebar.number_input("Min overlap alt disagree (0 = no minimum)", min_value=0, max_value=10000, value=0, step=1, key="min_overlap_disagree")
variant_called_sel = st.sidebar.selectbox("Variant called", ["All", "Yes", "No", "Unknown (no VCF/TSV)"], key="variant_called_sel")
_vf_has_data = _has_data("variant_filter")
if _vf_has_data:
    _vf_options = con.execute(
        f"SELECT DISTINCT variant_filter FROM {table_expr} WHERE variant_filter IS NOT NULL ORDER BY variant_filter"
    ).df()["variant_filter"].tolist()
    if "variant_filter_sel" not in st.session_state:
        st.session_state["variant_filter_sel"] = []
    variant_filter_sel = st.sidebar.multiselect(
        "Variant filter (blank = all)",
        _vf_options,
        key="variant_filter_sel",
        help="Filter values from the VCF FILTER field or variants TSV. Select one or more values to restrict to those loci.",
    )
else:
    variant_filter_sel = []
    st.sidebar.caption("Variant filter unavailable — run geac collect with --vcf or --variants-tsv to enable.")
on_target_sel = st.sidebar.selectbox("Target bases", ["All", "On target", "Off target"], key="on_target_sel")

_repeat_cols_present = _has_data("homopolymer_len")
if _repeat_cols_present:
    homopolymer_range = st.sidebar.slider("Homopolymer length range", 0, 20, (0, 20), step=1, key="homopolymer_range")
    str_len_range     = st.sidebar.slider("STR length range",         0, 50, (0, 50), step=1, key="str_len_range")
else:
    homopolymer_range = (0, 20)
    str_len_range     = (0, 50)
    st.sidebar.caption("Repeat filters unavailable — run geac collect with a newer build to enable.")
min_depth = st.sidebar.number_input("Min depth (0 = no minimum)", min_value=0, value=0, step=1, key="min_depth")
max_depth = st.sidebar.number_input("Max depth (0 = no maximum)", min_value=0, value=0, step=1, key="max_depth")


# ── Per-read filters (only when alt_reads table is present) ───────────────────
_reads_conditions = []
if _has_alt_reads:
    _reads_maxes = con.execute("""
        SELECT
            MAX(family_size),
            COALESCE(MAX(dist_from_read_end), 300),
            COALESCE(MAX(map_qual), 60),
            COUNT(insert_size) > 0
        FROM alt_reads
    """).fetchone()
    _fs_max_raw = _reads_maxes[0]   # None if all NULL
    _dfe_max    = int(_reads_maxes[1])
    _mq_max     = int(_reads_maxes[2])
    _fs_has_data = _fs_max_raw is not None
    _fs_max = int(_fs_max_raw) if _fs_has_data else 0
    _is_has_data = bool(_reads_maxes[3])

    st.sidebar.divider()
    st.sidebar.subheader("Per-read filters")
    st.sidebar.caption(
        "When active, alt_count and VAF are re-computed from the reads table "
        "using only reads that pass these filters. Loci with no passing reads are excluded."
    )

    if _fs_has_data:
        _fs_slider_col, _fs_toggle_col = st.sidebar.columns([3, 1])
        with _fs_slider_col:
            family_size_range = st.slider(
                "Family size range",
                min_value=0, max_value=_fs_max, value=(0, _fs_max), step=1,
                key="family_size_range",
                help="family_size = cD tag (total raw read count per molecule).",
            )
        with _fs_toggle_col:
            st.write("Mode")
            fs_exclude_mode = st.toggle(
                "Excl.",
                key="fs_exclude_mode",
                help="Off = include only reads within this range. "
                     "On = exclude reads within this range (keep reads outside it).",
            )
    else:
        family_size_range = (0, 0)
        fs_exclude_mode = False
        st.sidebar.caption("Family size unavailable — BAM has no fgbio cD tag.")

    _dfe_slider_col, _dfe_toggle_col = st.sidebar.columns([3, 1])
    with _dfe_slider_col:
        dist_from_end_range = st.slider(
            "Dist from read end",
            min_value=0, max_value=_dfe_max, value=(0, _dfe_max), step=1,
            key="dist_from_end_range",
            help="Raise the lower bound to exclude reads clustered at read ends (a common artefact).",
        )
    with _dfe_toggle_col:
        st.write("Mode")
        dfe_exclude_mode = st.toggle(
            "Excl.",
            key="dfe_exclude_mode",
            help="Off = include only reads within this range. "
                 "On = exclude reads within this range (keep reads outside it).",
        )

    _mq_slider_col, _mq_toggle_col = st.sidebar.columns([3, 1])
    with _mq_slider_col:
        map_qual_range = st.slider(
            "Mapping quality range",
            min_value=0, max_value=_mq_max, value=(0, _mq_max), step=1,
            key="map_qual_range",
            help="Filter alt-supporting reads by mapping quality (MAPQ).",
        )
    with _mq_toggle_col:
        st.write("Mode")
        mq_exclude_mode = st.toggle(
            "Excl.",
            key="mq_exclude_mode",
            help="Off = include only reads within this range. "
                 "On = exclude reads within this range (keep reads outside it).",
        )

    _IS_MIN, _IS_MAX = 20, 500
    if _is_has_data:
        insert_size_range = st.sidebar.slider(
            "Insert size range",
            min_value=_IS_MIN, max_value=_IS_MAX,
            value=(_IS_MIN, _IS_MAX), step=1,
            key="insert_size_range",
            help="Filter alt-supporting reads by template insert size (|TLEN|).",
        )
        _is_lo, _is_hi = insert_size_range
        if _is_lo == _IS_MIN and _is_hi == _IS_MAX:
            st.sidebar.caption(
                "Insert size: no filter active — reads with any insert size "
                "(including <20 and >500 bp) are accepted."
            )
        else:
            st.sidebar.caption(
                f"Insert size: keeping only reads with insert size "
                f"between {_is_lo} and {_is_hi} bp. "
                "Reads with very small or very large inserts (and unpaired reads) are excluded."
            )
    else:
        insert_size_range = (_IS_MIN, _IS_MAX)

    _fs_lo, _fs_hi = family_size_range
    _dfe_lo, _dfe_hi = dist_from_end_range
    _mq_lo, _mq_hi = map_qual_range
    _is_lo, _is_hi = insert_size_range

    if _fs_has_data and (_fs_lo > 0 or _fs_hi < _fs_max):
        if fs_exclude_mode:
            _reads_conditions.append(
                f"(family_size IS NULL OR family_size < {_fs_lo} OR family_size > {_fs_hi})"
            )
        else:
            _reads_conditions.append(f"family_size BETWEEN {_fs_lo} AND {_fs_hi}")

    if _dfe_lo > 0 or _dfe_hi < _dfe_max:
        if dfe_exclude_mode:
            _reads_conditions.append(
                f"(dist_from_read_end < {_dfe_lo} OR dist_from_read_end > {_dfe_hi})"
            )
        else:
            _reads_conditions.append(f"dist_from_read_end BETWEEN {_dfe_lo} AND {_dfe_hi}")

    if _mq_lo > 0 or _mq_hi < _mq_max:
        if mq_exclude_mode:
            _reads_conditions.append(
                f"(map_qual < {_mq_lo} OR map_qual > {_mq_hi})"
            )
        else:
            _reads_conditions.append(f"map_qual BETWEEN {_mq_lo} AND {_mq_hi}")

    if _is_has_data and (_is_lo > _IS_MIN or _is_hi < _IS_MAX):
        _reads_conditions.append(f"insert_size BETWEEN {_is_lo} AND {_is_hi}")
else:
    family_size_range = (0, 0)
    dist_from_end_range = (0, 0)
    map_qual_range = (0, 0)
    insert_size_range = (0, 0)
    fs_exclude_mode = False
    dfe_exclude_mode = False
    mq_exclude_mode = False
    _fs_lo = _fs_hi = _dfe_lo = _dfe_hi = _mq_lo = _mq_hi = _is_lo = _is_hi = 0

# When per-read filters are active, redefine table_expr as a JOIN subquery that
# replaces alt_count with the filtered fragment count. Every downstream query
# (counts, VAF, plots, IGV, drill-down) automatically sees the re-aggregated values.
_reads_active = bool(_reads_conditions)
if _reads_active:
    _reads_where = " AND ".join(_reads_conditions)
    table_expr = f"""(
        SELECT
            ab.* EXCLUDE (alt_count),
            ar_agg.filtered_alt_count AS alt_count,
            ROUND(ab.alt_count * 1.0 / ab.total_depth, 4) AS original_vaf
        FROM alt_bases ab
        INNER JOIN (
            SELECT sample_id, chrom, pos, alt_allele, COUNT(*) AS filtered_alt_count
            FROM alt_reads
            WHERE {_reads_where}
            GROUP BY sample_id, chrom, pos, alt_allele
        ) ar_agg ON ab.sample_id = ar_agg.sample_id
                 AND ab.chrom = ar_agg.chrom
                 AND ab.pos = ar_agg.pos
                 AND ab.alt_allele = ar_agg.alt_allele
    )"""

# ── IGV integration (sidebar) ─────────────────────────────────────────────────
st.sidebar.divider()
st.sidebar.header("IGV Integration")
auto_launch_igv = st.sidebar.checkbox(
    "Auto-launch IGV",
    value=_cfg.get("auto_launch_igv", False),
    help="Write session files to a temp directory and load them in IGV automatically. "
         "Uses IGV's REST API (port 60151) if IGV is already running, otherwise launches IGV.",
)

import os as _os
_default_manifest = (
    _cfg.get("manifest")
    or _os.path.join(_os.path.dirname(_os.path.abspath(path)), "manifest.tsv")
)

manifest_path = st.sidebar.text_input(
    "Manifest file (optional)",
    value=_default_manifest,
    help="Tab-separated file with columns: collaborator_sample_id, duplex_output_bam, duplex_output_bam_index, final_annotated_variants",
)
_genome_options = ["hg19", "hg38", "mm10", "mm39", "other"]
_cfg_genome = _cfg.get("genome_build") or _cfg.get("genome", "hg19")
if _cfg_genome in _genome_options:
    _genome_default_idx = _genome_options.index(_cfg_genome)
else:
    _genome_default_idx = _genome_options.index("other")
genome = st.sidebar.selectbox("Genome", _genome_options, index=_genome_default_idx)
if genome == "other":
    genome = st.sidebar.text_input("Genome ID", value=_cfg_genome if _cfg_genome not in _genome_options else "hg38")

@st.cache_data
def load_manifest(p: str) -> dict:
    mdf = pd.read_csv(p.strip(), sep="\t")
    result = {}
    for row in mdf.itertuples(index=False):
        bai = str(row.duplex_output_bam_index) if hasattr(row, "duplex_output_bam_index") and pd.notna(row.duplex_output_bam_index) else None
        variants = str(row.final_annotated_variants) if hasattr(row, "final_annotated_variants") and pd.notna(row.final_annotated_variants) else None
        result[str(row.collaborator_sample_id)] = {"bam": str(row.duplex_output_bam), "bai": bai, "variants_tsv": variants}
    return result

manifest = {}
if manifest_path and manifest_path.strip():
    try:
        manifest = load_manifest(manifest_path.strip())
        st.sidebar.success(f"{len(manifest):,} samples loaded from manifest")
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


def launch_igv_session(session_xml: str, bed: str) -> str:
    """Write session.xml and positions.bed to a temp dir and load in IGV.

    Tries the IGV REST API (localhost:60151) first — works if IGV is already
    running. If the connection is refused, launches IGV via subprocess.
    Returns a status message suitable for st.info / st.success / st.error.
    """
    import tempfile, subprocess, urllib.request, urllib.error

    tmp = tempfile.mkdtemp(prefix="geac_igv_")
    session_path = _os.path.join(tmp, "session.xml")
    bed_path     = _os.path.join(tmp, "positions.bed")
    with open(session_path, "w") as f:
        f.write(session_xml)
    with open(bed_path, "w") as f:
        f.write(bed)

    # Try REST API first
    url = f"http://localhost:60151/load?file={urllib.request.pathname2url(session_path)}&merge=false"
    try:
        urllib.request.urlopen(url, timeout=8)
        return f"Session loaded into running IGV instance."
    except (urllib.error.URLError, TimeoutError, OSError):
        pass

    # Fall back to launching IGV
    candidates = ["igv", "igv.sh"]
    if _os.path.exists("/Applications/IGV.app"):
        candidates = ["open", "-a", "IGV", session_path]
        try:
            subprocess.Popen(candidates)
            return f"Launched IGV with session at {session_path}"
        except FileNotFoundError:
            pass
    else:
        for cmd in ["igv", "igv.sh"]:
            try:
                subprocess.Popen([cmd, session_path])
                return f"Launched IGV with session at {session_path}"
            except FileNotFoundError:
                continue

    return (
        f"Could not launch IGV automatically. Session written to:\n{session_path}\n"
        "Open it manually in IGV, or install IGV to /Applications/IGV.app."
    )


def make_igv_session(df: pd.DataFrame, manifest: dict, genome: str) -> str:
    sample_ids = df["sample_id"].unique().tolist()
    first = df.sort_values(["chrom", "pos"]).iloc[0]
    locus = f"{first['chrom']}:{max(0, int(first['pos']) - 99)}-{int(first['pos']) + 101}"

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

    _extra_w = " AND ".join(conditions + extra_conditions)
    sample_ids = query_distinct_samples(con, table_expr, _extra_w)
    n = len(sample_ids)
    cap_samples = sample_ids[:IGV_CAP]

    missing = [sid for sid in sample_ids if str(sid) not in manifest]
    if missing:
        st.warning(
            f"Sample(s) not found in manifest — no BAM track will be added for: "
            f"{', '.join(str(s) for s in missing)}"
        )

    if n > IGV_CAP:
        _total_records = con.execute(
            f"SELECT COUNT(*) FROM {table_expr} WHERE {_extra_w}"
        ).fetchone()[0]
        st.warning(
            f"{n} samples in this selection. IGV session capped at {IGV_CAP}. "
            "Select specific samples below, or check the box to load all "
            "(may crash IGV — you're on your own)."
        )
        chosen = st.multiselect(
            "Samples to include in IGV session",
            options=sample_ids,
            default=cap_samples,
            key=f"{key}_sample_pick",
        )
        cap_samples = chosen if chosen else cap_samples
        _cap_list = ", ".join(f"'{s}'" for s in cap_samples)
        _chosen_records = con.execute(
            f"SELECT COUNT(*) FROM {table_expr} WHERE {_extra_w} "
            f"AND sample_id IN ({_cap_list})"
        ).fetchone()[0]
        st.caption(f"{_chosen_records:,} / {_total_records:,} records from selected samples")
        _override_label = f"Load all {n} samples instead"
        if n > 10:
            _override_label += "  :red[(Too many samples, it's a bad idea to click this!)]"
        if st.checkbox(_override_label, key=f"{key}_override"):
            cap_samples = sample_ids

    if st.button("Prepare IGV session", key=f"{key}_prepare"):
        _sample_clause = "sample_id IN ({})".format(
            ", ".join(f"'{s}'" for s in cap_samples)
        )
        w = " AND ".join(conditions + extra_conditions + [_sample_clause])
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

        igv_df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
        bed = make_bed(igv_df)
        session = make_igv_session(igv_df, manifest, genome)
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("session.xml", session)
            zf.writestr("positions.bed", bed)
        st.session_state[f"{key}_igv"]     = buf.getvalue()
        st.session_state[f"{key}_session"] = session
        st.session_state[f"{key}_bed"]     = bed

        if auto_launch_igv:
            msg = launch_igv_session(session, bed)
            st.info(msg)

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
if batch_sel:
    b_list = ", ".join(f"'{b}'" for b in batch_sel)
    conditions.append(f"batch IN ({b_list})")
if variant_sel:
    t_list = ", ".join(f"'{t}'" for t in variant_sel)
    conditions.append(f"variant_type IN ({t_list})")
if min_fwd_alt > 0:
    conditions.append(f"fwd_alt_count >= {min_fwd_alt}")
if min_rev_alt > 0:
    conditions.append(f"rev_alt_count >= {min_rev_alt}")
if min_overlap_agree > 0:
    conditions.append(f"overlap_alt_agree >= {min_overlap_agree}")
if min_overlap_disagree > 0:
    conditions.append(f"overlap_alt_disagree >= {min_overlap_disagree}")
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
if variant_filter_sel:
    vf_list = ", ".join(f"'{v.replace(chr(39), chr(39)*2)}'" for v in variant_filter_sel)
    conditions.append(f"variant_filter IN ({vf_list})")
if "on_target" in _schema_cols:
    if on_target_sel == "On target":
        conditions.append("on_target = true")
    elif on_target_sel == "Off target":
        conditions.append("on_target = false")
if gene_text.strip() and "gene" in _schema_cols:
    _gene_escaped = gene_text.strip().replace("'", "''")
    conditions.append(f"gene = '{_gene_escaped}'")

if _repeat_cols_present:
    conditions.append(f"homopolymer_len BETWEEN {homopolymer_range[0]} AND {homopolymer_range[1]}")
    conditions.append(f"str_len BETWEEN {str_len_range[0]} AND {str_len_range[1]}")

where = " AND ".join(conditions)

# ── Summary stats display ──────────────────────────────────────────────────────
fstats = con.execute(f"""
    SELECT
        COUNT(*)                                            AS n_records,
        COUNT(DISTINCT sample_id)                           AS n_samples,
        COALESCE(SUM(alt_count), 0)                         AS total_alt_bases,
        COALESCE(ROUND(AVG(alt_count * 1.0 / total_depth), 4), 0) AS mean_vaf,
        COALESCE(ROUND(AVG(total_depth), 1), 0)             AS mean_depth,
        COUNT(*) FILTER (WHERE variant_called IS NOT NULL)  AS n_annotated,
        COUNT(*) FILTER (WHERE variant_called = true)       AS n_called
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

def query_records(extra: list[str] = [], limit: int | None = None) -> pd.DataFrame:
    """Query records with current filters plus any extra conditions."""
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

st.info(f"**{total_count:,}** records match the current filters.")

if _reads_active:
    _active_parts = []
    if _fs_has_data and (_fs_lo > 0 or _fs_hi < _fs_max):
        _mode = "excluding" if fs_exclude_mode else "including only"
        _active_parts.append(f"family size: {_mode} {_fs_lo}–{_fs_hi}")
    if _dfe_lo > 0 or _dfe_hi < _dfe_max:
        _mode = "excluding" if dfe_exclude_mode else "including only"
        _active_parts.append(f"dist from read end: {_mode} {_dfe_lo}–{_dfe_hi}")
    if _mq_lo > 0 or _mq_hi < _mq_max:
        _mode = "excluding" if mq_exclude_mode else "including only"
        _active_parts.append(f"map qual: {_mode} {_mq_lo}–{_mq_hi}")
    st.warning(
        f"**Per-read filters active** ({'; '.join(_active_parts)}). "
        "alt_count and VAF are re-aggregated from reads passing the filter. "
        "original_vaf shows the unfiltered VAF for comparison. "
        "ref_count, total_depth, and strand/overlap columns still reflect unfiltered locus-level values."
    )

# ── Data table ────────────────────────────────────────────────────────────────
_tbl_limit_options = [100, 500, 1000, 5000, 10000, 50000, "All"]
_tbl_limit_sel = st.selectbox(
    "Table row limit",
    _tbl_limit_options,
    index=1,
    key="table_limit_sel",
    help="Limits rows shown in the data table below. Has no effect on plots — they always use the full filtered dataset.",
)
_tbl_limit = None if _tbl_limit_sel == "All" else int(_tbl_limit_sel)

df = query_records(limit=_tbl_limit)

_table_cols = [
    c for c in [
        "sample_id", "chrom", "pos", "ref_allele", "alt_allele",
        "variant_type", "vaf", *( ["original_vaf"] if _reads_active else []), "alt_count", "ref_count", "total_depth",
        "fwd_alt_count", "rev_alt_count", "overlap_alt_agree",
        "overlap_alt_disagree", "variant_called", "variant_filter", "on_target", "gene",
    ]
    if c in df.columns
]

with st.expander("Data table", expanded=True):
    _tbl_caption = (
        f"Showing {len(df):,} of {total_count:,} rows. Increase the table row limit above to see more."
        if len(df) < total_count else
        f"Showing all {total_count:,} rows."
    )
    st.caption(_tbl_caption)
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

    # ── Per-read detail (only when alt_reads table is present) ────────────────
    if _has_alt_reads:
        _reads_df = con.execute(f"""
            SELECT
                sample_id,
                alt_allele,
                dist_from_read_start,
                dist_from_read_end,
                read_length,
                ab_count,
                ba_count,
                family_size,
                base_qual,
                map_qual
            FROM alt_reads
            WHERE chrom = '{_chrom}' AND pos = {_pos}
            ORDER BY sample_id, alt_allele, family_size DESC NULLS LAST
        """).df()

        if _reads_df.empty:
            st.caption("No per-read detail available for this locus.")
        else:
            st.markdown(f"**Per-read detail** — {len(_reads_df):,} alt-supporting reads across {_reads_df['sample_id'].nunique()} sample(s)")

            _reads_summary = (
                _reads_df
                .groupby(["sample_id", "alt_allele"])
                .agg(
                    n_reads=("dist_from_read_end", "count"),
                    median_dist_from_end=("dist_from_read_end", "median"),
                    median_family_size=("family_size", "median"),
                    min_family_size=("family_size", "min"),
                    max_family_size=("family_size", "max"),
                    mean_base_qual=("base_qual", "mean"),
                )
                .reset_index()
                .round(1)
            )
            st.caption("Summary by sample / allele")
            st.dataframe(_reads_summary, use_container_width=True, hide_index=True)

            with st.expander("Individual reads"):
                st.dataframe(_reads_df, use_container_width=True, hide_index=True)

# ── Plots ─────────────────────────────────────────────────────────────────────
tab1, tab2, tab3, tab4, tab_cohort, tab_reads, tab_sig_cmp = st.tabs(["VAF distribution", "Error spectrum", "Strand bias", "Overlap agreement", "Cohort", "Reads", "Sig. Comparison (Exp.)"])

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

# ── SBS96 helpers (used by Error Spectrum tab and Cohort tab) ─────────────────
_COMP = str.maketrans('ACGT', 'TGCA')
_SBS_MUT_TYPES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
_SBS_COLORS    = {
    "C>A": "#1BBDEB", "C>G": "#808080", "C>T": "#E22926",
    "T>A": "#CBCACB", "T>C": "#97D54C", "T>G": "#ECC6C5",
}
_SBS_ORDER = [
    f"{b5}[{mt}]{b3}"
    for mt in _SBS_MUT_TYPES
    for b5 in "ACGT"
    for b3 in "ACGT"
]

def _sbs_label(trinuc_context, ref_allele, alt_allele):
    """Convert raw trinuc/ref/alt into a pyrimidine-normalised SBS96 label."""
    ctx, r, a = trinuc_context, ref_allele, alt_allele
    if not all(b in 'ACGT' for b in (ctx + r + a)):
        return None
    if r in ('A', 'G'):
        ctx = ctx[::-1].translate(_COMP)
        r = r.translate(_COMP)
        a = a.translate(_COMP)
    return f"{ctx[0]}[{r}>{a}]{ctx[2]}"


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


@st.cache_data
def _load_cosmic(p: str) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t", index_col=0)


with tab2:
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
            raw["sbs_label"] = raw.apply(
                lambda row: _sbs_label(row["trinuc_context"], row["ref_allele"], row["alt_allele"]),
                axis=1,
            )
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
            _total_snvs = spec96["count"].sum()
            spec96["fraction"] = spec96["count"] / _total_snvs if _total_snvs > 0 else 0.0

            _sbs_y_mode = st.radio(
                "Y axis", ["Fraction", "Count"],
                horizontal=True, key="sbs_y_mode",
            )
            _sbs_use_fraction = _sbs_y_mode == "Fraction"
            _sbs_y_field = "fraction" if _sbs_use_fraction else "count"
            _sbs_y_title = "Fraction of SNVs" if _sbs_use_fraction else "Count"
            _sbs_y_fmt   = ".3f" if _sbs_use_fraction else "d"

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
                        alt.Y(f"{_sbs_y_field}:Q", title=_sbs_y_title,
                              **({"axis": alt.Axis(format=".3f")} if _sbs_use_fraction else {})),
                        opacity=alt.condition(sel_param, alt.value(1.0), alt.value(0.4)),
                        tooltip=[
                            "sbs_label:N",
                            alt.Tooltip(f"{_sbs_y_field}:Q", title=_sbs_y_title, format=_sbs_y_fmt),
                        ],
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
                clicked_labels = [p.get("sbs_label") for p in pts if p.get("sbs_label")]
                if clicked_labels:
                    matching = raw[raw["sbs_label"].isin(clicked_labels)][["trinuc_context", "ref_allele", "alt_allele"]]
                    if not matching.empty:
                        or_clauses = " OR ".join(
                            f"(trinuc_context = '{r.trinuc_context}' AND ref_allele = '{r.ref_allele}' AND alt_allele = '{r.alt_allele}')"
                            for r in matching.itertuples(index=False)
                        )
                        extra_cond = f"variant_type = 'SNV' AND ({or_clauses})"
                        sel = query_records([extra_cond])
                        label_str = ", ".join(clicked_labels)
                        st.caption(f"{len(sel):,} records matching {len(clicked_labels)} selected context(s): {label_str}")
                        st.dataframe(sel[_table_cols], use_container_width=True)
                        igv_buttons([extra_cond], sel, key=f"sbs_{'_'.join(clicked_labels)}")

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
                value=_cfg.get("cosmic", ""),
                placeholder="/path/to/COSMIC_v3.4_SBS_GRCh37.txt",
                key="cosmic_path",
            )

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

                        # ── Observed vs reconstructed overlay ─────────────
                        recon_total = reconstructed.sum()
                        recon_df = spec96[["sbs_label", "mut_type"]].copy()
                        recon_df["recon_count"] = reconstructed
                        recon_df["recon_frac"] = (
                            reconstructed / recon_total if recon_total > 0 else 0.0
                        )
                        _recon_y = "recon_frac" if _sbs_use_fraction else "recon_count"
                        _recon_fmt = ".3f" if _sbs_use_fraction else "d"

                        _overlay_charts = []
                        for _mt in _SBS_MUT_TYPES:
                            _obs_sub   = spec96[spec96["mut_type"] == _mt].copy()
                            _recon_sub = recon_df[recon_df["mut_type"] == _mt].copy()
                            _order = [lbl for lbl in _SBS_ORDER if f"[{_mt}]" in lbl]
                            _bars = (
                                alt.Chart(_obs_sub)
                                .mark_bar(color=_SBS_COLORS[_mt], opacity=0.7)
                                .encode(
                                    alt.X("sbs_label:N", sort=_order, title=None,
                                          axis=alt.Axis(labelAngle=-90, labelFontSize=8)),
                                    alt.Y(f"{_sbs_y_field}:Q", title=_sbs_y_title,
                                          **({"axis": alt.Axis(format=".3f")} if _sbs_use_fraction else {})),
                                    tooltip=[
                                        "sbs_label:N",
                                        alt.Tooltip(f"{_sbs_y_field}:Q", title="Observed", format=_sbs_y_fmt),
                                    ],
                                )
                            )
                            _line = (
                                alt.Chart(_recon_sub)
                                .mark_point(color="black", size=15, filled=True, opacity=0.85)
                                .encode(
                                    alt.X("sbs_label:N", sort=_order),
                                    alt.Y(f"{_recon_y}:Q"),
                                    tooltip=[
                                        "sbs_label:N",
                                        alt.Tooltip(f"{_recon_y}:Q", title="Reconstructed", format=_recon_fmt),
                                    ],
                                )
                            )
                            _overlay_charts.append(
                                alt.layer(_bars, _line)
                                .properties(
                                    title=alt.TitleParams(_mt, color=_SBS_COLORS[_mt], fontSize=13, fontWeight="bold"),
                                    width=150, height=140,
                                )
                            )

                        overlay_chart = (
                            alt.concat(*_overlay_charts, columns=3)
                            .resolve_scale(y="shared")
                            .properties(title=alt.TitleParams(
                                "Observed vs Reconstructed Spectrum (bars = observed, dots = reconstructed)", fontSize=15,
                            ))
                        )
                        st.altair_chart(overlay_chart, use_container_width=True)
                        st.caption(
                            "Bars = observed spectrum. Black dots = COSMIC reconstruction (NNLS fit). "
                            "Contexts where the line deviates from the bars are poorly explained by the selected signatures."
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
    _sb_col1, _sb_col2 = st.columns(2)
    _sb_scale = _sb_col1.radio(
        "Axis scale", ["Linear", "log1p"], horizontal=True, key="sb_scale",
        help="log1p = log(1 + x), compresses large counts while keeping zeros visible.",
    )
    _use_log1p = _sb_scale == "log1p"

    _color_options = ["Variant type", "Sample"]
    if _has_data("batch"):
        _color_options.append("Batch")
    if _has_data("on_target"):
        _color_options.append("On target")
    if _has_data("variant_called"):
        _color_options.append("Called variant")
    _color_by = _sb_col2.radio(
        "Color by", _color_options, horizontal=True, key="sb_color_by",
    )

    _sb_opt_cols = ", ".join(filter(None, [
        "batch"          if _has_data("batch")          else None,
        "on_target"      if _has_data("on_target")      else None,
        "variant_called" if _has_data("variant_called") else None,
        "gene"           if _genes_available             else None,
    ]))
    _SB_SAMPLE_THRESHOLD = 5_000
    _sb_needs_sample = total_count > _SB_SAMPLE_THRESHOLD
    _sb_show_all = False
    if _sb_needs_sample:
        st.warning(
            f"**{total_count:,} loci** exceed the scatter plot threshold of "
            f"{_SB_SAMPLE_THRESHOLD:,}. Showing a random sample of "
            f"{_SB_SAMPLE_THRESHOLD:,} points.",
            icon="⚠️",
        )
        _sb_show_all = st.checkbox(
            "Show all loci (may be slow)",
            value=False,
            key="sb_show_all",
        )

    _sb_sql_base = f"""
        SELECT sample_id, chrom, pos, ref_allele, alt_allele, variant_type,
               fwd_alt_count, rev_alt_count,
               ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
               {(', ' + _sb_opt_cols) if _sb_opt_cols else ''}
        FROM {table_expr}
        WHERE {where}
    """
    if _sb_needs_sample and not _sb_show_all:
        sample_df = con.execute(
            f"{_sb_sql_base} USING SAMPLE reservoir({_SB_SAMPLE_THRESHOLD}) REPEATABLE(42)"
        ).df()
    else:
        if _sb_show_all:
            st.warning("Loading all loci — this may take a moment and slow down your browser.", icon="🐢")
        sample_df = con.execute(_sb_sql_base).df()

    # Round linear values used as tick labels when in log1p mode.
    _log1p_ticks_linear = [1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000]

    if _use_log1p:
        sample_df["fwd_plot"] = np.log1p(sample_df["fwd_alt_count"])
        sample_df["rev_plot"] = np.log1p(sample_df["rev_alt_count"])
        _x_title = "Forward alt reads"
        _y_title = "Reverse alt reads"
    else:
        sample_df["fwd_plot"] = sample_df["fwd_alt_count"]
        sample_df["rev_plot"] = sample_df["rev_alt_count"]
        _x_title = "Forward alt reads"
        _y_title = "Reverse alt reads"

    max_val = max(
        float(sample_df["fwd_plot"].max()) if len(sample_df) > 0 else 50,
        float(sample_df["rev_plot"].max()) if len(sample_df) > 0 else 50,
        1.0,
    )

    # 95% CI band — computed in original (linear) space, then optionally
    # transformed so the band correctly follows the data points.
    _x_lin = np.arange(0, int(sample_df["fwd_alt_count"].max()) + 1, dtype=float)
    _z = 1.96
    _s_lo = (-_z + np.sqrt(_z**2 + 8 * _x_lin)) / 2
    _s_hi = ( _z + np.sqrt(_z**2 + 8 * _x_lin)) / 2
    _rev_min_lin = np.maximum(_s_lo**2 - _x_lin, 0.0)
    _rev_max_lin = _s_hi**2 - _x_lin

    if _use_log1p:
        _ci_band = pd.DataFrame({
            "fwd":     np.log1p(_x_lin),
            "rev_min": np.log1p(_rev_min_lin),
            "rev_max": np.log1p(_rev_max_lin),
        })
        _diag = pd.DataFrame({"fwd": [0.0, max_val], "rev": [0.0, max_val]})
    else:
        _ci_band = pd.DataFrame({
            "fwd":     _x_lin,
            "rev_min": _rev_min_lin,
            "rev_max": _rev_max_lin,
        })
        _diag = pd.DataFrame({"fwd": [0.0, max_val], "rev": [0.0, max_val]})

    ci_lower = (
        alt.Chart(_ci_band)
        .mark_line(color="steelblue", opacity=0.6, strokeDash=[3, 3], tooltip=None)
        .encode(alt.X("fwd:Q"), alt.Y("rev_min:Q"))
    )
    ci_upper = (
        alt.Chart(_ci_band)
        .mark_line(color="steelblue", opacity=0.6, strokeDash=[3, 3], tooltip=None)
        .encode(alt.X("fwd:Q"), alt.Y("rev_max:Q"))
    )
    diag_line = (
        alt.Chart(_diag)
        .mark_line(strokeDash=[6, 4], color="gray", opacity=0.7)
        .encode(
            alt.X("fwd:Q"),
            alt.Y("rev:Q"),
        )
    )
    _sb_n_pts = len(sample_df)
    _sb_title = (
        f"Strand Bias (log1p scale) — solid: perfect balance; dashed: 95% CI under Binomial(n, 0.5) — showing {_sb_n_pts:,} points"
        if _use_log1p else
        f"Strand Bias — solid: perfect balance; dashed: 95% CI under Binomial(n, 0.5) — showing {_sb_n_pts:,} points"
    )

    if _use_log1p:
        _max_linear = max(
            int(sample_df["fwd_alt_count"].max()),
            int(sample_df["rev_alt_count"].max()),
        )
        _tick_vals = [np.log1p(v) for v in _log1p_ticks_linear if v <= _max_linear * 1.1]
        _log1p_axis = alt.Axis(
            values=_tick_vals,
            labelExpr="format(exp(datum.value) - 1, 'd')",
        )
        _enc_x = alt.X("fwd_plot:Q", title=_x_title, axis=_log1p_axis,
                        scale=alt.Scale(domain=[0, max_val]))
        _enc_y = alt.Y("rev_plot:Q", title=_y_title, axis=_log1p_axis,
                        scale=alt.Scale(domain=[0, max_val]))
    else:
        _enc_x = alt.X("fwd_plot:Q", title=_x_title)
        _enc_y = alt.Y("rev_plot:Q", title=_y_title)

    _color_field, _color_title, _color_scale = {
        "Variant type":   ("variant_type:N",  "Variant type",   alt.Scale()),
        "Sample":         ("sample_id:N",      "Sample",         alt.Scale()),
        "Batch":          ("batch:N",          "Batch",          alt.Scale()),
        "On target":      ("on_target:N",      "On target",
                           alt.Scale(domain=[True, False], range=["#2ca02c", "#d62728"])),
        "Called variant": ("variant_called:N", "Called variant",
                           alt.Scale(domain=[True, False], range=["#2ca02c", "#d62728"])),
    }[_color_by]

    _sb_point_sel = alt.selection_point(
        name="sb_select",
        fields=["sample_id", "chrom", "pos", "ref_allele", "alt_allele"],
        on="click",
        toggle="event.shiftKey",
    )

    scatter = (
        alt.Chart(sample_df)
        .mark_point(size=40)
        .encode(
            _enc_x,
            _enc_y,
            alt.Color(_color_field, title=_color_title, scale=_color_scale),
            opacity=alt.condition(_sb_point_sel, alt.value(1.0), alt.value(0.35)),
            size=alt.condition(_sb_point_sel, alt.value(80), alt.value(30)),
            tooltip=(
                ["sample_id", "chrom", "pos", "ref_allele", "alt_allele",
                 "variant_type", "fwd_alt_count", "rev_alt_count", "vaf"]
                + (["on_target"] if _has_data("on_target") else [])
                + (["variant_called"] if _has_data("variant_called") else [])
                + (["gene"] if _genes_available else [])
            ),
        )
        .add_params(_sb_point_sel)
        .properties(title=_sb_title, height=350)
    )
    sb_event = st.altair_chart(
        (ci_lower + ci_upper + diag_line + scatter).resolve_scale(color="independent"),
        use_container_width=True,
        on_select="rerun",
    )

    # ── Drill-down for selected points ────────────────────────────────────────
    sb_pts = (sb_event.selection or {}).get("sb_select", [])
    if sb_pts:
        # Build a dataframe of selected records from the full (non-sampled) data
        _sb_or_clauses = " OR ".join(
            f"(sample_id = '{p['sample_id']}' AND chrom = '{p['chrom']}' "
            f"AND pos = {int(p['pos'])} AND ref_allele = '{p['ref_allele']}' "
            f"AND alt_allele = '{p['alt_allele']}')"
            for p in sb_pts
            if all(k in p for k in ["sample_id", "chrom", "pos", "ref_allele", "alt_allele"])
        )
        if _sb_or_clauses:
            _sb_sel_df = con.execute(f"""
                SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
                FROM {table_expr}
                WHERE {_sb_or_clauses}
                ORDER BY sample_id, chrom, pos
            """).df()

            n_pts = len(_sb_sel_df)
            n_smp = _sb_sel_df["sample_id"].nunique()
            st.caption(
                f"{n_pts} selected loci across {n_smp} sample(s) — "
                "shift-click to select multiple points"
            )
            st.dataframe(_sb_sel_df[_table_cols], use_container_width=True, hide_index=True)
            igv_buttons(
                [f"({_sb_or_clauses})"],
                _sb_sel_df,
                key=f"sb_sel_{'_'.join(str(int(p['pos'])) for p in sb_pts[:5] if 'pos' in p)}",
            )
    else:
        st.caption("Click a point to select it; shift-click to select multiple.")

with tab4:
    _ov_df = con.execute(f"""
        SELECT
            ROUND(
                overlap_alt_agree * 1.0
                / NULLIF(overlap_alt_agree + overlap_alt_disagree, 0),
                3
            ) AS agree_frac,
            COUNT(*) AS n
        FROM {table_expr}
        WHERE {where}
          AND overlap_depth > 0
          AND overlap_alt_agree IS NOT NULL
        GROUP BY agree_frac
        ORDER BY agree_frac
    """).df()

    if _ov_df.empty:
        st.info("No overlapping fragments in current selection.")
    else:
        chart = (
            alt.Chart(_ov_df)
            .mark_bar(opacity=0.8)
            .encode(
                alt.X("agree_frac:Q", bin=alt.Bin(maxbins=20), title="Overlap agreement fraction"),
                alt.Y("n:Q", title="Count"),
                tooltip=[
                    alt.Tooltip("agree_frac:Q", format=".3f", title="Agreement fraction"),
                    alt.Tooltip("n:Q", title="Count"),
                ],
            )
            .properties(title="Overlap Agreement Fraction", height=350)
        )
        st.altair_chart(chart, use_container_width=True)

with tab_cohort:
    if not path.endswith(".duckdb"):
        st.info("Cohort view is available when loading a merged DuckDB file (`geac merge` output).")
    else:
        st.subheader("Per-sample summary")
        st.caption(
            "Applies all active sidebar filters except the sample selection. "
            "Click a row to focus all other tabs on that sample."
        )

        # Build a where clause without the sample_sel condition so all
        # samples always appear in the cohort summary table.
        _cohort_conditions = [c for c in conditions if not c.startswith("sample_id IN")]
        _cohort_where = " AND ".join(_cohort_conditions) if _cohort_conditions else "true"

        _has_overlap = _has_data("overlap_alt_agree")
        _overlap_col = """
            ROUND(
                SUM(overlap_alt_agree) * 1.0
                    / NULLIF(SUM(overlap_alt_agree + overlap_alt_disagree), 0),
                4
            ) AS overlap_concordance,
        """ if _has_overlap else "NULL AS overlap_concordance,"

        _cohort_stats = con.execute(f"""
            SELECT
                sample_id,
                {'ANY_VALUE(batch) AS batch,' if _has_data('batch') else ''}
                COUNT(*) FILTER (WHERE variant_type = 'SNV')       AS n_snv,
                COUNT(*) FILTER (WHERE variant_type = 'insertion')  AS n_insertion,
                COUNT(*) FILTER (WHERE variant_type = 'deletion')   AS n_deletion,
                ROUND(AVG(total_depth), 1)                          AS mean_depth,
                ROUND(AVG(alt_count * 1.0 / total_depth), 6)       AS mean_vaf,
                ROUND(
                    AVG(fwd_alt_count * 1.0
                        / NULLIF(fwd_alt_count + rev_alt_count, 0)),
                    4
                ) AS strand_balance,
                {_overlap_col}
                COUNT(*) AS n_total
            FROM {table_expr}
            WHERE {_cohort_where}
            GROUP BY sample_id
            ORDER BY sample_id
        """).df()

        if _cohort_stats.empty:
            st.warning("No records match the current filters.")
        else:
            _cohort_event = st.dataframe(
                _cohort_stats,
                use_container_width=True,
                on_select="rerun",
                selection_mode="single-row",
                hide_index=True,
            )

            _cohort_sel = (_cohort_event.selection or {}).get("rows", [])
            if _cohort_sel:
                _focused_sample = _cohort_stats.iloc[_cohort_sel[0]]["sample_id"]
                st.caption(
                    f"Focused on **{_focused_sample}** — "
                    "click button below to filter all tabs to this sample."
                )
                if st.button(f"Filter all tabs to {_focused_sample}"):
                    st.session_state["sample_sel"] = [_focused_sample]
                    st.rerun()

            # ── Step 2: VAF distribution overlay ──────────────────────────────
            st.divider()
            st.subheader("VAF distribution by sample")

            _vaf_raw = con.execute(f"""
                SELECT sample_id,
                       ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
                FROM {table_expr}
                WHERE {_cohort_where} AND variant_type = 'SNV'
            """).df()

            if _vaf_raw.empty:
                st.info("No SNVs in current selection.")
            else:
                _vaf_yscale = st.radio(
                    "Y-axis scale", ["Linear", "Symlog"], horizontal=True,
                    key="vaf_yscale",
                )
                _vaf_chart = (
                    alt.Chart(_vaf_raw)
                    .transform_density(
                        density="vaf",
                        groupby=["sample_id"],
                        extent=[0, 1],
                        steps=200,
                    )
                    .mark_line(opacity=0.8)
                    .encode(
                        alt.X("value:Q", title="VAF"),
                        alt.Y("density:Q", title="Density",
                              scale=alt.Scale(type="symlog" if _vaf_yscale == "Symlog" else "linear")),
                        alt.Color("sample_id:N", title="Sample"),
                        tooltip=["sample_id:N",
                                 alt.Tooltip("value:Q", format=".3f", title="VAF")],
                    )
                    .properties(
                        title="SNV VAF Distribution (all samples)",
                        height=350,
                    )
                )
                st.altair_chart(_vaf_chart, use_container_width=True)

            # ── Step 3: Strand balance scatter ────────────────────────────────
            st.divider()
            st.subheader("Strand balance by sample")
            st.caption(
                "Each dot is one sample. x = mean strand balance (0.5 = perfect), "
                "y = mean VAF. Outliers in either axis may indicate a problematic sample."
            )

            _strand_stats = con.execute(f"""
                SELECT
                    sample_id,
                    ROUND(AVG(alt_count * 1.0 / total_depth), 6) AS mean_vaf,
                    ROUND(
                        AVG(fwd_alt_count * 1.0
                            / NULLIF(fwd_alt_count + rev_alt_count, 0)),
                        4
                    ) AS mean_strand_balance,
                    COUNT(*) AS n_loci
                FROM {table_expr}
                WHERE {_cohort_where} AND variant_type = 'SNV'
                GROUP BY sample_id
            """).df()

            if _strand_stats.empty:
                st.info("No SNVs in current selection.")
            else:
                _strand_chart = (
                    alt.Chart(_strand_stats)
                    .mark_circle(size=80, opacity=0.85)
                    .encode(
                        alt.X("mean_strand_balance:Q",
                              title="Mean strand balance (0.5 = perfect)",
                              scale=alt.Scale(domain=[0, 1])),
                        alt.Y("mean_vaf:Q",
                              title="Mean VAF",
                              scale=alt.Scale(zero=False)),
                        alt.Color("sample_id:N", title="Sample"),
                        alt.Size("n_loci:Q", title="SNV loci",
                                 scale=alt.Scale(range=[40, 300])),
                        tooltip=[
                            "sample_id:N",
                            alt.Tooltip("mean_strand_balance:Q", format=".4f",
                                        title="Mean strand balance"),
                            alt.Tooltip("mean_vaf:Q", format=".6f", title="Mean VAF"),
                            alt.Tooltip("n_loci:Q", format=",", title="SNV loci"),
                        ],
                    )
                    .properties(
                        title="Strand Balance vs Mean VAF (per sample)",
                        height=350,
                    )
                )
                # Dashed reference line at x=0.5 (perfect strand balance)
                _ref_line = (
                    alt.Chart(pd.DataFrame({"x": [0.5]}))
                    .mark_rule(strokeDash=[4, 4], color="gray", opacity=0.6)
                    .encode(alt.X("x:Q"))
                )
                st.altair_chart(
                    (_strand_chart + _ref_line).properties(height=350),
                    use_container_width=True,
                )

            # ── Step 4: Alt loci count vs mean base quality ───────────────────
            st.divider()
            st.subheader("Alt loci count vs mean base quality (per sample)")
            st.caption(
                "Each dot is one sample. x = number of distinct alt loci, "
                "y = mean base quality across all alt-supporting reads. "
                "Samples with many loci but low base quality may be artefact-driven."
            )
            _bq_loci_df = con.execute(f"""
                SELECT
                    ab.sample_id,
                    COUNT(DISTINCT CONCAT(ab.chrom, ':', ab.pos, ':', ab.alt_allele)) AS n_alt_loci,
                    ROUND(AVG(ar.base_qual), 2) AS mean_base_qual,
                    COUNT(ar.rowid) AS n_reads
                FROM (
                    SELECT sample_id, chrom, pos, alt_allele
                    FROM {table_expr}
                    WHERE {_cohort_where}
                ) ab
                INNER JOIN alt_reads ar
                    ON  ab.sample_id  = ar.sample_id
                    AND ab.chrom      = ar.chrom
                    AND ab.pos        = ar.pos
                    AND ab.alt_allele = ar.alt_allele
                WHERE ar.base_qual IS NOT NULL
                GROUP BY ab.sample_id
            """).df()

            if _bq_loci_df.empty:
                st.info("No base quality data available (alt_reads table may be absent).")
            else:
                _bq_loci_chart = (
                    alt.Chart(_bq_loci_df)
                    .mark_circle(size=80, opacity=0.85)
                    .encode(
                        alt.X("n_alt_loci:Q", title="Number of alt loci",
                              scale=alt.Scale(zero=True)),
                        alt.Y("mean_base_qual:Q", title="Mean base quality (Phred)",
                              scale=alt.Scale(zero=False)),
                        alt.Color("sample_id:N", title="Sample"),
                        alt.Size("n_reads:Q", title="Alt-supporting reads",
                                 scale=alt.Scale(range=[40, 300])),
                        tooltip=[
                            "sample_id:N",
                            alt.Tooltip("n_alt_loci:Q", format=",", title="Alt loci"),
                            alt.Tooltip("mean_base_qual:Q", format=".1f", title="Mean base qual"),
                            alt.Tooltip("n_reads:Q", format=",", title="Alt-supporting reads"),
                        ],
                    )
                    .properties(height=350, title="Alt loci count vs mean base quality (per sample)")
                )
                st.altair_chart(_bq_loci_chart, use_container_width=True)

            # ── Step 5: SNV count bar chart stacked by SBS6 substitution ──────
            st.subheader("SNV Count by Sample (SBS6 breakdown)")
            _sbs6_df = con.execute(f"""
                SELECT
                    sample_id,
                    CASE
                        WHEN ref_allele IN ('C','G') AND alt_allele IN ('A','T') THEN 'C>A'
                        WHEN ref_allele IN ('C','G') AND alt_allele IN ('G','C') THEN 'C>G'
                        WHEN ref_allele IN ('C','G') AND alt_allele IN ('T','A') THEN 'C>T'
                        WHEN ref_allele IN ('T','A') AND alt_allele IN ('A','T') THEN 'T>A'
                        WHEN ref_allele IN ('T','A') AND alt_allele IN ('C','G') THEN 'T>C'
                        WHEN ref_allele IN ('T','A') AND alt_allele IN ('G','C') THEN 'T>G'
                        ELSE 'other'
                    END AS substitution,
                    COUNT(*) AS n_snv
                FROM {table_expr}
                WHERE {_cohort_where} AND variant_type = 'SNV'
                GROUP BY sample_id, substitution
                ORDER BY sample_id, substitution
            """).df()

            if _sbs6_df.empty:
                st.info("No SNVs in current selection.")
            else:
                _sbs6_color_scale = alt.Scale(
                    domain=["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"],
                    range=["#1BBDEB", "#808080", "#E22926", "#CBCACB", "#A1CE63", "#EDB5C0"],
                )
                _sbs6_chart = (
                    alt.Chart(_sbs6_df)
                    .mark_bar()
                    .encode(
                        alt.X("sample_id:N", title="Sample", sort="-y",
                              axis=alt.Axis(labelAngle=-45, labelLimit=200)),
                        alt.Y("n_snv:Q", title="SNV count", stack="zero"),
                        alt.Color("substitution:N", title="Substitution",
                                  scale=_sbs6_color_scale),
                        alt.Tooltip(["sample_id:N", "substitution:N", "n_snv:Q"]),
                    )
                    .properties(height=350, title="SNV count per sample colored by SBS6 substitution type")
                )
                st.altair_chart(_sbs6_chart, use_container_width=True)

            # ── Step 5: SBS96 heatmap ──────────────────────────────────────────
            st.subheader("SBS96 Heatmap (samples × trinucleotide contexts)")
            if not _has_data("trinuc_context"):
                st.info("Trinucleotide context unavailable — run geac collect with a reference FASTA.")
            else:
                _hm_raw = con.execute(f"""
                    SELECT sample_id, trinuc_context, ref_allele, alt_allele, COUNT(*) AS n
                    FROM {table_expr}
                    WHERE {_cohort_where} AND variant_type = 'SNV'
                      AND trinuc_context IS NOT NULL AND length(trinuc_context) = 3
                    GROUP BY sample_id, trinuc_context, ref_allele, alt_allele
                """).df()

                if _hm_raw.empty:
                    st.info("No SNVs with trinucleotide context in current selection.")
                else:
                    _hm_raw["sbs_label"] = _hm_raw.apply(
                        lambda row: _sbs_label(row["trinuc_context"], row["ref_allele"], row["alt_allele"]),
                        axis=1,
                    )
                    _hm_raw = _hm_raw.dropna(subset=["sbs_label"])
                    _hm_agg = _hm_raw.groupby(["sample_id", "sbs_label"], as_index=False)["n"].sum()

                    # Normalize per sample so colours reflect profile, not total count
                    _totals = _hm_agg.groupby("sample_id")["n"].transform("sum")
                    _hm_agg["fraction"] = _hm_agg["n"] / _totals

                    # Fill missing context/sample combinations with zero
                    _all_combos = pd.MultiIndex.from_product(
                        [_hm_agg["sample_id"].unique(), _SBS_ORDER],
                        names=["sample_id", "sbs_label"],
                    )
                    _hm_full = (
                        _hm_agg.set_index(["sample_id", "sbs_label"])
                        .reindex(_all_combos, fill_value=0)
                        .reset_index()
                    )
                    _hm_full["mut_type"] = _hm_full["sbs_label"].str.extract(r'\[([A-Z]>[A-Z])\]')[0]

                    _hm_chart = (
                        alt.Chart(_hm_full)
                        .mark_rect()
                        .encode(
                            alt.X("sbs_label:N", sort=_SBS_ORDER, title=None,
                                  axis=alt.Axis(labels=False, ticks=False)),
                            alt.Y("sample_id:N", title="Sample"),
                            alt.Color("fraction:Q", title="Fraction of SNVs",
                                      scale=alt.Scale(scheme="blues")),
                            alt.Tooltip(["sample_id:N", "sbs_label:N", "n:Q", "fraction:Q"]),
                        )
                        .properties(
                            height=max(200, 20 * _hm_full["sample_id"].nunique()),
                            title="Normalised SBS96 profile per sample (fraction of SNVs)",
                        )
                    )

                    # Label strip — one label per mutation type at the group midpoint
                    _hm_label_df = pd.DataFrame([
                        {"sbs_label": [l for l in _SBS_ORDER if f"[{mt}]" in l][8], "mut_type": mt}
                        for mt in _SBS_MUT_TYPES
                    ])
                    _hm_label_strip = (
                        alt.Chart(_hm_label_df)
                        .mark_text(align="center", fontSize=11, fontWeight="bold")
                        .encode(
                            alt.X("sbs_label:N", sort=_SBS_ORDER,
                                  axis=alt.Axis(labels=False, ticks=False, title=None)),
                            alt.Y(value=15),
                            alt.Color("mut_type:N", legend=None,
                                      scale=alt.Scale(
                                          domain=list(_SBS_COLORS.keys()),
                                          range=list(_SBS_COLORS.values()),
                                      )),
                            alt.Text("mut_type:N"),
                        )
                        .properties(height=30)
                    )
                    st.altair_chart(
                        alt.vconcat(_hm_chart, _hm_label_strip, spacing=2)
                        .resolve_scale(x="shared"),
                        use_container_width=True,
                    )
                    st.caption(
                        "Color = fraction of that sample's SNVs falling in each trinucleotide context. "
                        "Contexts ordered by mutation type (C>A, C>G, C>T, T>A, T>C, T>G) then flanking bases."
                    )

with tab_reads:
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

        # Build a subquery that joins alt_reads to the current filtered locus set.
        # This ensures all reads plots respect the active locus-level filters.
        _r_join = f"""
            alt_reads ar
            INNER JOIN (
                SELECT DISTINCT sample_id, chrom, pos, alt_allele
                FROM {table_expr}
                WHERE {where}
            ) _filt
            ON  ar.sample_id  = _filt.sample_id
            AND ar.chrom      = _filt.chrom
            AND ar.pos        = _filt.pos
            AND ar.alt_allele = _filt.alt_allele
        """

        # ── Row 1: Family size histogram + Read position bias ─────────────────
        _r_col1, _r_col2 = st.columns(2)

        with _r_col1:
            st.subheader("Family size distribution")
            _fs_ctrl_col1, _fs_ctrl_col2 = st.columns(2)
            _fs_color_options = ["All samples (aggregate)", "Sample"]
            if _has_data("batch"):
                _fs_color_options.append("Batch")
            _fs_color_by = _fs_ctrl_col1.radio(
                "Color by", _fs_color_options,
                horizontal=True, key="fs_color_by",
            )
            _fs_y_mode = _fs_ctrl_col2.radio(
                "Y axis", ["Fraction", "Count"],
                horizontal=True, key="fs_y_mode",
            )
            _fs_by_sample = _fs_color_by == "Sample"
            _fs_by_batch  = _fs_color_by == "Batch"
            _fs_normalize = _fs_y_mode == "Fraction"

            _fs_group_col = (
                "ar.sample_id" if _fs_by_sample else
                "ab.batch"     if _fs_by_batch  else
                None
            )
            _fs_select = f"{_fs_group_col}, " if _fs_group_col else ""
            _fs_group  = f"{_fs_group_col}, " if _fs_group_col else ""

            # For batch grouping we need to join against alt_bases to get the batch column
            _fs_source = (
                f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, batch FROM {table_expr} WHERE batch IS NOT NULL) ab ON ar.sample_id = ab.sample_id"
                if _fs_by_batch else _r_join
            )

            _fs_df = con.execute(f"""
                SELECT {_fs_select}ar.family_size, COUNT(*) AS n_reads
                FROM {_fs_source}
                WHERE ar.family_size IS NOT NULL
                GROUP BY {_fs_group}ar.family_size
                ORDER BY {_fs_group}ar.family_size
            """).df()

            if _fs_df.empty:
                st.info("No family size data (fgbio cD tag absent in this dataset).")
            else:
                _fs_label_col = (
                    "sample_id" if _fs_by_sample else
                    "batch"     if _fs_by_batch  else
                    None
                )
                if _fs_normalize:
                    if _fs_label_col:
                        _fs_df["y_val"] = _fs_df.groupby(_fs_label_col)["n_reads"].transform(
                            lambda x: x / x.sum()
                        )
                    else:
                        _fs_df["y_val"] = _fs_df["n_reads"] / _fs_df["n_reads"].sum()
                    _y_field = "y_val:Q"
                    _y_title = "Fraction of alt-supporting reads"
                    _y_fmt = ".3f"
                else:
                    _fs_df["y_val"] = _fs_df["n_reads"]
                    _y_field = "y_val:Q"
                    _y_title = "Alt-supporting reads"
                    _y_fmt = "d"

                _fs_enc = dict(
                    x=alt.X("family_size:Q", title="Family size (cD tag)", bin=False),
                    y=alt.Y(_y_field, title=_y_title),
                    tooltip=[
                        *([f"{_fs_label_col}:N"] if _fs_label_col else []),
                        "family_size:Q",
                        alt.Tooltip(_y_field, format=_y_fmt, title=_y_title),
                    ],
                )
                if _fs_label_col:
                    _fs_enc["color"] = alt.Color(f"{_fs_label_col}:N", scale=alt.Scale(scheme="tableau10"))
                    _fs_chart = (
                        alt.Chart(_fs_df)
                        .mark_line(point=True, opacity=0.8)
                        .encode(**_fs_enc)
                        .properties(height=300)
                    )
                else:
                    _fs_chart = (
                        alt.Chart(_fs_df)
                        .mark_line(point=True, color="#4c78a8")
                        .encode(**_fs_enc)
                        .properties(height=300)
                    )
                st.altair_chart(_fs_chart, use_container_width=True)
                _fs_norm_note = (
                    "Fraction mode normalizes each batch independently."
                    if _fs_by_batch else
                    "Fraction mode normalizes each sample independently, allowing shape comparison across samples with different read counts."
                )
                st.caption(f"Artefacts are enriched in singletons (family_size = 1). {_fs_norm_note}")

        with _r_col2:
            st.subheader("Read position bias")
            _dfe_ctrl1, _dfe_ctrl2 = st.columns(2)
            _dfe_color_options = ["All samples (aggregate)", "Sample"]
            if _has_data("batch"):
                _dfe_color_options.append("Batch")
            _dfe_color_by = _dfe_ctrl1.radio(
                "Color by", _dfe_color_options,
                horizontal=True, key="dfe_color_by",
            )
            _dfe_y_mode = _dfe_ctrl2.radio(
                "Y axis", ["Fraction", "Count"],
                horizontal=True, key="dfe_y_mode",
            )
            _dfe_by_sample = _dfe_color_by == "Sample"
            _dfe_by_batch  = _dfe_color_by == "Batch"
            _dfe_normalize = _dfe_y_mode == "Fraction"
            _dfe_label_col = (
                "sample_id" if _dfe_by_sample else
                "batch"     if _dfe_by_batch  else
                None
            )
            _dfe_select = f"ar.{_dfe_label_col}, " if _dfe_label_col else ""
            _dfe_group  = f"ar.{_dfe_label_col}, " if _dfe_label_col else ""
            _dfe_source = (
                f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, batch FROM {table_expr} WHERE batch IS NOT NULL) ab ON ar.sample_id = ab.sample_id"
                if _dfe_by_batch else _r_join
            )
            # For batch, group col is on the joined table
            _dfe_group_expr = (
                "ab.batch, " if _dfe_by_batch else
                f"ar.{_dfe_label_col}, " if _dfe_label_col else ""
            )
            _dfe_select_expr = (
                "ab.batch AS batch, " if _dfe_by_batch else
                f"ar.{_dfe_label_col}, " if _dfe_label_col else ""
            )

            _dfe_df = con.execute(f"""
                SELECT {_dfe_select_expr}ar.dist_from_read_start, COUNT(*) AS n_reads
                FROM {_dfe_source}
                GROUP BY {_dfe_group_expr}ar.dist_from_read_start
                ORDER BY {_dfe_group_expr}ar.dist_from_read_start
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
                    x=alt.X("dist_from_read_start:Q", title="Cycle (distance from read start)", bin=False),
                    y=alt.Y(_dfe_y_field, title=_dfe_y_title),
                    tooltip=[
                        *([f"{_dfe_label_col}:N"] if _dfe_label_col else []),
                        alt.Tooltip("dist_from_read_start:Q", title="Cycle"),
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
                st.altair_chart(_dfe_chart, use_container_width=True)
                st.caption(
                    "A spike at high cycle numbers indicates alt-supporting reads clustered at read ends — "
                    "a red flag for alignment artefacts or damaged bases."
                )

        # ── Row 2: Base qual vs dist from read end scatter ────────────────────
        st.subheader("Mean base quality by cycle")
        _bq_color_options = ["All samples (aggregate)", "Sample"]
        if _has_data("batch"):
            _bq_color_options.append("Batch")
        _bq_color_by = st.radio(
            "Color by", _bq_color_options,
            horizontal=True, key="bq_color_by",
        )
        _bq_by_sample = _bq_color_by == "Sample"
        _bq_by_batch  = _bq_color_by == "Batch"
        _bq_label_col = (
            "sample_id" if _bq_by_sample else
            "batch"     if _bq_by_batch  else
            None
        )
        _bq_source = (
            f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, batch FROM {table_expr} WHERE batch IS NOT NULL) ab ON ar.sample_id = ab.sample_id"
            if _bq_by_batch else _r_join
        )
        _bq_select_expr = (
            "ab.batch AS batch, " if _bq_by_batch else
            f"ar.{_bq_label_col}, " if _bq_label_col else ""
        )
        _bq_group_expr = (
            "ab.batch, " if _bq_by_batch else
            f"ar.{_bq_label_col}, " if _bq_label_col else ""
        )

        _bq_df = con.execute(f"""
            SELECT
                {_bq_select_expr}ar.dist_from_read_start,
                ROUND(AVG(ar.base_qual), 2) AS mean_base_qual,
                COUNT(*) AS n_reads
            FROM {_bq_source}
            GROUP BY {_bq_group_expr}ar.dist_from_read_start
            ORDER BY {_bq_group_expr}ar.dist_from_read_start
        """).df()

        if _bq_df.empty:
            st.info("No data.")
        else:
            _bq_enc = dict(
                x=alt.X("dist_from_read_start:Q", title="Cycle (distance from read start)"),
                y=alt.Y("mean_base_qual:Q", title="Mean base quality (Phred)",
                        scale=alt.Scale(zero=False)),
                tooltip=[
                    *([f"{_bq_label_col}:N"] if _bq_label_col else []),
                    alt.Tooltip("dist_from_read_start:Q", title="Cycle"),
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
            st.altair_chart(_bq_chart, use_container_width=True)
            st.caption(
                "A drop in mean base quality at high cycle numbers (late in the read) "
                "indicates that alt-supporting reads at those positions may be artefacts."
            )

        # ── Row 3: Insert size distribution ───────────────────────────────────
        if con.execute("SELECT COUNT(*) FROM alt_reads WHERE insert_size IS NOT NULL LIMIT 1").fetchone()[0] > 0:
            st.subheader("Insert size distribution")
            _ins_color_options = ["All samples (aggregate)", "Sample"]
            if _has_data("batch"):
                _ins_color_options.append("Batch")
            _ins_color_by = st.radio(
                "Color by", _ins_color_options,
                horizontal=True, key="ins_color_by",
            )
            _ins_by_sample = _ins_color_by == "Sample"
            _ins_by_batch  = _ins_color_by == "Batch"
            _ins_label_col = (
                "sample_id" if _ins_by_sample else
                "batch"     if _ins_by_batch  else
                None
            )
            _ins_source = (
                f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, batch FROM {table_expr} WHERE batch IS NOT NULL) ab ON ar.sample_id = ab.sample_id"
                if _ins_by_batch else _r_join
            )
            _ins_select_expr = (
                "ab.batch AS batch, " if _ins_by_batch else
                f"ar.{_ins_label_col}, " if _ins_label_col else ""
            )
            _ins_group_expr = (
                "ab.batch, " if _ins_by_batch else
                f"ar.{_ins_label_col}, " if _ins_label_col else ""
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
                _ins_use_freq = _ins_y_mode == "Frequency"
                if _ins_use_freq:
                    _ins_group_col = _ins_label_col or "__all__"
                    if _ins_label_col:
                        _ins_totals = _ins_df.groupby(_ins_label_col)["n_reads"].transform("sum")
                    else:
                        _ins_totals = _ins_df["n_reads"].sum()
                    _ins_df["frequency"] = _ins_df["n_reads"] / _ins_totals
                _ins_y_field = "frequency" if _ins_use_freq else "n_reads"
                _ins_y_title = "Frequency" if _ins_use_freq else "Alt-supporting reads"
                _ins_enc = dict(
                    x=alt.X("insert_size:Q", title="Insert size (bp)"),
                    y=alt.Y(f"{_ins_y_field}:Q", title=_ins_y_title,
                            **({"axis": alt.Axis(format=".3f")} if _ins_use_freq else {})),
                    tooltip=[
                        *([f"{_ins_label_col}:N"] if _ins_label_col else []),
                        alt.Tooltip("insert_size:Q", title="Insert size (bp)"),
                        alt.Tooltip(f"{_ins_y_field}:Q", title="Frequency" if _ins_use_freq else "Reads",
                                    **({"format": ".4f"} if _ins_use_freq else {})),
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
                st.altair_chart(_ins_chart, use_container_width=True)
                st.caption(
                    "Insert size distribution of alt-supporting reads. "
                    "A shift toward shorter inserts can indicate adapter contamination or artefacts."
                )

        # ── Row 3b: Insert size by AF class (germline vs somatic) ─────────────
        if con.execute("SELECT COUNT(*) FROM alt_reads WHERE insert_size IS NOT NULL LIMIT 1").fetchone()[0] > 0:
            st.subheader("Insert size by allele frequency class")
            _af_ins_color_options = ["All samples (aggregate)", "Sample"]
            if _has_data("batch"):
                _af_ins_color_options.append("Batch")
            _af_ins_ctrl1, _af_ins_ctrl2 = st.columns(2)
            _af_ins_color_by = _af_ins_ctrl1.radio(
                "Color by", _af_ins_color_options,
                horizontal=True, key="af_ins_color_by",
            )
            _af_ins_y_mode = _af_ins_ctrl2.radio(
                "Y axis", ["Frequency", "Count"],
                horizontal=True, key="af_ins_y_mode",
            )
            _af_ins_by_sample = _af_ins_color_by == "Sample"
            _af_ins_by_batch  = _af_ins_color_by == "Batch"
            _af_ins_group_col = (
                "sample_id" if _af_ins_by_sample else
                "batch"     if _af_ins_by_batch  else
                None
            )
            _af_ins_extra_select = (
                "ar.sample_id, " if _af_ins_by_sample else
                "_locus.batch, " if _af_ins_by_batch  else
                ""
            )
            _af_ins_extra_group = (
                "ar.sample_id, " if _af_ins_by_sample else
                "_locus.batch, " if _af_ins_by_batch  else
                ""
            )
            _af_ins_locus_extra = ", batch" if _af_ins_by_batch else ""

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
                    SELECT DISTINCT sample_id, chrom, pos, alt_allele, alt_count, total_depth{_af_ins_locus_extra}
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
                _af_ins_use_freq = _af_ins_y_mode == "Frequency"

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

                if _af_ins_use_freq:
                    _norm_key = "series" if _af_ins_group_col else "af_class"
                    _totals = _af_ins_df.groupby(_norm_key)["n_reads"].transform("sum")
                    _af_ins_df["frequency"] = _af_ins_df["n_reads"] / _totals

                _af_ins_y = (
                    alt.Y("frequency:Q", title="Frequency", axis=alt.Axis(format=".3f"))
                    if _af_ins_use_freq else
                    alt.Y("n_reads:Q", title="Alt-supporting reads")
                )
                _af_ins_y_field = "frequency" if _af_ins_use_freq else "n_reads"
                _af_ins_tooltip = [
                    *([f"{_af_ins_group_col}:N"] if _af_ins_group_col else []),
                    alt.Tooltip("insert_size:Q", title="Insert size (bp)"),
                    alt.Tooltip("af_class:N", title="AF class"),
                    alt.Tooltip(_af_ins_y_field + ":Q",
                                title="Frequency" if _af_ins_use_freq else "Reads",
                                **({"format": ".4f"} if _af_ins_use_freq else {})),
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
                st.altair_chart(_af_ins_chart, use_container_width=True)
                st.caption(
                    "Insert size distributions split by allele frequency. "
                    "Frequency is normalized independently per series so lines are directly comparable "
                    "regardless of how many reads fall in each group. "
                    "A shift in one class toward very short or very long inserts suggests artefacts in that group."
                )

        # ── Row 4: Family size vs VAF per locus ───────────────────────────────
        st.subheader("Family size vs VAF (per locus)")
        _fsvaf_df = con.execute(f"""
            SELECT
                ab.sample_id,
                ab.chrom,
                ab.pos,
                ab.alt_allele,
                ROUND(ab.alt_count * 1.0 / ab.total_depth, 4) AS vaf,
                AVG(ar.family_size) AS mean_family_size,
                COUNT(*)            AS n_reads
            FROM {table_expr} ab
            INNER JOIN alt_reads ar
                ON  ab.sample_id  = ar.sample_id
                AND ab.chrom      = ar.chrom
                AND ab.pos        = ar.pos
                AND ab.alt_allele = ar.alt_allele
            WHERE {where} AND ar.family_size IS NOT NULL
            GROUP BY ab.sample_id, ab.chrom, ab.pos, ab.alt_allele, ab.alt_count, ab.total_depth
        """).df()

        if _fsvaf_df.empty:
            st.info("No data with family size available for the current selection.")
        else:
            _fsvaf_plot_df = (
                _fsvaf_df
                .sort_values(["sample_id", "chrom", "pos", "alt_allele"])
                .head(3000)
                .reset_index(drop=True)
            )
            _fsvaf_chart = (
                alt.Chart(_fsvaf_plot_df)
                .mark_point(filled=True, size=60)
                .encode(
                    alt.X("mean_family_size:Q", title="Mean family size of alt reads"),
                    alt.Y("vaf:Q", title="VAF", scale=alt.Scale(domain=[0, 1])),
                    alt.Color("sample_id:N", title="Sample",
                              scale=alt.Scale(scheme="tableau10")),
                    tooltip=[
                        alt.Tooltip("sample_id:N"),
                        alt.Tooltip("chrom:N"),
                        alt.Tooltip("pos:Q"),
                        alt.Tooltip("alt_allele:N"),
                        alt.Tooltip("vaf:Q", format=".4f"),
                        alt.Tooltip("mean_family_size:Q", format=".1f", title="Mean family size"),
                        alt.Tooltip("n_reads:Q", title="Alt reads"),
                    ],
                )
                .properties(height=350)
            )
            st.altair_chart(_fsvaf_chart, use_container_width=True)
            st.caption(
                "True low-VAF variants should have reasonable mean family sizes. "
                "Artefacts at low VAF tend to cluster at low family size."
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
            st.altair_chart(_mq_chart, use_container_width=True)
            st.caption(
                "Stacked by locus type (repetitive = homopolymer ≥ 5 or STR length ≥ 6). "
                "Low MAPQ at repetitive loci indicates multi-mapping artefacts."
            )

        # ── Row 5: Cohort artefact family size comparison (DuckDB only) ───────
        if path.endswith(".duckdb"):
            st.subheader("Cohort artefact vs rare variant: family size comparison")
            _n_samples_total = con.execute(
                f"SELECT COUNT(DISTINCT sample_id) FROM {table_expr} WHERE {where}"
            ).fetchone()[0]

            if _n_samples_total < 2:
                st.info("Need at least 2 samples for cohort artefact comparison.")
            else:
                _cohort_fs_df = con.execute(f"""
                    WITH locus_counts AS (
                        SELECT chrom, pos, alt_allele,
                               COUNT(DISTINCT sample_id) AS n_samples_with_alt
                        FROM {table_expr}
                        WHERE {where}
                        GROUP BY chrom, pos, alt_allele
                    ),
                    labeled AS (
                        SELECT
                            CASE
                                WHEN lc.n_samples_with_alt = 1 THEN 'Seen in 1 sample'
                                WHEN lc.n_samples_with_alt <= 3 THEN 'Seen in 2–3 samples'
                                ELSE 'Seen in 4+ samples'
                            END AS cohort_freq,
                            ar.family_size
                        FROM alt_reads ar
                        INNER JOIN (
                            SELECT DISTINCT sample_id, chrom, pos, alt_allele
                            FROM {table_expr}
                            WHERE {where}
                        ) _filt ON  ar.sample_id  = _filt.sample_id
                                 AND ar.chrom      = _filt.chrom
                                 AND ar.pos        = _filt.pos
                                 AND ar.alt_allele = _filt.alt_allele
                        INNER JOIN locus_counts lc
                            ON  ar.chrom      = lc.chrom
                            AND ar.pos        = lc.pos
                            AND ar.alt_allele = lc.alt_allele
                        WHERE ar.family_size IS NOT NULL
                    )
                    SELECT
                        cohort_freq,
                        PERCENTILE_CONT(0.0)  WITHIN GROUP (ORDER BY family_size) AS min_val,
                        PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY family_size) AS q1,
                        PERCENTILE_CONT(0.5)  WITHIN GROUP (ORDER BY family_size) AS median,
                        PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY family_size) AS q3,
                        PERCENTILE_CONT(1.0)  WITHIN GROUP (ORDER BY family_size) AS max_val,
                        COUNT(*) AS n_reads
                    FROM labeled
                    GROUP BY cohort_freq
                """).df()

                if _cohort_fs_df.empty:
                    st.info("No family size data for cohort comparison.")
                else:
                    _cf_order = ["Seen in 1 sample", "Seen in 2–3 samples", "Seen in 4+ samples"]
                    # IQR box
                    _cf_box = (
                        alt.Chart(_cohort_fs_df)
                        .mark_bar(size=40)
                        .encode(
                            alt.X("cohort_freq:N", title="Cohort frequency", sort=_cf_order),
                            alt.Y("q1:Q", title="Family size", scale=alt.Scale(zero=True)),
                            alt.Y2("q3:Q"),
                            alt.Color("cohort_freq:N", legend=None,
                                      scale=alt.Scale(scheme="tableau10")),
                            tooltip=[
                                alt.Tooltip("cohort_freq:N", title="Group"),
                                alt.Tooltip("min_val:Q", title="Min"),
                                alt.Tooltip("q1:Q", title="Q1"),
                                alt.Tooltip("median:Q", title="Median"),
                                alt.Tooltip("q3:Q", title="Q3"),
                                alt.Tooltip("max_val:Q", title="Max"),
                                alt.Tooltip("n_reads:Q", title="Reads"),
                            ],
                        )
                    )
                    # Min–max whiskers
                    _cf_whisker = (
                        alt.Chart(_cohort_fs_df)
                        .mark_rule()
                        .encode(
                            alt.X("cohort_freq:N", sort=_cf_order),
                            alt.Y("min_val:Q"),
                            alt.Y2("max_val:Q"),
                            alt.Color("cohort_freq:N", legend=None,
                                      scale=alt.Scale(scheme="tableau10")),
                        )
                    )
                    # Median tick
                    _cf_median = (
                        alt.Chart(_cohort_fs_df)
                        .mark_tick(color="white", thickness=2, size=40)
                        .encode(
                            alt.X("cohort_freq:N", sort=_cf_order),
                            alt.Y("median:Q"),
                        )
                    )
                    _cf_chart = (_cf_whisker + _cf_box + _cf_median).properties(height=300)
                    st.altair_chart(_cf_chart, use_container_width=True)
                    st.caption(
                        "Cohort artefacts (seen in many samples) tend to have lower family sizes "
                        "than rare variants, confirming they are sequencing noise rather than "
                        "recurrent true variants."
                    )

# ── Sig. Comparison (Experimental) tab ───────────────────────────────────────
with tab_sig_cmp:
    st.subheader("Called vs Uncalled: Trinucleotide Spectrum & COSMIC Signatures")
    st.caption(
        "**Experimental.** Compares the SBS96 trinucleotide spectrum and COSMIC signature "
        "exposures between loci where a variant was called (variant_called = true) and "
        "where it was not (variant_called = false). "
        "Requires trinuc_context and variant_called columns."
    )

    if not _has_data("trinuc_context"):
        st.info("trinuc_context column not present — re-run geac collect with a gene annotation or reference FASTA.")
    elif not _has_data("variant_called"):
        st.info("variant_called column not present — provide a variants TSV or VCF when running geac collect.")
    else:
        def _build_spectrum(called_val):
            """Return (96-element count array, total_count) for one variant_called group."""
            _raw = con.execute(f"""
                SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
                FROM {table_expr}
                WHERE {where}
                  AND variant_type = 'SNV'
                  AND trinuc_context IS NOT NULL
                  AND length(trinuc_context) = 3
                  AND variant_called IS {'TRUE' if called_val else 'FALSE'}
                GROUP BY trinuc_context, ref_allele, alt_allele
            """).df()

            if _raw.empty:
                return np.zeros(96, dtype=float), 0

            _raw["sbs_label"] = _raw.apply(
                lambda row: _sbs_label(row["trinuc_context"], row["ref_allele"], row["alt_allele"]),
                axis=1,
            )
            _raw = _raw.dropna(subset=["sbs_label"])
            _agg = _raw.groupby("sbs_label", as_index=False)["count"].sum()
            _obs = (
                pd.Series(0, index=_SBS_ORDER, dtype=float)
                .add(_agg.set_index("sbs_label")["count"], fill_value=0)
                .reindex(_SBS_ORDER)
                .values.astype(float)
            )
            return _obs, int(_obs.sum())

        _obs_called,   _n_called   = _build_spectrum(True)
        _obs_uncalled, _n_uncalled = _build_spectrum(False)

        if _n_called == 0 and _n_uncalled == 0:
            st.info("No SNVs with trinucleotide context found in either group.")
        else:
            # ── Trinucleotide spectrum comparison ─────────────────────────────
            st.subheader("Trinucleotide mutation spectrum")

            def _make_spec96_df(obs_arr, label, n_total):
                """Convert a 96-element count array to a plottable DataFrame with fractions."""
                df_s = pd.DataFrame({
                    "sbs_label": _SBS_ORDER,
                    "mut_type":  [lbl[2:5] for lbl in _SBS_ORDER],
                    "count":     obs_arr.astype(int),
                    "group":     label,
                })
                df_s["fraction"] = df_s["count"] / n_total if n_total > 0 else 0.0
                return df_s

            _spec_called   = _make_spec96_df(_obs_called,   f"Called (n={_n_called:,})",   _n_called)
            _spec_uncalled = _make_spec96_df(_obs_uncalled, f"Uncalled (n={_n_uncalled:,})", _n_uncalled)

            def _sbs96_chart(df_s, title):
                _sub_charts = []
                for _mt in _SBS_MUT_TYPES:
                    _sub = df_s[df_s["mut_type"] == _mt].copy()
                    _order = [lbl for lbl in _SBS_ORDER if f"[{_mt}]" in lbl]
                    _c = (
                        alt.Chart(_sub)
                        .mark_bar(color=_SBS_COLORS[_mt])
                        .encode(
                            alt.X("sbs_label:N", sort=_order, title=None,
                                  axis=alt.Axis(labelAngle=-90, labelFontSize=7)),
                            alt.Y("fraction:Q", title="Fraction",
                                  axis=alt.Axis(format=".0%")),
                            tooltip=[
                                alt.Tooltip("sbs_label:N", title="Context"),
                                alt.Tooltip("count:Q", title="Count"),
                                alt.Tooltip("fraction:Q", format=".2%", title="Fraction"),
                            ],
                        )
                        .properties(
                            title=alt.TitleParams(_mt, color=_SBS_COLORS[_mt], fontSize=11, fontWeight="bold"),
                            width=130, height=110,
                        )
                    )
                    _sub_charts.append(_c)
                return (
                    alt.concat(*_sub_charts, columns=3)
                    .resolve_scale(y="shared")
                    .properties(title=alt.TitleParams(title, fontSize=13))
                )

            # Option 1 — stacked (original, kept for reference)
            with st.expander("Option 1: Stacked (original)", expanded=False):
                st.altair_chart(
                    _sbs96_chart(_spec_called, f"Called variants (n={_n_called:,})"),
                    use_container_width=True,
                )
                st.altair_chart(
                    _sbs96_chart(_spec_uncalled, f"Uncalled loci (n={_n_uncalled:,})"),
                    use_container_width=True,
                )

            st.divider()

            # Option 2 — mirrored (butterfly) chart
            st.markdown("**Option 2: Mirrored spectrum** — Called bars point up, Uncalled point down")
            _m_df = pd.concat([
                _spec_called.assign(y=_spec_called["fraction"],  group="Called"),
                _spec_uncalled.assign(y=-_spec_uncalled["fraction"], group="Uncalled"),
            ])
            _mirror_sub = []
            for _mt in _SBS_MUT_TYPES:
                _sub = _m_df[_m_df["mut_type"] == _mt]
                _order = [lbl for lbl in _SBS_ORDER if f"[{_mt}]" in lbl]
                _c = (
                    alt.Chart(_sub)
                    .mark_bar()
                    .encode(
                        alt.X("sbs_label:N", sort=_order, title=None,
                              axis=alt.Axis(labelAngle=-90, labelFontSize=7)),
                        alt.Y("y:Q", axis=alt.Axis(format=".0%"),
                              title="← Uncalled | Called →"),
                        alt.Color("group:N",
                                  scale=alt.Scale(
                                      domain=["Called", "Uncalled"],
                                      range=["#4c78a8", "#e45756"],
                                  ),
                                  legend=alt.Legend(orient="bottom")),
                        tooltip=[
                            alt.Tooltip("sbs_label:N", title="Context"),
                            alt.Tooltip("group:N"),
                            alt.Tooltip("fraction:Q", format=".2%", title="Fraction"),
                            alt.Tooltip("count:Q", title="Count"),
                        ],
                    )
                    .properties(
                        title=alt.TitleParams(_mt, color=_SBS_COLORS[_mt], fontSize=11, fontWeight="bold"),
                        width=130, height=150,
                    )
                )
                _mirror_sub.append(_c)
            st.altair_chart(
                alt.concat(*_mirror_sub, columns=3)
                .resolve_scale(y="shared")
                .properties(title=alt.TitleParams(
                    f"Called (n={_n_called:,}) vs Uncalled (n={_n_uncalled:,}) — mirrored",
                    fontSize=13,
                )),
                use_container_width=True,
            )

            st.divider()

            # ── COSMIC decomposition ──────────────────────────────────────────
            st.subheader("COSMIC signature decomposition")
            _sc_cosmic_path = st.text_input(
                "COSMIC SBS matrix path",
                value=_cfg.get("cosmic", ""),
                placeholder="/path/to/COSMIC_v3.4_SBS_GRCh37.txt",
                key="sc_cosmic_path",
            )

            if not _sc_cosmic_path or not _sc_cosmic_path.strip():
                st.info("Enter a COSMIC SBS matrix path above to run the decomposition.")
            else:
                try:
                    _sc_cosmic_df = _load_cosmic(_sc_cosmic_path.strip())
                    _sc_cosmic_aligned = _sc_cosmic_df.reindex(_SBS_ORDER)
                    _sc_missing = _sc_cosmic_aligned.isna().any(axis=1).sum()

                    if _sc_missing > 0:
                        st.warning(
                            f"{_sc_missing} context(s) not found in COSMIC matrix — "
                            "check that the file uses the standard A[C>A]A format."
                        )
                    else:
                        _sc_W = _sc_cosmic_aligned.values.astype(float)  # 96 × N

                        def _fit(obs):
                            h, _ = nnls(_sc_W, obs)
                            total = h.sum()
                            h_norm = h / total if total > 0 else h
                            reconstructed = _sc_W @ h
                            cos_sim = (
                                float(np.dot(obs, reconstructed))
                                / (np.linalg.norm(obs) * np.linalg.norm(reconstructed) + 1e-12)
                            )
                            residual_pct = (
                                float(np.linalg.norm(obs - reconstructed))
                                / (float(obs.sum()) + 1e-12) * 100
                            )
                            return h_norm, cos_sim, residual_pct

                        _h_called,   _cos_called,   _res_called   = _fit(_obs_called)
                        _h_uncalled, _cos_uncalled, _res_uncalled = _fit(_obs_uncalled)

                        # ── Goodness-of-fit metrics ────────────────────────
                        _gof_c1, _gof_c2, _gof_c3, _gof_c4 = st.columns(4)
                        _gof_c1.metric("Called SNVs", f"{_n_called:,}")
                        _gof_c2.metric("Cosine sim (called)", f"{_cos_called:.4f}")
                        _gof_c3.metric("Uncalled SNVs", f"{_n_uncalled:,}")
                        _gof_c4.metric("Cosine sim (uncalled)", f"{_cos_uncalled:.4f}")

                        _sc_top_n = st.slider(
                            "Top signatures to display", 3, min(30, len(_sc_cosmic_df.columns)), 8,
                            key="sc_top_n",
                        )

                        # Build combined DataFrame for grouped bar chart
                        _sig_names = _sc_cosmic_aligned.columns.tolist()
                        _cmp_df = pd.DataFrame({
                            "signature": _sig_names * 2,
                            "group": ["Called"] * len(_sig_names) + ["Uncalled"] * len(_sig_names),
                            "exposure": list(_h_called) + list(_h_uncalled),
                        })
                        _cmp_df["etiology"] = _cmp_df["signature"].map(
                            lambda s: _SBS_ETIOLOGY.get(s, "")
                        )

                        # Rank signatures by max exposure across both groups
                        _max_exp = (
                            _cmp_df.groupby("signature")["exposure"].max()
                            .sort_values(ascending=False)
                        )
                        _top_sigs = _max_exp.head(_sc_top_n).index.tolist()
                        _top_cmp_df = _cmp_df[_cmp_df["signature"].isin(_top_sigs)].copy()

                        _bars = (
                            alt.Chart(_top_cmp_df)
                            .mark_bar()
                            .encode(
                                alt.X("signature:N",
                                      sort=_top_sigs,
                                      title="Signature"),
                                alt.Y("exposure:Q",
                                      title="Exposure (proportion)",
                                      axis=alt.Axis(format=".0%")),
                                alt.Color("group:N",
                                          title="Variant group",
                                          scale=alt.Scale(
                                              domain=["Called", "Uncalled"],
                                              range=["#4c78a8", "#e45756"],
                                          )),
                                alt.XOffset("group:N"),
                                tooltip=[
                                    alt.Tooltip("signature:N"),
                                    alt.Tooltip("group:N"),
                                    alt.Tooltip("exposure:Q", format=".2%", title="Exposure"),
                                    alt.Tooltip("etiology:N", title="Etiology"),
                                ],
                            )
                        )
                        # Vertical separator lines between signature groups.
                        # Using the *next* signature at bandPosition=0 (left edge)
                        # places each rule exactly between the Uncalled bar of
                        # signature i and the Called bar of signature i+1.
                        _dividers = (
                            alt.Chart(pd.DataFrame({"signature": _top_sigs[1:]}))
                            .mark_rule(color="#888", strokeWidth=1, opacity=0.5)
                            .encode(
                                alt.X("signature:N", sort=_top_sigs, bandPosition=0),
                            )
                        )
                        _cmp_chart = (
                            alt.layer(_bars, _dividers)
                            .properties(
                                title=f"Top {_sc_top_n} COSMIC SBS Signatures — Called vs Uncalled",
                                height=350,
                            )
                        )
                        st.altair_chart(_cmp_chart, use_container_width=True)
                        st.caption(
                            "Blue = called variants; red = uncalled. "
                            "If called variants are enriched in known cancer signatures (e.g. SBS1, SBS5) "
                            "while uncalled are dominated by artefact signatures (e.g. SBS58), "
                            "this supports the quality of the variant calling."
                        )

                        # ── Detailed table ─────────────────────────────────
                        with st.expander("Full signature table"):
                            _pivot = (
                                _cmp_df[_cmp_df["exposure"] > 0]
                                .pivot(index="signature", columns="group", values="exposure")
                                .fillna(0)
                                .reset_index()
                            )
                            for col in ["Called", "Uncalled"]:
                                if col in _pivot.columns:
                                    _pivot[col] = _pivot[col].map("{:.2%}".format)
                            _pivot["etiology"] = _pivot["signature"].map(
                                lambda s: _SBS_ETIOLOGY.get(s, "")
                            )
                            st.dataframe(_pivot, use_container_width=True, hide_index=True)

                except Exception as exc:
                    st.error(f"Failed to load COSMIC matrix: {exc}")
