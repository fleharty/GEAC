import io
import json
import warnings
import zipfile
from pathlib import Path
import numpy as np
import streamlit as st
import duckdb
import altair as alt
import pandas as pd
from scipy.optimize import nnls
from signature_nmf import (
    build_signature_download_zip,
    build_signature_exposure_download_table,
    build_signature_download_table,
    compare_signatures_to_cosmic,
    fit_cosmic_augmented_nmf,
    fit_sbs_nmf,
)
from pipeline_compare_helpers import (
    build_unique_pipeline_characterization_df,
    summarize_unique_pipeline_groups,
)
from read_context_helpers import add_read_context_fraction_metrics

# Altair emits a spurious "Automatically deduplicated selection parameter" UserWarning
# when a shared cross-panel selection param is present in multiple sub-charts.  The
# deduplication it performs is intentionally correct behaviour; suppress the noise.
warnings.filterwarnings(
    "ignore",
    message="Automatically deduplicated",
    category=UserWarning,
)

from igv_helpers import (
    query_distinct_samples,
    per_read_warning_note,
    resolve_index_uri,
)
import geac_config
from explorer import GEAC_VERSION, MAIN_FILTER_KEYS, MAIN_FILTER_STATE, DataSource
from explorer.main_table import (
    render_position_drilldown,
    render_records_table,
    render_summary_metrics,
)
from explorer.tabs import TAB_MODULES

_IS_MIN, _IS_MAX = 20, 500  # insert size slider bounds


def _sql_str(value: str) -> str:
    """Escape a string value for safe interpolation into a SQL literal."""
    return value.replace("'", "''")


@st.cache_data
def _compute_recurrence_loci(
    path: str, sr_lo: int, sr_hi: int, scope_where: str = "TRUE"
) -> pd.DataFrame:
    """Return (chrom, pos, ref_allele, alt_allele) tuples seen in sr_lo..sr_hi samples.

    ``scope_where`` restricts which rows participate in the recurrence count
    (e.g. sample/batch/label/on-target filters) so the slider reflects only
    the currently visible cohort subset.

    Results are cached by (path, sr_lo, sr_hi, scope_where) so the expensive
    GROUP BY only re-executes when any of these values change.
    """
    if path.endswith(".duckdb"):
        _c = duckdb.connect(path, read_only=True)
        tbl = "alt_bases"
    else:
        _c = duckdb.connect()
        tbl = f"read_parquet('{path}', union_by_name=true)"
    result = _c.execute(f"""
        SELECT chrom, pos, ref_allele, alt_allele
        FROM {tbl}
        WHERE {scope_where}
        GROUP BY chrom, pos, ref_allele, alt_allele
        HAVING COUNT(DISTINCT sample_id) BETWEEN {sr_lo} AND {sr_hi}
    """).df()
    _c.close()
    return result

st.set_page_config(page_title="GEAC Explorer", layout="wide")

_LOGO = Path(__file__).parent / "geac-logo.svg"
_LOGO_COMPACT = Path(__file__).parent / "geac-logo-compact.svg"
if _LOGO.exists():
    st.image(str(_LOGO), width="stretch")
else:
    st.title("GEAC Explorer")
    st.markdown(
        "**Genomic Evidence Atlas of Cohorts** — inspect alt base metrics from "
        "per-sample Parquet files or a merged cohort DuckDB database."
    )

st.sidebar.caption(f"geac v{GEAC_VERSION}")

# ── Project config (geac.toml or --config flag) ───────────────────────────────
_cfg = geac_config.load()

# ── File input ────────────────────────────────────────────────────────────────
path = st.text_input(
    "Data file path",
    value=_cfg.get("data", ""),
    placeholder="/path/to/sample.parquet  or  cohort.duckdb",
)

if not path or not path.strip():
    st.info("Enter a Parquet or DuckDB file path above to begin.", icon="🗂️")
    st.stop()

path = path.strip()

try:
    data_source = DataSource.open_alt_bases(path)
except Exception as e:
    st.error(f"Could not open file: {e}")
    st.stop()

con = data_source.con
table_expr = data_source.table_expr

# Detect whether optional tables are available (only possible in DuckDB mode).
_has_alt_reads = data_source.has_optional_table("alt_reads")
_has_normal_evidence = data_source.has_optional_table("normal_evidence")
_has_pon_evidence = data_source.has_optional_table("pon_evidence")
_alt_reads_cols = (
    set(con.execute("SELECT * FROM alt_reads LIMIT 0").df().columns)
    if _has_alt_reads
    else set()
)


def _has_alt_reads_cols(*cols: str) -> bool:
    return all(col in _alt_reads_cols for col in cols)

# ── Version check ─────────────────────────────────────────────────────────────
if data_source.is_duckdb:
    _db_version = data_source.db_version
    _db_created = data_source.db_created
    if _db_version is None:
        st.warning(
            f"This database was created with a version of geac older than v0.3.14 "
            f"(no `geac_metadata` table found). The Explorer expects v{GEAC_VERSION}. "
            "Some columns or features may be missing.",
            icon="⚠️",
        )
    elif _db_version != GEAC_VERSION:
        st.warning(
            f"Version mismatch: database was created with geac v{_db_version}, "
            f"but this Explorer expects v{GEAC_VERSION}. "
            "Results may be incomplete or columns may differ.",
            icon="⚠️",
        )

if _db_created is not None:
    _created_str = str(_db_created)[:10]  # YYYY-MM-DD
    st.sidebar.caption(f"DB created {_created_str}")

_missing_required_cols = data_source.required_columns_missing()
if _missing_required_cols:
    st.warning(
        "This dataset is missing required `alt_bases` columns expected by the current Explorer: "
        + ", ".join(sorted(_missing_required_cols)),
        icon="⚠️",
    )

# ── Summary stats ─────────────────────────────────────────────────────────────
stats = data_source.summary_stats()

n_annotated = int(stats["n_annotated"][0])
n_called    = int(stats["n_called"][0])
pct_called  = f"{100 * n_called / n_annotated:.1f}%" if n_annotated > 0 else "N/A"

# ── Filters (sidebar) ─────────────────────────────────────────────────────────
_FILTER_KEYS = list(MAIN_FILTER_KEYS)

_sidebar_logo = _LOGO_COMPACT if _LOGO_COMPACT.exists() else _LOGO
if _sidebar_logo.exists():
    st.sidebar.image(str(_sidebar_logo), width="stretch")

chroms = data_source.distinct_values("chrom")
samples = data_source.distinct_values("sample_id")

_hdr_col, _btn_col = st.sidebar.columns([2, 1])
_hdr_col.header("🔧 Filters")
if _btn_col.button("Clear all", help="Reset all filters to defaults"):
    _r_fs_max  = st.session_state.get("_cached_fs_max", 0)
    _r_cycle_max = st.session_state.get("_cached_cycle_max", 300)
    _r_mq_max    = st.session_state.get("_cached_mq_max", 60)
    MAIN_FILTER_STATE.reset(
        st.session_state,
        overrides={
            "sample_recurrence": (1, len(samples)),
            "family_size_range": (0, _r_fs_max),
            "cycle_range": (1, _r_cycle_max),
            "map_qual_range": (0, _r_mq_max),
            "insert_size_range": (_IS_MIN, _IS_MAX),
        },
    )
    # Clear drill-down state — the underlying data is changing.
    st.session_state.pop("_drill_locus", None)
    st.rerun()

# Initialize all managed filter defaults into session state on first load.
# This avoids the Streamlit warning that fires when a widget has both key= and
# value= and the key is already in session_state (e.g. after "Clear all").
for _fk, _fv in MAIN_FILTER_STATE.defaults.items():
    if _fk not in st.session_state:
        st.session_state[_fk] = _fv

chrom_sel = st.sidebar.selectbox("Chromosome", ["All"] + chroms, key="chrom_sel")
if "sample_sel" not in st.session_state:
    st.session_state["sample_sel"] = []
sample_sel = st.sidebar.multiselect("Samples (blank = all)", samples, key="sample_sel")

_n_samples_total = len(samples)
if _n_samples_total > 1:
    # Clamp any stored session state to the valid range for the current dataset
    # (guards against loading a different dataset with fewer samples).
    if "sample_recurrence" not in st.session_state:
        st.session_state["sample_recurrence"] = (1, _n_samples_total)
    else:
        _sr = st.session_state["sample_recurrence"]
        st.session_state["sample_recurrence"] = (
            max(1, min(_sr[0], _n_samples_total)),
            max(1, min(_sr[1], _n_samples_total)),
        )
    sample_recurrence = st.sidebar.slider(
        "Sample recurrence (# samples with locus)",
        min_value=1,
        max_value=_n_samples_total,
        step=1,
        key="sample_recurrence",
        help="Filter loci by how many samples carry that alt allele. "
             "Set min > 1 to find recurrent sites; set max = 1 to find sample-unique sites.",
    )
else:
    sample_recurrence = (1, 1)

_schema_cols = set(data_source.schema_cols)

# Probe all optional columns in a single query rather than N+1 individual queries.
_optional_cols = data_source.has_non_null_batch([
    "batch", "label1", "label2", "label3", "gene",
    "variant_filter", "gnomad_af", "homopolymer_len",
    "trinuc_context", "on_target", "variant_called",
    "overlap_alt_agree",
])

def _has_data(col: str) -> bool:
    """True iff col exists in the schema AND has at least one non-null value."""
    return _optional_cols.get(col, data_source.has_non_null(col))

if _has_data("batch"):
    _batches = con.execute(f"SELECT DISTINCT batch FROM {table_expr} WHERE batch IS NOT NULL ORDER BY batch").df()["batch"].tolist()
    if "batch_sel" not in st.session_state:
        st.session_state["batch_sel"] = []
    batch_sel = st.sidebar.multiselect("Batch (blank = all)", _batches, key="batch_sel")
else:
    batch_sel = []

_pipelines = con.execute(f"SELECT DISTINCT pipeline FROM {table_expr} WHERE pipeline IS NOT NULL ORDER BY pipeline").df()["pipeline"].tolist()
if "pipeline_sel" not in st.session_state:
    st.session_state["pipeline_sel"] = []
pipeline_sel = st.sidebar.multiselect("Pipeline (blank = all)", _pipelines, key="pipeline_sel")

if _has_data("label1"):
    _label1_vals = con.execute(f"SELECT DISTINCT label1 FROM {table_expr} WHERE label1 IS NOT NULL ORDER BY label1").df()["label1"].tolist()
    if "label1_sel" not in st.session_state:
        st.session_state["label1_sel"] = []
    label1_sel = st.sidebar.multiselect("Label 1 (blank = all)", _label1_vals, key="label1_sel")
else:
    label1_sel = []

if _has_data("label2"):
    _label2_vals = con.execute(f"SELECT DISTINCT label2 FROM {table_expr} WHERE label2 IS NOT NULL ORDER BY label2").df()["label2"].tolist()
    if "label2_sel" not in st.session_state:
        st.session_state["label2_sel"] = []
    label2_sel = st.sidebar.multiselect("Label 2 (blank = all)", _label2_vals, key="label2_sel")
else:
    label2_sel = []

if _has_data("label3"):
    _label3_vals = con.execute(f"SELECT DISTINCT label3 FROM {table_expr} WHERE label3 IS NOT NULL ORDER BY label3").df()["label3"].tolist()
    if "label3_sel" not in st.session_state:
        st.session_state["label3_sel"] = []
    label3_sel = st.sidebar.multiselect("Label 3 (blank = all)", _label3_vals, key="label3_sel")
else:
    label3_sel = []

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
vaf_range = st.sidebar.slider("VAF range", 0.0, 1.0, step=0.01, key="vaf_range")
min_alt = st.sidebar.number_input("Min alt count", min_value=1, max_value=10000, step=1, key="min_alt")
max_alt = st.sidebar.number_input("Max alt count (0 = no maximum)", min_value=0, max_value=10000, step=1, key="max_alt")
min_fwd_alt = st.sidebar.number_input("Min fwd alt count (0 = no minimum)", min_value=0, max_value=10000, step=1, key="min_fwd_alt")
min_rev_alt = st.sidebar.number_input("Min rev alt count (0 = no minimum)", min_value=0, max_value=10000, step=1, key="min_rev_alt")
min_overlap_agree    = st.sidebar.number_input("Min overlap alt agree (0 = no minimum)",    min_value=0, max_value=10000, step=1, key="min_overlap_agree")
min_overlap_disagree = st.sidebar.number_input("Min overlap alt disagree (0 = no minimum)", min_value=0, max_value=10000, step=1, key="min_overlap_disagree")
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

_GNOMAD_AF_STEPS = ["0", "1e-6", "1e-5", "1e-4", "1e-3", "0.01", "0.1", "1.0"]
if _has_data("gnomad_af"):
    # Ensure session state holds a valid 2-tuple before the widget renders.
    # A plain string (e.g. "0") causes select_slider to render as single-value
    # and return a string, which breaks the [0]/[1] indexing below.
    # Streamlit can store the range as a list after user interaction, so
    # accept both list and tuple as valid; only reset on a bare string/None.
    # select_slider requires value= to be a tuple to render in range mode;
    # key= alone does not communicate range mode to Streamlit. Manage
    # session state manually so we can always pass an explicit tuple value.
    _af_ss = st.session_state.get("gnomad_af_range", ("0", "1.0"))
    if not (isinstance(_af_ss, (tuple, list)) and len(_af_ss) == 2):
        _af_ss = ("0", "1.0")
    _af_ss = tuple(_af_ss)
    gnomad_af_range = st.sidebar.select_slider(
        "gnomAD AF (log scale)",
        options=_GNOMAD_AF_STEPS,
        value=_af_ss,
        help="Filter by gnomAD allele frequency. Steps are logarithmic.",
    )
    st.session_state["gnomad_af_range"] = gnomad_af_range
    gnomad_include_null = st.sidebar.checkbox(
        "Include sites absent from gnomAD",
        key="gnomad_include_null",
    )
else:
    gnomad_af_range    = ("0", "1.0")
    gnomad_include_null = True

_repeat_cols_present = _has_data("homopolymer_len")
if _repeat_cols_present:
    homopolymer_range = st.sidebar.slider("Homopolymer length range", 0, 20, step=1, key="homopolymer_range")
    str_len_range     = st.sidebar.slider("STR length range",         0, 50, step=1, key="str_len_range")
else:
    homopolymer_range = (0, 20)
    str_len_range     = (0, 50)
    st.sidebar.caption("Repeat filters unavailable — run geac collect with a newer build to enable.")
min_depth = st.sidebar.number_input("Min depth (0 = no minimum)", min_value=0, step=1, key="min_depth")
max_depth = st.sidebar.number_input("Max depth (0 = no maximum)", min_value=0, step=1, key="max_depth")


# ── Per-read filters (only when alt_reads table is present) ───────────────────
_reads_conditions = []
if _has_alt_reads:
    if "_cached_fs_max" not in st.session_state:
        _reads_maxes = con.execute("""
            SELECT
                MAX(family_size),
                COALESCE(MAX(cycle), 300),
                COALESCE(MAX(map_qual), 60),
                COUNT(insert_size) > 0
            FROM alt_reads
        """).fetchone()
        st.session_state["_cached_fs_max"]     = _reads_maxes[0]   # None if all NULL
        st.session_state["_cached_cycle_max"]   = int(_reads_maxes[1])
        st.session_state["_cached_mq_max"]     = int(_reads_maxes[2])
        st.session_state["_cached_is_has_data"] = bool(_reads_maxes[3])
    _fs_max_raw  = st.session_state["_cached_fs_max"]
    _cycle_max   = st.session_state["_cached_cycle_max"]
    _mq_max      = st.session_state["_cached_mq_max"]
    _is_has_data = st.session_state["_cached_is_has_data"]
    _fs_has_data = _fs_max_raw is not None
    _fs_max = int(_fs_max_raw) if _fs_has_data else 0

    st.sidebar.divider()
    st.sidebar.subheader("Per-read filters")
    recompute_vaf = st.sidebar.checkbox(
        "Recompute alt count from filtered reads",
        key="recompute_vaf",
        help="When checked, alt_count is re-aggregated from the reads table using only "
             "reads that pass the filters below. VAF denominator (total_depth) is unchanged — "
             "so displayed VAF is a lower bound, not the true filtered VAF.\n\n"
             "When unchecked (default), filters control which loci are visible: loci with "
             "no passing reads are hidden, but alt_count and VAF reflect all reads.",
    )

    if _fs_has_data:
        family_size_range = st.sidebar.slider(
            "Family size range",
            min_value=0, max_value=_fs_max, value=(0, _fs_max), step=1,
            key="family_size_range",
            help="family_size = cD tag (total raw read count per molecule).",
        )
    else:
        family_size_range = (0, 0)
        st.sidebar.caption("Family size unavailable — BAM has no fgbio cD tag.")

    cycle_range = st.sidebar.slider(
        "Cycle number",
        min_value=1, max_value=_cycle_max, value=(1, _cycle_max), step=1,
        key="cycle_range",
        help="Filter alt-supporting reads by sequencing cycle (1-based position within the read). "
             "Lower the upper bound to exclude variants clustered at read ends (a common artefact).",
    )

    map_qual_range = st.sidebar.slider(
        "Mapping quality range",
        min_value=0, max_value=_mq_max, value=(0, _mq_max), step=1,
        key="map_qual_range",
        help="Filter alt-supporting reads by mapping quality (MAPQ).",
    )

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
                "(including unpaired reads with no insert size) are accepted."
            )
        else:
            st.sidebar.caption(
                f"Insert size: keeping only reads with insert size between {_is_lo} and {_is_hi} bp. "
                "Unpaired reads (no insert size) are excluded."
            )
    else:
        insert_size_range = (_IS_MIN, _IS_MAX)

    read_strand_sel = st.sidebar.radio(
        "Read",
        ["All", "R1 only", "R2 only"],
        horizontal=True,
        key="read_strand_sel",
        help="Filter to R1 reads (BAM flag 0x40 set), R2 reads, or all reads.",
    )

    _fs_lo, _fs_hi = family_size_range
    _cycle_lo, _cycle_hi = cycle_range
    _mq_lo, _mq_hi = map_qual_range
    _is_lo, _is_hi = insert_size_range

    if _fs_has_data and (_fs_lo > 0 or _fs_hi < _fs_max):
        _reads_conditions.append(f"family_size BETWEEN {_fs_lo} AND {_fs_hi}")

    if _cycle_lo > 1 or _cycle_hi < _cycle_max:
        _reads_conditions.append(f"cycle BETWEEN {_cycle_lo} AND {_cycle_hi}")

    if _mq_lo > 0 or _mq_hi < _mq_max:
        _reads_conditions.append(f"map_qual BETWEEN {_mq_lo} AND {_mq_hi}")

    if _is_has_data and (_is_lo > _IS_MIN or _is_hi < _IS_MAX):
        # Unpaired reads (insert_size IS NULL) are excluded when a range is active
        _reads_conditions.append(f"insert_size BETWEEN {_is_lo} AND {_is_hi}")

    if read_strand_sel == "R1 only":
        _reads_conditions.append("is_read1 = true")
    elif read_strand_sel == "R2 only":
        _reads_conditions.append("is_read1 = false")

    if _reads_conditions:
        if recompute_vaf:
            st.sidebar.caption(
                "Alt count re-aggregated from filtered reads. "
                "VAF = filtered alt count / total depth (all reads) — a lower bound, not true filtered VAF."
            )
        else:
            st.sidebar.caption(
                "Loci with no alt reads passing these filters are hidden. "
                "Alt count and VAF reflect all reads."
            )
else:
    family_size_range = (0, 0)
    cycle_range = (1, 1)
    map_qual_range = (0, 0)
    insert_size_range = (0, 0)
    read_strand_sel = "All"
    _fs_lo = _fs_hi = _cycle_lo = _cycle_hi = _mq_lo = _mq_hi = _is_lo = _is_hi = 0

_base_table_expr = table_expr  # pre-reads-filter; used for sample-recurrence counts

# When per-read filters are active, redefine table_expr as a JOIN subquery.
# Two modes controlled by the "Recompute alt count from filtered reads" checkbox:
#
# recompute_vaf=False (default, locus-inclusion mode):
#   INNER JOIN keeps only loci that have ≥1 passing read. alt_count and VAF
#   come from alt_bases unchanged — VAF is always meaningful and consistent.
#
# recompute_vaf=True (re-aggregation mode):
#   alt_count is replaced with the filtered read count. total_depth is unchanged
#   so displayed VAF = filtered_alt_count / total_depth (a lower bound).
_reads_active = bool(_reads_conditions)
if _reads_active:
    _reads_where = " AND ".join(_reads_conditions)
    if recompute_vaf:
        # LEFT JOIN so loci with no rows in alt_reads (e.g. indels) keep their
        # original alt_count rather than being dropped by an INNER JOIN.
        # A single grouped subquery computes both the filtered count and a
        # presence flag so we can distinguish two NULL cases:
        #   ar_agg.has_reads IS NULL  → locus has no alt_reads rows at all (indels)
        #                               → preserve original alt_count
        #   ar_agg.has_reads IS TRUE  → locus has reads but none passed the filter
        #                               → filtered_alt_count = 0
        table_expr = f"""(
            SELECT
                ab.* EXCLUDE (alt_count),
                CASE
                    WHEN ar_agg.has_reads IS NULL THEN ab.alt_count
                    ELSE COALESCE(ar_agg.filtered_alt_count, 0)
                END AS alt_count,
                ROUND(ab.alt_count * 1.0 / ab.total_depth, 4) AS original_vaf
            FROM alt_bases ab
            LEFT JOIN (
                SELECT sample_id, chrom, pos, alt_allele,
                       COUNT(*) FILTER (WHERE {_reads_where}) AS filtered_alt_count,
                       TRUE AS has_reads
                FROM alt_reads
                GROUP BY sample_id, chrom, pos, alt_allele
            ) ar_agg ON ab.sample_id = ar_agg.sample_id
                     AND ab.chrom = ar_agg.chrom
                     AND ab.pos = ar_agg.pos
                     AND ab.alt_allele = ar_agg.alt_allele
        )"""
    else:
        # Locus-inclusion mode: exclude loci that HAVE reads in alt_reads but
        # none pass the filter. Loci with NO reads in alt_reads (e.g. indels,
        # which are not yet written to alt_reads by geac collect) pass through.
        # Single LEFT JOIN with BOOL_OR: one pass over alt_reads, no correlated subqueries.
        # ar.sample_id IS NULL  → locus has no alt_reads rows at all (e.g. indels): pass through.
        # ar.any_passing        → locus has at least one read satisfying the filter: include.
        table_expr = f"""(
            SELECT ab.*
            FROM alt_bases ab
            LEFT JOIN (
                SELECT sample_id, chrom, pos, alt_allele,
                       BOOL_OR({_reads_where}) AS any_passing
                FROM alt_reads
                GROUP BY sample_id, chrom, pos, alt_allele
            ) ar ON ab.sample_id = ar.sample_id
                AND ab.chrom     = ar.chrom
                AND ab.pos       = ar.pos
                AND ab.alt_allele = ar.alt_allele
            WHERE ar.sample_id IS NULL OR ar.any_passing
        )"""

# ── IGV integration (sidebar) ─────────────────────────────────────────────────
st.sidebar.divider()
st.sidebar.header("🧭 IGV Integration")
auto_launch_igv = st.sidebar.checkbox(
    "Auto-launch IGV",
    value=_cfg.get("auto_launch_igv", True),
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
target_regions = st.sidebar.text_input(
    "Target regions (optional)",
    value=_cfg.get("target_regions", ""),
    help="Path to a BED or interval list file. When set, it is added as a track in every IGV session.",
)
if _has_data("gnomad_af"):
    gnomad_track = st.sidebar.text_input(
        "gnomAD track (optional)",
        value=_cfg.get("gnomad_track", ""),
        help="Path or URI to a gnomAD VCF/BCF. When set, it is added as a track in IGV sessions.",
    )
    gnomad_track_index = st.sidebar.text_input(
        "gnomAD track index (optional)",
        value=_cfg.get("gnomad_track_index", ""),
        help="Optional explicit .tbi / .csi path or URI for the gnomAD track. Leave blank to infer from the track path.",
    )
    _gnomad_index_path = resolve_index_uri(
        gnomad_track.strip(),
        gnomad_track_index.strip() or None,
    ) if gnomad_track.strip() else None
    if gnomad_track.strip() and _gnomad_index_path:
        _is_local_index = not (
            _gnomad_index_path.startswith("gs://")
            or _gnomad_index_path.startswith("http://")
            or _gnomad_index_path.startswith("https://")
        )
        if _is_local_index and not _os.path.exists(_gnomad_index_path):
            st.sidebar.warning(
                f"gnomAD index not found: {_gnomad_index_path}",
                icon="⚠️",
            )
        else:
            st.sidebar.caption(f"gnomAD index: {_gnomad_index_path}")
    elif gnomad_track.strip():
        st.sidebar.warning(
            "Could not infer a gnomAD index path from the track path. Set 'gnomAD track index' explicitly.",
            icon="⚠️",
        )
else:
    gnomad_track = ""
    gnomad_track_index = ""
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
    """Load a TSV manifest into a lookup dict.

    Keys are ``sample_id`` strings by default.  If the manifest has a
    ``pipeline`` column, keys are ``(sample_id, pipeline)`` tuples so that
    the same sample run on multiple pipelines can resolve to different BAMs.
    Both key styles may coexist: rows without a ``pipeline`` value (or from a
    manifest without the column) use the plain ``sample_id`` key.
    """
    mdf = pd.read_csv(p.strip(), sep="\t")
    _has_pipeline_col = "pipeline" in mdf.columns
    result = {}
    for row in mdf.itertuples(index=False):
        bai = str(row.duplex_output_bam_index) if hasattr(row, "duplex_output_bam_index") and pd.notna(row.duplex_output_bam_index) else None
        variants = str(row.final_annotated_variants) if hasattr(row, "final_annotated_variants") and pd.notna(row.final_annotated_variants) else None
        entry = {"bam": str(row.duplex_output_bam), "bai": bai, "variants_tsv": variants}
        sid = str(row.collaborator_sample_id)
        if _has_pipeline_col and pd.notna(row.pipeline) and str(row.pipeline).strip():
            result[(sid, str(row.pipeline).strip())] = entry
        else:
            result[sid] = entry
    return result

manifest = {}
if manifest_path and manifest_path.strip():
    try:
        manifest = load_manifest(manifest_path.strip())
        st.sidebar.success(f"{len(manifest):,} samples loaded from manifest")
    except Exception as e:
        st.sidebar.error(f"Could not load manifest: {e}")

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
                key="advanced_metadata_inputs",
            )

st.sidebar.divider()
with st.sidebar.expander("Filter state", expanded=False):
    st.caption("Save the current sidebar filters to a JSON file, or reload a previously saved state.")
    # ── Save ──────────────────────────────────────────────────────────────────
    _filter_json = MAIN_FILTER_STATE.to_json(st.session_state)
    st.download_button(
        "Save filter state",
        data=_filter_json,
        file_name="geac_filters.json",
        mime="application/json",
        help="Download all current sidebar filter values as a JSON file.",
        use_container_width=True,
    )
    # ── Load ──────────────────────────────────────────────────────────────────
    st.markdown("**Load filter state**")
    _uploaded_filter = st.file_uploader(
        "Upload filter JSON",
        type=["json"],
        key="filter_state_upload",
        label_visibility="collapsed",
    )
    if _uploaded_filter is not None:
        _raw = _uploaded_filter.read()
        _upload_hash = hash(_raw)
        # Only apply once per distinct uploaded file (avoid re-applying on every rerun).
        if st.session_state.get("_filter_upload_hash") != _upload_hash:
            st.session_state["_filter_upload_hash"] = _upload_hash
            _apply_warns = MAIN_FILTER_STATE.apply_json(_raw.decode(), st.session_state)
            for _w in _apply_warns:
                st.warning(_w)
            st.rerun()

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


def launch_igv_session(session_xml: str, bed: str, sort_locus: str = "") -> str:
    """Write session.xml and positions.bed to a temp dir and load in IGV.

    Tries the IGV REST API (localhost:60151) first — works if IGV is already
    running. If the connection is refused, launches IGV via subprocess.

    *sort_locus* (e.g. "chr1:12345") triggers a sort-by-base command via the
    REST API after loading.  IGV's XML session format does not support
    sort-on-load, so this only works when IGV is reachable via the REST API.

    Returns a status message suitable for st.info / st.success / st.error.
    """
    import tempfile, subprocess, urllib.request, urllib.error, time

    tmp = tempfile.mkdtemp(prefix="geac_igv_")
    session_path = _os.path.join(tmp, "session.xml")
    bed_path     = _os.path.join(tmp, "positions.bed")
    with open(session_path, "w") as f:
        f.write(session_xml)
    with open(bed_path, "w") as f:
        f.write(bed)

    def _igv_sort_by_base(locus: str) -> None:
        """Send a sort-by-base command to IGV via the socket command interface.

        Port 60151 accepts plain-text batch commands (not HTTP for sort).
        Format: ``sort base chr:pos``
        """
        if not locus:
            return
        import socket
        try:
            with socket.create_connection(("localhost", 60151), timeout=5) as sock:
                sock.sendall(f"sort base {locus}\n".encode())
                sock.recv(256)  # read "OK" response
        except (OSError, TimeoutError):
            pass  # best-effort; user can sort manually

    # Try REST API first
    url = f"http://localhost:60151/load?file={urllib.request.pathname2url(session_path)}&merge=false"
    try:
        urllib.request.urlopen(url, timeout=8)
        time.sleep(2)  # give IGV time to load BAM tracks before sorting
        _igv_sort_by_base(sort_locus)
        return "Session loaded into running IGV instance."
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


def make_igv_session(
    df: pd.DataFrame,
    manifest: dict,
    genome: str,
    target_regions: str = "",
    gnomad_track: str = "",
    gnomad_track_index: str = "",
) -> tuple[str, str]:
    """Build an IGV session XML and return (xml, sort_locus).

    *sort_locus* is ``"chrom:pos"`` (1-based) for the first locus — used by
    ``launch_igv_session`` to sort reads by base via the REST API.  IGV's XML
    session format does not support sort-on-load, so the XML itself contains
    no sort directives.
    """
    first = df.sort_values(["chrom", "pos"]).iloc[0]
    locus = f"{first['chrom']}:{max(0, int(first['pos']) - 99)}-{int(first['pos']) + 101}"

    resources, tracks = [], []

    if target_regions and target_regions.strip():
        tr = target_regions.strip()
        resources.append(f'        <Resource path="{tr}" name="Target regions"/>')
        tracks.append(f'        <Track id="{tr}" name="Target regions" color="0,100,200" height="40" featureVisibilityWindow="-1"/>')

    if gnomad_track and gnomad_track.strip():
        gt = gnomad_track.strip()
        gt_index = resolve_index_uri(gt, gnomad_track_index.strip() or None)
        index_attr = f' index="{gt_index}"' if gt_index else ""
        resources.append(f'        <Resource path="{gt}" name="gnomAD"{index_attr}/>')
        tracks.append('        <Track id="{0}" name="gnomAD" height="60"/>'.format(gt))

    sort_pos = int(first["pos"]) + 1  # IGV is 1-based; GEAC pos is 0-based
    sort_locus = f"{first['chrom']}:{sort_pos}"

    # Build (sample_id, pipeline, track_label) tuples.  When multiple pipelines
    # are present in df, emit one track per (sample_id, pipeline) pair so both
    # BAMs appear in the session.  Manifest lookup tries "{sid}__{pipeline}"
    # first (for manifests with pipeline-specific entries) then falls back to
    # "{sid}" for backwards compatibility.
    _has_pipeline = "pipeline" in df.columns and df["pipeline"].nunique() > 1
    if _has_pipeline:
        _sid_pipe_pairs = (
            df[["sample_id", "pipeline"]]
            .drop_duplicates()
            .itertuples(index=False, name=None)
        )
        _track_items = [
            (sid, pipe, f"{sid} ({pipe})")
            for sid, pipe in _sid_pipe_pairs
        ]
    else:
        _track_items = [
            (sid, None, str(sid))
            for sid in df["sample_id"].unique().tolist()
        ]

    _seen_bams: set[str] = set()
    for sid, pipe, label in _track_items:
        entry = (
            manifest.get((str(sid), str(pipe))) or manifest.get(str(sid))
            if pipe is not None
            else manifest.get(str(sid))
        )
        if entry:
            bam, bai = entry["bam"], entry["bai"]
            if bam in _seen_bams:
                continue
            _seen_bams.add(bam)
            index_attr = f' index="{bai}"' if bai else ""
            resources.append(f'        <Resource path="{bam}" name="{label}"{index_attr}/>')
            tracks.append(f'        <Track id="{bam}" name="{label}"/>')

    resources.append('        <Resource path="positions.bed" name="Selected positions"/>')
    tracks.append('        <Track id="positions.bed" name="Selected positions" color="255,0,0" height="40"/>')

    xml = (
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
    return xml, sort_locus


IGV_CAP = 5

_IGV_CHUNK = 10_000


def igv_buttons(
    extra_conditions: list[str],
    display_df: pd.DataFrame,
    key: str,
    use_global_filters: bool = True,
):
    """Render IGV Prepare + Download buttons with chunked progress.

    extra_conditions   — SQL WHERE fragments
    display_df         — already-fetched display DataFrame
    key                — unique widget key prefix
    use_global_filters — when True (default), prepend the sidebar ``conditions``
                         to *extra_conditions*; when False, use *extra_conditions*
                         alone (e.g. for the position drill-down, which shows all
                         samples at the locus regardless of sidebar filters).
    """
    if not manifest:
        st.caption("🧭 Add a manifest in the sidebar to enable IGV session download.")
        return

    _where_parts = (conditions + extra_conditions) if use_global_filters else extra_conditions
    _extra_w = " AND ".join(_where_parts) if _where_parts else "TRUE"
    sample_ids = query_distinct_samples(con, table_expr, _extra_w)
    n = len(sample_ids)
    cap_samples = sample_ids[:IGV_CAP]

    missing = [sid for sid in sample_ids if str(sid) not in manifest and not any(
        isinstance(k, tuple) and k[0] == str(sid) for k in manifest
    )]
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
        _cap_list = ", ".join(f"'{_sql_str(s)}'" for s in cap_samples)
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
            ", ".join(f"'{_sql_str(s)}'" for s in cap_samples)
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
        session, sort_locus = make_igv_session(
            igv_df,
            manifest,
            genome,
            target_regions,
            gnomad_track,
            gnomad_track_index,
        )
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("session.xml", session)
            zf.writestr("positions.bed", bed)
        st.session_state[f"{key}_igv"]        = buf.getvalue()
        st.session_state[f"{key}_session"]    = session
        st.session_state[f"{key}_bed"]        = bed
        st.session_state[f"{key}_sort_locus"] = sort_locus

        if auto_launch_igv:
            msg = launch_igv_session(session, bed, sort_locus)
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
conditions = []
if min_alt > 1:
    conditions.append(f"alt_count >= {min_alt}")
if vaf_range != (0.0, 1.0):
    conditions.append(f"alt_count * 1.0 / total_depth BETWEEN {vaf_range[0]} AND {vaf_range[1]}")
if max_alt > 0:
    conditions.append(f"alt_count <= {max_alt}")
if chrom_sel != "All":
    conditions.append(f"chrom = '{chrom_sel}'")
if sample_sel:
    s_list = ", ".join(f"'{_sql_str(s)}'" for s in sample_sel)
    conditions.append(f"sample_id IN ({s_list})")
_sr_lo, _sr_hi = sample_recurrence
if _n_samples_total > 1 and (_sr_lo > 1 or _sr_hi < _n_samples_total):
    # Build a scoping clause so recurrence is counted only among the
    # currently selected cohort subset (sample/batch/label/on-target).
    _scope_parts: list[str] = []
    if sample_sel:
        _scope_parts.append("sample_id IN ({})".format(
            ", ".join(f"'{_sql_str(s)}'" for s in sample_sel)))
    if batch_sel:
        _scope_parts.append("batch IN ({})".format(
            ", ".join(f"'{_sql_str(b)}'" for b in batch_sel)))
    if pipeline_sel:
        _scope_parts.append("pipeline IN ({})".format(
            ", ".join(f"'{_sql_str(p)}'" for p in pipeline_sel)))
    if label1_sel:
        _scope_parts.append("label1 IN ({})".format(
            ", ".join(f"'{_sql_str(v)}'" for v in label1_sel)))
    if label2_sel:
        _scope_parts.append("label2 IN ({})".format(
            ", ".join(f"'{_sql_str(v)}'" for v in label2_sel)))
    if label3_sel:
        _scope_parts.append("label3 IN ({})".format(
            ", ".join(f"'{_sql_str(v)}'" for v in label3_sel)))
    if "on_target" in _schema_cols:
        if on_target_sel == "On target":
            _scope_parts.append("on_target = true")
        elif on_target_sel == "Off target":
            _scope_parts.append("on_target = false")
    _scope_where = " AND ".join(_scope_parts) if _scope_parts else "TRUE"
    _rec_df = _compute_recurrence_loci(path, _sr_lo, _sr_hi, _scope_where)
    con.register("_recurrence_loci", _rec_df)
    conditions.append(
        "(chrom, pos, ref_allele, alt_allele) IN "
        "(SELECT chrom, pos, ref_allele, alt_allele FROM _recurrence_loci)"
    )
if batch_sel:
    b_list = ", ".join(f"'{_sql_str(b)}'" for b in batch_sel)
    conditions.append(f"batch IN ({b_list})")
if pipeline_sel:
    p_list = ", ".join(f"'{_sql_str(p)}'" for p in pipeline_sel)
    conditions.append(f"pipeline IN ({p_list})")
if label1_sel:
    l_list = ", ".join(f"'{_sql_str(v)}'" for v in label1_sel)
    conditions.append(f"label1 IN ({l_list})")
if label2_sel:
    l_list = ", ".join(f"'{_sql_str(v)}'" for v in label2_sel)
    conditions.append(f"label2 IN ({l_list})")
if label3_sel:
    l_list = ", ".join(f"'{_sql_str(v)}'" for v in label3_sel)
    conditions.append(f"label3 IN ({l_list})")
if variant_sel:
    t_list = ", ".join(f"'{_sql_str(t)}'" for t in variant_sel)
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
    vf_list = ", ".join(f"'{_sql_str(v)}'" for v in variant_filter_sel)
    conditions.append(f"variant_filter IN ({vf_list})")
if "on_target" in _schema_cols:
    if on_target_sel == "On target":
        conditions.append("on_target = true")
    elif on_target_sel == "Off target":
        conditions.append("on_target = false")
if "gnomad_af" in _schema_cols:
    _af_lo = float(gnomad_af_range[0])
    _af_hi = float(gnomad_af_range[1])
    _af_filtered = not (_af_lo == 0.0 and _af_hi == 1.0)
    if _af_filtered and gnomad_include_null:
        conditions.append(f"(gnomad_af BETWEEN {_af_lo} AND {_af_hi} OR gnomad_af IS NULL)")
    elif _af_filtered:
        conditions.append(f"gnomad_af BETWEEN {_af_lo} AND {_af_hi}")
    elif not gnomad_include_null:
        conditions.append("gnomad_af IS NOT NULL")
if gene_text.strip() and "gene" in _schema_cols:
    _gene_escaped = _sql_str(gene_text.strip())
    conditions.append(f"gene = '{_gene_escaped}'")

if _repeat_cols_present:
    if homopolymer_range != (0, 20):
        conditions.append(f"(homopolymer_len IS NULL OR homopolymer_len BETWEEN {homopolymer_range[0]} AND {homopolymer_range[1]})")
    if str_len_range != (0, 50):
        conditions.append(f"(str_len IS NULL OR str_len BETWEEN {str_len_range[0]} AND {str_len_range[1]})")

where = " AND ".join(conditions) if conditions else "TRUE"

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

render_summary_metrics(stats, fstats, pct_called, fpct_called)

def query_records(extra: list[str] = [], limit: int | None = None) -> pd.DataFrame:
    """Query records with current filters plus any extra conditions."""
    w = " AND ".join(conditions + extra)
    limit_clause = f"LIMIT {limit}" if limit is not None else ""
    return con.execute(f"""
        SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
        FROM {table_expr}
        WHERE {w}
        ORDER BY chrom, pos, alt_allele, sample_id
        {limit_clause}
    """).df()

total_count = con.execute(f"SELECT COUNT(*) FROM {table_expr} WHERE {where}").fetchone()[0]

if total_count == 0:
    st.warning("No records match the current filters.", icon="🔎")
    st.stop()

st.info(f"**{total_count:,}** records match the current filters.", icon="✅")

_reads_banner = st.empty()
if _reads_active:
    _active_parts = []
    if _fs_has_data and (_fs_lo > 0 or _fs_hi < _fs_max):
        _active_parts.append(f"family size: {_fs_lo}–{_fs_hi}")
    if _cycle_lo > 1 or _cycle_hi < _cycle_max:
        _active_parts.append(f"cycle number: {_cycle_lo}–{_cycle_hi}")
    if _mq_lo > 0 or _mq_hi < _mq_max:
        _active_parts.append(f"map qual: {_mq_lo}–{_mq_hi}")
    if _is_has_data and (_is_lo > _IS_MIN or _is_hi < _IS_MAX):
        _active_parts.append(f"insert size: {_is_lo}–{_is_hi}")
    if read_strand_sel != "All":
        _active_parts.append(read_strand_sel.lower())
    _reads_banner.warning(
        f"**Per-read filters active** ({'; '.join(_active_parts)}). "
        f"{per_read_warning_note(recompute_vaf)}"
    )


def _append_provenance_row(
    rows: list[dict[str, str]],
    section: str,
    name: str,
    value,
    *,
    active: bool = True,
) -> None:
    if not active:
        return
    if isinstance(value, bool):
        value_str = "true" if value else "false"
    elif isinstance(value, (list, tuple)):
        value_str = ", ".join(str(v) for v in value)
    else:
        value_str = str(value)
    rows.append({"section": section, "name": name, "value": value_str})


def _build_active_filter_provenance(
    *,
    discovery_mode: str,
    discovery_items: list[tuple[str, object]] | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    _append_provenance_row(rows, "data", "data_file", path)
    _append_provenance_row(rows, "query", "where_sql", where)
    _append_provenance_row(rows, "filters", "chromosome", chrom_sel, active=chrom_sel != "All")
    _append_provenance_row(rows, "filters", "samples", sample_sel, active=bool(sample_sel))
    _append_provenance_row(
        rows,
        "filters",
        "sample_recurrence",
        f"{_sr_lo}-{_sr_hi}",
        active=_n_samples_total > 1 and (_sr_lo > 1 or _sr_hi < _n_samples_total),
    )
    _append_provenance_row(rows, "filters", "batch", batch_sel, active=bool(batch_sel))
    _append_provenance_row(rows, "filters", "label1", label1_sel, active=bool(label1_sel))
    _append_provenance_row(rows, "filters", "label2", label2_sel, active=bool(label2_sel))
    _append_provenance_row(rows, "filters", "label3", label3_sel, active=bool(label3_sel))
    _append_provenance_row(
        rows,
        "filters",
        "gene",
        gene_text.strip(),
        active=bool(gene_text.strip()) and "gene" in _schema_cols,
    )
    _append_provenance_row(
        rows,
        "filters",
        "variant_type",
        variant_sel,
        active=bool(variant_sel) and set(variant_sel) != {"SNV", "insertion", "deletion"},
    )
    _append_provenance_row(rows, "filters", "vaf_range", f"{vaf_range[0]}-{vaf_range[1]}", active=vaf_range != (0.0, 1.0))
    _append_provenance_row(rows, "filters", "min_alt_count", min_alt, active=min_alt > 1)
    _append_provenance_row(rows, "filters", "max_alt_count", max_alt, active=max_alt > 0)
    _append_provenance_row(rows, "filters", "min_fwd_alt_count", min_fwd_alt, active=min_fwd_alt > 0)
    _append_provenance_row(rows, "filters", "min_rev_alt_count", min_rev_alt, active=min_rev_alt > 0)
    _append_provenance_row(rows, "filters", "min_overlap_alt_agree", min_overlap_agree, active=min_overlap_agree > 0)
    _append_provenance_row(rows, "filters", "min_overlap_alt_disagree", min_overlap_disagree, active=min_overlap_disagree > 0)
    _append_provenance_row(rows, "filters", "min_depth", min_depth, active=min_depth > 0)
    _append_provenance_row(rows, "filters", "max_depth", max_depth, active=max_depth > 0)
    _append_provenance_row(rows, "filters", "variant_called", variant_called_sel, active=variant_called_sel != "All")
    _append_provenance_row(rows, "filters", "variant_filter", variant_filter_sel, active=bool(variant_filter_sel))
    _append_provenance_row(rows, "filters", "target_bases", on_target_sel, active=on_target_sel != "All")
    _append_provenance_row(
        rows,
        "filters",
        "gnomad_af_range",
        f"{gnomad_af_range[0]}-{gnomad_af_range[1]}",
        active=("gnomad_af" in _schema_cols) and gnomad_af_range != ("0", "1.0"),
    )
    _append_provenance_row(
        rows,
        "filters",
        "include_sites_absent_from_gnomad",
        gnomad_include_null,
        active=("gnomad_af" in _schema_cols) and (gnomad_af_range != ("0", "1.0") or not gnomad_include_null),
    )
    _append_provenance_row(
        rows,
        "filters",
        "homopolymer_length_range",
        f"{homopolymer_range[0]}-{homopolymer_range[1]}",
        active=_repeat_cols_present and homopolymer_range != (0, 20),
    )
    _append_provenance_row(
        rows,
        "filters",
        "str_length_range",
        f"{str_len_range[0]}-{str_len_range[1]}",
        active=_repeat_cols_present and str_len_range != (0, 50),
    )
    _append_provenance_row(
        rows,
        "per_read_filters",
        "application_mode",
        "recompute alt count from filtered reads" if recompute_vaf else "hide loci with no passing reads",
        active=_reads_active,
    )
    _append_provenance_row(
        rows,
        "per_read_filters",
        "family_size",
        f"{_fs_lo}-{_fs_hi}",
        active=_has_alt_reads and _fs_has_data and (_fs_lo > 0 or _fs_hi < _fs_max),
    )
    _append_provenance_row(
        rows,
        "per_read_filters",
        "cycle_number",
        f"{_cycle_lo}-{_cycle_hi}",
        active=_has_alt_reads and (_cycle_lo > 1 or _cycle_hi < _cycle_max),
    )
    _append_provenance_row(
        rows,
        "per_read_filters",
        "mapping_quality",
        f"{_mq_lo}-{_mq_hi}",
        active=_has_alt_reads and (_mq_lo > 0 or _mq_hi < _mq_max),
    )
    _append_provenance_row(
        rows,
        "per_read_filters",
        "insert_size",
        f"{_is_lo}-{_is_hi}",
        active=_has_alt_reads and _is_has_data and (_is_lo > _IS_MIN or _is_hi < _IS_MAX),
    )
    _append_provenance_row(
        rows,
        "per_read_filters",
        "read",
        read_strand_sel,
        active=_has_alt_reads and read_strand_sel != "All",
    )
    _append_provenance_row(rows, "signature_discovery", "discovery_mode", discovery_mode)
    if discovery_items:
        for _name, _value in discovery_items:
            _append_provenance_row(
                rows,
                "signature_discovery",
                _name,
                _value,
                active=_value not in (None, "", [], ()),
            )

    return pd.DataFrame(rows, columns=["section", "name", "value"])

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
        "overlap_alt_disagree", "variant_called", "variant_filter", "on_target", "gene", "gnomad_af",
    ]
    if c in df.columns
]

_tbl_event = render_records_table(
    df,
    total_count,
    _table_cols,
    igv_buttons=igv_buttons,
    key="main",
)

# ── Position-level drill-down ──────────────────────────────────────────────────
_selected_rows = (_tbl_event.selection or {}).get("rows", [])
if _selected_rows:
    _row = df.iloc[_selected_rows[0]]
    _new_locus = (str(_row["chrom"]), int(_row["pos"]), str(_row["alt_allele"]))
    if st.session_state.get("_drill_locus") != _new_locus:
        st.session_state["_drill_locus"] = _new_locus

if "_drill_locus" in st.session_state:
    _chrom, _pos, _selected_alt = st.session_state["_drill_locus"]
    render_position_drilldown(
        con=con,
        table_expr=table_expr,
        schema_cols=_schema_cols,
        chrom=_chrom,
        pos=_pos,
        selected_alt=_selected_alt,
        has_alt_reads=_has_alt_reads,
        sql_str=_sql_str,
        igv_buttons=igv_buttons,
    )

# ── Shared alt_reads join subquery (used by Duplex/Simplex and Reads tabs) ────
# Joins alt_reads to the current filtered locus set so all reads plots respect
# the active sidebar filters.  Defined here so both tabs can share it.
_r_reads_filter = f"WHERE {_reads_where}" if _reads_active else ""
_r_join = f"""
    (SELECT * FROM alt_reads {_r_reads_filter}) ar
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

# ── Plots ─────────────────────────────────────────────────────────────────────
tab1, tab2, tab3, tab4, tab_cohort, tab_reads, tab_duplex, tab_tn, tab_pon, tab_pipeline, tab_read_type = st.tabs(
    [module.LABEL for module in TAB_MODULES],
    key="main_tabs",
)

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
              AND total_depth > 0
              AND alt_count <= total_depth
            GROUP BY vaf_bin, vaf_bin_end
            HAVING vaf_bin IS NOT NULL AND vaf_bin >= 0.0
            ORDER BY vaf_bin
        """).df()

        if counts.empty:
            st.info(f"No {vtype}s in current selection.")
        else:
            _sel_name = f"bar_click_{vtype}"
            sel_param = alt.selection_point(
                name=_sel_name,
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
            _vtype_count = int(counts["count"].sum())
            event = st.altair_chart(chart, width="stretch", on_select="rerun", key=f"vaf_chart_{vtype}")
            st.caption(
                f"{_vtype_count:,} alt-allele records (one per unique alt allele observed "
                f"at a locus in a sample). Click a bar to drill down."
            )

            pts = (event.selection or {}).get(_sel_name, [])
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
                    st.dataframe(sel[_table_cols], width="stretch")
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

    # ── Helpers for stratified SBS96 spectra (used in Error Spectrum and Reads tabs) ──
    def _to_spec96_strat(raw_df):
        if raw_df.empty:
            return None, 0
        df = raw_df.copy()
        df["sbs_label"] = df.apply(
            lambda r: _sbs_label(r["trinuc_context"], r["ref_allele"], r["alt_allele"]),
            axis=1,
        )
        df = df.dropna(subset=["sbs_label"])
        df["mut_type"] = df["sbs_label"].str.extract(r'\[([A-Z]>[A-Z])\]')[0]
        agg = df.groupby(["sbs_label", "mut_type"], as_index=False)["count"].sum()
        full = pd.DataFrame({
            "sbs_label": _SBS_ORDER,
            "mut_type":  [lbl[2:5] for lbl in _SBS_ORDER],
        })
        s96 = full.merge(agg, on=["sbs_label", "mut_type"], how="left")
        s96["count"] = s96["count"].fillna(0).astype(int)
        total = int(s96["count"].sum())
        s96["fraction"] = s96["count"] / total if total > 0 else 0.0
        return s96, total

    def _strat_sbs96_chart(spec_df, title, y_max=None, sel_name=None):
        _y_scale = alt.Scale(domain=[0, y_max]) if y_max is not None else alt.Undefined
        _sel = (
            alt.selection_point(name=sel_name, fields=["sbs_label"], on="click")
            if sel_name else None
        )
        panels = []
        for _mt in _SBS_MUT_TYPES:
            _s = spec_df[spec_df["mut_type"] == _mt]
            _order = [lbl for lbl in _SBS_ORDER if f"[{_mt}]" in lbl]
            _enc = dict(
                x=alt.X("sbs_label:N", sort=_order, title=None,
                         axis=alt.Axis(labelAngle=-90, labelFontSize=7)),
                y=alt.Y("fraction:Q", title="Fraction",
                         scale=_y_scale, axis=alt.Axis(format=".3f")),
                tooltip=["sbs_label:N",
                         alt.Tooltip("fraction:Q", format=".3f", title="Fraction")],
            )
            if _sel is not None:
                # Reference by name dict to avoid auto-embedding the param definition
                # in each sub-chart (which would cause Altair deduplication warnings).
                _enc["opacity"] = alt.condition({"param": sel_name}, alt.value(1.0), alt.value(0.4))
            panels.append(
                alt.Chart(_s).mark_bar(color=_SBS_COLORS[_mt]).encode(**_enc).properties(
                    title=alt.TitleParams(_mt, color=_SBS_COLORS[_mt],
                                          fontSize=11, fontWeight="bold"),
                    width=120, height=110,
                )
            )
        chart = (
            alt.concat(*panels, columns=3)
            .resolve_scale(y="shared")
            .properties(title=alt.TitleParams(title, fontSize=13))
        )
        if _sel is not None:
            chart = chart.add_params(_sel)
        return chart

    def _load_sample_sbs96_matrix():
        _has_batch = _has_data("batch")
        _id_sql = "sample_id || ' / ' || batch" if _has_batch else "sample_id"
        _group_by = "sample_id, batch" if _has_batch else "sample_id"
        _raw = con.execute(f"""
            SELECT {_id_sql} AS sample_label,
                   trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
            FROM (SELECT * FROM {table_expr}) _t
            WHERE {where} AND variant_type = 'SNV'
              AND trinuc_context IS NOT NULL AND length(trinuc_context) = 3
            GROUP BY {_group_by}, trinuc_context, ref_allele, alt_allele
        """).df()

        if _raw.empty:
            return pd.DataFrame(columns=_SBS_ORDER, dtype=float)

        _rows = []
        for _sid, _grp in _raw.groupby("sample_label"):
            _grp = _grp.copy()
            _grp["sbs_label"] = _grp.apply(
                lambda r: _sbs_label(r["trinuc_context"], r["ref_allele"], r["alt_allele"]),
                axis=1,
            )
            _grp = _grp.dropna(subset=["sbs_label"])
            _agg = _grp.groupby("sbs_label")["count"].sum()
            _counts = (
                pd.Series(0.0, index=_SBS_ORDER)
                .add(_agg, fill_value=0)
                .reindex(_SBS_ORDER)
            )
            if float(_counts.sum()) > 0:
                _rows.append(pd.Series(_counts.values, index=_SBS_ORDER, name=_sid))

        if not _rows:
            return pd.DataFrame(columns=_SBS_ORDER, dtype=float)

        return pd.DataFrame(_rows, dtype=float)

    def _profile_to_spec96(profile, value_name="fraction"):
        _df = pd.DataFrame({
            "sbs_label": _SBS_ORDER,
            "mut_type": [lbl[2:5] for lbl in _SBS_ORDER],
            value_name: pd.Series(profile, index=_SBS_ORDER).reindex(_SBS_ORDER).values,
        })
        return _df

    def _signature_profile_chart(spec_df, title, overlay_df=None, overlay_label=None):
        _sub_charts = []
        for _mt in _SBS_MUT_TYPES:
            _sub = spec_df[spec_df["mut_type"] == _mt]
            _order = [lbl for lbl in _SBS_ORDER if f"[{_mt}]" in lbl]
            _bars = (
                alt.Chart(_sub)
                .mark_bar(color=_SBS_COLORS[_mt])
                .encode(
                    alt.X("sbs_label:N", sort=_order, title=None,
                          axis=alt.Axis(labelAngle=-90, labelFontSize=7)),
                    alt.Y("fraction:Q", title="Fraction",
                          axis=alt.Axis(format=".3f")),
                    tooltip=[
                        "sbs_label:N",
                        alt.Tooltip("fraction:Q", format=".3f", title="Fraction"),
                    ],
                )
            )
            if overlay_df is not None:
                _overlay_sub = overlay_df[overlay_df["mut_type"] == _mt]
                _overlay = (
                    alt.Chart(_overlay_sub)
                    .mark_point(color="black", size=14, filled=True, opacity=0.85)
                    .encode(
                        alt.X("sbs_label:N", sort=_order),
                        alt.Y("fraction:Q"),
                        tooltip=[
                            "sbs_label:N",
                            alt.Tooltip("fraction:Q", format=".3f", title=overlay_label or "Overlay"),
                        ],
                    )
                )
                _panel = alt.layer(_bars, _overlay)
            else:
                _panel = _bars

            _sub_charts.append(
                _panel.properties(
                    title=alt.TitleParams(_mt, color=_SBS_COLORS[_mt], fontSize=11, fontWeight="bold"),
                    width=120,
                    height=110,
                )
            )

        return (
            alt.concat(*_sub_charts, columns=3)
            .resolve_scale(y="shared")
            .properties(title=alt.TitleParams(title, fontSize=13))
        )

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

            # ── COSMIC controls (above chart so reconstruction data is ready) ──
            _cos_col, _path_col = st.columns([1, 3])
            cosmic_path = _path_col.text_input(
                "COSMIC SBS matrix path (optional — enables reconstruction overlay)",
                value=_cfg.get("cosmic", ""),
                placeholder="/path/to/COSMIC_v3.4_SBS_GRCh37.txt",
                key="cosmic_path",
                label_visibility="visible",
            )

            # Shared with Called vs Uncalled section; set when COSMIC loads successfully
            _cosmic_W        = None
            _cosmic_aligned  = None
            reconstructed    = None
            recon_df         = None
            sig_df           = None
            top_df           = None
            top_n_sig        = None
            cos_sim          = None
            residual_pct     = None

            if cosmic_path and cosmic_path.strip():
                try:
                    cosmic_df      = _load_cosmic(cosmic_path.strip())
                    cosmic_aligned = cosmic_df.reindex(_SBS_ORDER)
                    missing        = cosmic_aligned.isna().any(axis=1).sum()
                    if missing > 0:
                        st.warning(
                            f"{missing} context(s) not found in COSMIC matrix — "
                            "check that the file uses the standard A[C>A]A format."
                        )
                    else:
                        W              = cosmic_aligned.values.astype(float)
                        _cosmic_W      = W
                        _cosmic_aligned = cosmic_aligned
                        obs = (
                            spec96.set_index("sbs_label")["count"]
                            .reindex(_SBS_ORDER).fillna(0).values.astype(float)
                        )

                        h, _ = nnls(W, obs)
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

                        top_n_sig = _cos_col.slider(
                            "Top signatures", 3, min(20, len(sig_df)), 4,
                            key="top_n_sig",
                            help="Number of top signatures used for the reconstruction overlay and exposure chart.",
                        )
                        top_df = sig_df.head(top_n_sig)

                        # Refit to top-N
                        _top_sig_names = top_df["signature"].tolist()
                        W_refit        = cosmic_aligned[_top_sig_names].values.astype(float)
                        h_refit, _     = nnls(W_refit, obs)
                        reconstructed  = W_refit @ h_refit

                        cos_sim = (
                            float(np.dot(obs, reconstructed))
                            / (np.linalg.norm(obs) * np.linalg.norm(reconstructed) + 1e-12)
                        )
                        residual_pct = (
                            float(np.linalg.norm(obs - reconstructed))
                            / (float(obs.sum()) + 1e-12) * 100
                        )

                        recon_df = spec96[["sbs_label", "mut_type"]].copy()
                        recon_df["recon_count"] = reconstructed
                        recon_df["recon_frac"]  = (
                            reconstructed / _total_snvs if _total_snvs > 0 else 0.0
                        )

                except Exception as exc:
                    st.error(f"Failed to load COSMIC matrix: {exc}")

            _recon_y   = "recon_frac"  if _sbs_use_fraction else "recon_count"
            _recon_fmt = ".3f"         if _sbs_use_fraction else "d"

            # ── Unified SBS96 chart: bars (observed) + optional dots (reconstruction) ──
            _chart_title = (
                f"SNV Trinucleotide Spectrum — bars = observed, dots = reconstruction (top {top_n_sig} sigs)"
                if recon_df is not None else
                "SNV Trinucleotide Spectrum (SBS96)"
            )
            sel_param = alt.selection_point(name="bar_click", fields=["sbs_label"], on="click")

            _sub_charts = []
            for _mt in _SBS_MUT_TYPES:
                _obs_sub = spec96[spec96["mut_type"] == _mt].copy()
                _order   = [lbl for lbl in _SBS_ORDER if f"[{_mt}]" in lbl]

                _bars = (
                    alt.Chart(_obs_sub)
                    .mark_bar(color=_SBS_COLORS[_mt])
                    .encode(
                        alt.X("sbs_label:N", sort=_order, title=None,
                              axis=alt.Axis(labelAngle=-90, labelFontSize=8)),
                        alt.Y(f"{_sbs_y_field}:Q", title=_sbs_y_title,
                              **({"axis": alt.Axis(format=".3f")} if _sbs_use_fraction else {})),
                        opacity=alt.condition({"param": "bar_click"}, alt.value(1.0), alt.value(0.4)),
                        tooltip=[
                            "sbs_label:N",
                            alt.Tooltip(f"{_sbs_y_field}:Q", title="Observed", format=_sbs_y_fmt),
                        ],
                    )
                )

                if recon_df is not None:
                    _recon_sub = recon_df[recon_df["mut_type"] == _mt].copy()
                    _dots = (
                        alt.Chart(_recon_sub)
                        .mark_point(color="black", size=15, filled=True, opacity=0.85)
                        .encode(
                            alt.X("sbs_label:N", sort=_order),
                            alt.Y(f"{_recon_y}:Q"),
                            tooltip=alt.value(None),
                        )
                    )
                    _panel = alt.layer(_bars, _dots)
                else:
                    _panel = _bars

                _sub_charts.append(
                    _panel.properties(
                        title=alt.TitleParams(_mt, color=_SBS_COLORS[_mt], fontSize=13, fontWeight="bold"),
                        width=150, height=140,
                    )
                )

            chart = (
                alt.concat(*_sub_charts, columns=3)
                .resolve_scale(y="shared")
                .properties(title=alt.TitleParams(_chart_title, fontSize=14))
                .add_params(sel_param)
            )
            event = st.altair_chart(chart, width="stretch", on_select="rerun", key="sbs96_spectrum")

            if recon_df is not None:
                st.caption(
                    f"{_total_snvs:,} SNV alt-allele loci. "
                    f"Black dots = COSMIC reconstruction refit to top {top_n_sig} signatures. "
                    "Click bars to drill down. "
                    "Contexts where dots deviate from bars are poorly explained by the selected signatures."
                )
            else:
                st.caption(
                    f"{_total_snvs:,} SNV alt-allele loci. "
                    "Click one or more bars to drill down and open in IGV. "
                    "Enter a COSMIC matrix path above to overlay the reconstruction."
                )

            # ── Click drill-down + IGV ─────────────────────────────────────────
            pts = (event.selection or {}).get("bar_click", [])
            if pts:
                clicked_labels = [p.get("sbs_label") for p in pts if p.get("sbs_label")]
                if clicked_labels:
                    matching = raw[raw["sbs_label"].isin(clicked_labels)][
                        ["trinuc_context", "ref_allele", "alt_allele"]
                    ]
                    if not matching.empty:
                        or_clauses = " OR ".join(
                            f"(trinuc_context = '{r.trinuc_context}' AND ref_allele = '{r.ref_allele}' AND alt_allele = '{r.alt_allele}')"
                            for r in matching.itertuples(index=False)
                        )
                        extra_cond = f"variant_type = 'SNV' AND ({or_clauses})"
                        sel = query_records([extra_cond])
                        label_str = ", ".join(clicked_labels)
                        st.caption(f"{len(sel):,} records matching {len(clicked_labels)} selected context(s): {label_str}")
                        st.dataframe(sel[_table_cols], width="stretch")
                        igv_buttons([extra_cond], sel, key=f"sbs_{'_'.join(clicked_labels)}")

            # ── COSMIC results (below chart) ───────────────────────────────────
            if cos_sim is not None:
                st.divider()
                st.subheader("🌌 COSMIC Signature Decomposition")
                fit_col1, fit_col2 = st.columns(2)
                fit_col1.metric(
                    "Cosine similarity",
                    f"{cos_sim:.4f}",
                    help="1.0 = perfect reconstruction using the top-N signatures. Values above 0.95 indicate a good fit.",
                )
                fit_col2.metric(
                    "Residual (% of counts)",
                    f"{residual_pct:.1f}%",
                    help="L2 norm of unexplained counts as a percentage of total SNV count (refit to top-N). Lower is better.",
                )

                _all_sigs = list(sig_df["signature"])
                sig_chart = (
                    alt.Chart(top_df)
                    .mark_bar()
                    .encode(
                        alt.X("signature:N", sort=list(top_df["signature"]), title="Signature"),
                        alt.Y("exposure:Q", title="Exposure (proportion)",
                              axis=alt.Axis(format=".0%")),
                        alt.Color("signature:N",
                                  scale=alt.Scale(domain=_all_sigs), legend=None),
                        tooltip=[
                            "signature:N",
                            alt.Tooltip("exposure:Q", format=".2%", title="Exposure"),
                            alt.Tooltip("etiology:N", title="Etiology"),
                        ],
                    )
                    .properties(title=f"Top {top_n_sig} COSMIC SBS Signatures (NNLS fit)", height=300)
                )
                st.altair_chart(sig_chart, width="stretch")

                display = top_df.copy()
                display["exposure"] = display["exposure"].map("{:.2%}".format)
                st.dataframe(display, width="stretch", hide_index=True)

                # ── Per-sample COSMIC decomposition (cohort stacked bar) ───────
                if data_source.is_duckdb:
                    st.divider()
                    st.subheader("Per-sample Signature Exposures")
                    st.caption(
                        "Each sample fitted independently against the full COSMIC matrix (NNLS). "
                        "Signatures with no exposure in any sample are hidden."
                    )

                    _sample_matrix = _load_sample_sbs96_matrix()
                    if _sample_matrix.empty:
                        st.info("No SNVs with trinucleotide context in current selection.")
                    else:
                        _ps_rows = []
                        for _sid, _row in _sample_matrix.iterrows():
                            _obs = _row.reindex(_SBS_ORDER).values.astype(float)
                            _n_snvs = int(_row.sum())
                            if _n_snvs == 0:
                                continue
                            _h, _ = nnls(_cosmic_W, _obs)
                            _total = _h.sum()
                            _h_norm = _h / _total if _total > 0 else _h
                            for _sig, _exp in zip(_cosmic_aligned.columns, _h_norm):
                                if _exp > 0:
                                    _ps_rows.append({
                                        "sample_label": _sid,
                                        "signature": _sig,
                                        "exposure":  float(_exp),
                                        "n_snvs":    _n_snvs,
                                        "etiology":  _SBS_ETIOLOGY.get(_sig, ""),
                                    })

                        if not _ps_rows:
                            st.info("NNLS returned no exposures for any sample.")
                        else:
                            _ps_df = pd.DataFrame(_ps_rows)

                            # Keep only signatures present in at least one sample
                            _ps_sigs = _ps_df["signature"].unique().tolist()

                            # Sort samples alphabetically
                            _ps_order = sorted(_ps_df["sample_label"].unique().tolist())

                            # Fill zeros for sample/signature combos with no exposure
                            _ps_full = (
                                pd.MultiIndex.from_product(
                                    [_ps_order, _ps_sigs],
                                    names=["sample_label", "signature"],
                                )
                                .to_frame(index=False)
                                .merge(_ps_df[["sample_label", "signature", "exposure", "etiology"]],
                                       on=["sample_label", "signature"], how="left")
                            )
                            _ps_full["exposure"] = _ps_full["exposure"].fillna(0.0)
                            _ps_full["etiology"] = _ps_full["etiology"].fillna("")

                            _ps_chart = (
                                alt.Chart(_ps_full)
                                .mark_rect()
                                .encode(
                                    alt.X("signature:N", sort=_ps_sigs, title="Signature",
                                          axis=alt.Axis(labelAngle=-45, labelLimit=200)),
                                    alt.Y("sample_label:N", sort=_ps_order, title="Sample"),
                                    alt.Color("exposure:Q", title="Exposure",
                                              scale=alt.Scale(scheme="blues"),
                                              legend=alt.Legend(format=".0%")),
                                    tooltip=[
                                        "sample_label:N",
                                        "signature:N",
                                        alt.Tooltip("exposure:Q", format=".2%", title="Exposure"),
                                        "etiology:N",
                                    ],
                                )
                                .properties(
                                    title="Per-sample COSMIC SBS signature exposures (NNLS)",
                                    height=max(150, 22 * len(_ps_order)),
                                )
                            )
                            st.altair_chart(_ps_chart, width="stretch")
                            st.caption(
                                f"{len(_ps_order)} samples · {len(_ps_sigs)} active signatures. "
                                "Color intensity = exposure proportion. "
                                "Signatures with no exposure in any sample are hidden."
                            )

            # ── NMF signature discovery (cohort / DuckDB only) ────────────────
            if data_source.is_duckdb:
                st.divider()
                st.subheader("NMF Signature Discovery")
                st.caption(
                    "Learn de novo SBS96 signatures across the currently selected samples. "
                    "This uses the current sidebar filters and requires at least two samples with SNVs."
                )

                _enable_nmf = st.checkbox(
                    "Run NMF signature discovery",
                    value=False,
                    key="run_nmf_signatures",
                )

                if _enable_nmf:
                    _sample_matrix = _load_sample_sbs96_matrix()
                    _n_samples_nmf = _sample_matrix.shape[0]

                    if _n_samples_nmf < 2:
                        st.info("NMF requires at least two samples with SNVs in the current selection.")
                    else:
                        _hybrid_available = (
                            _cosmic_aligned is not None
                            and sig_df is not None
                            and not sig_df.empty
                        )
                        if _hybrid_available:
                            _nmf_mode = st.radio(
                                "Discovery mode",
                                ["COSMIC-guided + one learned signature", "De novo NMF"],
                                horizontal=True,
                                key="nmf_mode",
                                help="Use top COSMIC signatures as fixed basis functions and learn one additional non-negative signature, or run fully de novo NMF.",
                            )
                        else:
                            _nmf_mode = "De novo NMF"

                        if _nmf_mode == "De novo NMF":
                            _nmf_max = min(8, _n_samples_nmf)
                            _nmf_default = min(3, _nmf_max)
                            _nmf_components = st.slider(
                                "Number of discovered signatures",
                                2,
                                _nmf_max,
                                _nmf_default,
                                key="nmf_components",
                            )

                            try:
                                _nmf_result = fit_sbs_nmf(_sample_matrix, _nmf_components)
                            except Exception as exc:
                                st.error(f"Could not run NMF signature discovery: {exc}")
                            else:
                                _nmf_profiles = _nmf_result["profiles"]
                                _nmf_exposure_fractions = _nmf_result["exposure_fractions"]
                                _nmf_matrix_cosine = _nmf_result["matrix_cosine"]
                                _nmf_relative_error_pct = _nmf_result["relative_error_pct"]
                                _nmf_iter = _nmf_result["n_iter"]

                                _m1, _m2, _m3, _m4 = st.columns(4)
                                _m1.metric("Samples", f"{_n_samples_nmf:,}")
                                _m2.metric("Signatures", str(_nmf_components))
                                _m3.metric("Matrix cosine", f"{_nmf_matrix_cosine:.4f}")
                                _m4.metric("Relative error", f"{_nmf_relative_error_pct:.1f}%")
                                st.caption(
                                    f"NMF converged in {_nmf_iter} iteration(s). "
                                    "Signatures are shown as normalized SBS96 profiles."
                                )

                                _nmf_long = (
                                    _nmf_exposure_fractions.reset_index(names="sample_label")
                                    .melt(id_vars="sample_label", var_name="signature", value_name="exposure")
                                )
                                _sig_order = _nmf_profiles.index.tolist()
                                _sample_order = sorted(_nmf_exposure_fractions.index.tolist())
                                _nmf_chart = (
                                    alt.Chart(_nmf_long)
                                    .mark_rect()
                                    .encode(
                                        alt.X("signature:N", sort=_sig_order, title="Discovered signature"),
                                        alt.Y("sample_label:N", sort=_sample_order, title="Sample"),
                                        alt.Color("exposure:Q", title="Exposure",
                                                  scale=alt.Scale(scheme="oranges"),
                                                  legend=alt.Legend(format=".0%")),
                                        tooltip=[
                                            "sample_label:N",
                                            "signature:N",
                                            alt.Tooltip("exposure:Q", format=".2%", title="Exposure"),
                                        ],
                                    )
                                    .properties(
                                        title="Per-sample NMF signature exposures",
                                        height=max(150, 22 * len(_sample_order)),
                                    )
                                )
                                st.altair_chart(_nmf_chart, width="stretch")

                                _nmf_match_df = None
                                if _cosmic_aligned is not None:
                                    try:
                                        _nmf_match_df, _ = compare_signatures_to_cosmic(
                                            _nmf_profiles,
                                            _cosmic_aligned,
                                        )
                                    except Exception as exc:
                                        st.warning(f"Could not compare NMF signatures to COSMIC: {exc}")
                                    else:
                                        _nmf_match_df["etiology"] = _nmf_match_df["most_similar_cosmic_signature"].map(
                                            lambda s: _SBS_ETIOLOGY.get(s, "")
                                        )
                                        _display = _nmf_match_df.copy()
                                        _display["most_similar_cosine_similarity"] = _display[
                                            "most_similar_cosine_similarity"
                                        ].map("{:.3f}".format)
                                        st.markdown("**Best COSMIC match for each discovered signature**")
                                        st.dataframe(_display, width="stretch", hide_index=True)
                                elif cosmic_path and cosmic_path.strip():
                                    st.info(
                                        "COSMIC comparison is unavailable because the matrix did not load successfully."
                                    )
                                else:
                                    st.info(
                                        "Provide a COSMIC SBS matrix path above to compare discovered signatures."
                                    )

                                _sig_tabs = st.tabs(_sig_order)
                                for _tab, _sig_name in zip(_sig_tabs, _sig_order):
                                    with _tab:
                                        _sig_spec = _profile_to_spec96(
                                            _nmf_profiles.loc[_sig_name],
                                            value_name="fraction",
                                        )
                                        _overlay_spec = None
                                        _overlay_label = None
                                        _caption = "De novo SBS96 signature learned by NMF."
                                        if _nmf_match_df is not None and not _nmf_match_df.empty:
                                            _match_row = (
                                                _nmf_match_df[_nmf_match_df["signature"] == _sig_name]
                                                .iloc[0]
                                            )
                                            _cosmic_sig = _match_row["most_similar_cosmic_signature"]
                                            _overlay_spec = _profile_to_spec96(
                                                _cosmic_aligned[_cosmic_sig],
                                                value_name="fraction",
                                            )
                                            _overlay_label = _cosmic_sig
                                            _caption = (
                                                f"Best COSMIC match: {_cosmic_sig} "
                                                f"(cosine {_match_row['most_similar_cosine_similarity']:.3f})."
                                            )
                                        st.altair_chart(
                                            _signature_profile_chart(
                                                _sig_spec,
                                                _sig_name,
                                                overlay_df=_overlay_spec,
                                                overlay_label=_overlay_label,
                                            ),
                                            width="stretch",
                                        )
                                        st.caption(_caption)
                        else:
                            _fixed_max = min(10, len(sig_df))
                            _fixed_default = min(4, _fixed_max)
                            _n_fixed = st.slider(
                                "Number of fixed COSMIC signatures",
                                1,
                                _fixed_max,
                                _fixed_default,
                                key="nmf_fixed_cosmic",
                                help="Use the top N COSMIC signatures from the current cohort-level fit as fixed basis signatures, then learn one additional non-negative residual signature.",
                            )
                            _fixed_signature_names = sig_df.head(_n_fixed)["signature"].tolist()
                            st.caption(
                                "Fixed COSMIC signatures: " + ", ".join(_fixed_signature_names)
                            )

                            try:
                                _guided_result = fit_cosmic_augmented_nmf(
                                    _sample_matrix,
                                    _cosmic_aligned,
                                    _fixed_signature_names,
                                )
                            except Exception as exc:
                                st.error(f"Could not run COSMIC-guided NMF discovery: {exc}")
                            else:
                                _guided_profiles = _guided_result["profiles"]
                                _guided_exposure_fractions = _guided_result["exposure_fractions"]
                                _guided_iter = _guided_result["n_iter"]
                                _guided_matrix_cosine = _guided_result["matrix_cosine"]
                                _guided_relative_error_pct = _guided_result["relative_error_pct"]
                                _guided_fixed_relative_error_pct = _guided_result[
                                    "fixed_only_relative_error_pct"
                                ]
                                _guided_improvement_pct = _guided_result[
                                    "relative_error_improvement_pct"
                                ]
                                _learned_names = list(_guided_result["learned_signature_names"])

                                _m1, _m2, _m3, _m4 = st.columns(4)
                                _m1.metric("Samples", f"{_n_samples_nmf:,}")
                                _m2.metric("Fixed COSMIC", str(_n_fixed))
                                _m3.metric("Matrix cosine", f"{_guided_matrix_cosine:.4f}")
                                _m4.metric("Residual improvement", f"{_guided_improvement_pct:.2f}%")
                                st.caption(
                                    f"Alternating constrained fit converged in {_guided_iter} iteration(s). "
                                    f"Fixed-only relative error was {_guided_fixed_relative_error_pct:.2f}%; "
                                    f"adding the learned signature reduced it to {_guided_relative_error_pct:.2f}%."
                                )

                                _guided_long = (
                                    _guided_exposure_fractions.reset_index(names="sample_label")
                                    .melt(id_vars="sample_label", var_name="signature", value_name="exposure")
                                )
                                _guided_order = _guided_profiles.index.tolist()
                                _sample_order = sorted(_guided_exposure_fractions.index.tolist())
                                _guided_chart = (
                                    alt.Chart(_guided_long)
                                    .mark_rect()
                                    .encode(
                                        alt.X("signature:N", sort=_guided_order, title="Signature component"),
                                        alt.Y("sample_label:N", sort=_sample_order, title="Sample"),
                                        alt.Color("exposure:Q", title="Exposure",
                                                  scale=alt.Scale(scheme="oranges"),
                                                  legend=alt.Legend(format=".0%")),
                                        tooltip=[
                                            "sample_label:N",
                                            "signature:N",
                                            alt.Tooltip("exposure:Q", format=".2%", title="Exposure"),
                                        ],
                                    )
                                    .properties(
                                        title="Per-sample COSMIC-guided signature exposures",
                                        height=max(150, 22 * len(_sample_order)),
                                    )
                                )
                                st.altair_chart(_guided_chart, width="stretch")

                                _guided_match_df = None
                                try:
                                    _guided_match_df, _ = compare_signatures_to_cosmic(
                                        _guided_profiles.loc[_learned_names],
                                        _cosmic_aligned,
                                    )
                                except Exception as exc:
                                    st.warning(f"Could not compare the learned signature to COSMIC: {exc}")
                                else:
                                    _guided_match_df["etiology"] = _guided_match_df["most_similar_cosmic_signature"].map(
                                        lambda s: _SBS_ETIOLOGY.get(s, "")
                                    )
                                    _display = _guided_match_df.copy()
                                    _display["most_similar_cosine_similarity"] = _display[
                                        "most_similar_cosine_similarity"
                                    ].map("{:.3f}".format)
                                    st.markdown("**Best COSMIC match for the learned residual signature**")
                                    st.dataframe(_display, width="stretch", hide_index=True)

                                _novel_sig_name = _learned_names[0]
                                _novel_match_row = None
                                if _guided_match_df is not None and not _guided_match_df.empty:
                                    _novel_match_row = (
                                        _guided_match_df[_guided_match_df["signature"] == _novel_sig_name]
                                        .iloc[0]
                                    )
                                _novel_download_df = build_signature_download_table(
                                    _profile_to_spec96(
                                        _guided_profiles.loc[_novel_sig_name],
                                        value_name="fraction",
                                    ),
                                    signature_name=_novel_sig_name,
                                    most_similar_cosmic_signature=(
                                        None
                                        if _novel_match_row is None
                                        else _novel_match_row["most_similar_cosmic_signature"]
                                    ),
                                    most_similar_cosine_similarity=(
                                        None
                                        if _novel_match_row is None
                                        else float(_novel_match_row["most_similar_cosine_similarity"])
                                    ),
                                    fixed_signature_names=_fixed_signature_names,
                                )
                                _novel_provenance_df = _build_active_filter_provenance(
                                    discovery_mode="COSMIC-guided + one learned signature",
                                    discovery_items=[
                                        ("cosmic_matrix_path", cosmic_path.strip()),
                                        ("fixed_cosmic_count", _n_fixed),
                                        ("fixed_cosmic_signatures", _fixed_signature_names),
                                        (
                                            "most_similar_cosmic_signature",
                                            None if _novel_match_row is None else _novel_match_row["most_similar_cosmic_signature"],
                                        ),
                                        (
                                            "most_similar_cosine_similarity",
                                            None
                                            if _novel_match_row is None
                                            else f"{float(_novel_match_row['most_similar_cosine_similarity']):.6f}",
                                        ),
                                    ],
                                )
                                _novel_exposures_df = build_signature_exposure_download_table(
                                    _guided_exposure_fractions
                                )
                                _novel_bundle = build_signature_download_zip(
                                    _novel_download_df,
                                    _novel_provenance_df,
                                    signature_name=_novel_sig_name,
                                    sample_exposures_df=_novel_exposures_df,
                                )
                                st.download_button(
                                    "Download novel signature bundle (.zip)",
                                    data=_novel_bundle,
                                    file_name=f"{_novel_sig_name.lower()}_signature_bundle.zip",
                                    mime="application/zip",
                                    key="nmf_download_novel_signature",
                                    help="Download a zip containing the learned SBS96 signature, the per-sample guided exposure table, and a provenance table of active filters and discovery settings.",
                                )

                                _sig_tabs = st.tabs(_guided_order)
                                for _tab, _sig_name in zip(_sig_tabs, _guided_order):
                                    with _tab:
                                        _sig_spec = _profile_to_spec96(
                                            _guided_profiles.loc[_sig_name],
                                            value_name="fraction",
                                        )
                                        if _sig_name in _fixed_signature_names:
                                            st.altair_chart(
                                                _signature_profile_chart(_sig_spec, _sig_name),
                                                width="stretch",
                                            )
                                            st.caption("Fixed COSMIC signature used in the constrained fit.")
                                        else:
                                            _overlay_spec = None
                                            _overlay_label = None
                                            _caption = "Learned non-negative residual signature."
                                            if _guided_match_df is not None and not _guided_match_df.empty:
                                                _match_row = (
                                                    _guided_match_df[_guided_match_df["signature"] == _sig_name]
                                                    .iloc[0]
                                                )
                                                _cosmic_sig = _match_row["most_similar_cosmic_signature"]
                                                _overlay_spec = _profile_to_spec96(
                                                    _cosmic_aligned[_cosmic_sig],
                                                    value_name="fraction",
                                                )
                                                _overlay_label = _cosmic_sig
                                                _caption = (
                                                    f"Best COSMIC match: {_cosmic_sig} "
                                                    f"(cosine {_match_row['most_similar_cosine_similarity']:.3f})."
                                                )
                                            st.altair_chart(
                                                _signature_profile_chart(
                                                    _sig_spec,
                                                    _sig_name,
                                                    overlay_df=_overlay_spec,
                                                    overlay_label=_overlay_label,
                                                ),
                                                width="stretch",
                                            )
                                            st.caption(_caption)

            # ── Called vs Uncalled Comparison ─────────────────────────────────
            if _has_data("variant_called"):
                st.divider()
                st.subheader("Called vs Uncalled Comparison")
                st.caption(
                    "Compares the SBS96 trinucleotide spectrum and COSMIC signature exposures "
                    "between loci where a variant was called and where it was not. "
                    "Requires the variant_called column (provide a VCF or variants TSV at collect time)."
                )

                def _build_spectrum(called_val):
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
                    _called_label   = f"Called (n={_n_called:,})"
                    _uncalled_label = f"Uncalled (n={_n_uncalled:,})"

                    def _make_spec96_df(obs_arr, n_total):
                        df_s = pd.DataFrame({
                            "sbs_label": _SBS_ORDER,
                            "mut_type":  [lbl[2:5] for lbl in _SBS_ORDER],
                            "count":     obs_arr.astype(int),
                        })
                        df_s["fraction"] = df_s["count"] / n_total if n_total > 0 else 0.0
                        return df_s

                    _spec_called   = _make_spec96_df(_obs_called,   _n_called)
                    _spec_uncalled = _make_spec96_df(_obs_uncalled, _n_uncalled)

                    # Mirrored (butterfly) chart
                    _m_df = pd.concat([
                        _spec_called.assign(y=_spec_called["fraction"],     group=_called_label),
                        _spec_uncalled.assign(y=-_spec_uncalled["fraction"], group=_uncalled_label),
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
                                          title=None,
                                          scale=alt.Scale(
                                              domain=[_called_label, _uncalled_label],
                                              range=["#4c78a8", "#e45756"],
                                          ),
                                          legend=alt.Legend(orient="bottom")),
                                tooltip=[
                                    alt.Tooltip("sbs_label:N", title="Context"),
                                    alt.Tooltip("group:N", title="Group"),
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
                            f"Called vs Uncalled — mirrored trinucleotide spectrum  ·  {_called_label} / {_uncalled_label}",
                            fontSize=13,
                        )),
                        width="stretch",
                    )
                    st.caption(
                        "Called variants point up (blue), uncalled loci point down (red). "
                        "Each group is normalised to its own fraction so differences in count don't dominate."
                    )

                    # COSMIC signature comparison (reuses matrix loaded in section above)
                    if _cosmic_W is not None:
                        st.markdown("**COSMIC signature comparison — Called vs Uncalled**")

                        def _fit_cmp(obs):
                            _h, _ = nnls(_cosmic_W, obs)
                            _total = _h.sum()
                            _h_norm = _h / _total if _total > 0 else _h
                            _recon  = _cosmic_W @ _h
                            _cos    = (
                                float(np.dot(obs, _recon))
                                / (np.linalg.norm(obs) * np.linalg.norm(_recon) + 1e-12)
                            )
                            return _h_norm, _cos

                        _h_called,   _cos_called   = _fit_cmp(_obs_called)
                        _h_uncalled, _cos_uncalled = _fit_cmp(_obs_uncalled)

                        _gof_c1, _gof_c2, _gof_c3, _gof_c4 = st.columns(4)
                        _gof_c1.metric("Called SNVs",           f"{_n_called:,}")
                        _gof_c2.metric("Cosine sim (called)",   f"{_cos_called:.4f}")
                        _gof_c3.metric("Uncalled SNVs",         f"{_n_uncalled:,}")
                        _gof_c4.metric("Cosine sim (uncalled)", f"{_cos_uncalled:.4f}")

                        _sig_names = _cosmic_aligned.columns.tolist()
                        _cmp_df = pd.DataFrame({
                            "signature": _sig_names * 2,
                            "group":     [_called_label]   * len(_sig_names) +
                                         [_uncalled_label] * len(_sig_names),
                            "exposure":  list(_h_called) + list(_h_uncalled),
                        })
                        _cmp_df["etiology"] = _cmp_df["signature"].map(
                            lambda s: _SBS_ETIOLOGY.get(s, "")
                        )

                        _cmp_top_n = st.slider(
                            "Top signatures to display (comparison)",
                            3, min(20, len(_sig_names)), 8,
                            key="cmp_top_n",
                        )
                        _max_exp = (
                            _cmp_df.groupby("signature")["exposure"].max()
                            .sort_values(ascending=False)
                        )
                        _top_sigs = _max_exp.head(_cmp_top_n).index.tolist()
                        _top_cmp_df = _cmp_df[_cmp_df["signature"].isin(_top_sigs)].copy()

                        _cmp_bars = (
                            alt.Chart(_top_cmp_df)
                            .mark_bar()
                            .encode(
                                alt.X("signature:N", sort=_top_sigs, title="Signature"),
                                alt.Y("exposure:Q", title="Exposure (proportion)",
                                      axis=alt.Axis(format=".0%")),
                                alt.Color("group:N",
                                          title=None,
                                          scale=alt.Scale(
                                              domain=[_called_label, _uncalled_label],
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
                        _cmp_dividers = (
                            alt.Chart(pd.DataFrame({"signature": _top_sigs[1:]}))
                            .mark_rule(color="#888", strokeWidth=1, opacity=0.5)
                            .encode(
                                alt.X("signature:N", sort=_top_sigs, bandPosition=0),
                            )
                        )
                        st.altair_chart(
                            alt.layer(_cmp_bars, _cmp_dividers)
                            .properties(
                                title=f"Top {_cmp_top_n} COSMIC SBS Signatures — Called vs Uncalled",
                                height=350,
                            ),
                            width="stretch",
                        )
                        st.caption(
                            "Blue = called variants; red = uncalled. "
                            "If called variants are enriched in known cancer signatures (e.g. SBS1, SBS5) "
                            "while uncalled are dominated by artefact signatures (e.g. SBS58), "
                            "this supports the quality of the variant calling."
                        )

                        with st.expander("Full signature table"):
                            _pivot = (
                                _cmp_df[_cmp_df["exposure"] > 0]
                                .pivot(index="signature", columns="group", values="exposure")
                                .fillna(0)
                                .reset_index()
                            )
                            for col in [_called_label, _uncalled_label]:
                                if col in _pivot.columns:
                                    _pivot[col] = _pivot[col].map("{:.2%}".format)
                            _pivot["etiology"] = _pivot["signature"].map(
                                lambda s: _SBS_ETIOLOGY.get(s, "")
                            )
                            st.dataframe(_pivot, width="stretch", hide_index=True)


            # ── VAF-stratified spectrum ────────────────────────────────────────
            st.divider()
            st.subheader("VAF-stratified Spectrum")
            st.caption(
                "Germline (VAF > 30%) vs somatic (VAF ≤ 30%) trinucleotide spectra. "
                "Differences between the two profiles reveal what drives low-VAF calls."
            )

            _vaf_germ_raw = con.execute(f"""
                SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
                FROM (SELECT * FROM {table_expr}) _t
                WHERE {where} AND variant_type = 'SNV'
                  AND trinuc_context IS NOT NULL AND length(trinuc_context) = 3
                  AND alt_count * 1.0 / total_depth > 0.3
                GROUP BY trinuc_context, ref_allele, alt_allele
            """).df()
            _vaf_som_raw = con.execute(f"""
                SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
                FROM (SELECT * FROM {table_expr}) _t
                WHERE {where} AND variant_type = 'SNV'
                  AND trinuc_context IS NOT NULL AND length(trinuc_context) = 3
                  AND alt_count * 1.0 / total_depth <= 0.3
                GROUP BY trinuc_context, ref_allele, alt_allele
            """).df()

            _vaf_germ_s96, _n_germ = _to_spec96_strat(_vaf_germ_raw)
            _vaf_som_s96,  _n_som  = _to_spec96_strat(_vaf_som_raw)

            _vc1, _vc2 = st.columns(2)
            with _vc1:
                if _vaf_germ_s96 is not None:
                    st.altair_chart(
                        _strat_sbs96_chart(_vaf_germ_s96, f"Germline VAF > 30% (n={_n_germ:,})"),
                        width="stretch",
                    )
                else:
                    st.info("No germline SNVs in current selection.")
            with _vc2:
                if _vaf_som_s96 is not None:
                    st.altair_chart(
                        _strat_sbs96_chart(_vaf_som_s96, f"Somatic VAF ≤ 30% (n={_n_som:,})"),
                        width="stretch",
                    )
                else:
                    st.info("No somatic SNVs in current selection.")

            # ── R1 / R2 stratified spectrum ────────────────────────────────────
            st.divider()
            st.subheader("R1 / R2 Stratified Spectrum")
            st.caption(
                "Trinucleotide spectra for Read 1 vs Read 2. "
                "Differences between the two profiles indicate read-level artefacts "
                "or strand-specific damage patterns (e.g. oxidative damage on R2)."
            )
            if not _has_alt_reads:
                st.info("R1/R2 spectrum requires the alt_reads table (DuckDB mode only).")
            else:
                _locus_snv = f"""
                    SELECT sample_id, chrom, pos, alt_allele,
                           trinuc_context, ref_allele
                    FROM {table_expr}
                    WHERE {where} AND variant_type = 'SNV'
                      AND trinuc_context IS NOT NULL AND length(trinuc_context) = 3
                """
                _r1_raw = con.execute(f"""
                    SELECT ab.trinuc_context, ab.ref_allele, ab.alt_allele,
                           COUNT(*) AS count
                    FROM (SELECT * FROM alt_reads {_r_reads_filter}) ar
                    INNER JOIN ({_locus_snv}) ab
                    ON  ar.sample_id  = ab.sample_id
                    AND ar.chrom      = ab.chrom
                    AND ar.pos        = ab.pos
                    AND ar.alt_allele = ab.alt_allele
                    WHERE ar.is_read1 = TRUE
                    GROUP BY ab.trinuc_context, ab.ref_allele, ab.alt_allele
                """).df()
                _r2_raw = con.execute(f"""
                    SELECT ab.trinuc_context, ab.ref_allele, ab.alt_allele,
                           COUNT(*) AS count
                    FROM (SELECT * FROM alt_reads {_r_reads_filter}) ar
                    INNER JOIN ({_locus_snv}) ab
                    ON  ar.sample_id  = ab.sample_id
                    AND ar.chrom      = ab.chrom
                    AND ar.pos        = ab.pos
                    AND ar.alt_allele = ab.alt_allele
                    WHERE ar.is_read1 = FALSE
                    GROUP BY ab.trinuc_context, ab.ref_allele, ab.alt_allele
                """).df()

                # Annotate raw dfs with sbs_label for drill-down
                for _rdf in [_r1_raw, _r2_raw]:
                    if not _rdf.empty:
                        _rdf["sbs_label"] = _rdf.apply(
                            lambda r: _sbs_label(r["trinuc_context"], r["ref_allele"], r["alt_allele"]),
                            axis=1,
                        )

                _r1_s96, _n_r1 = _to_spec96_strat(_r1_raw)
                _r2_s96, _n_r2 = _to_spec96_strat(_r2_raw)

                _r12_y_max = None
                if _r1_s96 is not None and _r2_s96 is not None:
                    _r12_y_max = max(_r1_s96["fraction"].max(), _r2_s96["fraction"].max())

                _ev_r1 = _ev_r2 = None
                _rc1, _rc2 = st.columns(2)
                with _rc1:
                    if _r1_s96 is not None:
                        _ev_r1 = st.altair_chart(
                            _strat_sbs96_chart(_r1_s96, f"Read 1 (n={_n_r1:,})", _r12_y_max, sel_name="r1_click"),
                            width="stretch",
                            on_select="rerun",
                            key="sbs96_r1",
                        )
                    else:
                        st.info("No R1 SNVs in current selection.")
                with _rc2:
                    if _r2_s96 is not None:
                        _ev_r2 = st.altair_chart(
                            _strat_sbs96_chart(_r2_s96, f"Read 2 (n={_n_r2:,})", _r12_y_max, sel_name="r2_click"),
                            width="stretch",
                            on_select="rerun",
                            key="sbs96_r2",
                        )
                    else:
                        st.info("No R2 SNVs in current selection.")

                st.caption("Click a bar to drill down to matching loci and open in IGV.")

                def _r12_drilldown(event, raw_df, read_label, sel_key_prefix):
                    if event is None:
                        return
                    pts = (event.selection or {}).get(f"{sel_key_prefix}_click", [])
                    if not pts:
                        return
                    clicked = [p.get("sbs_label") for p in pts if p.get("sbs_label")]
                    if not clicked or raw_df.empty:
                        return
                    matching = (
                        raw_df[raw_df["sbs_label"].isin(clicked)]
                        [["trinuc_context", "ref_allele", "alt_allele"]]
                        .drop_duplicates()
                    )
                    if matching.empty:
                        return
                    or_clauses = " OR ".join(
                        f"(trinuc_context = '{r.trinuc_context}' AND ref_allele = '{r.ref_allele}' AND alt_allele = '{r.alt_allele}')"
                        for r in matching.itertuples(index=False)
                    )
                    extra_cond = f"variant_type = 'SNV' AND ({or_clauses})"
                    sel = query_records([extra_cond])
                    label_str = ", ".join(clicked)
                    st.caption(f"**{read_label}** — {len(sel):,} loci · context(s): {label_str}")
                    st.dataframe(sel[_table_cols], width="stretch")
                    igv_buttons([extra_cond], sel, key=f"{sel_key_prefix}_sbs_{'_'.join(clicked)}")

                _r12_drilldown(_ev_r1, _r1_raw, "Read 1", "r1")
                _r12_drilldown(_ev_r2, _r2_raw, "Read 2", "r2")

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
            event = st.altair_chart(chart, width="stretch", on_select="rerun", key="snv_error_spectrum")

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
                    st.dataframe(sel[_table_cols], width="stretch")
                    igv_buttons([
                        "variant_type = 'SNV'",
                        f"ref_allele = '{ref}'",
                        f"alt_allele = '{alt_allele}'",
                    ], sel, key=f"spectrum_{sub}")

    # ── SBS96 heatmap (cohort / DuckDB only) ──────────────────────────────────
    if path.endswith(".duckdb"):
        st.divider()
        st.subheader("SBS96 Heatmap (samples × trinucleotide contexts)")
        if not _has_data("trinuc_context"):
            st.info("Trinucleotide context unavailable — run geac collect with a reference FASTA.")
        else:
            _hm_has_batch = _has_data("batch")
            _hm_id_sql    = "sample_id || ' / ' || batch" if _hm_has_batch else "sample_id"
            _hm_group_by  = "sample_id, batch" if _hm_has_batch else "sample_id"

            _hm_raw = con.execute(f"""
                SELECT {_hm_id_sql} AS sample_label,
                       trinuc_context, ref_allele, alt_allele, COUNT(*) AS n
                FROM {table_expr}
                WHERE {where} AND variant_type = 'SNV'
                  AND trinuc_context IS NOT NULL AND length(trinuc_context) = 3
                GROUP BY {_hm_group_by}, trinuc_context, ref_allele, alt_allele
            """).df()

            if _hm_raw.empty:
                st.info("No SNVs with trinucleotide context in current selection.")
            else:
                _hm_raw["sbs_label"] = _hm_raw.apply(
                    lambda row: _sbs_label(row["trinuc_context"], row["ref_allele"], row["alt_allele"]),
                    axis=1,
                )
                _hm_raw = _hm_raw.dropna(subset=["sbs_label"])
                _hm_agg = _hm_raw.groupby(["sample_label", "sbs_label"], as_index=False)["n"].sum()

                _totals = _hm_agg.groupby("sample_label")["n"].transform("sum")
                _hm_agg["fraction"] = _hm_agg["n"] / _totals

                _all_combos = pd.MultiIndex.from_product(
                    [_hm_agg["sample_label"].unique(), _SBS_ORDER],
                    names=["sample_label", "sbs_label"],
                )
                _hm_full = (
                    _hm_agg.set_index(["sample_label", "sbs_label"])
                    .reindex(_all_combos, fill_value=0)
                    .reset_index()
                )
                _n_hm_samples = _hm_agg["sample_label"].nunique()
                _n_hm_loci    = int(_hm_agg["n"].sum())

                _hm_full["mut_type"] = _hm_full["sbs_label"].str.extract(r'\[([A-Z]>[A-Z])\]')[0]

                _hm_chart = (
                    alt.Chart(_hm_full)
                    .mark_rect()
                    .encode(
                        alt.X("sbs_label:N", sort=_SBS_ORDER, title=None,
                              axis=alt.Axis(labels=False, ticks=False)),
                        alt.Y("sample_label:N", title="Sample"),
                        alt.Color("fraction:Q", title="Fraction of SNVs",
                                  scale=alt.Scale(scheme="blues")),
                        tooltip=[
                            alt.Tooltip("sample_label:N", title="Sample"),
                            alt.Tooltip("sbs_label:N",    title="Context"),
                            alt.Tooltip("n:Q",            title="Alt loci"),
                            alt.Tooltip("fraction:Q",     title="Fraction"),
                        ],
                    )
                    .properties(
                        height=max(200, 20 * _hm_full["sample_label"].nunique()),
                        title="Normalised SBS96 profile per sample (fraction of SNVs)",
                    )
                )

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
                    width="stretch",
                )
                st.caption(
                    f"{_n_hm_samples:,} samples · {_n_hm_loci:,} alt loci. "
                    "Color = fraction of that sample's SNVs falling in each trinucleotide context. "
                    "Contexts ordered by mutation type (C>A, C>G, C>T, T>A, T>C, T>G) then flanking bases."
                )

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
    if _has_data("label1"):
        _color_options.append("Label 1")
    if _has_data("label2"):
        _color_options.append("Label 2")
    if _has_data("label3"):
        _color_options.append("Label 3")
    if _has_data("on_target"):
        _color_options.append("On target")
    if _has_data("variant_called"):
        _color_options.append("Called variant")
    _color_by = _sb_col2.radio(
        "Color by", _color_options, horizontal=True, key="sb_color_by",
    )

    _sb_opt_cols = ", ".join(filter(None, [
        "batch"          if _has_data("batch")          else None,
        "label1"         if _has_data("label1")         else None,
        "label2"         if _has_data("label2")         else None,
        "label3"         if _has_data("label3")         else None,
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
        f"Strand Bias (log1p scale) — solid: perfect balance; dashed: 95% CI under Binomial(n, 0.5) — showing {_sb_n_pts:,} alt loci"
        if _use_log1p else
        f"Strand Bias — solid: perfect balance; dashed: 95% CI under Binomial(n, 0.5) — showing {_sb_n_pts:,} alt loci"
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
        "Label 1":        ("label1:N",         "Label 1",        alt.Scale()),
        "Label 2":        ("label2:N",         "Label 2",        alt.Scale()),
        "Label 3":        ("label3:N",         "Label 3",        alt.Scale()),
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
        width="stretch",
        on_select="rerun",
        key="strand_bias_scatter",
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
            st.dataframe(_sb_sel_df[_table_cols], width="stretch", hide_index=True)
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
                alt.Y("n:Q", title="Alt loci"),
                tooltip=[
                    alt.Tooltip("agree_frac:Q", format=".3f", title="Agreement fraction"),
                    alt.Tooltip("n:Q", title="Alt loci"),
                ],
            )
            .properties(title="Overlap Agreement Fraction", height=350)
        )
        st.altair_chart(chart, width="stretch")
        st.caption(
            "For alt loci where both reads of a pair overlap the variant site, the agreement "
            "fraction is the proportion of overlapping pairs where R1 and R2 both support the "
            "same alt base (agree / (agree + disagree)). A value of 1.0 means every overlapping "
            "pair agreed on the alt; lower values indicate that one read in the pair called ref "
            "or a different base. True somatic variants typically cluster near 1.0; strand "
            "artefacts — where only one read of a pair carries the alt — tend to show lower "
            "agreement or are absent from this plot entirely (they contribute to overlap_alt_disagree "
            "but not overlap_alt_agree)."
        )

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

        # When batch is present each row is (sample_id, batch); combine into a
        # single display label so every chart has one series per unique unit.
        _has_batch = _has_data("batch")
        _cohort_id_sql = (
            "sample_id || ' / ' || batch" if _has_batch else "sample_id"
        )
        _cohort_group_by = "sample_id, batch" if _has_batch else "sample_id"

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
                {_cohort_id_sql}                                     AS sample_label,
                sample_id,
                {'batch,' if _has_batch else ''}
                COUNT(*) FILTER (WHERE variant_type = 'SNV')        AS n_snv,
                COUNT(*) FILTER (WHERE variant_type = 'insertion')   AS n_insertion,
                COUNT(*) FILTER (WHERE variant_type = 'deletion')    AS n_deletion,
                ROUND(AVG(total_depth), 1)                           AS mean_depth,
                ROUND(AVG(alt_count * 1.0 / total_depth), 6)        AS mean_vaf,
                ROUND(
                    AVG(fwd_alt_count * 1.0
                        / NULLIF(fwd_alt_count + rev_alt_count, 0)),
                    4
                ) AS strand_balance,
                {_overlap_col}
                COUNT(*) AS n_total
            FROM {table_expr}
            WHERE {_cohort_where}
            GROUP BY {_cohort_group_by}
            ORDER BY {_cohort_group_by}
        """).df()

        if _cohort_stats.empty:
            st.warning("No records match the current filters.", icon="🔎")
        else:
            _cohort_event = st.dataframe(
                _cohort_stats,
                width="stretch",
                on_select="rerun",
                selection_mode="single-row",
                hide_index=True,
                key="cohort_data_table",
            )

            _cohort_sel = (_cohort_event.selection or {}).get("rows", [])
            if _cohort_sel:
                _focused_row    = _cohort_stats.iloc[_cohort_sel[0]]
                _focused_sample = _focused_row["sample_id"]
                _focused_label  = _focused_row["sample_label"]
                st.caption(
                    f"Focused on **{_focused_label}** — "
                    "click button below to filter all tabs to this sample."
                )
                if st.button(f"Filter all tabs to {_focused_label}"):
                    st.session_state["sample_sel"] = [_focused_sample]
                    st.rerun()

            # ── Step 2: Strand balance scatter ────────────────────────────────
            st.divider()
            st.subheader("Strand balance by sample")

            _strand_stats = con.execute(f"""
                SELECT
                    {_cohort_id_sql} AS sample_label,
                    ROUND(AVG(alt_count * 1.0 / total_depth), 6) AS mean_vaf,
                    ROUND(
                        AVG(fwd_alt_count * 1.0
                            / NULLIF(fwd_alt_count + rev_alt_count, 0)),
                        4
                    ) AS mean_strand_balance,
                    COUNT(*) AS n_loci
                FROM {table_expr}
                WHERE {_cohort_where} AND variant_type = 'SNV'
                GROUP BY {_cohort_group_by}
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
                        alt.Color("sample_label:N", title="Sample"),
                        alt.Size("n_loci:Q", title="SNV loci",
                                 scale=alt.Scale(range=[40, 300])),
                        tooltip=[
                            "sample_label:N",
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
                    width="stretch",
                )
                st.caption(
                    "Each dot is one sample. x = mean strand balance (0.5 = perfect), "
                    "y = mean VAF. Outliers in either axis may indicate a problematic sample."
                )

            # ── Step 4: Alt loci count vs mean base quality ───────────────────
            st.divider()
            st.subheader("Alt loci count vs mean base quality (per sample)")
            _bq_label_sql  = "ab.sample_id || ' / ' || ab.batch" if _has_batch else "ab.sample_id"
            _bq_group_sql  = "ab.sample_id, ab.batch"            if _has_batch else "ab.sample_id"
            _bq_batch_sel  = "batch,"                            if _has_batch else ""
            _bq_loci_df = con.execute(f"""
                SELECT
                    {_bq_label_sql} AS sample_label,
                    COUNT(DISTINCT CONCAT(ab.chrom, ':', ab.pos, ':', ab.alt_allele)) AS n_alt_loci,
                    ROUND(AVG(ar.base_qual), 2) AS mean_base_qual,
                    COUNT(ar.rowid) AS n_reads
                FROM (
                    SELECT sample_id, {_bq_batch_sel} chrom, pos, alt_allele
                    FROM {table_expr}
                    WHERE {_cohort_where}
                ) ab
                INNER JOIN alt_reads ar
                    ON  ab.sample_id  = ar.sample_id
                    AND ab.chrom      = ar.chrom
                    AND ab.pos        = ar.pos
                    AND ab.alt_allele = ar.alt_allele
                WHERE ar.base_qual IS NOT NULL
                GROUP BY {_bq_group_sql}
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
                        alt.Color("sample_label:N", title="Sample"),
                        alt.Size("n_reads:Q", title="Alt-supporting reads",
                                 scale=alt.Scale(range=[40, 300])),
                        tooltip=[
                            "sample_label:N",
                            alt.Tooltip("n_alt_loci:Q", format=",", title="Alt loci"),
                            alt.Tooltip("mean_base_qual:Q", format=".1f", title="Mean base qual"),
                            alt.Tooltip("n_reads:Q", format=",", title="Alt-supporting reads"),
                        ],
                    )
                    .properties(height=350, title="Alt loci count vs mean base quality (per sample)")
                )
                st.altair_chart(_bq_loci_chart, width="stretch")
                st.caption(
                    "Each dot is one sample. x = number of distinct alt loci, "
                    "y = mean base quality across all alt-supporting reads. "
                    "Samples with many loci but low base quality may be artefact-driven."
                )

            # ── Step 5: SNV count bar chart stacked by SBS6 substitution ──────
            st.subheader("SNV Count by Sample (SBS6 breakdown)")
            _sbs6_df = con.execute(f"""
                SELECT
                    {_cohort_id_sql} AS sample_label,
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
                GROUP BY {_cohort_group_by}, substitution
                ORDER BY {_cohort_group_by}, substitution
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
                        alt.X("sample_label:N", title="Sample", sort="-y",
                              axis=alt.Axis(labelAngle=-45, labelLimit=200)),
                        alt.Y("n_snv:Q", title="SNV count", stack="zero"),
                        alt.Color("substitution:N", title="Substitution",
                                  scale=_sbs6_color_scale),
                        alt.Tooltip(["sample_label:N", "substitution:N", "n_snv:Q"]),
                    )
                    .properties(height=350, title="SNV count per sample colored by SBS6 substitution type")
                )
                _sbs6_total = int(_sbs6_df["n_snv"].sum())
                st.altair_chart(_sbs6_chart, width="stretch")
                st.caption(
                    f"{_sbs6_total:,} SNVs across {_sbs6_df['sample_label'].nunique():,} samples. "
                    "Each bar shows the total SNV count for one sample, stacked by the six SBS "
                    "substitution classes (C>A, C>G, C>T, T>A, T>C, T>G). Samples are sorted by "
                    "total SNV count. A dominant C>T signal can indicate UV or FFPE damage; elevated "
                    "C>A is associated with oxidative damage or smoking; a relatively flat distribution "
                    "is typical of background noise."
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

        # ── Row 1: Read position bias ──────────────────────────────────────────
        st.subheader("Read position bias")
        _dfe_ctrl1, _dfe_ctrl2, _dfe_ctrl3 = st.columns([3, 2, 1])
        _dfe_color_options = ["All samples (aggregate)", "Sample"]
        if _has_data("batch"):
            _dfe_color_options.append("Batch")
        if _has_data("label1"):
            _dfe_color_options.append("Label 1")
        if _has_data("label2"):
            _dfe_color_options.append("Label 2")
        if _has_data("label3"):
            _dfe_color_options.append("Label 3")
        _dfe_color_by = _dfe_ctrl1.radio(
            "Color by", _dfe_color_options,
            horizontal=True, key="dfe_color_by",
        )
        _dfe_y_mode = _dfe_ctrl2.radio(
            "Y axis", ["Fraction", "Count"],
            horizontal=True, key="dfe_y_mode",
        )
        _dfe_by_read   = _dfe_ctrl3.checkbox("Show R1/R2", value=False, key="dfe_show_r1r2")
        _dfe_by_sample = _dfe_color_by == "Sample"
        _dfe_by_batch  = _dfe_color_by == "Batch"
        _dfe_by_label  = _dfe_color_by in ("Label 1", "Label 2", "Label 3")
        _dfe_lbl_col   = {"Label 1": "label1", "Label 2": "label2", "Label 3": "label3"}.get(_dfe_color_by)
        _dfe_normalize = _dfe_y_mode == "Fraction"
        _DFE_READ_EXPR = "CASE WHEN ar.is_read1 THEN 'R1' ELSE 'R2' END"
        _dfe_batch_src = (
            f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, batch FROM {table_expr}"
            f" WHERE batch IS NOT NULL) ab ON ar.sample_id = ab.sample_id"
        )
        _dfe_label_src = (
            f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, {_dfe_lbl_col} FROM {table_expr}"
            f" WHERE {_dfe_lbl_col} IS NOT NULL) _lbl ON ar.sample_id = _lbl.sample_id"
        ) if _dfe_by_label else _r_join
        if _dfe_by_batch and _dfe_by_read:
            _dfe_source      = _dfe_batch_src
            _dfe_select_expr = f"ab.batch || ' ' || {_DFE_READ_EXPR} AS label, "
            _dfe_group_expr  = "label, "
            _dfe_label_col   = "label"
        elif _dfe_by_label and _dfe_by_read:
            _dfe_source      = _dfe_label_src
            _dfe_select_expr = f"_lbl.{_dfe_lbl_col} || ' ' || {_DFE_READ_EXPR} AS label, "
            _dfe_group_expr  = "label, "
            _dfe_label_col   = "label"
        elif _dfe_by_sample and _dfe_by_read:
            _dfe_source      = _r_join
            _dfe_select_expr = f"ar.sample_id || ' ' || {_DFE_READ_EXPR} AS label, "
            _dfe_group_expr  = "label, "
            _dfe_label_col   = "label"
        elif _dfe_by_read:
            _dfe_source      = _r_join
            _dfe_select_expr = f"{_DFE_READ_EXPR} AS read, "
            _dfe_group_expr  = "read, "
            _dfe_label_col   = "read"
        elif _dfe_by_batch:
            _dfe_source      = _dfe_batch_src
            _dfe_select_expr = "ab.batch AS batch, "
            _dfe_group_expr  = "ab.batch, "
            _dfe_label_col   = "batch"
        elif _dfe_by_label:
            _dfe_source      = _dfe_label_src
            _dfe_select_expr = f"_lbl.{_dfe_lbl_col} AS {_dfe_lbl_col}, "
            _dfe_group_expr  = f"_lbl.{_dfe_lbl_col}, "
            _dfe_label_col   = _dfe_lbl_col
        elif _dfe_by_sample:
            _dfe_source      = _r_join
            _dfe_select_expr = "ar.sample_id, "
            _dfe_group_expr  = "ar.sample_id, "
            _dfe_label_col   = "sample_id"
        else:
            _dfe_source      = _r_join
            _dfe_select_expr = ""
            _dfe_group_expr  = ""
            _dfe_label_col   = None

        _dfe_df = con.execute(f"""
            SELECT {_dfe_select_expr}ar.cycle, COUNT(*) AS n_reads
            FROM {_dfe_source}
            GROUP BY {_dfe_group_expr}ar.cycle
            ORDER BY {_dfe_group_expr}ar.cycle
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
                x=alt.X("cycle:Q", title="Cycle (1-based)", bin=False),
                y=alt.Y(_dfe_y_field, title=_dfe_y_title),
                tooltip=[
                    *([f"{_dfe_label_col}:N"] if _dfe_label_col else []),
                    alt.Tooltip("cycle:Q", title="Cycle"),
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
            st.altair_chart(_dfe_chart, width="stretch")
            st.caption(
                "A spike at high cycle numbers indicates alt-supporting reads clustered at read ends — "
                "a red flag for alignment artefacts or damaged bases."
            )

        # ── Row 2: Base qual vs dist from read end scatter ────────────────────
        st.subheader("Mean base quality by cycle")
        _bq_ctrl1, _bq_ctrl2 = st.columns([4, 1])
        _bq_color_options = ["All samples (aggregate)", "Sample"]
        if _has_data("batch"):
            _bq_color_options.append("Batch")
        if _has_data("label1"):
            _bq_color_options.append("Label 1")
        if _has_data("label2"):
            _bq_color_options.append("Label 2")
        if _has_data("label3"):
            _bq_color_options.append("Label 3")
        _bq_color_by = _bq_ctrl1.radio(
            "Color by", _bq_color_options,
            horizontal=True, key="bq_color_by",
        )
        _bq_by_read   = _bq_ctrl2.checkbox("Show R1/R2", value=False, key="bq_show_r1r2")
        _bq_by_sample = _bq_color_by == "Sample"
        _bq_by_batch  = _bq_color_by == "Batch"
        _bq_by_label  = _bq_color_by in ("Label 1", "Label 2", "Label 3")
        _bq_lbl_col   = {"Label 1": "label1", "Label 2": "label2", "Label 3": "label3"}.get(_bq_color_by)
        _BQ_READ_EXPR = "CASE WHEN ar.is_read1 THEN 'R1' ELSE 'R2' END"
        _bq_batch_src = (
            f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, batch FROM {table_expr}"
            f" WHERE batch IS NOT NULL) ab ON ar.sample_id = ab.sample_id"
        )
        _bq_label_src = (
            f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, {_bq_lbl_col} FROM {table_expr}"
            f" WHERE {_bq_lbl_col} IS NOT NULL) _lbl ON ar.sample_id = _lbl.sample_id"
        ) if _bq_by_label else _r_join
        if _bq_by_batch and _bq_by_read:
            _bq_source      = _bq_batch_src
            _bq_select_expr = f"ab.batch || ' ' || {_BQ_READ_EXPR} AS label, "
            _bq_group_expr  = "label, "
            _bq_label_col   = "label"
        elif _bq_by_label and _bq_by_read:
            _bq_source      = _bq_label_src
            _bq_select_expr = f"_lbl.{_bq_lbl_col} || ' ' || {_BQ_READ_EXPR} AS label, "
            _bq_group_expr  = "label, "
            _bq_label_col   = "label"
        elif _bq_by_sample and _bq_by_read:
            _bq_source      = _r_join
            _bq_select_expr = f"ar.sample_id || ' ' || {_BQ_READ_EXPR} AS label, "
            _bq_group_expr  = "label, "
            _bq_label_col   = "label"
        elif _bq_by_read:
            _bq_source      = _r_join
            _bq_select_expr = f"{_BQ_READ_EXPR} AS read, "
            _bq_group_expr  = "read, "
            _bq_label_col   = "read"
        elif _bq_by_batch:
            _bq_source      = _bq_batch_src
            _bq_select_expr = "ab.batch AS batch, "
            _bq_group_expr  = "ab.batch, "
            _bq_label_col   = "batch"
        elif _bq_by_label:
            _bq_source      = _bq_label_src
            _bq_select_expr = f"_lbl.{_bq_lbl_col} AS {_bq_lbl_col}, "
            _bq_group_expr  = f"_lbl.{_bq_lbl_col}, "
            _bq_label_col   = _bq_lbl_col
        elif _bq_by_sample:
            _bq_source      = _r_join
            _bq_select_expr = "ar.sample_id, "
            _bq_group_expr  = "ar.sample_id, "
            _bq_label_col   = "sample_id"
        else:
            _bq_source      = _r_join
            _bq_select_expr = ""
            _bq_group_expr  = ""
            _bq_label_col   = None

        _bq_df = con.execute(f"""
            SELECT
                {_bq_select_expr}ar.cycle,
                ROUND(AVG(ar.base_qual), 2) AS mean_base_qual,
                COUNT(*) AS n_reads
            FROM {_bq_source}
            GROUP BY {_bq_group_expr}ar.cycle
            ORDER BY {_bq_group_expr}ar.cycle
        """).df()

        if _bq_df.empty:
            st.info("No data.")
        else:
            _bq_enc = dict(
                x=alt.X("cycle:Q", title="Cycle (1-based)"),
                y=alt.Y("mean_base_qual:Q", title="Mean base quality (Phred)",
                        scale=alt.Scale(zero=False)),
                tooltip=[
                    *([f"{_bq_label_col}:N"] if _bq_label_col else []),
                    alt.Tooltip("cycle:Q", title="Cycle"),
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
            st.altair_chart(_bq_chart, width="stretch")
            st.caption(
                "A drop in mean base quality at high cycle numbers (late in the read) "
                "indicates that alt-supporting reads at those positions may be artefacts."
            )

        # ── Row 3: N context around alt-supporting reads ──────────────────────
        st.subheader("N context around alt-supporting reads")
        if not _has_alt_reads_cols(
            "n_before_alt",
            "n_after_alt",
            "n_n_before_alt",
            "n_n_after_alt",
            "leading_n_run_len",
            "trailing_n_run_len",
        ):
            st.info(
                "N-context analysis requires a newer `alt_reads` schema. "
                "Re-run `geac collect --reads-output` and `geac merge` to enable these plots."
            )
        else:
            _nctx_ctrl1, _nctx_ctrl2 = st.columns([4, 1])
            _nctx_color_mode = _nctx_ctrl1.radio(
                "Group by",
                ["All reads", "R1/R2"],
                horizontal=True,
                key="nctx_group_by",
            )
            _nctx_norm_mode = _nctx_ctrl2.radio(
                "Y axis",
                ["Fraction", "Count"],
                horizontal=True,
                key="nctx_y_mode",
            )
            _nctx_df = con.execute(f"""
                SELECT
                    ar.is_read1,
                    ar.n_before_alt,
                    ar.n_after_alt,
                    ar.n_n_before_alt,
                    ar.n_n_after_alt,
                    ar.leading_n_run_len,
                    ar.trailing_n_run_len
                FROM {_r_join}
            """).df()

            if _nctx_df.empty:
                st.info("No read-context data available under current filters.")
            else:
                _nctx_df = add_read_context_fraction_metrics(_nctx_df)
                _nctx_df["read_group"] = np.where(_nctx_df["is_read1"], "R1", "R2")
                _nctx_group_col = "read_group" if _nctx_color_mode == "R1/R2" else None

                def _nctx_dist_chart(
                    df: pd.DataFrame,
                    value_col: str,
                    title: str,
                    x_title: str,
                    *,
                    domain=None,
                    integer_bins: bool = False,
                    fraction_format: bool = False,
                ):
                    chart_df = df.dropna(subset=[value_col]).copy()
                    if chart_df.empty:
                        return None
                    y_enc = (
                        alt.Y("proportion:Q", title="Fraction of alt-supporting reads")
                        if _nctx_norm_mode == "Fraction"
                        else alt.Y("count():Q", title="Alt-supporting reads")
                    )
                    if _nctx_norm_mode == "Fraction":
                        if _nctx_group_col:
                            chart_df["proportion"] = chart_df.groupby(_nctx_group_col)[value_col].transform(
                                lambda x: 1.0 / len(x)
                            )
                        else:
                            chart_df["proportion"] = 1.0 / len(chart_df)
                    bin_def = alt.Bin(step=1) if integer_bins else alt.Bin(maxbins=40)
                    enc = dict(
                        x=alt.X(
                            f"{value_col}:Q",
                            title=x_title,
                            bin=bin_def,
                            scale=alt.Scale(domain=domain) if domain is not None else alt.Scale(),
                        ),
                        y=y_enc,
                        tooltip=[
                            *([f"{_nctx_group_col}:N"] if _nctx_group_col else []),
                            alt.Tooltip("count():Q", title="Count"),
                        ],
                    )
                    if _nctx_group_col:
                        enc["color"] = alt.Color(f"{_nctx_group_col}:N", title="Read group")
                    return (
                        alt.Chart(chart_df)
                        .mark_bar(opacity=0.7)
                        .encode(**enc)
                        .properties(title=title, height=260)
                    )

                _nctx_row1 = st.columns(2)
                with _nctx_row1[0]:
                    _trailing_chart = _nctx_dist_chart(
                        _nctx_df,
                        "trailing_n_run_len",
                        "Trailing N run length",
                        "Contiguous N run immediately after alt",
                        integer_bins=True,
                    )
                    if _trailing_chart is not None:
                        st.altair_chart(_trailing_chart, width="stretch")
                with _nctx_row1[1]:
                    _frac_after_df = _nctx_df[_nctx_df["n_after_alt"] > 0].copy()
                    _frac_after_chart = _nctx_dist_chart(
                        _frac_after_df,
                        "frac_n_after_alt",
                        "Fraction N after alt",
                        "Fraction of bases after alt that are N",
                        domain=[0, 1],
                        fraction_format=True,
                    )
                    if _frac_after_chart is not None:
                        st.altair_chart(_frac_after_chart, width="stretch")

                _nctx_compare_rows = []
                for _side, _col in [
                    ("Before alt", "frac_n_before_alt"),
                    ("After alt", "frac_n_after_alt"),
                ]:
                    _sub = _nctx_df.dropna(subset=[_col]).copy()
                    if not _sub.empty:
                        _sub["side"] = _side
                        _sub["frac_n"] = _sub[_col]
                        _nctx_compare_rows.append(_sub[["frac_n", "side"]])
                _nctx_row2 = st.columns(2)
                with _nctx_row2[0]:
                    if _nctx_compare_rows:
                        _nctx_compare_df = pd.concat(_nctx_compare_rows, ignore_index=True)
                        _nctx_compare_chart = (
                            alt.Chart(_nctx_compare_df)
                            .mark_boxplot(extent="min-max")
                            .encode(
                                x=alt.X("side:N", title=None),
                                y=alt.Y(
                                    "frac_n:Q",
                                    title="Fraction N",
                                    scale=alt.Scale(domain=[0, 1]),
                                ),
                                color=alt.Color("side:N", title=None),
                                tooltip=[
                                    "side:N",
                                    alt.Tooltip("median(frac_n):Q", title="Median", format=".3f"),
                                ],
                            )
                            .properties(title="Leading vs trailing N burden", height=260)
                        )
                        st.altair_chart(_nctx_compare_chart, width="stretch")
                with _nctx_row2[1]:
                    _delta_chart = _nctx_dist_chart(
                        _nctx_df.dropna(subset=["delta_n_fraction"]),
                        "delta_n_fraction",
                        "Before vs after asymmetry",
                        "delta_n_fraction = frac_N_after - frac_N_before",
                    )
                    if _delta_chart is not None:
                        st.altair_chart(_delta_chart, width="stretch")

                _nctx_trailing_frac = (_nctx_df["trailing_n_run_len"] > 0).mean()
                _nctx_median_delta = _nctx_df["delta_n_fraction"].median()
                st.caption(
                    f"{_nctx_trailing_frac:.1%} of alt-supporting reads have trailing N run length > 0. "
                    f"Median delta_n_fraction = {_nctx_median_delta:.3f} "
                    "(`frac_N_after - frac_N_before`). Positive values indicate greater N burden after the alt."
                )

        # ── Row 4: Insert size distribution ───────────────────────────────────
        if con.execute("SELECT COUNT(*) FROM alt_reads WHERE insert_size IS NOT NULL LIMIT 1").fetchone()[0] > 0:
            # Median read length — used for gap correction in both insert size plots.
            # Fragments longer than 2×read_len have a coverage gap; the correction
            # weights each count by 1 / min(1, 2R/L) to recover the unbiased distribution.
            _read_len_median = con.execute(
                "SELECT MEDIAN(read_length) FROM alt_reads WHERE read_length IS NOT NULL"
            ).fetchone()[0] or 0
            _gap_threshold = 2 * _read_len_median  # bp where gap effect begins

            if _gap_threshold > 0:
                st.info(
                    f"**Gap correction**: for paired-end reads of length {int(_read_len_median)} bp, "
                    f"fragments longer than {int(_gap_threshold)} bp have a gap between R1 and R2 where neither "
                    "read covers the position. Longer fragments are therefore underrepresented in the raw count, "
                    "producing a kink at that threshold. **Frequency** mode divides each bin by its capture "
                    "probability \u2014 min(1, 2R\u2009/\u2009L) \u2014 to recover the unbiased fragment size distribution. "
                    "Switch to **Count** to see the raw read counts."
                )

            st.subheader("Insert size distribution")
            _ins_color_options = ["All samples (aggregate)", "Sample"]
            if _has_data("batch"):
                _ins_color_options.append("Batch")
            if _has_data("label1"):
                _ins_color_options.append("Label 1")
            if _has_data("label2"):
                _ins_color_options.append("Label 2")
            if _has_data("label3"):
                _ins_color_options.append("Label 3")
            _ins_color_by = st.radio(
                "Color by", _ins_color_options,
                horizontal=True, key="ins_color_by",
            )
            _ins_by_sample = _ins_color_by == "Sample"
            _ins_by_batch  = _ins_color_by == "Batch"
            _ins_by_label  = _ins_color_by in ("Label 1", "Label 2", "Label 3")
            _ins_lbl_col   = {"Label 1": "label1", "Label 2": "label2", "Label 3": "label3"}.get(_ins_color_by)
            _ins_label_col = (
                "sample_id"   if _ins_by_sample else
                "batch"       if _ins_by_batch  else
                _ins_lbl_col  if _ins_by_label  else
                None
            )
            _ins_source = (
                f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, batch FROM {table_expr} WHERE batch IS NOT NULL) ab ON ar.sample_id = ab.sample_id"
                if _ins_by_batch else
                f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, {_ins_lbl_col} FROM {table_expr} WHERE {_ins_lbl_col} IS NOT NULL) _lbl ON ar.sample_id = _lbl.sample_id"
                if _ins_by_label else _r_join
            )
            _ins_select_expr = (
                "ab.batch AS batch, "                        if _ins_by_batch  else
                f"_lbl.{_ins_lbl_col} AS {_ins_lbl_col}, "  if _ins_by_label  else
                f"ar.{_ins_label_col}, "                     if _ins_label_col else ""
            )
            _ins_group_expr = (
                "ab.batch, "                  if _ins_by_batch  else
                f"_lbl.{_ins_lbl_col}, "      if _ins_by_label  else
                f"ar.{_ins_label_col}, "      if _ins_label_col else ""
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
                _ins_use_corrected = _ins_y_mode == "Frequency"

                # Apply gap correction: weight = min(1, 2R/L); corrected = raw / weight
                if _ins_use_corrected and _read_len_median > 0:
                    _ins_df["n_eff"] = _ins_df["n_reads"] / _ins_df["insert_size"].apply(
                        lambda x: min(1.0, 2.0 * _read_len_median / x) if x > 0 else 1.0
                    )
                else:
                    _ins_df["n_eff"] = _ins_df["n_reads"].astype(float)

                if _ins_use_corrected:
                    if _ins_label_col:
                        _ins_totals = _ins_df.groupby(_ins_label_col)["n_eff"].transform("sum")
                    else:
                        _ins_totals = _ins_df["n_eff"].sum()
                    _ins_df["frequency"] = _ins_df["n_eff"] / _ins_totals

                _ins_y_field = "frequency" if _ins_use_corrected else "n_reads"
                _ins_y_title = "Frequency" if _ins_use_corrected else "Alt-supporting reads"
                _ins_enc = dict(
                    x=alt.X("insert_size:Q", title="Insert size (bp)"),
                    y=alt.Y(f"{_ins_y_field}:Q", title=_ins_y_title,
                            **({"axis": alt.Axis(format=".3f")} if _ins_use_corrected else {})),
                    tooltip=[
                        *([f"{_ins_label_col}:N"] if _ins_label_col else []),
                        alt.Tooltip("insert_size:Q", title="Insert size (bp)"),
                        alt.Tooltip(f"{_ins_y_field}:Q",
                                    title=_ins_y_title,
                                    **({"format": ".4f"} if _ins_use_corrected else {})),
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
                st.altair_chart(_ins_chart, width="stretch")
                _ins_caption = (
                    "Insert size distribution of alt-supporting reads. "
                    "A shift toward shorter inserts can indicate adapter contamination or artefacts."
                )
                if _gap_threshold > 0:
                    _ins_caption += (
                        f" A kink near {int(_gap_threshold)} bp (2× read length) is expected: "
                        "fragments longer than this have a coverage gap between R1 and R2, "
                        "so not every fragment is captured at every locus. "
                        "Frequency mode divides each bin by its capture probability "
                        "min(1, 2R/L), recovering the unbiased fragment size distribution."
                    )
                st.caption(_ins_caption)

        # ── Row 3b: Insert size by AF class (germline vs somatic) ─────────────
        if con.execute("SELECT COUNT(*) FROM alt_reads WHERE insert_size IS NOT NULL LIMIT 1").fetchone()[0] > 0:
            st.subheader("Insert size by allele frequency class")
            _af_ins_color_options = ["All samples (aggregate)", "Sample"]
            if _has_data("batch"):
                _af_ins_color_options.append("Batch")
            if _has_data("label1"):
                _af_ins_color_options.append("Label 1")
            if _has_data("label2"):
                _af_ins_color_options.append("Label 2")
            if _has_data("label3"):
                _af_ins_color_options.append("Label 3")
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
            _af_ins_by_label  = _af_ins_color_by in ("Label 1", "Label 2", "Label 3")
            _af_ins_lbl_col   = {"Label 1": "label1", "Label 2": "label2", "Label 3": "label3"}.get(_af_ins_color_by)
            _af_ins_group_col = (
                "sample_id"       if _af_ins_by_sample else
                "batch"           if _af_ins_by_batch  else
                _af_ins_lbl_col   if _af_ins_by_label  else
                None
            )
            _af_ins_extra_select = (
                "ar.sample_id, "                           if _af_ins_by_sample else
                "_locus.batch, "                           if _af_ins_by_batch  else
                f"_locus.{_af_ins_lbl_col}, "              if _af_ins_by_label  else
                ""
            )
            _af_ins_extra_group = (
                "ar.sample_id, "                           if _af_ins_by_sample else
                "_locus.batch, "                           if _af_ins_by_batch  else
                f"_locus.{_af_ins_lbl_col}, "              if _af_ins_by_label  else
                ""
            )
            _af_ins_locus_extra = (
                ", batch"                   if _af_ins_by_batch else
                f", {_af_ins_lbl_col}"      if _af_ins_by_label else
                ""
            )

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
                _af_ins_use_corrected = _af_ins_y_mode == "Frequency"

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

                # Apply gap correction
                if _af_ins_use_corrected and _read_len_median > 0:
                    _af_ins_df["n_eff"] = _af_ins_df["n_reads"] / _af_ins_df["insert_size"].apply(
                        lambda x: min(1.0, 2.0 * _read_len_median / x) if x > 0 else 1.0
                    )
                else:
                    _af_ins_df["n_eff"] = _af_ins_df["n_reads"].astype(float)

                if _af_ins_use_corrected:
                    _norm_key = "series" if _af_ins_group_col else "af_class"
                    _af_totals = _af_ins_df.groupby(_norm_key)["n_eff"].transform("sum")
                    _af_ins_df["frequency"] = _af_ins_df["n_eff"] / _af_totals

                _af_ins_y_title = "Frequency" if _af_ins_use_corrected else "Alt-supporting reads"
                _af_ins_y_field = "frequency" if _af_ins_use_corrected else "n_reads"
                _af_ins_y = (
                    alt.Y(f"{_af_ins_y_field}:Q", title=_af_ins_y_title,
                          axis=alt.Axis(format=".3f"))
                    if _af_ins_use_corrected else
                    alt.Y("n_reads:Q", title="Alt-supporting reads")
                )
                _af_ins_tooltip = [
                    *([f"{_af_ins_group_col}:N"] if _af_ins_group_col else []),
                    alt.Tooltip("insert_size:Q", title="Insert size (bp)"),
                    alt.Tooltip("af_class:N", title="AF class"),
                    alt.Tooltip(_af_ins_y_field + ":Q",
                                title=_af_ins_y_title,
                                **({"format": ".4f"} if _af_ins_use_corrected else {})),
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
                st.altair_chart(_af_ins_chart, width="stretch")
                _af_ins_caption = (
                    "Insert size distributions split by allele frequency. "
                    "Each series is normalised independently so lines are directly comparable "
                    "regardless of how many reads fall in each group. "
                    "A shift in one class toward very short or very long inserts suggests artefacts in that group."
                )
                if _gap_threshold > 0:
                    _af_ins_caption += (
                        f" Select 'Frequency' to remove the coverage-gap bias near {int(_gap_threshold)} bp."
                    )
                st.caption(_af_ins_caption)

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
            st.altair_chart(_mq_chart, width="stretch")
            st.caption(
                "Stacked by locus type (repetitive = homopolymer ≥ 5 or STR length ≥ 6). "
                "Low MAPQ at repetitive loci indicates multi-mapping artefacts."
            )

with tab_duplex:
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

        # ── Family size distribution ───────────────────────────────────────────
        st.subheader("Family size distribution")
        _fs_ctrl_col1, _fs_ctrl_col2 = st.columns(2)
        _fs_color_options = ["All samples (aggregate)", "Sample"]
        if _has_data("batch"):
            _fs_color_options.append("Batch")
        if _has_data("label1"):
            _fs_color_options.append("Label 1")
        if _has_data("label2"):
            _fs_color_options.append("Label 2")
        if _has_data("label3"):
            _fs_color_options.append("Label 3")
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
        _fs_by_label  = _fs_color_by in ("Label 1", "Label 2", "Label 3")
        _fs_lbl_col   = {"Label 1": "label1", "Label 2": "label2", "Label 3": "label3"}.get(_fs_color_by)
        _fs_normalize = _fs_y_mode == "Fraction"

        _fs_group_col = (
            "ar.sample_id"          if _fs_by_sample else
            "ab.batch"              if _fs_by_batch  else
            f"_lbl.{_fs_lbl_col}"   if _fs_by_label  else
            None
        )
        _fs_select = f"{_fs_group_col}, " if _fs_group_col else ""
        _fs_group  = f"{_fs_group_col}, " if _fs_group_col else ""

        _fs_source = (
            f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, batch FROM {table_expr} WHERE batch IS NOT NULL) ab ON ar.sample_id = ab.sample_id"
            if _fs_by_batch else
            f"{_r_join} INNER JOIN (SELECT DISTINCT sample_id, {_fs_lbl_col} FROM {table_expr} WHERE {_fs_lbl_col} IS NOT NULL) _lbl ON ar.sample_id = _lbl.sample_id"
            if _fs_by_label else _r_join
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
            _fs_x_min = int(_fs_df["family_size"].min())
            _fs_x_max = int(_fs_df["family_size"].max())
            if _fs_x_min < _fs_x_max:
                _fs_x_range = st.slider(
                    "X-axis range (family size)",
                    min_value=_fs_x_min,
                    max_value=_fs_x_max,
                    value=(_fs_x_min, _fs_x_max),
                    key="fs_x_range",
                )
                _fs_df = _fs_df[
                    (_fs_df["family_size"] >= _fs_x_range[0])
                    & (_fs_df["family_size"] <= _fs_x_range[1])
                ]

            _fs_label_col = (
                "sample_id"  if _fs_by_sample else
                "batch"      if _fs_by_batch  else
                _fs_lbl_col  if _fs_by_label  else
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
            st.altair_chart(_fs_chart, width="stretch")
            _fs_norm_note = (
                "Fraction mode normalizes each batch independently."
                if _fs_by_batch else
                "Fraction mode normalizes each sample independently, allowing shape comparison across samples with different read counts."
            )
            st.caption(f"Artefacts are enriched in singletons (family_size = 1). {_fs_norm_note}")

        # ── Family size vs VAF (per locus) ────────────────────────────────────
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
            FROM (SELECT * FROM {table_expr} WHERE {where}) ab
            INNER JOIN alt_reads ar
                ON  ab.sample_id  = ar.sample_id
                AND ab.chrom      = ar.chrom
                AND ab.pos        = ar.pos
                AND ab.alt_allele = ar.alt_allele
            WHERE ar.family_size IS NOT NULL
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
            st.altair_chart(_fsvaf_chart, width="stretch")
            st.caption(
                "True low-VAF variants should have reasonable mean family sizes. "
                "Artefacts at low VAF tend to cluster at low family size."
            )

        # ── Family-size stratified spectrum ───────────────────────────────────
        if _has_alt_reads and _fs_has_data:
            st.divider()
            st.subheader("Family-size stratified Spectrum")
            st.caption(
                "Singleton reads (family_size = 1) vs multi-member families (family_size > 1). "
                "Singletons are enriched for sequencing errors; differences in profile shape "
                "reveal the true variant signal from the error process."
            )

            _locus_fs_filter = f"AND {_reads_where}" if _reads_active else ""
            _fs_strat_raw = con.execute(f"""
                WITH locus_fs AS (
                    SELECT sample_id, chrom, pos, alt_allele,
                           MEDIAN(family_size) AS median_fs
                    FROM alt_reads
                    WHERE family_size IS NOT NULL {_locus_fs_filter}
                    GROUP BY sample_id, chrom, pos, alt_allele
                )
                SELECT
                    _t.trinuc_context, _t.ref_allele, _t.alt_allele,
                    CASE WHEN COALESCE(lfs.median_fs, 1) <= 1
                         THEN 'singleton' ELSE 'multi' END AS fs_group,
                    COUNT(*) AS count
                FROM (SELECT * FROM {table_expr} WHERE {where}) _t
                LEFT JOIN locus_fs lfs
                    ON  lfs.sample_id  = _t.sample_id
                    AND lfs.chrom      = _t.chrom
                    AND lfs.pos        = _t.pos
                    AND lfs.alt_allele = _t.alt_allele
                WHERE _t.variant_type = 'SNV'
                  AND _t.trinuc_context IS NOT NULL
                  AND length(_t.trinuc_context) = 3
                GROUP BY _t.trinuc_context, _t.ref_allele, _t.alt_allele, fs_group
            """).df()

            _fs_sing_raw  = _fs_strat_raw[_fs_strat_raw["fs_group"] == "singleton"].drop(columns="fs_group")
            _fs_multi_raw = _fs_strat_raw[_fs_strat_raw["fs_group"] == "multi"].drop(columns="fs_group")

            _fs_sing_s96,  _n_sing  = _to_spec96_strat(_fs_sing_raw)
            _fs_multi_s96, _n_multi = _to_spec96_strat(_fs_multi_raw)

            _fc1, _fc2 = st.columns(2)
            with _fc1:
                if _fs_sing_s96 is not None:
                    st.altair_chart(
                        _strat_sbs96_chart(_fs_sing_s96, f"Singleton (family_size = 1, n={_n_sing:,})"),
                        width="stretch",
                    )
                else:
                    st.info("No singleton loci in current selection.")
            with _fc2:
                if _fs_multi_s96 is not None:
                    st.altair_chart(
                        _strat_sbs96_chart(_fs_multi_s96, f"Multi-member (family_size > 1, n={_n_multi:,})"),
                        width="stretch",
                    )
                else:
                    st.info("No multi-member loci in current selection.")

        # ── Cohort artefact vs rare variant: family size comparison ───────────
        if path.endswith(".duckdb"):
            st.divider()
            st.subheader("Cohort artefact vs rare variant: family size comparison")
            _n_samples_total = con.execute(
                f"SELECT COUNT(DISTINCT sample_id) FROM {table_expr} WHERE {where}"
            ).fetchone()[0]

            if _n_samples_total < 2:
                st.info("Need at least 2 samples for cohort artefact comparison.")
            else:
                _cohort_fs_df = con.execute(f"""
                    WITH _base AS (
                        SELECT sample_id,
                               chrom,
                               CAST(pos AS BIGINT) AS pos,
                               alt_allele
                        FROM {table_expr}
                        WHERE {where}
                    ),
                    locus_counts AS (
                        SELECT chrom, pos, alt_allele,
                               COUNT(DISTINCT sample_id) AS n_samples_with_alt
                        FROM _base
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
                        INNER JOIN _base _filt
                            ON  ar.sample_id  = _filt.sample_id
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
                    _cf_median = (
                        alt.Chart(_cohort_fs_df)
                        .mark_tick(color="white", thickness=2, size=40)
                        .encode(
                            alt.X("cohort_freq:N", sort=_cf_order),
                            alt.Y("median:Q"),
                        )
                    )
                    _cf_chart = (_cf_whisker + _cf_box + _cf_median).properties(height=300)
                    st.altair_chart(_cf_chart, width="stretch")
                    st.caption(
                        "Cohort artefacts (seen in many samples) tend to have lower family sizes "
                        "than rare variants, confirming they are sequencing noise rather than "
                        "recurrent true variants."
                    )

        # ── AB vs BA strand count scatter ─────────────────────────────────────
        _ab_has_data = con.execute(
            "SELECT COUNT(*) FROM alt_reads WHERE ab_count IS NOT NULL LIMIT 1"
        ).fetchone()[0] > 0

        if _ab_has_data:
            st.divider()
            st.subheader("AB vs BA strand counts")
            st.caption(
                "`ab_count` (aD tag) and `ba_count` (bD tag) are the number of raw reads from "
                "the AB (top) and BA (bottom) strands that contributed to each consensus read. "
                "Points on the diagonal indicate balanced duplex support from both strands. "
                "Points on the axes (ba_count = 0 or ab_count = 0) came from only one strand."
            )

            _ab_heat_df = con.execute(f"""
                SELECT
                    ar.ab_count,
                    ar.ba_count,
                    COUNT(*) AS n_reads
                FROM {_r_join}
                WHERE ar.ab_count IS NOT NULL AND ar.ba_count IS NOT NULL
                GROUP BY ar.ab_count, ar.ba_count
            """).df()

            if _ab_heat_df.empty:
                st.info("No AB/BA data in current selection.")
            else:
                _ab_sel = alt.selection_point(
                    name="ab_ba_click", fields=["ab_count", "ba_count"], on="click"
                )
                _ab_heat_chart = (
                    alt.Chart(_ab_heat_df)
                    .mark_rect()
                    .encode(
                        alt.X("ab_count:O", title="AB strand count (aD tag)"),
                        alt.Y("ba_count:O", title="BA strand count (bD tag)"),
                        alt.Color("n_reads:Q", title="Alt-supporting reads",
                                  scale=alt.Scale(scheme="blues")),
                        opacity=alt.condition(_ab_sel, alt.value(1.0), alt.value(0.4)),
                        tooltip=[
                            alt.Tooltip("ab_count:O", title="AB count"),
                            alt.Tooltip("ba_count:O", title="BA count"),
                            alt.Tooltip("n_reads:Q", title="Reads"),
                        ],
                    )
                    .add_params(_ab_sel)
                    .properties(height=400)
                )
                # key= is required: without a stable key Streamlit cannot
                # associate the selection state with the widget across reruns.
                _ev_ab = st.altair_chart(
                    _ab_heat_chart, width="stretch", on_select="rerun",
                    key="ab_ba_heatmap",
                )
                st.caption(
                    "Colour intensity = number of alt-supporting reads with each (ab_count, ba_count) combination. "
                    "Reads on the diagonal have balanced strand support; reads on the axes came from one strand only. "
                    "Click a cell to drill down."
                )

                _pts_ab = (_ev_ab.selection or {}).get("ab_ba_click", [])
                _valid_ab = [p for p in _pts_ab if "ab_count" in p and "ba_count" in p]

                if _pts_ab:
                    _ab_or = " OR ".join(
                        f"(ab_count = {int(p['ab_count'])} AND ba_count = {int(p['ba_count'])})"
                        for p in _pts_ab
                    )
                    _ab_cond = (
                        f"(sample_id, chrom, pos, alt_allele) IN ("
                        f"SELECT sample_id, chrom, pos, alt_allele "
                        f"FROM alt_reads WHERE {_ab_or})"
                    )
                    _ab_sel_df = query_records([_ab_cond])
                    _pairs_str = ", ".join(
                        f"({int(p['ab_count'])},{int(p['ba_count'])})"
                        for p in _pts_ab
                    )
                    st.caption(f"{len(_ab_sel_df):,} loci · selected (AB, BA): {_pairs_str}")
                    st.dataframe(_ab_sel_df[_table_cols], width="stretch")
                    _ab_key = "ab_ba_" + "_".join(
                        f"{int(p['ab_count'])}x{int(p['ba_count'])}"
                        for p in _pts_ab
                    )
                    igv_buttons([_ab_cond], _ab_sel_df, key=_ab_key)

# ── Tumor/Normal tab ──────────────────────────────────────────────────────────
with tab_tn:
    if not _has_normal_evidence:
        st.info(
            "No `normal_evidence` table found in this database. "
            "Run `geac annotate-normal` on each tumor/normal pair and include the "
            "`.normal_evidence.parquet` files when running `geac merge`."
        )
    else:
        st.subheader("Normal evidence at tumor loci")
        st.caption(
            "Joins `alt_bases` to `normal_evidence` to show, for every tumor alt locus, "
            "how much support for the same allele was observed in the paired normal. "
            "All active sidebar filters are applied to the tumor loci."
        )

        # ── Build the summary JOIN ─────────────────────────────────────────────
        # For each tumor locus (from the filtered alt_bases), look up:
        #   normal_depth       — from the anchor row (normal_alt_allele IS NULL)
        #   normal_alt_count   — count of reads in normal supporting the same allele
        # Two LEFT JOINs: one for the anchor row, one for the matching-allele row.
        _tn_df = con.execute(f"""
            WITH ne_anchor AS (
                SELECT tumor_sample_id, chrom, pos, tumor_alt_allele,
                       normal_sample_id, normal_depth
                FROM normal_evidence
                WHERE normal_alt_allele IS NULL
            ),
            ne_match AS (
                SELECT tumor_sample_id, chrom, pos, tumor_alt_allele,
                       SUM(normal_alt_count) AS normal_alt_count
                FROM normal_evidence
                WHERE normal_alt_allele = tumor_alt_allele
                GROUP BY tumor_sample_id, chrom, pos, tumor_alt_allele
            )
            SELECT
                ab.sample_id                                       AS tumor_sample_id,
                ab.chrom,
                ab.pos,
                ab.alt_allele                                      AS tumor_alt_allele,
                ab.variant_type,
                ROUND(ab.alt_count * 1.0 / ab.total_depth, 4)    AS tumor_vaf,
                ab.alt_count                                       AS tumor_alt_count,
                ab.total_depth                                     AS tumor_depth,
                a.normal_sample_id,
                COALESCE(a.normal_depth, 0)                       AS normal_depth,
                COALESCE(m.normal_alt_count, 0)                   AS normal_alt_count,
                CASE
                    WHEN a.normal_depth IS NULL OR a.normal_depth = 0 THEN 0.0
                    ELSE ROUND(COALESCE(m.normal_alt_count, 0) * 1.0 / a.normal_depth, 4)
                END AS normal_vaf,
                CASE
                    WHEN a.normal_depth IS NULL            THEN 'No normal data'
                    WHEN a.normal_depth = 0                THEN 'No normal coverage'
                    WHEN COALESCE(m.normal_alt_count, 0) = 0 THEN 'Somatic candidate'
                    WHEN COALESCE(m.normal_alt_count, 0) * 1.0 / a.normal_depth >= 0.2
                                                           THEN 'Germline-like'
                    ELSE                                        'Artifact-like'
                END AS classification
            FROM (SELECT * FROM {table_expr} WHERE {where}) ab
            LEFT JOIN ne_anchor a ON ab.sample_id = a.tumor_sample_id
                AND ab.chrom = a.chrom AND ab.pos = a.pos
                AND ab.alt_allele = a.tumor_alt_allele
            LEFT JOIN ne_match m ON ab.sample_id = m.tumor_sample_id
                AND ab.chrom = m.chrom AND ab.pos = m.pos
                AND ab.alt_allele = m.tumor_alt_allele
        """).df()

        if _tn_df.empty:
            st.info("No records match the current filters.")
        else:
            # ── Classification summary ─────────────────────────────────────────
            _cls_counts = (
                _tn_df["classification"]
                .value_counts()
                .rename_axis("classification")
                .reset_index(name="n_loci")
            )
            _cls_total = len(_tn_df)
            _cls_counts["pct"] = (_cls_counts["n_loci"] / _cls_total * 100).round(1)

            _cls_color = {
                "Somatic candidate":   "#2ca02c",
                "Germline-like":       "#d62728",
                "Artifact-like":       "#ff7f0e",
                "No normal coverage":  "#aec7e8",
                "No normal data":      "#c7c7c7",
            }
            _cls_order = ["Somatic candidate", "Artifact-like", "Germline-like",
                          "No normal coverage", "No normal data"]

            _cls_chart = (
                alt.Chart(_cls_counts)
                .mark_bar()
                .encode(
                    alt.X("n_loci:Q", title="Loci"),
                    alt.Y("classification:N", sort=_cls_order, title=None),
                    alt.Color("classification:N",
                              scale=alt.Scale(
                                  domain=list(_cls_color.keys()),
                                  range=list(_cls_color.values()),
                              ),
                              legend=None),
                    tooltip=[
                        alt.Tooltip("classification:N"),
                        alt.Tooltip("n_loci:Q", title="Loci"),
                        alt.Tooltip("pct:Q", title="%", format=".1f"),
                    ],
                )
                .properties(height=200, title="Loci by normal-evidence classification")
            )
            st.altair_chart(_cls_chart, width="stretch")
            st.caption(
                "**Somatic candidate**: no alt allele detected in the normal (normal_depth > 0).  "
                "**Germline-like**: same allele seen in normal at ≥ 20 % frequency.  "
                "**Artifact-like**: same allele seen in normal at < 20 % frequency (possible sequencing noise in normal).  "
                "**No normal coverage**: locus has no reads in the normal BAM.  "
                "**No normal data**: locus was not found in the normal_evidence table."
            )

            st.divider()

            # ── Tumor VAF vs Normal VAF scatter ───────────────────────────────
            st.subheader("Tumor VAF vs Normal VAF")

            _scatter_df = _tn_df[_tn_df["classification"] != "No normal data"].copy()

            if _scatter_df.empty:
                st.info("No records with normal evidence data.")
            else:
                _tn_scatter = (
                    alt.Chart(_scatter_df)
                    .mark_circle(opacity=0.6, size=40)
                    .encode(
                        alt.X("tumor_vaf:Q", title="Tumor VAF",
                              scale=alt.Scale(domain=[0, 1])),
                        alt.Y("normal_vaf:Q", title="Normal VAF",
                              scale=alt.Scale(domain=[0, 1])),
                        alt.Color("classification:N",
                                  scale=alt.Scale(
                                      domain=list(_cls_color.keys()),
                                      range=list(_cls_color.values()),
                                  )),
                        tooltip=[
                            alt.Tooltip("tumor_sample_id:N", title="Tumor sample"),
                            alt.Tooltip("normal_sample_id:N", title="Normal sample"),
                            alt.Tooltip("chrom:N"),
                            alt.Tooltip("pos:Q"),
                            alt.Tooltip("tumor_alt_allele:N", title="Alt allele"),
                            alt.Tooltip("tumor_vaf:Q", title="Tumor VAF", format=".4f"),
                            alt.Tooltip("normal_vaf:Q", title="Normal VAF", format=".4f"),
                            alt.Tooltip("tumor_alt_count:Q", title="Tumor alt count"),
                            alt.Tooltip("normal_alt_count:Q", title="Normal alt count"),
                            alt.Tooltip("normal_depth:Q", title="Normal depth"),
                            alt.Tooltip("classification:N"),
                        ],
                    )
                    .properties(height=500)
                )
                # Diagonal reference line y = x
                _diag_df = pd.DataFrame({"x": [0, 1], "y": [0, 1]})
                _diag_line = (
                    alt.Chart(_diag_df)
                    .mark_line(strokeDash=[4, 4], color="grey", opacity=0.5)
                    .encode(x="x:Q", y="y:Q")
                )
                st.altair_chart((_diag_line + _tn_scatter).properties(height=500), width="stretch")
                st.caption(
                    "Each point is one tumor alt locus. Points near the diagonal suggest germline variants. "
                    "Points on the left edge (low tumor VAF, no normal support) are somatic candidates."
                )

            st.divider()

            # ── Normal depth distribution at tumor loci ────────────────────────
            st.subheader("Normal depth at tumor loci")

            _depth_df = _tn_df[_tn_df["classification"] != "No normal data"][["normal_depth"]].copy()

            if not _depth_df.empty:
                _nd_chart = (
                    alt.Chart(_depth_df)
                    .mark_bar(color="#4c78a8")
                    .encode(
                        alt.X("normal_depth:Q", bin=alt.Bin(maxbins=50), title="Normal depth"),
                        alt.Y("count():Q", title="Loci"),
                        tooltip=[
                            alt.Tooltip("normal_depth:Q", title="Normal depth", bin=True),
                            alt.Tooltip("count():Q", title="Loci"),
                        ],
                    )
                    .properties(height=300)
                )
                st.altair_chart(_nd_chart, width="stretch")
                st.caption(
                    "Distribution of normal read depth at tumor alt loci. "
                    "Low-depth positions may have unreliable normal evidence."
                )

            st.divider()

            # ── Normal evidence data table ─────────────────────────────────────
            with st.expander("Normal evidence data table"):
                _tbl_cols = [c for c in [
                    "tumor_sample_id", "normal_sample_id", "chrom", "pos",
                    "tumor_alt_allele", "variant_type",
                    "tumor_vaf", "tumor_alt_count", "tumor_depth",
                    "normal_vaf", "normal_alt_count", "normal_depth",
                    "classification",
                ] if c in _tn_df.columns]
                st.dataframe(
                    _tn_df[_tbl_cols].sort_values(
                        ["classification", "tumor_vaf"], ascending=[True, False]
                    ),
                    width="stretch",
                    hide_index=True,
                )

# ── Panel of Normals tab ──────────────────────────────────────────────────────
with tab_pon:
    if not _has_pon_evidence:
        st.info(
            "No `pon_evidence` table found in this database. "
            "To build one: run `geac collect` on each normal sample, "
            "`geac merge` the results into a PoN DuckDB, then run "
            "`geac annotate-pon --tumor-parquet <tumor.parquet> --pon-db <pon.duckdb>` "
            "and include the `.pon_evidence.parquet` output when running `geac merge` "
            "for the cohort."
        )
    else:
        st.subheader("Panel of Normals evidence at tumor loci")
        st.caption(
            "Joins `alt_bases` to `pon_evidence` to show, for every tumor alt locus, "
            "how frequently the same allele was observed across the Panel of Normals. "
            "Loci common in the PoN are likely sequencing artefacts or germline variants. "
            "All active sidebar filters are applied to the tumor loci."
        )

        # ── Main JOIN query ────────────────────────────────────────────────────
        _pon_df = con.execute(f"""
            SELECT
                ab.sample_id                                      AS tumor_sample_id,
                ab.chrom,
                ab.pos,
                ab.alt_allele                                     AS tumor_alt_allele,
                ab.variant_type,
                ROUND(ab.alt_count * 1.0 / ab.total_depth, 4)   AS tumor_vaf,
                COALESCE(pe.n_pon_samples, 0)                    AS n_pon_samples,
                pe.pon_total_samples,
                CASE
                    WHEN pe.pon_total_samples IS NULL OR pe.pon_total_samples = 0 THEN NULL
                    ELSE ROUND(COALESCE(pe.n_pon_samples, 0) * 1.0 / pe.pon_total_samples, 4)
                END AS pon_sample_fraction,
                pe.max_pon_vaf,
                pe.mean_pon_vaf,
                CASE
                    WHEN pe.pon_total_samples IS NULL THEN 'No PoN data'
                    WHEN COALESCE(pe.n_pon_samples, 0) = 0 THEN 'PoN clean'
                    WHEN COALESCE(pe.n_pon_samples, 0) * 1.0 / pe.pon_total_samples >= 0.1
                        THEN 'Common in PoN'
                    ELSE 'Rare in PoN'
                END AS pon_classification
            FROM (SELECT * FROM {table_expr} WHERE {where}) ab
            LEFT JOIN pon_evidence pe
                   ON pe.tumor_sample_id = ab.sample_id
                  AND pe.chrom           = ab.chrom
                  AND pe.pos             = ab.pos
                  AND pe.tumor_alt_allele = ab.alt_allele
        """).df()

        if _pon_df.empty:
            st.info("No records match the current filters.")
        else:
            # ── Summary metrics ────────────────────────────────────────────────
            _pon_total = int(_pon_df["pon_total_samples"].dropna().iloc[0]) if not _pon_df["pon_total_samples"].dropna().empty else 0
            _n_clean   = int((_pon_df["pon_classification"] == "PoN clean").sum())
            _n_rare    = int((_pon_df["pon_classification"] == "Rare in PoN").sum())
            _n_common  = int((_pon_df["pon_classification"] == "Common in PoN").sum())
            _n_nodata  = int((_pon_df["pon_classification"] == "No PoN data").sum())

            _pc1, _pc2, _pc3, _pc4, _pc5 = st.columns(5)
            _pc1.metric("PoN samples",    f"{_pon_total:,}")
            _pc2.metric("PoN clean",      f"{_n_clean:,}", help="No alt allele seen in any PoN sample")
            _pc3.metric("Rare in PoN",    f"{_n_rare:,}",  help="Alt seen in < 10 % of PoN samples")
            _pc4.metric("Common in PoN",  f"{_n_common:,}", help="Alt seen in ≥ 10 % of PoN samples")
            _pc5.metric("No PoN data",    f"{_n_nodata:,}", help="Locus not in pon_evidence table")

            st.divider()

            # ── Classification bar chart ───────────────────────────────────────
            _pon_cls_counts = (
                _pon_df["pon_classification"]
                .value_counts()
                .rename_axis("classification")
                .reset_index(name="n_loci")
            )
            _pon_cls_total = len(_pon_df)
            _pon_cls_counts["pct"] = (_pon_cls_counts["n_loci"] / _pon_cls_total * 100).round(1)

            _pon_cls_color = {
                "PoN clean":       "#2ca02c",
                "Rare in PoN":     "#ff7f0e",
                "Common in PoN":   "#d62728",
                "No PoN data":     "#c7c7c7",
            }
            _pon_cls_order = ["PoN clean", "Rare in PoN", "Common in PoN", "No PoN data"]

            _pon_cls_chart = (
                alt.Chart(_pon_cls_counts)
                .mark_bar()
                .encode(
                    alt.X("n_loci:Q", title="Loci"),
                    alt.Y("classification:N", sort=_pon_cls_order, title=None),
                    alt.Color("classification:N",
                              scale=alt.Scale(
                                  domain=list(_pon_cls_color.keys()),
                                  range=list(_pon_cls_color.values()),
                              ),
                              legend=None),
                    tooltip=[
                        alt.Tooltip("classification:N"),
                        alt.Tooltip("n_loci:Q", title="Loci"),
                        alt.Tooltip("pct:Q", title="%", format=".1f"),
                    ],
                )
                .properties(height=180, title="Loci by PoN classification")
            )
            st.altair_chart(_pon_cls_chart, width="stretch")

            st.divider()

            # ── Tumor VAF vs PoN sample fraction scatter ───────────────────────
            st.subheader("Tumor VAF vs PoN sample fraction")

            _pon_scatter_df = _pon_df[_pon_df["pon_classification"] != "No PoN data"].copy()

            if not _pon_scatter_df.empty:
                _pon_scatter = (
                    alt.Chart(_pon_scatter_df)
                    .mark_circle(opacity=0.6, size=40)
                    .encode(
                        alt.X("tumor_vaf:Q", title="Tumor VAF",
                              scale=alt.Scale(domain=[0, 1])),
                        alt.Y("pon_sample_fraction:Q",
                              title=f"PoN sample fraction (N={_pon_total})",
                              scale=alt.Scale(domain=[0, 1])),
                        alt.Color("pon_classification:N",
                                  scale=alt.Scale(
                                      domain=list(_pon_cls_color.keys()),
                                      range=list(_pon_cls_color.values()),
                                  )),
                        tooltip=[
                            alt.Tooltip("tumor_sample_id:N", title="Tumor sample"),
                            alt.Tooltip("chrom:N"),
                            alt.Tooltip("pos:Q"),
                            alt.Tooltip("tumor_alt_allele:N", title="Alt allele"),
                            alt.Tooltip("variant_type:N", title="Variant type"),
                            alt.Tooltip("tumor_vaf:Q", title="Tumor VAF", format=".4f"),
                            alt.Tooltip("n_pon_samples:Q", title="PoN samples with alt"),
                            alt.Tooltip("pon_sample_fraction:Q", title="PoN fraction", format=".3f"),
                            alt.Tooltip("max_pon_vaf:Q", title="Max PoN VAF", format=".4f"),
                            alt.Tooltip("mean_pon_vaf:Q", title="Mean PoN VAF", format=".4f"),
                        ],
                    )
                    .properties(height=500)
                )
                st.altair_chart(_pon_scatter, width="stretch")
                st.caption(
                    "Each point is one tumor alt locus. "
                    "Somatic candidates cluster near Y = 0 (absent from PoN). "
                    "Recurrent artefacts and germline variants appear at higher Y values. "
                    f"PoN fraction threshold for 'Common in PoN' is ≥ 10 % ({max(1, round(_pon_total * 0.1))} / {_pon_total} samples)."
                )

            st.divider()

            # ── Max PoN VAF distribution (loci present in PoN only) ────────────
            st.subheader("Max PoN VAF distribution")

            _pon_vaf_df = _pon_df[_pon_df["max_pon_vaf"].notna()][["max_pon_vaf"]].copy()

            if not _pon_vaf_df.empty:
                _pon_vaf_chart = (
                    alt.Chart(_pon_vaf_df)
                    .mark_bar(color="#d62728")
                    .encode(
                        alt.X("max_pon_vaf:Q", bin=alt.Bin(maxbins=50),
                              title="Max PoN VAF (highest alt frequency seen in any PoN sample)"),
                        alt.Y("count():Q", title="Loci"),
                        tooltip=[
                            alt.Tooltip("max_pon_vaf:Q", title="Max PoN VAF", bin=True),
                            alt.Tooltip("count():Q", title="Loci"),
                        ],
                    )
                    .properties(height=280)
                )
                st.altair_chart(_pon_vaf_chart, width="stretch")
                st.caption(
                    "Distribution of the highest alt VAF seen across all PoN samples at each locus "
                    "(restricted to loci with ≥ 1 PoN sample showing the alt). "
                    "Low max VAF suggests sequencing noise; high max VAF suggests a germline polymorphism."
                )
            else:
                st.info("No loci with PoN alt evidence in current selection.")

            st.divider()

            # ── Data table ─────────────────────────────────────────────────────
            with st.expander("PoN evidence data table"):
                _pon_tbl_cols = [c for c in [
                    "tumor_sample_id", "chrom", "pos", "tumor_alt_allele",
                    "variant_type", "tumor_vaf",
                    "n_pon_samples", "pon_total_samples", "pon_sample_fraction",
                    "max_pon_vaf", "mean_pon_vaf", "pon_classification",
                ] if c in _pon_df.columns]
                st.dataframe(
                    _pon_df[_pon_tbl_cols].sort_values(
                        "pon_sample_fraction", ascending=False
                    ),
                    width="stretch",
                    hide_index=True,
                )

# ──────────────────────────────────────────────────────────────────────────────
# Tab 10 — Pipeline comparison
# ──────────────────────────────────────────────────────────────────────────────
with tab_pipeline:
    if not data_source.is_duckdb:
        st.info(
            "Pipeline comparison requires a merged DuckDB file. "
            "Run `geac collect` twice with different `--pipeline` values and the same "
            "`--sample-id`, then `geac merge` both Parquets into one DuckDB."
        )
    else:
        _pc_pipelines = data_source.distinct_values("pipeline")
        if len(_pc_pipelines) < 2:
            st.info(
                "Pipeline comparison requires at least 2 distinct `pipeline` values in the database. "
                f"This database contains only: {', '.join(str(p) for p in _pc_pipelines)}. "
                "Re-run `geac collect` with a different `--pipeline` value for the same sample "
                "and merge both Parquets into one DuckDB."
            )
        else:
            st.caption(
                "Compare the same sample processed through two different pipelines. "
                "All active sidebar filters (chromosome, sample, VAF, etc.) are applied to both pipelines."
            )

            # ── Pipeline selectors ─────────────────────────────────────────────
            _pc_col1, _pc_col2 = st.columns(2)
            _pc_pipe_a = _pc_col1.selectbox(
                "Pipeline A", _pc_pipelines, index=0, key="pipe_cmp_a"
            )
            _pc_pipe_b_opts = [p for p in _pc_pipelines if p != _pc_pipe_a]
            _pc_pipe_b = _pc_col2.selectbox(
                "Pipeline B", _pc_pipe_b_opts, index=0, key="pipe_cmp_b"
            )

            # Strip any pipeline_sel filter from conditions — this tab manages
            # pipeline selection internally via its own A/B selectors.
            _pc_base_conds = [c for c in conditions if not c.startswith("pipeline IN")]
            _pc_wa = " AND ".join(
                _pc_base_conds + [f"pipeline = '{_sql_str(str(_pc_pipe_a))}'"]
            ) if _pc_base_conds else f"pipeline = '{_sql_str(str(_pc_pipe_a))}'"
            _pc_wb = " AND ".join(
                _pc_base_conds + [f"pipeline = '{_sql_str(str(_pc_pipe_b))}'"]
            ) if _pc_base_conds else f"pipeline = '{_sql_str(str(_pc_pipe_b))}'"

            # ── Build FULL OUTER JOIN dataset once, reused by all views ───────
            _pc_df = con.execute(f"""
                WITH a AS (
                    SELECT sample_id, chrom, pos, alt_allele, variant_type,
                           ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                           total_depth, alt_count,
                           trinuc_context, ref_allele
                    FROM {table_expr}
                    WHERE {_pc_wa}
                ),
                b AS (
                    SELECT sample_id, chrom, pos, alt_allele, variant_type,
                           ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                           total_depth, alt_count,
                           trinuc_context, ref_allele
                    FROM {table_expr}
                    WHERE {_pc_wb}
                )
                SELECT
                    COALESCE(a.sample_id,    b.sample_id)    AS sample_id,
                    COALESCE(a.chrom,        b.chrom)        AS chrom,
                    COALESCE(a.pos,          b.pos)          AS pos,
                    COALESCE(a.alt_allele,   b.alt_allele)   AS alt_allele,
                    COALESCE(a.variant_type, b.variant_type) AS variant_type,
                    a.vaf         AS vaf_a,
                    b.vaf         AS vaf_b,
                    a.total_depth AS depth_a,
                    b.total_depth AS depth_b,
                    a.alt_count   AS alt_count_a,
                    b.alt_count   AS alt_count_b,
                    COALESCE(a.trinuc_context, b.trinuc_context) AS trinuc_context,
                    COALESCE(a.ref_allele, b.ref_allele) AS ref_allele,
                    CASE
                        WHEN a.chrom IS NOT NULL AND b.chrom IS NOT NULL THEN 'shared'
                        WHEN a.chrom IS NOT NULL                         THEN 'only_a'
                        ELSE                                                  'only_b'
                    END AS concordance
                FROM a
                FULL OUTER JOIN b
                    ON  a.sample_id  = b.sample_id
                    AND a.chrom      = b.chrom
                    AND a.pos        = b.pos
                    AND a.alt_allele = b.alt_allele
            """).df()

            if _pc_df.empty:
                st.info("No records for either pipeline under the current filters.")
            else:
                _pc_label_map = {
                    "shared": "Shared",
                    "only_a": f"Only {_pc_pipe_a}",
                    "only_b": f"Only {_pc_pipe_b}",
                }
                _pc_df["concordance_label"] = _pc_df["concordance"].map(_pc_label_map)
                _pc_shared = _pc_df[_pc_df["concordance"] == "shared"]

                # ── View 1: Concordance summary ────────────────────────────────
                st.subheader("Locus concordance")

                _pc_n_shared = int((_pc_df["concordance"] == "shared").sum())
                _pc_n_only_a = int((_pc_df["concordance"] == "only_a").sum())
                _pc_n_only_b = int((_pc_df["concordance"] == "only_b").sum())
                _pc_m1, _pc_m2, _pc_m3 = st.columns(3)
                _pc_m1.metric("Shared loci", f"{_pc_n_shared:,}")
                _pc_m2.metric(f"Only {_pc_pipe_a}", f"{_pc_n_only_a:,}")
                _pc_m3.metric(f"Only {_pc_pipe_b}", f"{_pc_n_only_b:,}")

                _pc_conc_counts = (
                    _pc_df.groupby(["concordance_label", "variant_type"])
                    .size()
                    .reset_index(name="n_loci")
                )
                _pc_conc_order = ["Shared", f"Only {_pc_pipe_a}", f"Only {_pc_pipe_b}"]
                _pc_conc_chart = (
                    alt.Chart(_pc_conc_counts)
                    .mark_bar()
                    .encode(
                        x=alt.X("concordance_label:N", title=None, sort=_pc_conc_order),
                        y=alt.Y("n_loci:Q", title="Loci"),
                        color=alt.Color("variant_type:N", title="Variant type"),
                        tooltip=["concordance_label", "variant_type", "n_loci"],
                    )
                    .properties(height=280)
                )
                st.altair_chart(_pc_conc_chart, width="stretch")

                st.divider()

                # ── View 2: VAF correlation ────────────────────────────────────
                st.subheader("VAF correlation (shared loci)")

                if _pc_shared.empty:
                    st.info("No shared loci for current filters.")
                else:
                    _pc_diag = (
                        alt.Chart(
                            pd.DataFrame({"x": [0.0, 1.0], "y": [0.0, 1.0]})
                        )
                        .mark_line(color="gray", strokeDash=[4, 4], opacity=0.6)
                        .encode(x="x:Q", y="y:Q")
                    )
                    _pc_vaf_sel = alt.selection_point(
                        name="pc_vaf_select",
                        fields=["sample_id", "chrom", "pos", "alt_allele"],
                        on="click",
                        toggle="event.shiftKey",
                    )
                    _pc_vaf_scatter = (
                        alt.Chart(_pc_shared)
                        .mark_circle(size=40)
                        .encode(
                            x=alt.X(
                                "vaf_a:Q",
                                title=f"VAF — {_pc_pipe_a}",
                                scale=alt.Scale(domain=[0.0, 1.0]),
                            ),
                            y=alt.Y(
                                "vaf_b:Q",
                                title=f"VAF — {_pc_pipe_b}",
                                scale=alt.Scale(domain=[0.0, 1.0]),
                            ),
                            color=alt.Color("variant_type:N", title="Variant type"),
                            opacity=alt.condition(_pc_vaf_sel, alt.value(1.0), alt.value(0.4)),
                            size=alt.condition(_pc_vaf_sel, alt.value(120), alt.value(40)),
                            tooltip=[
                                "sample_id", "chrom", "pos", "alt_allele", "variant_type",
                                alt.Tooltip("vaf_a:Q", title=f"VAF ({_pc_pipe_a})", format=".4f"),
                                alt.Tooltip("vaf_b:Q", title=f"VAF ({_pc_pipe_b})", format=".4f"),
                            ],
                        )
                        .add_params(_pc_vaf_sel)
                        .properties(width=450, height=400)
                        .interactive()
                    )
                    _pc_vaf_event = st.altair_chart(
                        _pc_diag + _pc_vaf_scatter,
                        width="stretch",
                        on_select="rerun",
                        key="pc_vaf_scatter",
                    )

                    _pc_vaf_pts = (_pc_vaf_event.selection or {}).get("pc_vaf_select", [])
                    if _pc_vaf_pts:
                        _pc_vaf_or = " OR ".join(
                            f"(sample_id = '{_sql_str(p['sample_id'])}' AND chrom = '{_sql_str(p['chrom'])}' "
                            f"AND pos = {int(p['pos'])} AND alt_allele = '{_sql_str(p['alt_allele'])}')"
                            for p in _pc_vaf_pts
                            if all(k in p for k in ["sample_id", "chrom", "pos", "alt_allele"])
                        )
                        if _pc_vaf_or:
                            _pc_vaf_sel_df = con.execute(f"""
                                SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
                                FROM {table_expr}
                                WHERE ({_pc_vaf_or})
                                ORDER BY sample_id, chrom, pos, pipeline
                            """).df()
                            _pc_vaf_show_cols = [c for c in _table_cols + ["pipeline"] if c in _pc_vaf_sel_df.columns]
                            st.caption(
                                f"{len(_pc_vaf_sel_df)} records for {len(_pc_vaf_pts)} selected "
                                f"loci (both pipelines shown) — shift-click to select multiple"
                            )
                            st.dataframe(_pc_vaf_sel_df[_pc_vaf_show_cols], width="stretch", hide_index=True)
                            igv_buttons(
                                [f"({_pc_vaf_or})"],
                                _pc_vaf_sel_df,
                                key=f"pc_vaf_{'_'.join(str(int(p['pos'])) for p in _pc_vaf_pts[:5] if 'pos' in p)}",
                                use_global_filters=False,
                            )
                    else:
                        st.caption("Click a point to select it; shift-click to select multiple.")

                    if len(_pc_shared) >= 2:
                        _pc_r = _pc_shared["vaf_a"].corr(_pc_shared["vaf_b"])
                        st.caption(
                            f"Pearson r = {_pc_r:.4f} across {len(_pc_shared):,} shared loci. "
                            "Points off the diagonal are discordant calls."
                        )

                st.divider()

                # ── View 3: Unique-to-pipeline loci table ──────────────────────
                st.subheader("Loci unique to one pipeline")
                st.caption(
                    "One pipeline calls these loci; the other does not. "
                    "Highest-priority candidates for manual review."
                )

                _pc_uniq = _pc_df[_pc_df["concordance"] != "shared"].copy()
                _pc_uniq["pipeline"] = _pc_uniq["concordance"].map(
                    {"only_a": _pc_pipe_a, "only_b": _pc_pipe_b}
                )
                _pc_uniq_filter = st.radio(
                    "Show",
                    ["Both", f"Only {_pc_pipe_a}", f"Only {_pc_pipe_b}"],
                    horizontal=True,
                    key="pipe_cmp_uniq_filter",
                )
                if _pc_uniq_filter == f"Only {_pc_pipe_a}":
                    _pc_uniq_show = _pc_uniq[_pc_uniq["concordance"] == "only_a"]
                elif _pc_uniq_filter == f"Only {_pc_pipe_b}":
                    _pc_uniq_show = _pc_uniq[_pc_uniq["concordance"] == "only_b"]
                else:
                    _pc_uniq_show = _pc_uniq

                _pc_uniq_cols = [c for c in [
                    "pipeline", "sample_id", "chrom", "pos", "alt_allele",
                    "variant_type", "vaf_a", "vaf_b", "depth_a", "depth_b",
                ] if c in _pc_uniq_show.columns]

                if _pc_uniq_show.empty:
                    st.success("No unique loci for current filter selection.")
                else:
                    st.dataframe(
                        _pc_uniq_show[_pc_uniq_cols]
                        .rename(columns={
                            "vaf_a":   f"vaf_{_pc_pipe_a}",
                            "vaf_b":   f"vaf_{_pc_pipe_b}",
                            "depth_a": f"depth_{_pc_pipe_a}",
                            "depth_b": f"depth_{_pc_pipe_b}",
                        })
                        .sort_values(["pipeline", "chrom", "pos"]),
                        width="stretch",
                        hide_index=True,
                        key="pipe_cmp_uniq_table",
                    )

                st.divider()

                # ── View 4: Unique-loci characterization ──────────────────────
                st.subheader("Unique loci characterization")
                st.caption(
                    "These plots compare loci called by only one pipeline under the current filters. "
                    "They are intended to help explain whether unique calls differ in VAF, support, "
                    "or mutational context between the two pipelines."
                )

                _pc_uniq_cmp = build_unique_pipeline_characterization_df(
                    _pc_uniq, str(_pc_pipe_a), str(_pc_pipe_b)
                )
                _pc_uniq_summary = summarize_unique_pipeline_groups(
                    _pc_uniq_cmp, str(_pc_pipe_a), str(_pc_pipe_b)
                )
                _pc_uniq_group_a = f"Only {_pc_pipe_a}"
                _pc_uniq_group_b = f"Only {_pc_pipe_b}"
                _pc_uniq_a = _pc_uniq_cmp[_pc_uniq_cmp["group"] == _pc_uniq_group_a]
                _pc_uniq_b = _pc_uniq_cmp[_pc_uniq_cmp["group"] == _pc_uniq_group_b]

                if _pc_uniq_a.empty or _pc_uniq_b.empty:
                    st.info(
                        "Unique-loci characterization needs at least one unique locus in each pipeline "
                        "under the current filters."
                    )
                else:
                    def _pc_fmt_metric(v, fmt):
                        return "NA" if v is None or pd.isna(v) else format(v, fmt)

                    _pc_um = _pc_uniq_summary["metrics"]
                    _pc_mcols = st.columns(6)
                    _pc_mcols[0].metric(f"{_pc_pipe_a} unique loci", f"{_pc_um[_pc_uniq_group_a]['count']:,}")
                    _pc_mcols[1].metric(
                        f"{_pc_pipe_a} median VAF",
                        _pc_fmt_metric(_pc_um[_pc_uniq_group_a]["median_vaf"], ".4f"),
                    )
                    _pc_mcols[2].metric(
                        f"{_pc_pipe_a} median depth",
                        _pc_fmt_metric(_pc_um[_pc_uniq_group_a]["median_depth"], ".1f"),
                    )
                    _pc_mcols[3].metric(f"{_pc_pipe_b} unique loci", f"{_pc_um[_pc_uniq_group_b]['count']:,}")
                    _pc_mcols[4].metric(
                        f"{_pc_pipe_b} median VAF",
                        _pc_fmt_metric(_pc_um[_pc_uniq_group_b]["median_vaf"], ".4f"),
                    )
                    _pc_mcols[5].metric(
                        f"{_pc_pipe_b} median depth",
                        _pc_fmt_metric(_pc_um[_pc_uniq_group_b]["median_depth"], ".1f"),
                    )
                    st.caption(_pc_uniq_summary["summary"])

                    def _pc_overlay_dist_chart(
                        df: pd.DataFrame,
                        value_col: str,
                        title: str,
                        x_title: str,
                        key_col: str = "group",
                        domain=None,
                    ):
                        _counts = df.groupby(key_col).size()
                        if not _counts.empty and int(_counts.min()) >= 5:
                            _density_kwargs = {
                                "as_": [value_col, "density"],
                                "groupby": [key_col],
                                "steps": 100,
                            }
                            if domain is not None:
                                _density_kwargs["extent"] = domain
                            return (
                                alt.Chart(df)
                                .transform_density(value_col, **_density_kwargs)
                                .mark_line(strokeWidth=2.5)
                                .encode(
                                    x=alt.X(
                                        f"{value_col}:Q",
                                        title=x_title,
                                        scale=alt.Scale(domain=domain) if domain else alt.Scale(),
                                    ),
                                    y=alt.Y("density:Q", title="Density", stack=None),
                                    color=alt.Color(f"{key_col}:N", title="Unique set"),
                                    tooltip=[
                                        alt.Tooltip(f"{value_col}:Q", format=".4f" if value_col == "vaf" else ".2f"),
                                        f"{key_col}:N",
                                        alt.Tooltip("density:Q", format=".3f"),
                                    ],
                                )
                                .properties(title=title, height=260)
                            )
                        return (
                            alt.Chart(df)
                            .mark_bar(opacity=0.45)
                            .encode(
                                x=alt.X(
                                    f"{value_col}:Q",
                                    title=x_title,
                                    bin=alt.Bin(maxbins=40),
                                    scale=alt.Scale(domain=domain) if domain else alt.Scale(),
                                ),
                                y=alt.Y("count():Q", title="Count", stack=None),
                                color=alt.Color(f"{key_col}:N", title="Unique set"),
                                tooltip=[
                                    f"{key_col}:N",
                                    alt.Tooltip("count():Q", title="Count"),
                                ],
                            )
                            .properties(title=title, height=260)
                        )

                    _pc_dist_cols = st.columns(2)
                    with _pc_dist_cols[0]:
                        st.altair_chart(
                            _pc_overlay_dist_chart(
                                _pc_uniq_cmp.dropna(subset=["vaf"]),
                                "vaf",
                                "VAF distribution",
                                "VAF",
                                domain=[0, 1],
                            ),
                            width="stretch",
                        )
                    with _pc_dist_cols[1]:
                        st.altair_chart(
                            _pc_overlay_dist_chart(
                                _pc_uniq_cmp.dropna(subset=["depth"]),
                                "depth",
                                "Total depth distribution",
                                "Total depth",
                            ),
                            width="stretch",
                        )

                    st.altair_chart(
                        _pc_overlay_dist_chart(
                            _pc_uniq_cmp.dropna(subset=["alt_count"]),
                            "alt_count",
                            "Alt-count distribution",
                            "Alt count",
                        ),
                        width="stretch",
                    )

                    st.markdown("**Unique-only SBS96 spectrum**")
                    if not _has_data("trinuc_context"):
                        st.info(
                            "SBS96 spectrum requires the `trinuc_context` column. "
                            "Re-run `geac collect` with a reference FASTA."
                        )
                    else:
                        def _pc_uniq_sbs96(group_name: str):
                            _raw = (
                                _pc_uniq_cmp[
                                    (_pc_uniq_cmp["group"] == group_name)
                                    & (_pc_uniq_cmp["variant_type"] == "SNV")
                                    & (_pc_uniq_cmp["trinuc_context"].notna())
                                ][["trinuc_context", "ref_allele", "alt_allele"]]
                                .groupby(["trinuc_context", "ref_allele", "alt_allele"])
                                .size()
                                .reset_index(name="count")
                            )
                            return _to_spec96_strat(_raw)

                        _pc_u96_a, _pc_u_total_a = _pc_uniq_sbs96(_pc_uniq_group_a)
                        _pc_u96_b, _pc_u_total_b = _pc_uniq_sbs96(_pc_uniq_group_b)
                        if _pc_u96_a is None or _pc_u96_b is None:
                            st.info("Insufficient unique SNV data for unique-only SBS96 comparison.")
                        else:
                            _pc_u_sc1, _pc_u_sc2 = st.columns(2)
                            with _pc_u_sc1:
                                st.markdown(f"**{_pc_pipe_a} unique loci** ({_pc_u_total_a:,} SNVs)")
                                st.altair_chart(
                                    _strat_sbs96_chart(_pc_u96_a, f"Only {_pc_pipe_a}"),
                                    width="stretch",
                                )
                            with _pc_u_sc2:
                                st.markdown(f"**{_pc_pipe_b} unique loci** ({_pc_u_total_b:,} SNVs)")
                                st.altair_chart(
                                    _strat_sbs96_chart(_pc_u96_b, f"Only {_pc_pipe_b}"),
                                    width="stretch",
                                )

                st.divider()

                # ── View 5: SBS96 spectrum side-by-side ────────────────────────
                st.subheader("SBS96 error spectrum")

                if not _has_data("trinuc_context"):
                    st.info(
                        "SBS96 spectrum requires the `trinuc_context` column. "
                        "Re-run `geac collect` with a reference FASTA."
                    )
                else:
                    def _pc_sbs96(pipe_where: str):
                        _raw = con.execute(f"""
                            SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
                            FROM {table_expr}
                            WHERE {pipe_where}
                              AND variant_type = 'SNV'
                              AND trinuc_context IS NOT NULL
                              AND length(trinuc_context) = 3
                            GROUP BY trinuc_context, ref_allele, alt_allele
                        """).df()
                        return _to_spec96_strat(_raw)

                    _pc_s96_a, _pc_total_a = _pc_sbs96(_pc_wa)
                    _pc_s96_b, _pc_total_b = _pc_sbs96(_pc_wb)

                    if _pc_s96_a is None or _pc_s96_b is None:
                        st.info("Insufficient SNV data for spectrum comparison.")
                    else:
                        _pc_sc1, _pc_sc2 = st.columns(2)
                        with _pc_sc1:
                            st.markdown(f"**{_pc_pipe_a}** ({_pc_total_a:,} SNVs)")
                            st.altair_chart(
                                _strat_sbs96_chart(_pc_s96_a, _pc_pipe_a),
                                width="stretch",
                            )
                        with _pc_sc2:
                            st.markdown(f"**{_pc_pipe_b}** ({_pc_total_b:,} SNVs)")
                            st.altair_chart(
                                _strat_sbs96_chart(_pc_s96_b, _pc_pipe_b),
                                width="stretch",
                            )

                st.divider()

                # ── View 6: Depth comparison ───────────────────────────────────
                st.subheader("Depth comparison (shared loci)")

                if _pc_shared.empty:
                    st.info("No shared loci for current filters.")
                else:
                    _pc_depth_max = int(
                        max(_pc_shared["depth_a"].max(), _pc_shared["depth_b"].max())
                    )
                    _pc_depth_diag = (
                        alt.Chart(
                            pd.DataFrame({"x": [0, _pc_depth_max], "y": [0, _pc_depth_max]})
                        )
                        .mark_line(color="gray", strokeDash=[4, 4], opacity=0.6)
                        .encode(x="x:Q", y="y:Q")
                    )
                    _pc_depth_sel = alt.selection_point(
                        name="pc_depth_select",
                        fields=["sample_id", "chrom", "pos", "alt_allele"],
                        on="click",
                        toggle="event.shiftKey",
                    )
                    _pc_depth_scatter = (
                        alt.Chart(_pc_shared)
                        .mark_circle(size=40)
                        .encode(
                            x=alt.X("depth_a:Q", title=f"Total depth — {_pc_pipe_a}"),
                            y=alt.Y("depth_b:Q", title=f"Total depth — {_pc_pipe_b}"),
                            color=alt.Color("variant_type:N", title="Variant type"),
                            opacity=alt.condition(_pc_depth_sel, alt.value(1.0), alt.value(0.4)),
                            size=alt.condition(_pc_depth_sel, alt.value(120), alt.value(40)),
                            tooltip=[
                                "sample_id", "chrom", "pos", "alt_allele", "variant_type",
                                alt.Tooltip("depth_a:Q", title=f"Depth ({_pc_pipe_a})"),
                                alt.Tooltip("depth_b:Q", title=f"Depth ({_pc_pipe_b})"),
                            ],
                        )
                        .add_params(_pc_depth_sel)
                        .properties(width=450, height=350)
                        .interactive()
                    )
                    _pc_depth_event = st.altair_chart(
                        _pc_depth_diag + _pc_depth_scatter,
                        width="stretch",
                        on_select="rerun",
                        key="pc_depth_scatter",
                    )

                    _pc_depth_pts = (_pc_depth_event.selection or {}).get("pc_depth_select", [])
                    if _pc_depth_pts:
                        _pc_depth_or = " OR ".join(
                            f"(sample_id = '{_sql_str(p['sample_id'])}' AND chrom = '{_sql_str(p['chrom'])}' "
                            f"AND pos = {int(p['pos'])} AND alt_allele = '{_sql_str(p['alt_allele'])}')"
                            for p in _pc_depth_pts
                            if all(k in p for k in ["sample_id", "chrom", "pos", "alt_allele"])
                        )
                        if _pc_depth_or:
                            _pc_depth_sel_df = con.execute(f"""
                                SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf
                                FROM {table_expr}
                                WHERE ({_pc_depth_or})
                                ORDER BY sample_id, chrom, pos, pipeline
                            """).df()
                            _pc_depth_show_cols = [
                                c for c in _table_cols + ["pipeline"] if c in _pc_depth_sel_df.columns
                            ]
                            st.caption(
                                f"{len(_pc_depth_sel_df)} records for {len(_pc_depth_pts)} selected "
                                f"loci (both pipelines shown) — shift-click to select multiple"
                            )
                            st.dataframe(
                                _pc_depth_sel_df[_pc_depth_show_cols],
                                width="stretch",
                                hide_index=True,
                            )
                            igv_buttons(
                                [f"({_pc_depth_or})"],
                                _pc_depth_sel_df,
                                key=f"pc_depth_{'_'.join(str(int(p['pos'])) for p in _pc_depth_pts[:5] if 'pos' in p)}",
                                use_global_filters=False,
                            )
                    else:
                        st.caption("Click a point to select it; shift-click to select multiple.")

                    st.caption(
                        f"Points above the diagonal: higher depth in {_pc_pipe_b}. "
                        f"Points below: higher depth in {_pc_pipe_a}. "
                        "Systematic offset suggests different duplicate-collapsing or overlap behaviour."
                    )

with tab_read_type:
    if not data_source.is_duckdb:
        st.info(
            "Read-type comparison requires a merged DuckDB file. "
            "Run `geac collect` separately for each BAM (raw, simplex, duplex) with the same "
            "`--sample-id`, then `geac merge` all Parquets into one DuckDB."
        )
    else:
        _rt_read_types = data_source.distinct_values("read_type")
        if len(_rt_read_types) < 2:
            st.info(
                "Read-type comparison requires at least 2 distinct `read_type` values in the database. "
                f"This database contains only: {', '.join(str(r) for r in _rt_read_types)}. "
                "Re-run `geac collect` with a different `--read-type` value for the same sample "
                "and merge both Parquets into one DuckDB."
            )
        else:
            st.caption(
                "Compare the same sample processed at two different read types (e.g. raw vs duplex). "
                "All active sidebar filters (chromosome, sample, VAF, etc.) are applied to both. "
                "The goal is to quantify what duplex consensus processing removes vs retains."
            )

            # ── Read-type selectors ────────────────────────────────────────────
            _rt_col1, _rt_col2 = st.columns(2)
            _rt_a = _rt_col1.selectbox(
                "Read type A", _rt_read_types, index=0, key="rt_cmp_a"
            )
            _rt_b_opts = [r for r in _rt_read_types if r != _rt_a]
            _rt_b = _rt_col2.selectbox(
                "Read type B", _rt_b_opts, index=0, key="rt_cmp_b"
            )

            _rt_wa = " AND ".join(
                conditions + [f"read_type = '{_sql_str(str(_rt_a))}'"]
            ) if conditions else f"read_type = '{_sql_str(str(_rt_a))}'"
            _rt_wb = " AND ".join(
                conditions + [f"read_type = '{_sql_str(str(_rt_b))}'"]
            ) if conditions else f"read_type = '{_sql_str(str(_rt_b))}'"

            # ── Build FULL OUTER JOIN dataset once, reused by all views ───────
            _rt_df = con.execute(f"""
                WITH a AS (
                    SELECT sample_id, chrom, pos, alt_allele, variant_type,
                           ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                           total_depth,
                           fwd_alt_count, rev_alt_count,
                           trinuc_context, ref_allele
                    FROM {table_expr}
                    WHERE {_rt_wa}
                ),
                b AS (
                    SELECT sample_id, chrom, pos, alt_allele, variant_type,
                           ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
                           total_depth,
                           fwd_alt_count, rev_alt_count
                    FROM {table_expr}
                    WHERE {_rt_wb}
                )
                SELECT
                    COALESCE(a.sample_id,    b.sample_id)    AS sample_id,
                    COALESCE(a.chrom,        b.chrom)        AS chrom,
                    COALESCE(a.pos,          b.pos)          AS pos,
                    COALESCE(a.alt_allele,   b.alt_allele)   AS alt_allele,
                    COALESCE(a.variant_type, b.variant_type) AS variant_type,
                    a.vaf            AS vaf_a,
                    b.vaf            AS vaf_b,
                    a.total_depth    AS depth_a,
                    b.total_depth    AS depth_b,
                    a.fwd_alt_count  AS fwd_alt_a,
                    a.rev_alt_count  AS rev_alt_a,
                    b.fwd_alt_count  AS fwd_alt_b,
                    b.rev_alt_count  AS rev_alt_b,
                    a.trinuc_context,
                    a.ref_allele,
                    CASE
                        WHEN a.chrom IS NOT NULL AND b.chrom IS NOT NULL THEN 'shared'
                        WHEN a.chrom IS NOT NULL                         THEN 'only_a'
                        ELSE                                                  'only_b'
                    END AS concordance
                FROM a
                FULL OUTER JOIN b
                    ON  a.sample_id  = b.sample_id
                    AND a.chrom      = b.chrom
                    AND a.pos        = b.pos
                    AND a.alt_allele = b.alt_allele
            """).df()

            if _rt_df.empty:
                st.info("No records for either read type under the current filters.")
            else:
                _rt_label_map = {
                    "shared": "Shared",
                    "only_a": f"Only {_rt_a}",
                    "only_b": f"Only {_rt_b}",
                }
                _rt_df["concordance_label"] = _rt_df["concordance"].map(_rt_label_map)
                _rt_shared = _rt_df[_rt_df["concordance"] == "shared"]

                # ── View 1: Concordance summary ────────────────────────────────
                st.subheader("Locus concordance")
                st.caption(
                    "Loci unique to the less-filtered read type were likely removed by consensus "
                    "processing. Shared loci are retained across both read types."
                )

                _rt_n_shared = int((_rt_df["concordance"] == "shared").sum())
                _rt_n_only_a = int((_rt_df["concordance"] == "only_a").sum())
                _rt_n_only_b = int((_rt_df["concordance"] == "only_b").sum())
                _rt_m1, _rt_m2, _rt_m3 = st.columns(3)
                _rt_m1.metric("Shared loci", f"{_rt_n_shared:,}")
                _rt_m2.metric(f"Only {_rt_a}", f"{_rt_n_only_a:,}")
                _rt_m3.metric(f"Only {_rt_b}", f"{_rt_n_only_b:,}")

                _rt_conc_counts = (
                    _rt_df.groupby(["concordance_label", "variant_type"])
                    .size()
                    .reset_index(name="n_loci")
                )
                _rt_conc_order = ["Shared", f"Only {_rt_a}", f"Only {_rt_b}"]
                _rt_conc_chart = (
                    alt.Chart(_rt_conc_counts)
                    .mark_bar()
                    .encode(
                        x=alt.X("concordance_label:N", title=None, sort=_rt_conc_order),
                        y=alt.Y("n_loci:Q", title="Loci"),
                        color=alt.Color("variant_type:N", title="Variant type"),
                        tooltip=["concordance_label", "variant_type", "n_loci"],
                    )
                    .properties(height=280)
                )
                st.altair_chart(_rt_conc_chart, width="stretch")

                st.divider()

                # ── View 2: VAF distribution overlay ──────────────────────────
                st.subheader("VAF distribution")
                st.caption(
                    "Duplex processing typically shifts the VAF distribution by removing "
                    "low-VAF artefacts that arise from PCR or sequencing errors."
                )

                _rt_vaf_a = _rt_df[_rt_df["concordance"].isin(["shared", "only_a"])][["vaf_a", "variant_type"]].dropna(subset=["vaf_a"]).rename(columns={"vaf_a": "vaf"})
                _rt_vaf_a["read_type"] = str(_rt_a)
                _rt_vaf_b = _rt_df[_rt_df["concordance"].isin(["shared", "only_b"])][["vaf_b", "variant_type"]].dropna(subset=["vaf_b"]).rename(columns={"vaf_b": "vaf"})
                _rt_vaf_b["read_type"] = str(_rt_b)
                _rt_vaf_both = pd.concat([_rt_vaf_a, _rt_vaf_b], ignore_index=True)

                if not _rt_vaf_both.empty:
                    _rt_vaf_chart = (
                        alt.Chart(_rt_vaf_both)
                        .transform_density(
                            "vaf",
                            as_=["vaf", "density"],
                            groupby=["read_type"],
                            extent=[0, 1],
                            steps=100,
                        )
                        .mark_area(opacity=0.4)
                        .encode(
                            x=alt.X("vaf:Q", title="VAF", scale=alt.Scale(domain=[0, 1])),
                            y=alt.Y("density:Q", title="Density", stack=None),
                            color=alt.Color("read_type:N", title="Read type"),
                            tooltip=[
                                alt.Tooltip("vaf:Q", format=".3f"),
                                "read_type:N",
                                alt.Tooltip("density:Q", format=".3f"),
                            ],
                        )
                        .properties(height=280)
                    )
                    st.altair_chart(_rt_vaf_chart, width="stretch")

                st.divider()

                # ── View 3: VAF correlation (shared loci) ──────────────────────
                st.subheader("VAF correlation (shared loci)")

                if _rt_shared.empty:
                    st.info("No shared loci for current filters.")
                else:
                    _rt_vaf_diag = (
                        alt.Chart(
                            pd.DataFrame({"x": [0.0, 1.0], "y": [0.0, 1.0]})
                        )
                        .mark_line(color="gray", strokeDash=[4, 4], opacity=0.6)
                        .encode(x="x:Q", y="y:Q")
                    )
                    _rt_vaf_scatter = (
                        alt.Chart(_rt_shared)
                        .mark_circle(opacity=0.5, size=40)
                        .encode(
                            x=alt.X(
                                "vaf_a:Q",
                                title=f"VAF — {_rt_a}",
                                scale=alt.Scale(domain=[0.0, 1.0]),
                            ),
                            y=alt.Y(
                                "vaf_b:Q",
                                title=f"VAF — {_rt_b}",
                                scale=alt.Scale(domain=[0.0, 1.0]),
                            ),
                            color=alt.Color("variant_type:N", title="Variant type"),
                            tooltip=[
                                "sample_id", "chrom", "pos", "alt_allele", "variant_type",
                                alt.Tooltip("vaf_a:Q", title=f"VAF ({_rt_a})", format=".4f"),
                                alt.Tooltip("vaf_b:Q", title=f"VAF ({_rt_b})", format=".4f"),
                            ],
                        )
                        .properties(width=450, height=380)
                        .interactive()
                    )
                    st.altair_chart(_rt_vaf_diag + _rt_vaf_scatter, width="stretch")

                    if len(_rt_shared) >= 2:
                        _rt_r = _rt_shared["vaf_a"].corr(_rt_shared["vaf_b"])
                        st.caption(
                            f"Pearson r = {_rt_r:.4f} across {len(_rt_shared):,} shared loci. "
                            "Points below the diagonal indicate higher VAF in read type B."
                        )

                st.divider()

                # ── View 4: Strand balance comparison ─────────────────────────
                st.subheader("Strand balance")
                st.caption(
                    "Fraction of alt-supporting reads on the forward strand. "
                    "Values near 0.5 indicate balanced strand support; "
                    "values near 0 or 1 suggest strand artefacts."
                )

                _rt_strand_rows = []
                for _rt_lbl, _rt_fw_col, _rt_rv_col in [
                    (str(_rt_a), "fwd_alt_a", "rev_alt_a"),
                    (str(_rt_b), "fwd_alt_b", "rev_alt_b"),
                ]:
                    _rt_sub = _rt_df[[_rt_fw_col, _rt_rv_col, "variant_type"]].dropna()
                    _rt_total_alt = _rt_sub[_rt_fw_col] + _rt_sub[_rt_rv_col]
                    _rt_sub = _rt_sub[_rt_total_alt > 0].copy()
                    _rt_total_alt = _rt_total_alt[_rt_total_alt > 0]
                    if not _rt_sub.empty:
                        _rt_sub["frac_fwd"] = _rt_sub[_rt_fw_col] / _rt_total_alt
                        _rt_sub["read_type"] = _rt_lbl
                        _rt_strand_rows.append(_rt_sub[["frac_fwd", "variant_type", "read_type"]])

                if _rt_strand_rows:
                    _rt_strand_df = pd.concat(_rt_strand_rows, ignore_index=True)
                    _rt_strand_chart = (
                        alt.Chart(_rt_strand_df)
                        .transform_density(
                            "frac_fwd",
                            as_=["frac_fwd", "density"],
                            groupby=["read_type"],
                            extent=[0, 1],
                            steps=80,
                        )
                        .mark_area(opacity=0.4)
                        .encode(
                            x=alt.X("frac_fwd:Q", title="Fraction of alt reads on forward strand",
                                    scale=alt.Scale(domain=[0, 1])),
                            y=alt.Y("density:Q", title="Density", stack=None),
                            color=alt.Color("read_type:N", title="Read type"),
                            tooltip=[
                                alt.Tooltip("frac_fwd:Q", format=".3f", title="Frac forward"),
                                "read_type:N",
                            ],
                        )
                        .properties(height=250)
                    )
                    st.altair_chart(_rt_strand_chart, width="stretch")
                else:
                    st.info("No strand balance data available under current filters.")

                st.divider()

                # ── View 5: SBS96 spectrum side-by-side ────────────────────────
                st.subheader("SBS96 error spectrum")

                if not _has_data("trinuc_context"):
                    st.info(
                        "SBS96 spectrum requires the `trinuc_context` column. "
                        "Re-run `geac collect` with a reference FASTA."
                    )
                else:
                    def _rt_sbs96(rt_where: str):
                        _raw = con.execute(f"""
                            SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
                            FROM {table_expr}
                            WHERE {rt_where}
                              AND variant_type = 'SNV'
                              AND trinuc_context IS NOT NULL
                              AND length(trinuc_context) = 3
                            GROUP BY trinuc_context, ref_allele, alt_allele
                        """).df()
                        return _to_spec96_strat(_raw)

                    _rt_s96_a, _rt_total_a = _rt_sbs96(_rt_wa)
                    _rt_s96_b, _rt_total_b = _rt_sbs96(_rt_wb)

                    if _rt_s96_a is None or _rt_s96_b is None:
                        st.info("Insufficient SNV data for spectrum comparison.")
                    else:
                        _rt_sc1, _rt_sc2 = st.columns(2)
                        with _rt_sc1:
                            st.markdown(f"**{_rt_a}** ({_rt_total_a:,} SNVs)")
                            st.altair_chart(
                                _strat_sbs96_chart(_rt_s96_a, str(_rt_a)),
                                width="stretch",
                            )
                        with _rt_sc2:
                            st.markdown(f"**{_rt_b}** ({_rt_total_b:,} SNVs)")
                            st.altair_chart(
                                _strat_sbs96_chart(_rt_s96_b, str(_rt_b)),
                                width="stretch",
                            )

                st.divider()

                # ── View 6: Loci unique to one read type ───────────────────────
                st.subheader("Loci unique to one read type")
                st.caption(
                    "These loci are called in one read type but absent in the other. "
                    "Loci only in the less-filtered read type are candidates for "
                    "consensus-removed artefacts."
                )

                _rt_uniq = _rt_df[_rt_df["concordance"] != "shared"].copy()
                _rt_uniq["read_type"] = _rt_uniq["concordance"].map(
                    {"only_a": str(_rt_a), "only_b": str(_rt_b)}
                )
                _rt_uniq_filter = st.radio(
                    "Show",
                    ["Both", f"Only {_rt_a}", f"Only {_rt_b}"],
                    horizontal=True,
                    key="rt_cmp_uniq_filter",
                )
                if _rt_uniq_filter == f"Only {_rt_a}":
                    _rt_uniq_show = _rt_uniq[_rt_uniq["concordance"] == "only_a"]
                elif _rt_uniq_filter == f"Only {_rt_b}":
                    _rt_uniq_show = _rt_uniq[_rt_uniq["concordance"] == "only_b"]
                else:
                    _rt_uniq_show = _rt_uniq

                _rt_uniq_cols = [c for c in [
                    "read_type", "sample_id", "chrom", "pos", "alt_allele",
                    "variant_type", "vaf_a", "vaf_b", "depth_a", "depth_b",
                ] if c in _rt_uniq_show.columns]

                if _rt_uniq_show.empty:
                    st.success("No unique loci for current filter selection.")
                else:
                    st.dataframe(
                        _rt_uniq_show[_rt_uniq_cols]
                        .rename(columns={
                            "vaf_a":   f"vaf_{_rt_a}",
                            "vaf_b":   f"vaf_{_rt_b}",
                            "depth_a": f"depth_{_rt_a}",
                            "depth_b": f"depth_{_rt_b}",
                        })
                        .sort_values(["read_type", "chrom", "pos"]),
                        width="stretch",
                        hide_index=True,
                        key="rt_cmp_uniq_table",
                    )
