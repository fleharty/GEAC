import io
import json
import os
import textwrap
import time
import traceback
import warnings
import zipfile
from contextlib import contextmanager
from pathlib import Path

try:
    import anthropic as _anthropic
    _HAS_ANTHROPIC = True
except ImportError:
    _HAS_ANTHROPIC = False
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
    fit_cosmic_cohort_greedy,
    fit_cosmic_single_sample_greedy,
    fit_sbs_nmf,
)
from pipeline_compare_helpers import (
    build_unique_pipeline_characterization_df,
    summarize_unique_pipeline_groups,
)
from read_context_helpers import (
    add_read_context_fraction_metrics,
    compute_locus_n_asymmetry,
    compute_locus_n_asymmetry_precomputed,
)

# Altair emits a spurious "Automatically deduplicated selection parameter" UserWarning
# when a shared cross-panel selection param is present in multiple sub-charts.  The
# deduplication it performs is intentionally correct behaviour; suppress the noise.
warnings.filterwarnings(
    "ignore",
    message="Automatically deduplicated",
    category=UserWarning,
)

from igv_helpers import (
    make_vcf,
    query_distinct_samples,
    per_read_warning_note,
    resolve_index_uri,
)
import geac_config
from explorer import (
    GEAC_VERSION,
    MAIN_FILTER_KEYS,
    MAIN_FILTER_STATE,
    DataSource,
    PerReadFilters,
    LociThresholds,
    TabContext,
)
from explorer.data_source import cached_distinct_values, sort_chroms, parse_region, _sql_str
from explorer.main_table import (
    render_position_drilldown,
    render_records_table,
    render_summary_metrics,
)
from explorer.tabs import (
    TAB_MODULES,
    TAB_SUMMARY,
    TAB_VAF_DISTRIBUTION,
    TAB_ERROR_SPECTRUM,
    TAB_STRAND_BIAS,
    TAB_COHORT,
    TAB_READS,
    TAB_DUPLEX_SIMPLEX,
    TAB_TUMOR_NORMAL,
    TAB_PANEL_OF_NORMALS,
    TAB_PIPELINE_COMPARISON,
    TAB_READ_TYPE_COMPARISON,
    TAB_OVERLAP_AGREEMENT,
    TAB_MNV_CANDIDATES,
    TAB_AI_PLOT_BUILDER,
    TAB_FRAGMENTOMICS,
)

_IS_MIN, _IS_MAX = 20, 500  # insert size slider bounds


def _init_session_state(key: str, default):
    if key not in st.session_state:
        st.session_state[key] = default
    return st.session_state[key]


def _optional_distinct_multiselect(
    *,
    data_key: str,
    cache_key: str,
    label: str,
    key: str,
    caption: str | None = None,
) -> list[str]:
    if not _has_data(data_key):
        if caption:
            st.sidebar.caption(caption)
        return []
    values = _init_session_state(
        cache_key,
        cached_distinct_values(con, table_expr, data_key, cache_key=cache_key),
    )
    _init_session_state(key, [])
    return st.sidebar.multiselect(label, values, key=key)


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

# ── Debug timing ──────────────────────────────────────────────────────────────
# Enabled via the sidebar checkbox "Show query timings". All @_timed("label")
# calls accumulate into _TIMINGS for the current rerun; the table is rendered
# at the very end of the script when the flag is on.
if "_debug_timings" not in st.session_state:
    st.session_state["_debug_timings"] = []
_TIMINGS: list[tuple[str, float]] = st.session_state["_debug_timings"]
_TIMINGS.clear()  # reset at the start of every rerun


@contextmanager
def _timed(label: str):
    """Context manager that records wall-clock milliseconds for *label*."""
    t0 = time.perf_counter()
    try:
        yield
    finally:
        _TIMINGS.append((label, (time.perf_counter() - t0) * 1000))


st.set_page_config(page_title="GEAC Explorer", layout="wide")

# Preserve scroll position across reruns (sidebar filter changes reset scroll otherwise).
import streamlit.components.v1 as _stc
_stc.html(
    """
    <script>
    (function () {
        var KEY = 'geac_scroll_y';
        try {
            var win = window.parent;
            var doc = win.document;
            var main = doc.querySelector('section.main') ||
                       doc.querySelector('[data-testid="stMain"]') ||
                       doc.querySelector('.main');
            if (!main) return;

            // Save scroll position on every scroll event (attach once per element).
            if (!main._geac_scroll_attached) {
                main._geac_scroll_attached = true;
                main.addEventListener('scroll', function () {
                    win.localStorage.setItem(KEY, main.scrollTop);
                }, { passive: true });
            }

            // Restore after Streamlit finishes updating the DOM.
            var saved = parseInt(win.localStorage.getItem(KEY) || '0', 10);
            if (saved > 0) {
                win.requestAnimationFrame(function () {
                    win.requestAnimationFrame(function () {
                        main.scrollTop = saved;
                    });
                });
            }
        } catch (e) {}
    })();
    </script>
    """,
    height=0,
)

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
_has_sample_metrics = data_source.has_optional_table("sample_metrics")
_has_fragments = data_source.has_optional_table("fragments")

# Extract per-sample resource paths and cohort-level gnomAD path embedded in the DuckDB.
# These are populated by `geac collect` when the new bam_path/gnomad_path columns are present.
# Falls back gracefully to empty results for Parquet sources or older DuckDB files.
_db_resources = data_source.sample_resources()
_db_gnomad_paths = data_source.embedded_gnomad_paths()
_embedded_gnomad = _db_gnomad_paths[0] if len(_db_gnomad_paths) == 1 else ""
if "_alt_reads_cols_cached" not in st.session_state:
    st.session_state["_alt_reads_cols_cached"] = (
        set(con.execute("SELECT * FROM alt_reads LIMIT 0").df().columns)
        if _has_alt_reads
        else set()
    )
_alt_reads_cols = st.session_state["_alt_reads_cols_cached"]


def _has_alt_reads_cols(*cols: str) -> bool:
    return all(col in _alt_reads_cols for col in cols)


_alt_reads_has_pipeline  = _has_alt_reads_cols("pipeline")
_alt_reads_has_read_type = _has_alt_reads_cols("read_type")
_alt_reads_has_batch     = _has_alt_reads_cols("batch")
_alt_reads_has_label1    = _has_alt_reads_cols("label1")
_alt_reads_has_label2    = _has_alt_reads_cols("label2")
_alt_reads_has_label3    = _has_alt_reads_cols("label3")
_alt_reads_has_timepoint = _has_alt_reads_cols("timepoint")
_alt_reads_has_frag_gc   = _has_alt_reads_cols("frag_gc")


def _alt_reads_meta_join(alias: str) -> str:
    parts = [f"ar.sample_id = {alias}.sample_id"]
    if _alt_reads_has_pipeline:
        parts.append(f"ar.pipeline = {alias}.pipeline")
    if _alt_reads_has_read_type:
        parts.append(f"ar.read_type = {alias}.read_type")
    return " AND ".join(parts)

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

_debug_mode = st.sidebar.checkbox(
    "Show query timings",
    value=False,
    key="debug_mode",
    help="Display a per-query wall-clock timing table directly below this checkbox after each rerun.",
)
# Placeholder reserved here so the timing table appears immediately below the
# checkbox rather than at the very bottom of the sidebar after all filters.
_timing_placeholder = st.sidebar.empty()

if "_cached_missing_required_cols" not in st.session_state:
    st.session_state["_cached_missing_required_cols"] = data_source.required_columns_missing()
_missing_required_cols = st.session_state["_cached_missing_required_cols"]
if _missing_required_cols:
    st.warning(
        "This dataset is missing required `alt_bases` columns expected by the current Explorer: "
        + ", ".join(sorted(_missing_required_cols)),
        icon="⚠️",
    )

# ── Summary stats ─────────────────────────────────────────────────────────────
if "_cached_summary_stats" not in st.session_state:
    st.session_state["_cached_summary_stats"] = data_source.summary_stats()
stats = st.session_state["_cached_summary_stats"]

n_annotated = int(stats["n_annotated"][0])
n_called    = int(stats["n_called"][0])
pct_called  = f"{100 * n_called / n_annotated:.1f}%" if n_annotated > 0 else "N/A"

# ── Filters (sidebar) ─────────────────────────────────────────────────────────
_FILTER_KEYS = list(MAIN_FILTER_KEYS)

_sidebar_logo = _LOGO_COMPACT if _LOGO_COMPACT.exists() else _LOGO
if _sidebar_logo.exists():
    st.sidebar.image(str(_sidebar_logo), width="stretch")

if "_cached_chroms" not in st.session_state:
    st.session_state["_cached_chroms"] = sort_chroms(data_source.distinct_values("chrom"))
chroms = st.session_state["_cached_chroms"]

if "_cached_samples" not in st.session_state:
    st.session_state["_cached_samples"] = data_source.distinct_values("sample_id")
samples = st.session_state["_cached_samples"]

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

chrom_sel = st.sidebar.text_input(
    "Region or gene",
    placeholder="chr1 · 1 · chr1:100000 · 1:100000-200000 · TP53",
    key="chrom_sel",
)
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
# Cache result — column presence doesn't change during a session.
if "_cached_optional_cols" not in st.session_state:
    st.session_state["_cached_optional_cols"] = data_source.has_non_null_batch([
        "batch", "label1", "label2", "label3", "timepoint", "gene",
        "variant_filter", "gnomad_af", "homopolymer_len",
        "trinuc_context", "on_target", "variant_called",
        "overlap_alt_agree",
    ])
_optional_cols = st.session_state["_cached_optional_cols"]

def _has_data(col: str) -> bool:
    """True iff col exists in the schema AND has at least one non-null value."""
    return _optional_cols.get(col, data_source.has_non_null(col))


def _render_gnomad_filters() -> tuple[tuple[str, str], bool]:
    if not _has_data("gnomad_af"):
        return ("0", "1.0"), True

    # Streamlit needs an explicit tuple value to render select_slider in range mode.
    _af_ss = st.session_state.get("gnomad_af_range", ("0", "1.0"))
    if not (isinstance(_af_ss, (tuple, list)) and len(_af_ss) == 2):
        _af_ss = ("0", "1.0")
    gnomad_af_range = st.sidebar.select_slider(
        "gnomAD AF (log scale)",
        options=_GNOMAD_AF_STEPS,
        value=tuple(_af_ss),
        help="Filter by gnomAD allele frequency. Steps are logarithmic.",
    )
    st.session_state["gnomad_af_range"] = gnomad_af_range
    gnomad_include_null = st.sidebar.checkbox(
        "Include sites absent from gnomAD",
        key="gnomad_include_null",
    )
    return gnomad_af_range, gnomad_include_null


def _render_repeat_filters() -> tuple[tuple[int, int], tuple[int, int]]:
    if not _has_data("homopolymer_len"):
        st.sidebar.caption("Repeat filters unavailable — run geac collect with a newer build to enable.")
        return (0, 20), (0, 50)
    return (
        st.sidebar.slider("Homopolymer length range", 0, 20, step=1, key="homopolymer_range"),
        st.sidebar.slider("STR length range", 0, 50, step=1, key="str_len_range"),
    )


def _per_read_filter_limits() -> tuple[object, int, int, bool, int, bool, bool, bool, int]:
    if "_cached_fs_max" not in st.session_state:
        _n_total_expr = (
            "COALESCE(MAX(n_n_before_alt + n_n_after_alt), 0)"
            if _has_alt_reads_cols("n_n_before_alt", "n_n_after_alt")
            else "0"
        )
        _gc_has_data_expr = (
            "COUNT(frag_gc) > 0" if _alt_reads_has_frag_gc else "false"
        )
        _reads_maxes = con.execute(f"""
            SELECT
                MAX(family_size),
                COALESCE(MAX(cycle), 300),
                COALESCE(MAX(map_qual), 60),
                COUNT(insert_size) > 0,
                {_n_total_expr},
                {_gc_has_data_expr}
            FROM alt_reads
        """).fetchone()
        st.session_state["_cached_fs_max"]     = _reads_maxes[0]   # None if all NULL
        st.session_state["_cached_cycle_max"]   = int(_reads_maxes[1])
        st.session_state["_cached_mq_max"]     = int(_reads_maxes[2])
        st.session_state["_cached_is_has_data"] = bool(_reads_maxes[3])
        st.session_state["_cached_n_total_max"] = int(_reads_maxes[4])
        st.session_state["_cached_gc_has_data"] = bool(_reads_maxes[5])

    _fs_max_raw = st.session_state["_cached_fs_max"]
    _cycle_max = st.session_state["_cached_cycle_max"]
    _mq_max = st.session_state["_cached_mq_max"]
    _is_has_data = st.session_state["_cached_is_has_data"]
    _n_total_max = st.session_state.get("_cached_n_total_max", 0)
    _gc_has_data = st.session_state.get("_cached_gc_has_data", False)
    _n_total_has_data = _has_alt_reads_cols("n_n_before_alt", "n_n_after_alt")
    _fs_has_data = _fs_max_raw is not None
    _fs_max = int(_fs_max_raw) if _fs_has_data else 0
    return (
        _fs_max_raw,
        _cycle_max,
        _mq_max,
        _is_has_data,
        _n_total_max,
        _gc_has_data,
        _n_total_has_data,
        _fs_has_data,
        _fs_max,
    )


def _append_active_range_filter(
    conditions: list[str],
    *,
    enabled: bool,
    sql_expr: str,
    lo,
    hi,
    default_lo,
    default_hi,
) -> None:
    if enabled and (lo > default_lo or hi < default_hi):
        conditions.append(f"{sql_expr} BETWEEN {lo} AND {hi}")


batch_sel = _optional_distinct_multiselect(
    data_key="batch",
    cache_key="_cached_batches",
    label="Batch (blank = all)",
    key="batch_sel",
)

_pipelines = cached_distinct_values(con, table_expr, "pipeline", cache_key="_cached_pipelines")
_init_session_state("pipeline_sel", [])
pipeline_sel = st.sidebar.multiselect("Pipeline (blank = all)", _pipelines, key="pipeline_sel")

label1_sel = _optional_distinct_multiselect(
    data_key="label1",
    cache_key="_cached_label1_vals",
    label="Label 1 (blank = all)",
    key="label1_sel",
)
label2_sel = _optional_distinct_multiselect(
    data_key="label2",
    cache_key="_cached_label2_vals",
    label="Label 2 (blank = all)",
    key="label2_sel",
)
label3_sel = _optional_distinct_multiselect(
    data_key="label3",
    cache_key="_cached_label3_vals",
    label="Label 3 (blank = all)",
    key="label3_sel",
)
timepoint_sel = _optional_distinct_multiselect(
    data_key="timepoint",
    cache_key="_cached_timepoint_vals",
    label="Timepoint (blank = all)",
    key="timepoint_sel",
)

_genes_available = _has_data("gene")
if "variant_sel" not in st.session_state:
    st.session_state["variant_sel"] = ["SNV", "insertion", "deletion"]
variant_sel = st.sidebar.multiselect(
    "Variant type",
    ["SNV", "insertion", "deletion"],
    key="variant_sel",
)
indel_len_range = st.sidebar.slider(
    "Indel length range (bp)",
    min_value=1, max_value=500, step=1,
    key="indel_len_range",
    help="Filter insertions/deletions by length (bp). Length is derived from `alt_allele` (`+`/`-` prefix excluded). SNVs are unaffected — use Variant type to exclude them.",
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
    _vf_options = _init_session_state(
        "_cached_vf_options",
        cached_distinct_values(con, table_expr, "variant_filter", cache_key="_cached_vf_options"),
    )
    _init_session_state("variant_filter_sel", [])
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
gnomad_af_range, gnomad_include_null = _render_gnomad_filters()

_repeat_cols_present = _has_data("homopolymer_len")
homopolymer_range, str_len_range = _render_repeat_filters()
min_depth = st.sidebar.number_input("Min depth (0 = no minimum)", min_value=0, step=1, key="min_depth")
max_depth = st.sidebar.number_input("Max depth (0 = no maximum)", min_value=0, step=1, key="max_depth")


# ── Per-read filters (only when alt_reads table is present) ───────────────────
_reads_conditions = []
if _has_alt_reads:
    (
        _fs_max_raw,
        _cycle_max,
        _mq_max,
        _is_has_data,
        _n_total_max,
        _gc_has_data,
        _n_total_has_data,
        _fs_has_data,
        _fs_max,
    ) = _per_read_filter_limits()

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

    if _n_total_has_data:
        n_in_read_range = st.sidebar.slider(
            "N in read range",
            min_value=0,
            max_value=max(0, _n_total_max),
            value=(0, max(0, _n_total_max)),
            step=1,
            key="n_in_read_range",
            help="Filter alt-supporting reads by total tracked `N` bases in the stored read sequence, "
                 "defined as `n_n_before_alt + n_n_after_alt`.",
        )
    else:
        n_in_read_range = (0, 0)

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

    if _gc_has_data:
        frag_gc_range = st.sidebar.slider(
            "Fragment GC range",
            min_value=0.0, max_value=1.0,
            value=(0.0, 1.0), step=0.01,
            key="frag_gc_range",
            help="Filter alt-supporting reads by fragment GC fraction "
                 "(GC of the reference span from fragment start to fragment start + |TLEN|).",
        )
        _gc_lo, _gc_hi = frag_gc_range
        if _gc_lo == 0.0 and _gc_hi == 1.0:
            st.sidebar.caption(
                "Fragment GC: no filter active — reads with any GC (including those "
                "with no fragment GC available) are accepted."
            )
        else:
            st.sidebar.caption(
                f"Fragment GC: keeping only reads with GC between {_gc_lo:.2f} and {_gc_hi:.2f}. "
                "Reads with no fragment GC (unpaired) are excluded."
            )
    else:
        frag_gc_range = (0.0, 1.0)

    read_strand_sel = st.sidebar.radio(
        "Read",
        ["All", "R1 only", "R2 only"],
        horizontal=True,
        key="read_strand_sel",
        help="Filter to R1 reads (BAM flag 0x40 set), R2 reads, or all reads.",
    )

    if _n_total_has_data:
        min_delta_n_frac = st.sidebar.slider(
            "Min N-asymmetry score (locus discovery)",
            min_value=0.0,
            max_value=1.0,
            value=st.session_state.get("min_delta_n_frac", 0.0),
            step=0.05,
            key="min_delta_n_frac",
            help="In the Read Context tab, show only loci where mean(frac_N_after) − mean(frac_N_before) "
                 "is at least this value. Does not affect other tabs.",
        )
    else:
        min_delta_n_frac = 0.0

    _fs_lo, _fs_hi = family_size_range
    _cycle_lo, _cycle_hi = cycle_range
    _mq_lo, _mq_hi = map_qual_range
    _n_total_lo, _n_total_hi = n_in_read_range
    _is_lo, _is_hi = insert_size_range
    _gc_lo, _gc_hi = frag_gc_range

    _append_active_range_filter(
        _reads_conditions,
        enabled=_fs_has_data,
        sql_expr="family_size",
        lo=_fs_lo,
        hi=_fs_hi,
        default_lo=0,
        default_hi=_fs_max,
    )
    _append_active_range_filter(
        _reads_conditions,
        enabled=True,
        sql_expr="cycle",
        lo=_cycle_lo,
        hi=_cycle_hi,
        default_lo=1,
        default_hi=_cycle_max,
    )
    _append_active_range_filter(
        _reads_conditions,
        enabled=True,
        sql_expr="map_qual",
        lo=_mq_lo,
        hi=_mq_hi,
        default_lo=0,
        default_hi=_mq_max,
    )
    _append_active_range_filter(
        _reads_conditions,
        enabled=_n_total_has_data,
        sql_expr="(COALESCE(n_n_before_alt, 0) + COALESCE(n_n_after_alt, 0))",
        lo=_n_total_lo,
        hi=_n_total_hi,
        default_lo=0,
        default_hi=_n_total_max,
    )
    _append_active_range_filter(
        _reads_conditions,
        enabled=_is_has_data,
        sql_expr="insert_size",
        lo=_is_lo,
        hi=_is_hi,
        default_lo=_IS_MIN,
        default_hi=_IS_MAX,
    )
    _append_active_range_filter(
        _reads_conditions,
        enabled=_gc_has_data,
        sql_expr="frag_gc",
        lo=_gc_lo,
        hi=_gc_hi,
        default_lo=0.0,
        default_hi=1.0,
    )

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
    n_in_read_range = (0, 0)
    insert_size_range = (0, 0)
    frag_gc_range = (0.0, 1.0)
    read_strand_sel = "All"
    recompute_vaf = False
    _fs_has_data = False
    _is_has_data = False
    _gc_has_data = False
    _n_total_has_data = False
    _fs_max = 0
    _cycle_max = 1
    _mq_max = 0
    _n_total_max = 0
    _fs_lo = _fs_hi = _cycle_lo = _cycle_hi = _mq_lo = _mq_hi = _n_total_lo = _n_total_hi = _is_lo = _is_hi = 0
    _gc_lo, _gc_hi = 0.0, 1.0

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
    if len(_db_gnomad_paths) > 1:
        st.sidebar.warning(
            "Multiple gnomAD paths found in database; set 'gnomAD track' explicitly or via geac.toml.",
            icon="⚠️",
        )
    gnomad_track = st.sidebar.text_input(
        "gnomAD track (optional)",
        value=_cfg.get("gnomad_track") or _embedded_gnomad,
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

    Relative paths in the BAM/index columns are resolved against the manifest
    file's directory so that IGV session files (written to a temp dir) receive
    absolute paths that IGV can actually open.
    """
    manifest_dir = _os.path.dirname(_os.path.abspath(p.strip()))

    def _abs(val: str | None) -> str | None:
        if val is None:
            return None
        if not _os.path.isabs(val) and not val.startswith(("gs://", "s3://", "https://", "http://")):
            return _os.path.normpath(_os.path.join(manifest_dir, val))
        return val

    mdf = pd.read_csv(p.strip(), sep="\t")
    _has_pipeline_col = "pipeline" in mdf.columns
    result = {}
    for row in mdf.itertuples(index=False):
        bai = str(row.duplex_output_bam_index) if hasattr(row, "duplex_output_bam_index") and pd.notna(row.duplex_output_bam_index) else None
        variants = str(row.final_annotated_variants) if hasattr(row, "final_annotated_variants") and pd.notna(row.final_annotated_variants) else None
        entry = {"bam": _abs(str(row.duplex_output_bam)), "bai": _abs(bai), "variants_tsv": variants}
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

# Merge: start from embedded DB paths, then overlay the user-provided manifest.
# Per-sample entries from manifest.tsv override entries from DuckDB.
if _db_resources:
    manifest = {**_db_resources, **manifest}
    if not manifest_path or not manifest_path.strip():
        st.sidebar.caption(f"{len(_db_resources):,} samples resolved from database")

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
def launch_igv_session(session_xml: str, vcf: str, sort_locus: str = "") -> str:
    """Write session.xml and positions.vcf to a temp dir and load in IGV.

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
    vcf_path     = _os.path.join(tmp, "positions.vcf")
    with open(session_path, "w") as f:
        f.write(session_xml)
    with open(vcf_path, "w") as f:
        f.write(vcf)

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
    force_pipelines: list[str] | None = None,
) -> tuple[str, str]:
    """Build an IGV session XML and return (xml, sort_locus).

    *sort_locus* is ``"chrom:pos"`` (1-based) for the first locus — used by
    ``launch_igv_session`` to sort reads by base via the REST API.  IGV's XML
    session format does not support sort-on-load, so the XML itself contains
    no sort directives.

    *force_pipelines* — when provided, emit one BAM track per (sample_id, pipeline)
    for every pipeline in the list, regardless of which pipelines appear in *df*.
    Use this in pipeline-comparison views where unique loci only exist in one
    pipeline's rows but you still want both BAMs loaded for comparison.
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

    # Build (sample_id, pipeline, track_label) tuples.
    # force_pipelines: always emit one track per (sample, pipeline) for the given
    # pipelines — used when igv_df only contains rows from one pipeline (e.g.
    # unique loci) but the user wants both BAMs loaded for side-by-side comparison.
    if force_pipelines is not None:
        _sample_ids = sorted(df["sample_id"].dropna().unique().tolist(), key=str)
        _track_items = [
            (sid, pipe, f"{sid} ({pipe})")
            for sid in _sample_ids
            for pipe in force_pipelines
        ]
    elif "pipeline" in df.columns and df["pipeline"].nunique() > 1:
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
            # For local paths, only emit index= if the file actually exists.
            # A wrong explicit index overrides IGV's auto-detection and silently
            # produces an empty track; omitting it lets IGV find the index itself.
            _bai_is_cloud = bai and bai.startswith(("gs://", "s3://", "https://", "http://"))
            _bai_ok = bai and (_bai_is_cloud or _os.path.isfile(bai))
            index_attr = f' index="{bai}"' if _bai_ok else ""
            resources.append(f'        <Resource path="{bam}" name="{label}"{index_attr}/>')
            tracks.append(f'        <Track id="{bam}" name="{label}"/>')

    resources.append('        <Resource path="positions.vcf" name="Selected positions"/>')
    tracks.append('        <Track id="positions.vcf" name="Selected positions" color="255,0,0" height="40"/>')

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
    available_samples: list[str] | None = None,
    force_pipelines: list[str] | None = None,
):
    """Render IGV Prepare + Download buttons with chunked progress.

    extra_conditions   — SQL WHERE fragments
    display_df         — already-fetched display DataFrame
    key                — unique widget key prefix
    use_global_filters — when True (default), prepend the sidebar ``conditions``
                         to *extra_conditions*; when False, use *extra_conditions*
                         alone (e.g. for the position drill-down, which shows all
                         samples at the locus regardless of sidebar filters).
    available_samples  — when provided, used as the multiselect options instead of
                         querying the DB filtered to extra_conditions (useful when
                         extra_conditions encodes only a subset of loci, e.g. the
                         first 500 of many, but the full sample list is known).
    force_pipelines    — when provided, passed to make_igv_session so that one BAM
                         track is emitted per (sample, pipeline) for every listed
                         pipeline, even if igv_df only contains rows from one
                         pipeline (e.g. unique loci in a pipeline-comparison view).
    """
    if not manifest:
        st.caption("🧭 Add a manifest in the sidebar to enable IGV session download.")
        return

    _where_parts = (conditions + extra_conditions) if use_global_filters else extra_conditions
    _extra_w = " AND ".join(_where_parts) if _where_parts else "TRUE"
    sample_ids = available_samples if available_samples is not None else query_distinct_samples(con, table_expr, _extra_w)
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

    # Always show the multiselect so the user can see and control which samples
    # are included regardless of whether n is above or below the cap.  When the
    # multiselect was inside the n > IGV_CAP block, filter changes that reduced n
    # below the cap silently reset cap_samples to all current samples, causing
    # positions from unintended samples to appear in the BED.
    chosen = st.multiselect(
        "Samples to include in IGV session",
        options=sample_ids,
        default=cap_samples,
        key=f"{key}_sample_pick",
    )
    cap_samples = chosen if chosen else cap_samples

    if n > IGV_CAP:
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
        vcf = make_vcf(igv_df)
        session, sort_locus = make_igv_session(
            igv_df,
            manifest,
            genome,
            target_regions,
            gnomad_track,
            gnomad_track_index,
            force_pipelines=force_pipelines,
        )
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("session.xml", session)
            zf.writestr("positions.vcf", vcf)
        st.session_state[f"{key}_igv"]        = buf.getvalue()
        st.session_state[f"{key}_session"]    = session
        st.session_state[f"{key}_vcf"]        = vcf
        st.session_state[f"{key}_sort_locus"] = sort_locus

        if auto_launch_igv:
            msg = launch_igv_session(session, vcf, sort_locus)
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
# threshold_conditions: subset of conditions[] that compare per-pipeline count
# columns (alt_count, vaf, depth, fwd/rev alt, overlap).  The pipeline comparison
# tab strips these from its CTEs and applies them post-join so that a locus with
# sub-threshold evidence in one pipeline is not misclassified as absent.
threshold_conditions: list[str] = []
conditions = []
if min_alt > 1:
    _c = f"alt_count >= {min_alt}"
    conditions.append(_c)
    threshold_conditions.append(_c)
if vaf_range != (0.0, 1.0):
    _c = f"alt_count * 1.0 / total_depth BETWEEN {vaf_range[0]} AND {vaf_range[1]}"
    conditions.append(_c)
    threshold_conditions.append(_c)
if max_alt > 0:
    _c = f"alt_count <= {max_alt}"
    conditions.append(_c)
    threshold_conditions.append(_c)
_region = parse_region(chrom_sel, chroms)
if _region.chrom is not None:
    conditions.append(f"chrom = '{_sql_str(_region.chrom)}'")
if _region.start is not None:
    if _region.end is not None:
        conditions.append(f"pos BETWEEN {_region.start} AND {_region.end}")
    else:
        conditions.append(f"pos = {_region.start}")
if _region.gene is not None and _genes_available:
    conditions.append(f"gene = '{_sql_str(_region.gene)}'")
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
    if timepoint_sel:
        _scope_parts.append("timepoint IN ({})".format(
            ", ".join(f"'{_sql_str(v)}'" for v in timepoint_sel)))
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
if timepoint_sel:
    t_list = ", ".join(f"'{_sql_str(v)}'" for v in timepoint_sel)
    conditions.append(f"timepoint IN ({t_list})")
if variant_sel:
    t_list = ", ".join(f"'{_sql_str(t)}'" for t in variant_sel)
    conditions.append(f"variant_type IN ({t_list})")
if indel_len_range != (1, 500):
    _ilo, _ihi = indel_len_range
    conditions.append(
        "(variant_type = 'SNV' "
        f"OR LENGTH(alt_allele) - 1 BETWEEN {_ilo} AND {_ihi})"
    )
if min_fwd_alt > 0:
    _c = f"fwd_alt_count >= {min_fwd_alt}"
    conditions.append(_c)
    threshold_conditions.append(_c)
if min_rev_alt > 0:
    _c = f"rev_alt_count >= {min_rev_alt}"
    conditions.append(_c)
    threshold_conditions.append(_c)
if min_overlap_agree > 0:
    _c = f"overlap_alt_agree >= {min_overlap_agree}"
    conditions.append(_c)
    threshold_conditions.append(_c)
if min_overlap_disagree > 0:
    _c = f"overlap_alt_disagree >= {min_overlap_disagree}"
    conditions.append(_c)
    threshold_conditions.append(_c)
if min_depth > 0:
    _c = f"total_depth >= {min_depth}"
    conditions.append(_c)
    threshold_conditions.append(_c)
if max_depth > 0:
    _c = f"total_depth <= {max_depth}"
    conditions.append(_c)
    threshold_conditions.append(_c)
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

if _repeat_cols_present:
    if homopolymer_range != (0, 20):
        conditions.append(f"(homopolymer_len IS NULL OR homopolymer_len BETWEEN {homopolymer_range[0]} AND {homopolymer_range[1]})")
    if str_len_range != (0, 50):
        conditions.append(f"(str_len IS NULL OR str_len BETWEEN {str_len_range[0]} AND {str_len_range[1]})")

where = " AND ".join(conditions) if conditions else "TRUE"

# ── Shared alt_reads join subquery (used by Duplex/Simplex and Reads tabs) ────
# Joins alt_reads to the current filtered locus set so all reads plots respect
# the active sidebar filters.  Defined here so both tabs can share it.
_r_reads_filter = f"WHERE {_reads_where}" if _reads_active else ""
_r_filt_extra_cols = []
_r_join_extra_parts = []
if _alt_reads_has_pipeline:
    _r_filt_extra_cols.append("pipeline")
    _r_join_extra_parts.append("AND ar.pipeline = _filt.pipeline")
if _alt_reads_has_read_type:
    _r_filt_extra_cols.append("read_type")
    _r_join_extra_parts.append("AND ar.read_type = _filt.read_type")
_r_filt_extra_sql = "".join(f", {col}" for col in _r_filt_extra_cols)
_r_join_extra_sql = "\n    ".join(_r_join_extra_parts)
_r_join = f"""
    (SELECT * FROM alt_reads {_r_reads_filter}) ar
    INNER JOIN (
        SELECT DISTINCT sample_id, chrom, pos, alt_allele{_r_filt_extra_sql}
        FROM {table_expr}
        WHERE {where}
    ) _filt
    ON  ar.sample_id  = _filt.sample_id
    AND ar.chrom      = _filt.chrom
    AND ar.pos        = _filt.pos
    AND ar.alt_allele = _filt.alt_allele
    {_r_join_extra_sql}
"""


# ── Plots ─────────────────────────────────────────────────────────────────────
_MAIN_TAB_LABELS = [module.LABEL for module in TAB_MODULES]


def query_records(extra: list[str] = [], limit: int | None = None) -> pd.DataFrame:
    """Query records with current filters plus any extra conditions."""
    w = " AND ".join(conditions + extra)
    limit_clause = f"LIMIT {limit}" if limit is not None else ""
    return con.execute(f"""
        SELECT *, ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
               pos + 1 AS pos_display
        FROM {table_expr}
        WHERE {w}
        ORDER BY chrom, pos, alt_allele, sample_id
        {limit_clause}
    """).df()


def _display_table_cols(_df: pd.DataFrame) -> list[str]:
    return [
        c for c in [
            "sample_id", "chrom", "pos_display", "ref_allele", "alt_allele",
            "variant_type", "vaf", *( ["original_vaf"] if _reads_active else []), "alt_count", "ref_count", "total_depth",
            "fwd_alt_count", "rev_alt_count", "overlap_alt_agree",
            "overlap_alt_disagree", "variant_called", "variant_filter", "on_target", "gene", "gnomad_af",
        ]
        if c in _df.columns
    ]


_table_cols = _display_table_cols(query_records(limit=1))

with _timed("total_count"):
    total_count = con.execute(f"SELECT COUNT(*) FROM {table_expr} WHERE {where}").fetchone()[0]

if total_count == 0:
    st.warning("No records match the current filters.", icon="🔎")
    st.stop()

_active_main_tab = st.radio(
    "Main section",
    _MAIN_TAB_LABELS,
    key="main_tab_label",
    horizontal=True,
    label_visibility="collapsed",
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
    _region_active = bool(chrom_sel and chrom_sel.strip() and chrom_sel.strip().lower() != "all")
    _append_provenance_row(rows, "filters", "region_or_gene", chrom_sel.strip() or "All", active=_region_active)
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
    _append_provenance_row(rows, "filters", "timepoint", timepoint_sel, active=bool(timepoint_sel))
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
        "filters",
        "indel_length_range",
        f"{indel_len_range[0]}-{indel_len_range[1]}",
        active=indel_len_range != (1, 500),
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
        "n_in_read",
        f"{_n_total_lo}-{_n_total_hi}",
        active=_has_alt_reads and _n_total_has_data and (_n_total_lo > 0 or _n_total_hi < _n_total_max),
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
        "frag_gc",
        f"{_gc_lo:.2f}-{_gc_hi:.2f}",
        active=_has_alt_reads and _gc_has_data and (_gc_lo > 0.0 or _gc_hi < 1.0),
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


# ── Sidebar-filter SQL fragments against the fragments table ──────────────────
# Pre-built here so the Fragmentomics tab can append its own local clauses
# (motif, insert_size, sample selection) without reaching back into sidebar
# state. Mirrors `conditions` but uses `midpoint` instead of `pos`.
_fragments_where_parts: list[str] = []
if _region.chrom is not None:
    _fragments_where_parts.append(f"chrom = '{_sql_str(_region.chrom)}'")
if _region.start is not None and _region.end is not None:
    _fragments_where_parts.append(f"midpoint BETWEEN {_region.start} AND {_region.end}")
elif _region.start is not None:
    _fragments_where_parts.append(f"midpoint = {_region.start}")
if batch_sel:
    _fragments_where_parts.append(
        "batch IN ({})".format(", ".join(f"'{_sql_str(b)}'" for b in batch_sel))
    )
if label1_sel:
    _fragments_where_parts.append(
        "label1 IN ({})".format(", ".join(f"'{_sql_str(v)}'" for v in label1_sel))
    )
if label2_sel:
    _fragments_where_parts.append(
        "label2 IN ({})".format(", ".join(f"'{_sql_str(v)}'" for v in label2_sel))
    )
if label3_sel:
    _fragments_where_parts.append(
        "label3 IN ({})".format(", ".join(f"'{_sql_str(v)}'" for v in label3_sel))
    )
if timepoint_sel:
    _fragments_where_parts.append(
        "timepoint IN ({})".format(", ".join(f"'{_sql_str(v)}'" for v in timepoint_sel))
    )


# ── Tab dispatch context ──────────────────────────────────────────────────────
# Built once per render. Each tab's render(ctx) reads from this rather than
# reaching back into module globals — that's what makes tabs extractable.
ctx = TabContext(
    con=con,
    data_source=data_source,
    table_expr=table_expr,
    where=where,
    conditions=conditions,
    threshold_conditions=threshold_conditions,
    loci_thresholds=LociThresholds(
        min_alt=min_alt,
        max_alt=max_alt,
        vaf_lo=float(vaf_range[0]),
        vaf_hi=float(vaf_range[1]),
        min_depth=min_depth,
        max_depth=max_depth,
        min_fwd_alt=min_fwd_alt,
        min_rev_alt=min_rev_alt,
        min_overlap_agree=min_overlap_agree,
        min_overlap_disagree=min_overlap_disagree,
    ),
    schema_cols=_schema_cols,
    has_alt_reads=_has_alt_reads,
    has_normal_evidence=_has_normal_evidence,
    has_pon_evidence=_has_pon_evidence,
    stats=stats,
    pct_called=pct_called,
    total_count=total_count,
    table_cols=_table_cols,
    alt_reads_cols=_alt_reads_cols,
    cfg=_cfg,
    reads=PerReadFilters(
        active=_reads_active,
        recompute_vaf=recompute_vaf,
        read_strand_sel=read_strand_sel,
        fs_has_data=_fs_has_data,
        fs_lo=_fs_lo,
        fs_hi=_fs_hi,
        fs_max=_fs_max,
        cycle_lo=_cycle_lo,
        cycle_hi=_cycle_hi,
        cycle_max=_cycle_max,
        mq_lo=_mq_lo,
        mq_hi=_mq_hi,
        mq_max=_mq_max,
        n_total_has_data=_n_total_has_data,
        n_total_lo=_n_total_lo,
        n_total_hi=_n_total_hi,
        n_total_max=_n_total_max,
        is_has_data=_is_has_data,
        is_lo=_is_lo,
        is_hi=_is_hi,
        is_min=_IS_MIN,
        is_max=_IS_MAX,
        gc_has_data=_gc_has_data,
        gc_lo=_gc_lo,
        gc_hi=_gc_hi,
    ),
    timed=_timed,
    query_records=query_records,
    igv_buttons=igv_buttons,
    sql_str=_sql_str,
    has_data=_has_data,
    alt_reads_has=_has_alt_reads_cols,
    build_provenance=_build_active_filter_provenance,
    path=path,
    r_join=_r_join,
    reads_where=_reads_where if _reads_active else "",
    has_fragments=_has_fragments,
    fragments_where_parts=_fragments_where_parts,
)


if _active_main_tab == TAB_SUMMARY.LABEL:
    TAB_SUMMARY.render(ctx)

if _active_main_tab == TAB_VAF_DISTRIBUTION.LABEL:
    TAB_VAF_DISTRIBUTION.render(ctx)

if _active_main_tab == TAB_ERROR_SPECTRUM.LABEL:
    TAB_ERROR_SPECTRUM.render(ctx)

if _active_main_tab == TAB_STRAND_BIAS.LABEL:
    TAB_STRAND_BIAS.render(ctx)

if _active_main_tab == TAB_COHORT.LABEL:
    TAB_COHORT.render(ctx)


if _active_main_tab == TAB_READS.LABEL:
    TAB_READS.render(ctx)

if _active_main_tab == TAB_DUPLEX_SIMPLEX.LABEL:
    TAB_DUPLEX_SIMPLEX.render(ctx)

# ── Tumor/Normal tab ──────────────────────────────────────────────────────────
if _active_main_tab == TAB_TUMOR_NORMAL.LABEL:
    TAB_TUMOR_NORMAL.render(ctx)

# ── Panel of Normals tab ──────────────────────────────────────────────────────
if _active_main_tab == TAB_PANEL_OF_NORMALS.LABEL:
    TAB_PANEL_OF_NORMALS.render(ctx)

# ──────────────────────────────────────────────────────────────────────────────
# Tab 10 — Pipeline comparison
# ──────────────────────────────────────────────────────────────────────────────
if _active_main_tab == TAB_PIPELINE_COMPARISON.LABEL:
    TAB_PIPELINE_COMPARISON.render(ctx)

# ──────────────────────────────────────────────────────────────────────────────
# Tab 11 — Read-type comparison
# ──────────────────────────────────────────────────────────────────────────────
if _active_main_tab == TAB_READ_TYPE_COMPARISON.LABEL:
    TAB_READ_TYPE_COMPARISON.render(ctx)

# ── Overlap agreement ─────────────────────────────────────────────────────────
if _active_main_tab == TAB_OVERLAP_AGREEMENT.LABEL:
    TAB_OVERLAP_AGREEMENT.render(ctx)

# ── MNV candidates ─────────────────────────────────────────────────────────────
if _active_main_tab == TAB_MNV_CANDIDATES.LABEL:
    TAB_MNV_CANDIDATES.render(ctx)

# ── AI Plot Builder ────────────────────────────────────────────────────────────
if _active_main_tab == TAB_AI_PLOT_BUILDER.LABEL:
    TAB_AI_PLOT_BUILDER.render(ctx)

# ── Fragmentomics ───────────────────────────────────────────────────────────────────────
if _active_main_tab == TAB_FRAGMENTOMICS.LABEL:
    TAB_FRAGMENTOMICS.render(ctx)

# ── Query timing debug panel ───────────────────────────────────────────────────
# Write into the placeholder reserved directly below the checkbox so the table
# is always visible without scrolling the sidebar.
if _debug_mode and _TIMINGS:
    _timing_df = pd.DataFrame(_TIMINGS, columns=["query", "ms"])
    _timing_df = _timing_df.sort_values("ms", ascending=False).reset_index(drop=True)
    _total_ms = _timing_df["ms"].sum()
    with _timing_placeholder.container():
        st.caption(f"**Query timings — {_total_ms:,.0f} ms instrumented**")
        st.dataframe(
            _timing_df.style.format({"ms": "{:.1f}"}).bar(subset=["ms"], color="#f58518"),
            hide_index=True,
            use_container_width=True,
        )
else:
    _timing_placeholder.empty()
