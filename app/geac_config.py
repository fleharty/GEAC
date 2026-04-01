"""Load GEAC project configuration from a TOML file.

Discovery order:
1. Path passed via `streamlit run geac_explorer.py -- --config /path/to/geac.toml`
2. `geac.toml` in the current working directory
3. Empty config (all fields blank — user fills in manually)

Supported keys (all optional):

    data             = "/path/to/cohort.duckdb"   # or a .parquet file
    manifest         = "/path/to/manifest.tsv"
    cosmic           = "/path/to/COSMIC_v3.4_SBS_GRCh37.txt"
    genome_build     = "hg19"                     # hg19 | hg38 | mm10 | mm39 | other  (preferred)
    genome           = "hg19"                     # alias for genome_build (kept for compatibility)
    auto_launch_igv  = false                      # auto-launch IGV when a session is created
    target_regions   = "/path/to/targets.bed"     # BED or interval list added as a track in IGV sessions
    gnomad_track     = "/path/to/gnomad.vcf.gz"   # optional VCF/BCF added as a track in IGV sessions
    gnomad_track_index = "/path/to/gnomad.vcf.gz.tbi"  # optional explicit index path for the gnomAD track
"""

from __future__ import annotations

import os
import sys

# tomllib is stdlib in Python 3.11+; fall back to tomli for earlier versions.
try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib  # type: ignore[no-redef]
    except ImportError:
        tomllib = None  # type: ignore[assignment]


_STRING_KEYS = {
    "data",
    "manifest",
    "cosmic",
    "genome",
    "genome_build",
    "target_regions",
    "gnomad_track",
    "gnomad_track_index",
}
_PATH_KEYS = {
    "data",
    "manifest",
    "cosmic",
    "target_regions",
    "gnomad_track",
    "gnomad_track_index",
}
_BOOL_KEYS   = {"auto_launch_igv"}
_KNOWN_KEYS  = _STRING_KEYS | _BOOL_KEYS
_URI_PREFIXES = ("gs://", "http://", "https://")


def load() -> dict:
    """Return config dict for known keys.

    String keys are returned as str; boolean keys are returned as bool.
    Missing keys are absent from the returned dict (not None), so callers
    can use .get("data", "") or .get("auto_launch_igv", False) for safe defaults.
    """
    path = _find_config_path()
    if path is None:
        return {}

    if tomllib is None:
        # Can't parse TOML without the library — surface a clear message.
        import streamlit as st
        st.warning(
            f"Found {path} but could not parse it: install `tomli` "
            "(`pip install tomli`) or upgrade to Python 3.11+."
        )
        return {}

    try:
        with open(path, "rb") as fh:
            raw = tomllib.load(fh)
    except Exception as exc:
        import streamlit as st
        st.warning(f"Could not read config file {path}: {exc}")
        return {}

    unknown = set(raw) - _KNOWN_KEYS
    if unknown:
        import streamlit as st
        st.warning(f"geac.toml: unknown key(s) ignored: {', '.join(sorted(unknown))}")

    result = {}
    config_dir = os.path.dirname(path)
    for k, v in raw.items():
        if k in _STRING_KEYS:
            value = str(v)
            if k in _PATH_KEYS:
                value = _normalize_local_path(value, config_dir)
            result[k] = value
        elif k in _BOOL_KEYS:
            result[k] = bool(v)
    return result


def _normalize_local_path(value: str, config_dir: str) -> str:
    """Return *value* resolved relative to *config_dir* when it is local and relative."""
    lower = value.lower()
    if lower.startswith(_URI_PREFIXES):
        return value

    expanded = os.path.expanduser(value)
    if os.path.isabs(expanded):
        return expanded

    return os.path.abspath(os.path.join(config_dir, expanded))


def _find_config_path() -> str | None:
    """Return the config file path to use, or None if none found."""
    # Check for --config flag in sys.argv (passed after `--` to streamlit)
    args = sys.argv[1:]
    for i, arg in enumerate(args):
        if arg == "--config" and i + 1 < len(args):
            p = args[i + 1]
            if os.path.isfile(p):
                return os.path.abspath(p)
            # File specified but not found — warn but don't crash
            import streamlit as st
            st.warning(f"--config path not found: {p}")
            return None

    # Fall back to geac.toml in cwd
    cwd_toml = os.path.join(os.getcwd(), "geac.toml")
    if os.path.isfile(cwd_toml):
        return cwd_toml

    return None
