"""Load GEAC project configuration from a TOML file.

Discovery order:
1. Path passed via `streamlit run geac_explorer.py -- --config /path/to/geac.toml`
2. `geac.toml` in the current working directory
3. Empty config (all fields blank — user fills in manually)

Supported keys (all optional):

    data     = "/path/to/cohort.duckdb"   # or a .parquet file
    manifest = "/path/to/manifest.tsv"
    cosmic   = "/path/to/COSMIC_v3.4_SBS_GRCh37.txt"
    genome   = "hg38"                     # hg19 | hg38 | mm10 | mm39 | other
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


_KNOWN_KEYS = {"data", "manifest", "cosmic", "genome"}


def load() -> dict[str, str]:
    """Return config dict with string values for known keys.

    Missing keys are absent from the returned dict (not None), so callers
    can use .get("data", "") to obtain a safe default.
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

    return {k: str(v) for k, v in raw.items() if k in _KNOWN_KEYS}


def _find_config_path() -> str | None:
    """Return the config file path to use, or None if none found."""
    # Check for --config flag in sys.argv (passed after `--` to streamlit)
    args = sys.argv[1:]
    for i, arg in enumerate(args):
        if arg == "--config" and i + 1 < len(args):
            p = args[i + 1]
            if os.path.isfile(p):
                return p
            # File specified but not found — warn but don't crash
            import streamlit as st
            st.warning(f"--config path not found: {p}")
            return None

    # Fall back to geac.toml in cwd
    cwd_toml = os.path.join(os.getcwd(), "geac.toml")
    if os.path.isfile(cwd_toml):
        return cwd_toml

    return None
