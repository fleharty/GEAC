"""Pytest configuration for the Explorer Python helper tests.

These tests cover pure-Python logic (schema parsing, region parsing, chrom
sorting, DataSource DB queries, tab metadata) and are documented to run without
a Streamlit runtime (see CLAUDE.md). The UI tab modules under ``explorer.tabs``
import ``streamlit`` at module top but only call into it from within functions,
so when Streamlit isn't installed we register a stub module that satisfies the
import without pulling in the full UI runtime. When Streamlit *is* installed,
the real package is used unchanged.
"""

import importlib.util
import sys
from unittest.mock import MagicMock

# UI-presentation packages that the tab modules import at module top but only
# call into from within render functions. They are not needed to exercise the
# pure-Python logic under test, so stub any that aren't installed. Core data
# packages (pandas/numpy/scipy/duckdb) are intentionally *not* stubbed — the
# tests genuinely depend on them and a stub would mask real failures.
_OPTIONAL_UI_MODULES = ("streamlit", "altair")

for _mod in _OPTIONAL_UI_MODULES:
    if importlib.util.find_spec(_mod) is None:
        # MagicMock transparently supports attribute access, calls, and
        # decorator usage (e.g. ``@st.cache_data``), which is all the tab
        # modules need at import time.
        sys.modules[_mod] = MagicMock()
