"""AI Plot Builder tab — natural-language → executable visualization code.

Sends a description of the schema and the user's request to Claude, streams
back Python code, and `exec`s it in a namespace pre-populated with the
filtered DuckDB connection, altair, pandas, numpy, and Streamlit. Gated by:
(1) `anthropic` package installed; (2) API key from config or env or input;
(3) explicit acknowledgement of data-privacy risk.
"""
from __future__ import annotations

import os
import textwrap
import traceback

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from explorer.tab_context import TabContext

try:
    import anthropic as _anthropic
    _HAS_ANTHROPIC = True
except ImportError:
    _anthropic = None
    _HAS_ANTHROPIC = False

LABEL = "🤖 AI Plot Builder"

_AI_ACK_PHRASE = "I understand the data privacy risks with using the AI Plot Builder."


def render(ctx: TabContext) -> None:
    st.header("AI Plot Builder")
    st.caption(
        "Describe a plot in plain English and Claude will write the code to render it. "
        "The generated code has access to the current filtered dataset, so all sidebar "
        "filters apply automatically."
    )
    st.error(
        "⚠️ **Data privacy warning:** Generating a plot using the AI Plot Builder sends a description of your data "
        "schema and query results to the Anthropic API. "
        "**Do not use this feature with clinical, patient-identifiable, or proprietary data.**",
        icon="🔴",
    )

    if not _HAS_ANTHROPIC:
        st.error(
            "`anthropic` package is not installed. "
            "Run `pip install anthropic` in your environment then restart the explorer."
        )
        return

    _ai_api_key = (
        ctx.cfg.get("anthropic_api_key")
        or os.environ.get("ANTHROPIC_API_KEY", "")
    )
    if not _ai_api_key:
        _ai_api_key = st.text_input(
            "Anthropic API key",
            type="password",
            key="ai_api_key_input",
            help=(
                "Set ANTHROPIC_API_KEY in your environment or add anthropic_api_key "
                "to config.toml to avoid entering it here."
            ),
        )

    if not _ai_api_key:
        st.info("Enter your Anthropic API key above to enable the AI Plot Builder.")
        return

    _ai_ack = st.text_input(
        f'Type exactly: "{_AI_ACK_PHRASE}"',
        key="ai_ack_input",
        help="Required acknowledgement before the AI Plot Builder can be used.",
    )
    if _ai_ack.strip() != _AI_ACK_PHRASE:
        st.info(
            "Type the acknowledgement phrase above (exactly as shown) to enable "
            "the AI Plot Builder."
        )
        return

    _system = _build_system_prompt(ctx)
    _request = st.text_area(
        "Describe the plot or analysis you want",
        placeholder=(
            "e.g. 'Show a scatter plot of VAF vs. sequencing cycle for alt reads, "
            "colored by variant type' or 'Plot the distribution of base quality scores "
            "split by R1 and R2'"
        ),
        height=100,
        key="ai_plot_request",
    )

    _col1, _col2 = st.columns([1, 5])
    _generate = _col1.button("Generate", key="ai_generate_btn", type="primary")
    if _col2.button("Clear", key="ai_clear_btn"):
        st.session_state.pop("ai_last_code", None)
        st.session_state.pop("ai_last_request", None)
        st.rerun()

    if _generate and _request.strip():
        st.session_state["ai_last_request"] = _request.strip()
        with st.spinner("Generating…"):
            try:
                _client = _anthropic.Anthropic(api_key=_ai_api_key)
                _parts: list[str] = []
                _placeholder = st.empty()
                with _client.messages.stream(
                    model="claude-opus-4-6",
                    max_tokens=4096,
                    system=_system,
                    messages=[{"role": "user", "content": _request.strip()}],
                ) as _stream:
                    for _chunk in _stream.text_stream:
                        _parts.append(_chunk)
                        _placeholder.code("".join(_parts), language="python")
                st.session_state["ai_last_code"] = "".join(_parts).strip()
            except Exception as _err:
                st.error(f"Claude API error: {_err}")

    if "ai_last_code" in st.session_state:
        _code = st.session_state["ai_last_code"]
        with st.expander("Generated code", expanded=False):
            st.code(_code, language="python")

        st.divider()
        _namespace = {
            "con": ctx.con,
            "table_expr": ctx.table_expr,
            "where": ctx.where,
            "alt": alt,
            "pd": pd,
            "np": np,
            "st": st,
            "_r_join": ctx.r_join,
        }
        try:
            exec(_code, _namespace)  # noqa: S102
        except Exception:
            st.error("The generated code raised an error:")
            st.code(traceback.format_exc(), language="text")
            st.caption("Try rephrasing your request or click Generate again.")


def _build_system_prompt(ctx: TabContext) -> str:
    _locus_cols = ", ".join(ctx.table_cols) if ctx.table_cols else (
        "sample_id, chrom, pos, alt_allele, variant_type, total_depth, alt_count, "
        "ref_allele, fwd_alt_count, rev_alt_count, pipeline, batch, label1, "
        "trinuc_context, gnomad_af, on_target, gene, homopolymer_len, str_len, "
        "variant_called, variant_filter"
    )
    _reads_cols = ", ".join(sorted(ctx.alt_reads_cols)) if ctx.alt_reads_cols else (
        "sample_id, chrom, pos, alt_allele, cycle, read_length, is_read1, base_qual, "
        "map_qual, family_size, insert_size, n_before_alt, n_after_alt, "
        "n_n_before_alt, n_n_after_alt, leading_n_run_len, trailing_n_run_len"
    )

    return textwrap.dedent(f"""
        You are an expert data visualization assistant for the GEAC genomic cohort explorer.
        The user will describe a plot or analysis. You must return ONLY executable Python code — no
        explanation, no markdown code fences, no comments unless they clarify non-obvious logic.

        ## Available variables (already in scope when the code runs)

        | Variable | Description |
        |----------|-------------|
        | `con` | DuckDB connection to the cohort database |
        | `table_expr` | SQL table/view expression for the current locus selection (respects all sidebar filters) |
        | `where` | SQL WHERE clause string encoding the current sidebar filters |
        | `alt` | altair module |
        | `pd` | pandas module |
        | `np` | numpy module |
        | `st` | streamlit module — use `st.altair_chart`, `st.dataframe`, `st.metric`, etc. to display output |
        | `_r_join` | SQL FROM expression that joins alt_reads to the current filtered loci (use for read-level queries) |

        ## Schema

        ### Locus table — `{{table_expr}}`
        Columns: {_locus_cols}
        - `vaf` is not stored; compute as `ROUND(alt_count * 1.0 / total_depth, 4)` in SQL or `alt_count / total_depth` in pandas.
        - Query pattern: `con.execute(f"SELECT ... FROM {{table_expr}} WHERE {{where}}").df()`

        ### Alt-reads table — `alt_reads`
        Columns: {_reads_cols}
        - Query pattern: `con.execute(f"SELECT ... FROM {{_r_join}}").df()`
        - `_r_join` already filters alt_reads to only the loci matching the current sidebar filters.

        ## Display conventions
        - Use `st.altair_chart(chart, use_container_width=True)` to render Altair charts.
        - Use `st.dataframe(df, use_container_width=True, hide_index=True)` for tables.
        - If the query returns no rows, show `st.info("No data available under current filters.")`.
        - Prefer Altair for charts. Only use matplotlib if the user explicitly asks for it.
        - Do not call `st.set_page_config` or `st.sidebar` anything.
    """).strip()
