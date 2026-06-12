"""Sample Identity / Duplicates tab.

Detects samples that are genetically the same individual (sample swaps, unknown
duplicates, technical replicates) and flags possible cross-contamination, using
germline-SNP fingerprinting over the merged cohort's ``alt_bases`` table.

DuckDB-mode only (needs the full cohort). The heavy pairwise self-join is cached
in ``st.session_state`` keyed on the control values — see CHALLENGES.md on why
uncached tab-body queries make the Explorer slow.
"""
from __future__ import annotations

import altair as alt
import pandas as pd
import streamlit as st

from explorer.sample_identity_helpers import (
    IdentityParams,
    classify_flags,
    compute_all_loci_jaccard,
    compute_contamination,
    compute_pairwise_identity,
    subject_map,
)
from explorer.tab_context import TabContext

LABEL = "🧬 Sample Identity"

_FLAG_HELP = {
    "UNKNOWN_DUPLICATE": "Genetically identical but different/absent subject_id — a likely unrecognized duplicate.",
    "POSSIBLE_SWAP": "Same subject_id but genetically divergent — a likely sample swap or mislabel.",
    "EXPECTED_MATCH": "Same subject_id and high similarity — expected; shown as a sanity check.",
}


def render(ctx: TabContext) -> None:
    if not ctx.path.endswith(".duckdb"):
        st.info("Sample identity is available when loading a merged DuckDB file (`geac merge` output).")
        return

    st.subheader("Sample identity & duplicate detection")
    st.caption(
        "Germline-SNP fingerprinting over `alt_bases`. The primary signal is the "
        "**Jaccard overlap** of each sample's germline-variant set; **concordance** "
        "(matching het/hom genotypes at shared sites) confirms identity. Same-individual "
        "pairs score near 1.0 on both; unrelated pairs share only common SNPs."
    )

    has_gnomad = ctx.has_data("gnomad_af")
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        min_depth = st.number_input(
            "Min depth", min_value=1, max_value=10_000, step=5,
            key="si_min_depth",
            help="Minimum total_depth for an SNV call to be usable as a marker.",
        )
    with c2:
        het_lo = st.slider(
            "Germline VAF floor", min_value=0.05, max_value=0.45, step=0.01,
            key="si_het_lo",
            help="Minimum VAF to treat an alt record as a germline call.",
        )
    with c3:
        min_recurrence = st.number_input(
            "Min marker recurrence", min_value=2, max_value=1000, step=1,
            key="si_min_recurrence",
            help="A marker must appear in at least this many samples (selects common population SNPs).",
        )
    with c4:
        use_gnomad = st.checkbox(
            "Common gnomAD SNPs only", key="si_use_gnomad", disabled=not has_gnomad,
            help="Restrict markers to gnomAD AF 0.05–0.95." + ("" if has_gnomad else " (gnomad_af not in this cohort.)"),
        )

    t1, t2 = st.columns(2)
    with t1:
        t_high = st.slider(
            "Duplicate Jaccard threshold", min_value=0.3, max_value=0.99, step=0.01,
            key="si_t_high",
            help="Pairs at/above this Jaccard (with high concordance) and a different subject_id are flagged as unknown duplicates.",
        )
    with t2:
        t_low = st.slider(
            "Swap Jaccard threshold", min_value=0.0, max_value=0.7, step=0.01,
            key="si_t_low",
            help="Pairs sharing a subject_id but below this Jaccard (or low concordance) are flagged as possible swaps.",
        )

    p = IdentityParams(
        min_depth=int(min_depth),
        het_lo=float(het_lo),
        min_recurrence=int(min_recurrence),
        use_gnomad=bool(use_gnomad and has_gnomad),
        t_high=float(t_high),
        t_low=float(t_low),
    )

    cache_key = (
        ctx.path, p.min_depth, p.het_lo, p.min_recurrence, p.use_gnomad,
        p.max_markers, p.t_low,
    )
    if st.session_state.get("_si_cache_key") != cache_key:
        with ctx.timed("sample identity pairwise [cache miss]"):
            pairs = compute_pairwise_identity(ctx.con, p)
            subjects = subject_map(ctx.con)
            # All-loci Jaccard only for candidate pairs (cheap second pass).
            candidates = pairs[pairs["jaccard"] >= max(p.t_low, 0.3)]
            cand_samples = sorted(
                set(candidates["sample_a"]) | set(candidates["sample_b"])
            )
            all_loci = compute_all_loci_jaccard(ctx.con, cand_samples, p)
            contamination = compute_contamination(ctx.con, p)
        st.session_state["_si_cache_key"] = cache_key
        st.session_state["_si_pairs"] = pairs
        st.session_state["_si_subjects"] = subjects
        st.session_state["_si_all_loci"] = all_loci
        st.session_state["_si_contamination"] = contamination

    pairs = st.session_state["_si_pairs"]
    subjects = st.session_state["_si_subjects"]
    all_loci = st.session_state["_si_all_loci"]
    contamination = st.session_state["_si_contamination"]

    if pairs.empty:
        st.warning(
            "No sample pairs share any germline marker under these settings. "
            "Try lowering min depth, the VAF floor, or marker recurrence.",
            icon="🔎",
        )
        return

    flagged = classify_flags(pairs, subjects, p)
    if not all_loci.empty:
        flagged = flagged.merge(all_loci, on=["sample_a", "sample_b"], how="left")
    else:
        flagged["all_loci_jaccard"] = pd.NA

    _render_flagged_pairs(flagged, p)
    st.divider()
    _render_heatmap(flagged)
    st.divider()
    _render_contamination(contamination)


def _render_flagged_pairs(flagged: pd.DataFrame, p: IdentityParams) -> None:
    flagged_only = flagged[flagged["flag"] != ""].copy()
    st.markdown("#### Flagged pairs")
    if flagged_only.empty:
        st.success(
            "No suspicious pairs under the current thresholds — no unknown duplicates "
            "or swaps detected.",
            icon="✅",
        )
    else:
        for flag, desc in _FLAG_HELP.items():
            n = int((flagged_only["flag"] == flag).sum())
            if n:
                st.caption(f"**{flag}** ({n}): {desc}")
        _cols = [
            "flag", "sample_a", "subject_a", "sample_b", "subject_b",
            "jaccard", "concordance", "all_loci_jaccard",
            "n_shared", "n_nonref_a", "n_nonref_b",
        ]
        _disp = flagged_only.sort_values(
            ["flag", "jaccard"], ascending=[True, False]
        )[_cols]
        st.dataframe(
            _disp.style.format(
                {"jaccard": "{:.3f}", "concordance": "{:.3f}", "all_loci_jaccard": "{:.3f}"}
            ),
            width="stretch",
            hide_index=True,
        )
        st.caption(
            "A flagged same-individual pair whose **all_loci_jaccard** is also very "
            "high (matching even low-VAF/somatic loci) is likely a *technical "
            "replicate* of one sample rather than two distinct biological samples."
        )

    csv = flagged.sort_values("jaccard", ascending=False).to_csv(index=False)
    st.download_button(
        "Download all pairwise scores (CSV)",
        data=csv,
        file_name="sample_identity_pairs.csv",
        mime="text/csv",
    )


def _render_heatmap(flagged: pd.DataFrame) -> None:
    st.markdown("#### Pairwise similarity heatmap")
    st.caption("Color = germline Jaccard. Off-diagonal hot blocks indicate duplicate clusters.")

    # Symmetrize for a full matrix; diagonal = 1.0 (self-identity).
    a = flagged[["sample_a", "sample_b", "jaccard"]].rename(
        columns={"sample_a": "x", "sample_b": "y"}
    )
    b = flagged[["sample_a", "sample_b", "jaccard"]].rename(
        columns={"sample_a": "y", "sample_b": "x"}
    )
    samples = sorted(set(flagged["sample_a"]) | set(flagged["sample_b"]))
    diag = pd.DataFrame({"x": samples, "y": samples, "jaccard": 1.0})
    mat = pd.concat([a, b, diag], ignore_index=True)

    if len(samples) > 60:
        st.info(
            f"{len(samples)} samples share markers — heatmap omitted for legibility. "
            "Use the flagged-pairs table and CSV download instead."
        )
        return

    chart = (
        alt.Chart(mat)
        .mark_rect()
        .encode(
            alt.X("x:N", title=None, axis=alt.Axis(labelAngle=-45, labelLimit=120)),
            alt.Y("y:N", title=None, axis=alt.Axis(labelLimit=120)),
            alt.Color("jaccard:Q", title="Jaccard", scale=alt.Scale(scheme="magma", domain=[0, 1])),
            tooltip=[
                "x:N", "y:N",
                alt.Tooltip("jaccard:Q", format=".3f", title="Jaccard"),
            ],
        )
        .properties(height=min(600, 24 * len(samples) + 80))
    )
    st.altair_chart(chart, width="stretch")


def _render_contamination(contamination: pd.DataFrame) -> None:
    st.markdown("#### Per-sample contamination indicators")
    st.caption(
        "Heuristic only (not a calibrated estimate like VerifyBamID). "
        "**het VAF dispersion** = mean |VAF − 0.5| over germline het markers; it "
        "drifts up when foreign reads skew heterozygous balance. **low-VAF SNV "
        "count** is the burden of sub-germline alt SNVs. High values on both axes "
        "warrant a closer look."
    )
    if contamination.empty:
        st.info("No germline het markers available to assess contamination.")
        return

    chart = (
        alt.Chart(contamination.dropna(subset=["het_vaf_dispersion"]))
        .mark_circle(size=90, opacity=0.85)
        .encode(
            alt.X("het_vaf_dispersion:Q", title="Het VAF dispersion (mean |VAF − 0.5|)"),
            alt.Y("low_vaf_snv_count:Q", title="Low-VAF SNV count"),
            alt.Color("sample_id:N", title="Sample", legend=None),
            tooltip=[
                "sample_id:N",
                alt.Tooltip("het_vaf_dispersion:Q", format=".4f", title="Het VAF dispersion"),
                alt.Tooltip("n_het_markers:Q", format=",", title="Het markers"),
                alt.Tooltip("low_vaf_snv_count:Q", format=",", title="Low-VAF SNVs"),
            ],
        )
        .properties(height=320)
    )
    st.altair_chart(chart, width="stretch")
    st.dataframe(
        contamination.style.format({"het_vaf_dispersion": "{:.4f}"}),
        width="stretch",
        hide_index=True,
    )
