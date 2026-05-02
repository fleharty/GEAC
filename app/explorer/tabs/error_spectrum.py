"""Error Spectrum tab — SBS96 trinucleotide spectrum and COSMIC/NMF signature analysis.

Sections (rendered top to bottom):
  1. SBS96 spectrum (96 trinucleotide contexts, optional COSMIC reconstruction overlay)
  2. Per-sample COSMIC signature exposures (DuckDB cohort only)
  3. NMF signature discovery (de novo, or COSMIC-guided + one learned residual)
  4. Called vs Uncalled mirrored spectrum + COSMIC compare (when variant_called present)
  5. VAF-stratified spectrum (germline VAF > 30% vs somatic VAF ≤ 30%)
  6. R1 / R2 stratified spectrum (DuckDB / alt_reads only)
  7. Cohort SBS96 heatmap (samples × contexts; DuckDB only)

Falls back to a simple ref>alt bar chart when trinuc_context is unavailable.
"""
from __future__ import annotations

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
from scipy.optimize import nnls

from explorer.sbs96 import (
    SBS_COLORS,
    SBS_ETIOLOGY,
    SBS_MUT_TYPES,
    SBS_ORDER,
    load_cosmic,
    sbs_label,
    strat_sbs96_chart,
    to_spec96_strat,
)
from explorer.tab_context import TabContext
from signature_nmf import (
    build_signature_download_table,
    build_signature_download_zip,
    build_signature_exposure_download_table,
    compare_signatures_to_cosmic,
    fit_cosmic_augmented_nmf,
    fit_cosmic_cohort_greedy,
    fit_cosmic_single_sample_greedy,
    fit_sbs_nmf,
)

LABEL = "🧬 Error spectrum"


def render(ctx: TabContext) -> None:
    if ctx.has_data("trinuc_context"):
        _render_full_sbs96(ctx)
    else:
        _render_simple_fallback(ctx)

    if ctx.path.endswith(".duckdb"):
        _render_cohort_heatmap(ctx)


# ─────────────────────────────────────────────────────────────────────────────
# SBS96 helpers
# ─────────────────────────────────────────────────────────────────────────────


def _load_sample_sbs96_matrix(ctx: TabContext) -> pd.DataFrame:
    has_batch = ctx.has_data("batch")
    id_sql = "sample_id || ' / ' || batch" if has_batch else "sample_id"
    group_by = "sample_id, batch" if has_batch else "sample_id"
    raw = ctx.con.execute(f"""
        SELECT {id_sql} AS sample_label,
               trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
        FROM (SELECT * FROM {ctx.table_expr}) _t
        WHERE {ctx.where} AND variant_type = 'SNV'
          AND trinuc_context IS NOT NULL AND length(trinuc_context) = 3
        GROUP BY {group_by}, trinuc_context, ref_allele, alt_allele
    """).df()

    if raw.empty:
        return pd.DataFrame(columns=SBS_ORDER, dtype=float)

    rows = []
    for sid, grp in raw.groupby("sample_label"):
        grp = grp.copy()
        grp["sbs_label"] = grp.apply(
            lambda r: sbs_label(r["trinuc_context"], r["ref_allele"], r["alt_allele"]),
            axis=1,
        )
        grp = grp.dropna(subset=["sbs_label"])
        agg = grp.groupby("sbs_label")["count"].sum()
        counts = (
            pd.Series(0.0, index=SBS_ORDER)
            .add(agg, fill_value=0)
            .reindex(SBS_ORDER)
        )
        if float(counts.sum()) > 0:
            rows.append(pd.Series(counts.values, index=SBS_ORDER, name=sid))

    if not rows:
        return pd.DataFrame(columns=SBS_ORDER, dtype=float)

    return pd.DataFrame(rows, dtype=float)


def _profile_to_spec96(profile, value_name: str = "fraction") -> pd.DataFrame:
    return pd.DataFrame({
        "sbs_label": SBS_ORDER,
        "mut_type": [lbl[2:5] for lbl in SBS_ORDER],
        value_name: pd.Series(profile, index=SBS_ORDER).reindex(SBS_ORDER).values,
    })


def _signature_profile_chart(spec_df, title, overlay_df=None, overlay_label=None):
    sub_charts = []
    for mt in SBS_MUT_TYPES:
        sub = spec_df[spec_df["mut_type"] == mt]
        order = [lbl for lbl in SBS_ORDER if f"[{mt}]" in lbl]
        bars = (
            alt.Chart(sub)
            .mark_bar(color=SBS_COLORS[mt])
            .encode(
                alt.X("sbs_label:N", sort=order, title=None,
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
            overlay_sub = overlay_df[overlay_df["mut_type"] == mt]
            overlay = (
                alt.Chart(overlay_sub)
                .mark_point(color="black", size=14, filled=True, opacity=0.85)
                .encode(
                    alt.X("sbs_label:N", sort=order),
                    alt.Y("fraction:Q"),
                    tooltip=[
                        "sbs_label:N",
                        alt.Tooltip("fraction:Q", format=".3f", title=overlay_label or "Overlay"),
                    ],
                )
            )
            panel = alt.layer(bars, overlay)
        else:
            panel = bars

        sub_charts.append(
            panel.properties(
                title=alt.TitleParams(mt, color=SBS_COLORS[mt], fontSize=11, fontWeight="bold"),
                width=120,
                height=110,
            )
        )

    return (
        alt.concat(*sub_charts, columns=3)
        .resolve_scale(y="shared")
        .properties(title=alt.TitleParams(title, fontSize=13))
    )


# ─────────────────────────────────────────────────────────────────────────────
# Main SBS96 path
# ─────────────────────────────────────────────────────────────────────────────


def _render_full_sbs96(ctx: TabContext) -> None:
    raw = ctx.con.execute(f"""
        SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
        FROM {ctx.table_expr}
        WHERE {ctx.where} AND variant_type = 'SNV' AND trinuc_context IS NOT NULL
          AND length(trinuc_context) = 3
        GROUP BY trinuc_context, ref_allele, alt_allele
    """).df()

    if raw.empty:
        st.info("No SNVs with trinucleotide context in current selection.")
        return

    raw["sbs_label"] = raw.apply(
        lambda row: sbs_label(row["trinuc_context"], row["ref_allele"], row["alt_allele"]),
        axis=1,
    )
    raw = raw.dropna(subset=["sbs_label"])
    raw["mut_type"] = raw["sbs_label"].str.extract(r'\[([A-Z]>[A-Z])\]')[0]

    agg = raw.groupby(["sbs_label", "mut_type"], as_index=False)["count"].sum()
    full = pd.DataFrame({
        "sbs_label": SBS_ORDER,
        "mut_type":  [lbl[2:5] for lbl in SBS_ORDER],
    })
    spec96 = full.merge(agg, on=["sbs_label", "mut_type"], how="left")
    spec96["count"] = spec96["count"].fillna(0).astype(int)
    total_snvs = spec96["count"].sum()
    spec96["fraction"] = spec96["count"] / total_snvs if total_snvs > 0 else 0.0

    sbs_y_mode = st.radio("Y axis", ["Fraction", "Count"], horizontal=True, key="sbs_y_mode")
    sbs_use_fraction = sbs_y_mode == "Fraction"
    sbs_y_field = "fraction" if sbs_use_fraction else "count"
    sbs_y_title = "Fraction of SNVs" if sbs_use_fraction else "Count"
    sbs_y_fmt   = ".3f" if sbs_use_fraction else "d"

    cos_col, path_col = st.columns([1, 3])
    cosmic_path = path_col.text_input(
        "COSMIC SBS matrix path (optional — enables reconstruction overlay)",
        value=ctx.cfg.get("cosmic", ""),
        placeholder="/path/to/COSMIC_v3.4_SBS_GRCh37.txt",
        key="cosmic_path",
        label_visibility="visible",
    )

    cosmic_state = _maybe_fit_cosmic(
        ctx, cos_col, cosmic_path, spec96, total_snvs,
    )

    _render_spectrum_chart(
        ctx, spec96, raw, total_snvs,
        sbs_use_fraction, sbs_y_field, sbs_y_title, sbs_y_fmt,
        cosmic_state,
    )

    if cosmic_state is not None and cosmic_state["cos_sim"] is not None:
        _render_cosmic_results(ctx, cosmic_path, cosmic_state)

    if ctx.data_source.is_duckdb:
        _render_nmf_discovery(ctx, cosmic_path, cosmic_state)

    if ctx.has_data("variant_called"):
        _render_called_vs_uncalled(ctx, cosmic_state)

    _render_vaf_stratified(ctx)

    _render_r1r2_stratified(ctx)


def _maybe_fit_cosmic(ctx, cos_col, cosmic_path, spec96, total_snvs):
    """Load COSMIC matrix and fit signatures. Returns dict of intermediate state, or None."""
    if not (cosmic_path and cosmic_path.strip()):
        return None

    state = {
        "cosmic_W": None,
        "cosmic_aligned": None,
        "reconstructed": None,
        "recon_df": None,
        "sig_df": None,
        "top_df": None,
        "top_n_sig": None,
        "cos_sim": None,
        "residual_pct": None,
        "fit_method": None,
        "add_penalty": None,
        "remove_penalty": None,
        "use_connected": None,
    }

    try:
        cosmic_df      = load_cosmic(cosmic_path.strip())
        cosmic_aligned = cosmic_df.reindex(SBS_ORDER)
        missing        = cosmic_aligned.isna().any(axis=1).sum()
        if missing > 0:
            st.warning(
                f"{missing} context(s) not found in COSMIC matrix — "
                "check that the file uses the standard A[C>A]A format."
            )
            return state

        W = cosmic_aligned.values.astype(float)
        state["cosmic_W"] = W
        state["cosmic_aligned"] = cosmic_aligned
        obs = (
            spec96.set_index("sbs_label")["count"]
            .reindex(SBS_ORDER).fillna(0).values.astype(float)
        )

        fit_method = cos_col.radio(
            "Fitting method",
            ["NNLS (top-N)", "Greedy add/remove (SPA-style)"],
            horizontal=False,
            key="cosmic_fit_method",
            help=(
                "**NNLS**: fit all signatures at once, then restrict to the top-N "
                "by exposure. Simple and fast.\n\n"
                "**Greedy add/remove**: iteratively add signatures that reduce "
                "normalized L2 by more than *add penalty*, then remove signatures "
                "whose removal degrades L2 by less than *remove penalty*. "
                "Mirrors SigProfilerAssignment and produces sparser, more "
                "interpretable results."
            ),
        )
        state["fit_method"] = fit_method

        if fit_method == "NNLS (top-N)":
            h, _ = nnls(W, obs)
            total = h.sum()
            h_norm = h / total if total > 0 else h

            sig_df = pd.DataFrame({
                "signature": cosmic_aligned.columns.tolist(),
                "exposure":  h_norm,
            })
            sig_df["etiology"] = sig_df["signature"].map(lambda s: SBS_ETIOLOGY.get(s, ""))
            sig_df = sig_df[sig_df["exposure"] > 0].sort_values(
                "exposure", ascending=False
            ).reset_index(drop=True)

            top_n_sig = cos_col.slider(
                "Top signatures", 3, min(20, len(sig_df)), 4,
                key="top_n_sig",
                help="Number of top signatures used for the reconstruction overlay and exposure chart.",
            )
            top_df = sig_df.head(top_n_sig)

            top_sig_names = top_df["signature"].tolist()
            W_refit    = cosmic_aligned[top_sig_names].values.astype(float)
            h_refit, _ = nnls(W_refit, obs)
            reconstructed = W_refit @ h_refit
        else:
            add_penalty = cos_col.slider(
                "Add penalty",
                min_value=0.001, max_value=0.20, value=0.05, step=0.005,
                key="cosmic_add_penalty", format="%.3f",
                help="Minimum fractional L2 improvement required to add a signature. "
                     "Higher = sparser. SPA default: 0.05.",
            )
            remove_penalty = cos_col.slider(
                "Remove penalty",
                min_value=0.001, max_value=0.10, value=0.01, step=0.002,
                key="cosmic_remove_penalty", format="%.3f",
                help="Maximum fractional L2 degradation allowed when removing a signature. "
                     "Higher = more aggressive pruning. SPA default: 0.01.",
            )
            use_connected = cos_col.checkbox(
                "Expand connected groups (SBS2/13, SBS7a–d, SBS10a/b, SBS17a/b)",
                value=True, key="cosmic_connected_sigs",
                help="If any member of a biologically linked group is selected, "
                     "all group members are fitted together.",
            )

            greedy_result = fit_cosmic_single_sample_greedy(
                obs, W, cosmic_aligned.columns.tolist(),
                add_penalty=add_penalty,
                remove_penalty=remove_penalty,
                connected_sigs=use_connected,
            )
            h_norm = greedy_result["exposure_fractions"]
            sig_df = pd.DataFrame({
                "signature": cosmic_aligned.columns.tolist(),
                "exposure":  h_norm,
            })
            sig_df["etiology"] = sig_df["signature"].map(lambda s: SBS_ETIOLOGY.get(s, ""))
            sig_df = sig_df[sig_df["exposure"] > 0].sort_values(
                "exposure", ascending=False
            ).reset_index(drop=True)
            top_df = sig_df
            top_n_sig = len(top_df)

            top_sig_names = top_df["signature"].tolist()
            if top_sig_names:
                W_refit    = cosmic_aligned[top_sig_names].values.astype(float)
                h_refit    = greedy_result["exposures"][
                    [cosmic_aligned.columns.tolist().index(n) for n in top_sig_names]
                ]
                reconstructed = W_refit @ h_refit
            else:
                reconstructed = np.zeros_like(obs)

            state["add_penalty"] = add_penalty
            state["remove_penalty"] = remove_penalty
            state["use_connected"] = use_connected

        state["reconstructed"] = reconstructed
        state["sig_df"] = sig_df
        state["top_df"] = top_df
        state["top_n_sig"] = top_n_sig
        state["cos_sim"] = (
            float(np.dot(obs, reconstructed))
            / (np.linalg.norm(obs) * np.linalg.norm(reconstructed) + 1e-12)
        )
        state["residual_pct"] = (
            float(np.linalg.norm(obs - reconstructed))
            / (float(obs.sum()) + 1e-12) * 100
        )

        recon_df = spec96[["sbs_label", "mut_type"]].copy()
        recon_df["recon_count"] = reconstructed
        recon_df["recon_frac"]  = reconstructed / total_snvs if total_snvs > 0 else 0.0
        state["recon_df"] = recon_df

    except Exception as exc:
        st.error(f"Failed to load COSMIC matrix: {exc}")

    return state


def _render_spectrum_chart(
    ctx, spec96, raw, total_snvs,
    sbs_use_fraction, sbs_y_field, sbs_y_title, sbs_y_fmt,
    cosmic_state,
):
    recon_df  = cosmic_state["recon_df"]   if cosmic_state else None
    top_n_sig = cosmic_state["top_n_sig"]  if cosmic_state else None
    recon_y   = "recon_frac"  if sbs_use_fraction else "recon_count"

    chart_title = (
        f"SNV Trinucleotide Spectrum — bars = observed, dots = reconstruction (top {top_n_sig} sigs)"
        if recon_df is not None else
        "SNV Trinucleotide Spectrum (SBS96)"
    )
    sel_param = alt.selection_point(name="bar_click", fields=["sbs_label"], on="click")

    sub_charts = []
    for mt in SBS_MUT_TYPES:
        obs_sub = spec96[spec96["mut_type"] == mt].copy()
        order   = [lbl for lbl in SBS_ORDER if f"[{mt}]" in lbl]

        bars = (
            alt.Chart(obs_sub)
            .mark_bar(color=SBS_COLORS[mt])
            .encode(
                alt.X("sbs_label:N", sort=order, title=None,
                      axis=alt.Axis(labelAngle=-90, labelFontSize=8)),
                alt.Y(f"{sbs_y_field}:Q", title=sbs_y_title,
                      **({"axis": alt.Axis(format=".3f")} if sbs_use_fraction else {})),
                opacity=alt.condition({"param": "bar_click"}, alt.value(1.0), alt.value(0.4)),
                tooltip=[
                    "sbs_label:N",
                    alt.Tooltip(f"{sbs_y_field}:Q", title="Observed", format=sbs_y_fmt),
                ],
            )
        )

        if recon_df is not None:
            recon_sub = recon_df[recon_df["mut_type"] == mt].copy()
            dots = (
                alt.Chart(recon_sub)
                .mark_point(color="black", size=15, filled=True, opacity=0.85)
                .encode(
                    alt.X("sbs_label:N", sort=order),
                    alt.Y(f"{recon_y}:Q"),
                    tooltip=alt.value(None),
                )
            )
            panel = alt.layer(bars, dots)
        else:
            panel = bars

        sub_charts.append(
            panel.properties(
                title=alt.TitleParams(mt, color=SBS_COLORS[mt], fontSize=13, fontWeight="bold"),
                width=150, height=140,
            )
        )

    chart = (
        alt.concat(*sub_charts, columns=3)
        .resolve_scale(y="shared")
        .properties(title=alt.TitleParams(chart_title, fontSize=14))
        .add_params(sel_param)
    )
    event = st.altair_chart(chart, width="stretch", on_select="rerun", key="sbs96_spectrum")

    if recon_df is not None:
        st.caption(
            f"{total_snvs:,} SNV alt-allele loci. "
            f"Black dots = COSMIC reconstruction refit to top {top_n_sig} signatures. "
            "Click bars to drill down. "
            "Contexts where dots deviate from bars are poorly explained by the selected signatures."
        )
    else:
        st.caption(
            f"{total_snvs:,} SNV alt-allele loci. "
            "Click one or more bars to drill down and open in IGV. "
            "Enter a COSMIC matrix path above to overlay the reconstruction."
        )

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
                sel = ctx.query_records([extra_cond])
                label_str = ", ".join(clicked_labels)
                st.caption(f"{len(sel):,} records matching {len(clicked_labels)} selected context(s): {label_str}")
                st.dataframe(sel[ctx.table_cols], width="stretch")
                ctx.igv_buttons([extra_cond], sel, key=f"sbs_{'_'.join(clicked_labels)}")


def _render_cosmic_results(ctx, cosmic_path, cosmic_state):
    sig_df       = cosmic_state["sig_df"]
    top_df       = cosmic_state["top_df"]
    top_n_sig    = cosmic_state["top_n_sig"]
    cos_sim      = cosmic_state["cos_sim"]
    residual_pct = cosmic_state["residual_pct"]
    fit_method   = cosmic_state["fit_method"]

    st.divider()
    st.subheader("🌌 COSMIC Signature Decomposition")
    fit_col1, fit_col2 = st.columns(2)
    fit_col1.metric(
        "Cosine similarity", f"{cos_sim:.4f}",
        help="1.0 = perfect reconstruction using the top-N signatures. Values above 0.95 indicate a good fit.",
    )
    fit_col2.metric(
        "Residual (% of counts)", f"{residual_pct:.1f}%",
        help="L2 norm of unexplained counts as a percentage of total SNV count (refit to top-N). Lower is better.",
    )

    all_sigs = list(sig_df["signature"])
    sig_chart = (
        alt.Chart(top_df)
        .mark_bar()
        .encode(
            alt.X("signature:N", sort=list(top_df["signature"]), title="Signature"),
            alt.Y("exposure:Q", title="Exposure (proportion)",
                  axis=alt.Axis(format=".0%")),
            alt.Color("signature:N", scale=alt.Scale(domain=all_sigs), legend=None),
            tooltip=[
                "signature:N",
                alt.Tooltip("exposure:Q", format=".2%", title="Exposure"),
                alt.Tooltip("etiology:N", title="Etiology"),
            ],
        )
        .properties(
            title=(
                f"Top {top_n_sig} COSMIC SBS Signatures (NNLS fit)"
                if fit_method == "NNLS (top-N)" else
                f"COSMIC SBS Signatures — Greedy add/remove ({top_n_sig} active)"
            ),
            height=300,
        )
    )
    st.altair_chart(sig_chart, width="stretch")

    display = top_df.copy()
    display["exposure"] = display["exposure"].map("{:.2%}".format)
    st.dataframe(display, width="stretch", hide_index=True)

    if ctx.data_source.is_duckdb:
        _render_per_sample_exposures(ctx, cosmic_state)


def _render_per_sample_exposures(ctx, cosmic_state):
    cosmic_W       = cosmic_state["cosmic_W"]
    cosmic_aligned = cosmic_state["cosmic_aligned"]
    fit_method     = cosmic_state["fit_method"]
    add_penalty    = cosmic_state["add_penalty"]
    remove_penalty = cosmic_state["remove_penalty"]
    use_connected  = cosmic_state["use_connected"]

    st.divider()
    st.subheader("Per-sample Signature Exposures")
    if fit_method == "NNLS (top-N)":
        ps_caption = (
            "Each sample fitted independently against the full COSMIC matrix "
            "with NNLS (no top-N restriction — all non-zero exposures shown)."
        )
    else:
        ps_caption = (
            "Each sample fitted independently against the full COSMIC matrix "
            "with greedy add/remove (SPA-style)."
        )
    st.caption(f"{ps_caption} Signatures with no exposure in any sample are hidden.")

    sample_matrix = _load_sample_sbs96_matrix(ctx)
    if sample_matrix.empty:
        st.info("No SNVs with trinucleotide context in current selection.")
        return

    ps_rows = []
    cohort_result = None
    if fit_method == "Greedy add/remove (SPA-style)":
        cohort_result = fit_cosmic_cohort_greedy(
            sample_matrix, cosmic_aligned,
            add_penalty=add_penalty,
            remove_penalty=remove_penalty,
            connected_sigs=use_connected,
        )
        for sid, frac_row in cohort_result["exposure_fractions"].iterrows():
            n_snvs = int(sample_matrix.loc[sid].sum())
            if n_snvs == 0:
                continue
            for sig, exp in frac_row.items():
                if exp > 0:
                    ps_rows.append({
                        "sample_label": sid,
                        "signature": sig,
                        "exposure":  float(exp),
                        "n_snvs":    n_snvs,
                        "etiology":  SBS_ETIOLOGY.get(sig, ""),
                    })
    else:
        for sid, row in sample_matrix.iterrows():
            obs = row.reindex(SBS_ORDER).values.astype(float)
            n_snvs = int(row.sum())
            if n_snvs == 0:
                continue
            h, _ = nnls(cosmic_W, obs)
            total = h.sum()
            h_norm = h / total if total > 0 else h
            for sig, exp in zip(cosmic_aligned.columns, h_norm):
                if exp > 0:
                    ps_rows.append({
                        "sample_label": sid,
                        "signature": sig,
                        "exposure":  float(exp),
                        "n_snvs":    n_snvs,
                        "etiology":  SBS_ETIOLOGY.get(sig, ""),
                    })

    if not ps_rows:
        st.info("No exposures returned for any sample.")
        return

    ps_df = pd.DataFrame(ps_rows)
    ps_sigs  = ps_df["signature"].unique().tolist()
    ps_order = sorted(ps_df["sample_label"].unique().tolist())

    ps_full = (
        pd.MultiIndex.from_product(
            [ps_order, ps_sigs],
            names=["sample_label", "signature"],
        )
        .to_frame(index=False)
        .merge(ps_df[["sample_label", "signature", "exposure", "etiology"]],
               on=["sample_label", "signature"], how="left")
    )
    ps_full["exposure"] = ps_full["exposure"].fillna(0.0)
    ps_full["etiology"] = ps_full["etiology"].fillna("")

    ps_chart = (
        alt.Chart(ps_full)
        .mark_rect()
        .encode(
            alt.X("signature:N", sort=ps_sigs, title="Signature",
                  axis=alt.Axis(labelAngle=-45, labelLimit=200)),
            alt.Y("sample_label:N", sort=ps_order, title="Sample"),
            alt.Color("exposure:Q", title="Exposure",
                      scale=alt.Scale(scheme="blues"),
                      legend=alt.Legend(format=".0%")),
            tooltip=[
                "sample_label:N", "signature:N",
                alt.Tooltip("exposure:Q", format=".2%", title="Exposure"),
                "etiology:N",
            ],
        )
        .properties(
            title=(
                "Per-sample COSMIC SBS signature exposures (NNLS)"
                if fit_method == "NNLS (top-N)" else
                "Per-sample COSMIC SBS signature exposures (greedy add/remove)"
            ),
            height=max(150, 22 * len(ps_order)),
        )
    )
    st.altair_chart(ps_chart, width="stretch")
    st.caption(
        f"{len(ps_order)} samples · {len(ps_sigs)} active signatures. "
        "Color intensity = exposure proportion. "
        "Signatures with no exposure in any sample are hidden."
    )
    if fit_method == "Greedy add/remove (SPA-style)" and cohort_result is not None:
        metrics_df = cohort_result["per_sample_metrics"]
        st.caption("**Per-sample fit quality**")
        st.dataframe(
            metrics_df.style.format({
                "l2_distance": "{:.3f}",
                "cosine_similarity": "{:.4f}",
            }),
            width="stretch", hide_index=True,
        )


# ─────────────────────────────────────────────────────────────────────────────
# NMF discovery
# ─────────────────────────────────────────────────────────────────────────────


def _render_nmf_discovery(ctx, cosmic_path, cosmic_state):
    st.divider()
    st.subheader("NMF Signature Discovery")
    st.caption(
        "Learn de novo SBS96 signatures across the currently selected samples. "
        "This uses the current sidebar filters and requires at least two samples with SNVs."
    )

    enable_nmf = st.checkbox(
        "Run NMF signature discovery",
        value=False, key="run_nmf_signatures",
    )
    if not enable_nmf:
        return

    sample_matrix = _load_sample_sbs96_matrix(ctx)
    n_samples_nmf = sample_matrix.shape[0]

    if n_samples_nmf < 2:
        st.info("NMF requires at least two samples with SNVs in the current selection.")
        return

    cosmic_aligned = cosmic_state["cosmic_aligned"] if cosmic_state else None
    sig_df         = cosmic_state["sig_df"]         if cosmic_state else None

    hybrid_available = (
        cosmic_aligned is not None
        and sig_df is not None
        and not sig_df.empty
    )
    if hybrid_available:
        nmf_mode = st.radio(
            "Discovery mode",
            ["COSMIC-guided + one learned signature", "De novo NMF"],
            horizontal=True, key="nmf_mode",
            help="Use top COSMIC signatures as fixed basis functions and learn one additional non-negative signature, or run fully de novo NMF.",
        )
    else:
        nmf_mode = "De novo NMF"

    if nmf_mode == "De novo NMF":
        _render_de_novo_nmf(ctx, sample_matrix, n_samples_nmf, cosmic_aligned, cosmic_path)
    else:
        _render_cosmic_guided_nmf(
            ctx, sample_matrix, n_samples_nmf, cosmic_aligned, sig_df, cosmic_path,
        )


def _render_de_novo_nmf(ctx, sample_matrix, n_samples_nmf, cosmic_aligned, cosmic_path):
    nmf_max     = min(8, n_samples_nmf)
    nmf_default = min(3, nmf_max)
    nmf_components = st.slider(
        "Number of discovered signatures",
        2, nmf_max, nmf_default, key="nmf_components",
    )

    try:
        nmf_result = fit_sbs_nmf(sample_matrix, nmf_components)
    except Exception as exc:
        st.error(f"Could not run NMF signature discovery: {exc}")
        return

    nmf_profiles            = nmf_result["profiles"]
    nmf_exposure_fractions  = nmf_result["exposure_fractions"]
    nmf_matrix_cosine       = nmf_result["matrix_cosine"]
    nmf_relative_error_pct  = nmf_result["relative_error_pct"]
    nmf_iter                = nmf_result["n_iter"]

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Samples", f"{n_samples_nmf:,}")
    m2.metric("Signatures", str(nmf_components))
    m3.metric("Matrix cosine", f"{nmf_matrix_cosine:.4f}")
    m4.metric("Relative error", f"{nmf_relative_error_pct:.1f}%")
    st.caption(
        f"NMF converged in {nmf_iter} iteration(s). "
        "Signatures are shown as normalized SBS96 profiles."
    )

    nmf_long = (
        nmf_exposure_fractions.reset_index(names="sample_label")
        .melt(id_vars="sample_label", var_name="signature", value_name="exposure")
    )
    sig_order    = nmf_profiles.index.tolist()
    sample_order = sorted(nmf_exposure_fractions.index.tolist())
    nmf_chart = (
        alt.Chart(nmf_long)
        .mark_rect()
        .encode(
            alt.X("signature:N", sort=sig_order, title="Discovered signature"),
            alt.Y("sample_label:N", sort=sample_order, title="Sample"),
            alt.Color("exposure:Q", title="Exposure",
                      scale=alt.Scale(scheme="oranges"),
                      legend=alt.Legend(format=".0%")),
            tooltip=[
                "sample_label:N", "signature:N",
                alt.Tooltip("exposure:Q", format=".2%", title="Exposure"),
            ],
        )
        .properties(
            title="Per-sample NMF signature exposures",
            height=max(150, 22 * len(sample_order)),
        )
    )
    st.altair_chart(nmf_chart, width="stretch")

    nmf_match_df = None
    if cosmic_aligned is not None:
        try:
            nmf_match_df, _ = compare_signatures_to_cosmic(nmf_profiles, cosmic_aligned)
        except Exception as exc:
            st.warning(f"Could not compare NMF signatures to COSMIC: {exc}")
        else:
            nmf_match_df["etiology"] = nmf_match_df["most_similar_cosmic_signature"].map(
                lambda s: SBS_ETIOLOGY.get(s, "")
            )
            display = nmf_match_df.copy()
            display["most_similar_cosine_similarity"] = display[
                "most_similar_cosine_similarity"
            ].map("{:.3f}".format)
            st.markdown("**Best COSMIC match for each discovered signature**")
            st.dataframe(display, width="stretch", hide_index=True)
    elif cosmic_path and cosmic_path.strip():
        st.info("COSMIC comparison is unavailable because the matrix did not load successfully.")
    else:
        st.info("Provide a COSMIC SBS matrix path above to compare discovered signatures.")

    sig_tabs = st.tabs(sig_order)
    for tab, sig_name in zip(sig_tabs, sig_order):
        with tab:
            sig_spec = _profile_to_spec96(
                nmf_profiles.loc[sig_name], value_name="fraction",
            )
            overlay_spec  = None
            overlay_label = None
            caption = "De novo SBS96 signature learned by NMF."
            if nmf_match_df is not None and not nmf_match_df.empty:
                match_row = (
                    nmf_match_df[nmf_match_df["signature"] == sig_name].iloc[0]
                )
                cosmic_sig = match_row["most_similar_cosmic_signature"]
                overlay_spec = _profile_to_spec96(
                    cosmic_aligned[cosmic_sig], value_name="fraction",
                )
                overlay_label = cosmic_sig
                caption = (
                    f"Best COSMIC match: {cosmic_sig} "
                    f"(cosine {match_row['most_similar_cosine_similarity']:.3f})."
                )
            st.altair_chart(
                _signature_profile_chart(
                    sig_spec, sig_name,
                    overlay_df=overlay_spec, overlay_label=overlay_label,
                ),
                width="stretch",
            )
            st.caption(caption)


def _render_cosmic_guided_nmf(
    ctx, sample_matrix, n_samples_nmf, cosmic_aligned, sig_df, cosmic_path,
):
    fixed_max     = min(10, len(sig_df))
    fixed_default = min(4, fixed_max)
    n_fixed = st.slider(
        "Number of fixed COSMIC signatures",
        1, fixed_max, fixed_default, key="nmf_fixed_cosmic",
        help="Use the top N COSMIC signatures from the current cohort-level fit as fixed basis signatures, then learn one additional non-negative residual signature.",
    )
    fixed_signature_names = sig_df.head(n_fixed)["signature"].tolist()
    st.caption("Fixed COSMIC signatures: " + ", ".join(fixed_signature_names))

    try:
        guided_result = fit_cosmic_augmented_nmf(
            sample_matrix, cosmic_aligned, fixed_signature_names,
        )
    except Exception as exc:
        st.error(f"Could not run COSMIC-guided NMF discovery: {exc}")
        return

    guided_profiles            = guided_result["profiles"]
    guided_exposure_fractions  = guided_result["exposure_fractions"]
    guided_iter                = guided_result["n_iter"]
    guided_matrix_cosine       = guided_result["matrix_cosine"]
    guided_relative_error_pct  = guided_result["relative_error_pct"]
    guided_fixed_relative_error_pct = guided_result["fixed_only_relative_error_pct"]
    guided_improvement_pct     = guided_result["relative_error_improvement_pct"]
    learned_names              = list(guided_result["learned_signature_names"])

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Samples", f"{n_samples_nmf:,}")
    m2.metric("Fixed COSMIC", str(n_fixed))
    m3.metric("Matrix cosine", f"{guided_matrix_cosine:.4f}")
    m4.metric("Residual improvement", f"{guided_improvement_pct:.2f}%")
    st.caption(
        f"Alternating constrained fit converged in {guided_iter} iteration(s). "
        f"Fixed-only relative error was {guided_fixed_relative_error_pct:.2f}%; "
        f"adding the learned signature reduced it to {guided_relative_error_pct:.2f}%."
    )

    guided_long = (
        guided_exposure_fractions.reset_index(names="sample_label")
        .melt(id_vars="sample_label", var_name="signature", value_name="exposure")
    )
    guided_order = guided_profiles.index.tolist()
    sample_order = sorted(guided_exposure_fractions.index.tolist())
    guided_chart = (
        alt.Chart(guided_long)
        .mark_rect()
        .encode(
            alt.X("signature:N", sort=guided_order, title="Signature component"),
            alt.Y("sample_label:N", sort=sample_order, title="Sample"),
            alt.Color("exposure:Q", title="Exposure",
                      scale=alt.Scale(scheme="oranges"),
                      legend=alt.Legend(format=".0%")),
            tooltip=[
                "sample_label:N", "signature:N",
                alt.Tooltip("exposure:Q", format=".2%", title="Exposure"),
            ],
        )
        .properties(
            title="Per-sample COSMIC-guided signature exposures",
            height=max(150, 22 * len(sample_order)),
        )
    )
    st.altair_chart(guided_chart, width="stretch")

    guided_match_df = None
    try:
        guided_match_df, _ = compare_signatures_to_cosmic(
            guided_profiles.loc[learned_names], cosmic_aligned,
        )
    except Exception as exc:
        st.warning(f"Could not compare the learned signature to COSMIC: {exc}")
    else:
        guided_match_df["etiology"] = guided_match_df["most_similar_cosmic_signature"].map(
            lambda s: SBS_ETIOLOGY.get(s, "")
        )
        display = guided_match_df.copy()
        display["most_similar_cosine_similarity"] = display[
            "most_similar_cosine_similarity"
        ].map("{:.3f}".format)
        st.markdown("**Best COSMIC match for the learned residual signature**")
        st.dataframe(display, width="stretch", hide_index=True)

    novel_sig_name = learned_names[0]
    novel_match_row = None
    if guided_match_df is not None and not guided_match_df.empty:
        novel_match_row = (
            guided_match_df[guided_match_df["signature"] == novel_sig_name].iloc[0]
        )
    novel_download_df = build_signature_download_table(
        _profile_to_spec96(guided_profiles.loc[novel_sig_name], value_name="fraction"),
        signature_name=novel_sig_name,
        most_similar_cosmic_signature=(
            None if novel_match_row is None else novel_match_row["most_similar_cosmic_signature"]
        ),
        most_similar_cosine_similarity=(
            None if novel_match_row is None
            else float(novel_match_row["most_similar_cosine_similarity"])
        ),
        fixed_signature_names=fixed_signature_names,
    )
    novel_provenance_df = ctx.build_provenance(
        discovery_mode="COSMIC-guided + one learned signature",
        discovery_items=[
            ("cosmic_matrix_path", cosmic_path.strip()),
            ("fixed_cosmic_count", n_fixed),
            ("fixed_cosmic_signatures", fixed_signature_names),
            (
                "most_similar_cosmic_signature",
                None if novel_match_row is None else novel_match_row["most_similar_cosmic_signature"],
            ),
            (
                "most_similar_cosine_similarity",
                None if novel_match_row is None
                else f"{float(novel_match_row['most_similar_cosine_similarity']):.6f}",
            ),
        ],
    )
    novel_exposures_df = build_signature_exposure_download_table(guided_exposure_fractions)
    novel_bundle = build_signature_download_zip(
        novel_download_df, novel_provenance_df,
        signature_name=novel_sig_name,
        sample_exposures_df=novel_exposures_df,
    )
    st.download_button(
        "Download novel signature bundle (.zip)",
        data=novel_bundle,
        file_name=f"{novel_sig_name.lower()}_signature_bundle.zip",
        mime="application/zip",
        key="nmf_download_novel_signature",
        help="Download a zip containing the learned SBS96 signature, the per-sample guided exposure table, and a provenance table of active filters and discovery settings.",
    )

    sig_tabs = st.tabs(guided_order)
    for tab, sig_name in zip(sig_tabs, guided_order):
        with tab:
            sig_spec = _profile_to_spec96(
                guided_profiles.loc[sig_name], value_name="fraction",
            )
            if sig_name in fixed_signature_names:
                st.altair_chart(
                    _signature_profile_chart(sig_spec, sig_name),
                    width="stretch",
                )
                st.caption("Fixed COSMIC signature used in the constrained fit.")
            else:
                overlay_spec  = None
                overlay_label = None
                caption = "Learned non-negative residual signature."
                if guided_match_df is not None and not guided_match_df.empty:
                    match_row = (
                        guided_match_df[guided_match_df["signature"] == sig_name].iloc[0]
                    )
                    cosmic_sig = match_row["most_similar_cosmic_signature"]
                    overlay_spec = _profile_to_spec96(
                        cosmic_aligned[cosmic_sig], value_name="fraction",
                    )
                    overlay_label = cosmic_sig
                    caption = (
                        f"Best COSMIC match: {cosmic_sig} "
                        f"(cosine {match_row['most_similar_cosine_similarity']:.3f})."
                    )
                st.altair_chart(
                    _signature_profile_chart(
                        sig_spec, sig_name,
                        overlay_df=overlay_spec, overlay_label=overlay_label,
                    ),
                    width="stretch",
                )
                st.caption(caption)


# ─────────────────────────────────────────────────────────────────────────────
# Called vs Uncalled
# ─────────────────────────────────────────────────────────────────────────────


def _render_called_vs_uncalled(ctx, cosmic_state):
    st.divider()
    st.subheader("Called vs Uncalled Comparison")
    st.caption(
        "Compares the SBS96 trinucleotide spectrum and COSMIC signature exposures "
        "between loci where a variant was called and where it was not. "
        "Requires the variant_called column (provide a VCF or variants TSV at collect time)."
    )

    def _build_spectrum(called_val):
        raw = ctx.con.execute(f"""
            SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
            FROM {ctx.table_expr}
            WHERE {ctx.where}
              AND variant_type = 'SNV'
              AND trinuc_context IS NOT NULL
              AND length(trinuc_context) = 3
              AND variant_called IS {'TRUE' if called_val else 'FALSE'}
            GROUP BY trinuc_context, ref_allele, alt_allele
        """).df()
        if raw.empty:
            return np.zeros(96, dtype=float), 0
        raw["sbs_label"] = raw.apply(
            lambda row: sbs_label(row["trinuc_context"], row["ref_allele"], row["alt_allele"]),
            axis=1,
        )
        raw = raw.dropna(subset=["sbs_label"])
        agg = raw.groupby("sbs_label", as_index=False)["count"].sum()
        obs = (
            pd.Series(0, index=SBS_ORDER, dtype=float)
            .add(agg.set_index("sbs_label")["count"], fill_value=0)
            .reindex(SBS_ORDER)
            .values.astype(float)
        )
        return obs, int(obs.sum())

    obs_called,   n_called   = _build_spectrum(True)
    obs_uncalled, n_uncalled = _build_spectrum(False)

    if n_called == 0 and n_uncalled == 0:
        st.info("No SNVs with trinucleotide context found in either group.")
        return

    called_label   = f"Called (n={n_called:,})"
    uncalled_label = f"Uncalled (n={n_uncalled:,})"

    def _make_spec96_df(obs_arr, n_total):
        df_s = pd.DataFrame({
            "sbs_label": SBS_ORDER,
            "mut_type":  [lbl[2:5] for lbl in SBS_ORDER],
            "count":     obs_arr.astype(int),
        })
        df_s["fraction"] = df_s["count"] / n_total if n_total > 0 else 0.0
        return df_s

    spec_called   = _make_spec96_df(obs_called,   n_called)
    spec_uncalled = _make_spec96_df(obs_uncalled, n_uncalled)

    m_df = pd.concat([
        spec_called.assign(y=spec_called["fraction"],     group=called_label),
        spec_uncalled.assign(y=-spec_uncalled["fraction"], group=uncalled_label),
    ])
    mirror_sub = []
    for mt in SBS_MUT_TYPES:
        sub = m_df[m_df["mut_type"] == mt]
        order = [lbl for lbl in SBS_ORDER if f"[{mt}]" in lbl]
        c = (
            alt.Chart(sub)
            .mark_bar()
            .encode(
                alt.X("sbs_label:N", sort=order, title=None,
                      axis=alt.Axis(labelAngle=-90, labelFontSize=7)),
                alt.Y("y:Q", axis=alt.Axis(format=".0%"),
                      title="← Uncalled | Called →"),
                alt.Color("group:N", title=None,
                          scale=alt.Scale(
                              domain=[called_label, uncalled_label],
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
                title=alt.TitleParams(mt, color=SBS_COLORS[mt], fontSize=11, fontWeight="bold"),
                width=130, height=150,
            )
        )
        mirror_sub.append(c)
    st.altair_chart(
        alt.concat(*mirror_sub, columns=3)
        .resolve_scale(y="shared")
        .properties(title=alt.TitleParams(
            f"Called vs Uncalled — mirrored trinucleotide spectrum  ·  {called_label} / {uncalled_label}",
            fontSize=13,
        )),
        width="stretch",
    )
    st.caption(
        "Called variants point up (blue), uncalled loci point down (red). "
        "Each group is normalised to its own fraction so differences in count don't dominate."
    )

    cosmic_W = cosmic_state["cosmic_W"] if cosmic_state else None
    if cosmic_W is None:
        return

    cosmic_aligned = cosmic_state["cosmic_aligned"]
    st.markdown("**COSMIC signature comparison — Called vs Uncalled**")

    def _fit_cmp(obs):
        h, _ = nnls(cosmic_W, obs)
        total = h.sum()
        h_norm = h / total if total > 0 else h
        recon  = cosmic_W @ h
        cos    = (
            float(np.dot(obs, recon))
            / (np.linalg.norm(obs) * np.linalg.norm(recon) + 1e-12)
        )
        return h_norm, cos

    h_called,   cos_called   = _fit_cmp(obs_called)
    h_uncalled, cos_uncalled = _fit_cmp(obs_uncalled)

    gof_c1, gof_c2, gof_c3, gof_c4 = st.columns(4)
    gof_c1.metric("Called SNVs",           f"{n_called:,}")
    gof_c2.metric("Cosine sim (called)",   f"{cos_called:.4f}")
    gof_c3.metric("Uncalled SNVs",         f"{n_uncalled:,}")
    gof_c4.metric("Cosine sim (uncalled)", f"{cos_uncalled:.4f}")

    sig_names = cosmic_aligned.columns.tolist()
    cmp_df = pd.DataFrame({
        "signature": sig_names * 2,
        "group":     [called_label] * len(sig_names) + [uncalled_label] * len(sig_names),
        "exposure":  list(h_called) + list(h_uncalled),
    })
    cmp_df["etiology"] = cmp_df["signature"].map(lambda s: SBS_ETIOLOGY.get(s, ""))

    cmp_top_n = st.slider(
        "Top signatures to display (comparison)",
        3, min(20, len(sig_names)), 8, key="cmp_top_n",
    )
    max_exp = (
        cmp_df.groupby("signature")["exposure"].max()
        .sort_values(ascending=False)
    )
    top_sigs = max_exp.head(cmp_top_n).index.tolist()
    top_cmp_df = cmp_df[cmp_df["signature"].isin(top_sigs)].copy()

    cmp_bars = (
        alt.Chart(top_cmp_df)
        .mark_bar()
        .encode(
            alt.X("signature:N", sort=top_sigs, title="Signature"),
            alt.Y("exposure:Q", title="Exposure (proportion)",
                  axis=alt.Axis(format=".0%")),
            alt.Color("group:N", title=None,
                      scale=alt.Scale(
                          domain=[called_label, uncalled_label],
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
    cmp_dividers = (
        alt.Chart(pd.DataFrame({"signature": top_sigs[1:]}))
        .mark_rule(color="#888", strokeWidth=1, opacity=0.5)
        .encode(alt.X("signature:N", sort=top_sigs, bandPosition=0))
    )
    st.altair_chart(
        alt.layer(cmp_bars, cmp_dividers)
        .properties(
            title=f"Top {cmp_top_n} COSMIC SBS Signatures — Called vs Uncalled",
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
        pivot = (
            cmp_df[cmp_df["exposure"] > 0]
            .pivot(index="signature", columns="group", values="exposure")
            .fillna(0)
            .reset_index()
        )
        for col in [called_label, uncalled_label]:
            if col in pivot.columns:
                pivot[col] = pivot[col].map("{:.2%}".format)
        pivot["etiology"] = pivot["signature"].map(lambda s: SBS_ETIOLOGY.get(s, ""))
        st.dataframe(pivot, width="stretch", hide_index=True)


# ─────────────────────────────────────────────────────────────────────────────
# VAF and R1/R2 stratified spectra
# ─────────────────────────────────────────────────────────────────────────────


def _render_vaf_stratified(ctx):
    st.divider()
    st.subheader("VAF-stratified Spectrum")
    st.caption(
        "Germline (VAF > 30%) vs somatic (VAF ≤ 30%) trinucleotide spectra. "
        "Differences between the two profiles reveal what drives low-VAF calls."
    )

    germ_raw = ctx.con.execute(f"""
        SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
        FROM (SELECT * FROM {ctx.table_expr}) _t
        WHERE {ctx.where} AND variant_type = 'SNV'
          AND trinuc_context IS NOT NULL AND length(trinuc_context) = 3
          AND alt_count * 1.0 / total_depth > 0.3
        GROUP BY trinuc_context, ref_allele, alt_allele
    """).df()
    som_raw = ctx.con.execute(f"""
        SELECT trinuc_context, ref_allele, alt_allele, COUNT(*) AS count
        FROM (SELECT * FROM {ctx.table_expr}) _t
        WHERE {ctx.where} AND variant_type = 'SNV'
          AND trinuc_context IS NOT NULL AND length(trinuc_context) = 3
          AND alt_count * 1.0 / total_depth <= 0.3
        GROUP BY trinuc_context, ref_allele, alt_allele
    """).df()

    germ_s96, n_germ = to_spec96_strat(germ_raw)
    som_s96,  n_som  = to_spec96_strat(som_raw)

    vc1, vc2 = st.columns(2)
    with vc1:
        if germ_s96 is not None:
            st.altair_chart(
                strat_sbs96_chart(germ_s96, f"Germline VAF > 30% (n={n_germ:,})"),
                width="stretch",
            )
        else:
            st.info("No germline SNVs in current selection.")
    with vc2:
        if som_s96 is not None:
            st.altair_chart(
                strat_sbs96_chart(som_s96, f"Somatic VAF ≤ 30% (n={n_som:,})"),
                width="stretch",
            )
        else:
            st.info("No somatic SNVs in current selection.")


def _render_r1r2_stratified(ctx):
    st.divider()
    st.subheader("R1 / R2 Stratified Spectrum")
    st.caption(
        "Trinucleotide spectra for Read 1 vs Read 2. "
        "Differences between the two profiles indicate read-level artefacts "
        "or strand-specific damage patterns (e.g. oxidative damage on R2)."
    )
    if not ctx.has_alt_reads:
        st.info("R1/R2 spectrum requires the alt_reads table (DuckDB mode only).")
        return

    r_reads_filter = f"WHERE {ctx.reads_where}" if ctx.reads_where else ""
    locus_snv = f"""
        SELECT sample_id, chrom, pos, alt_allele,
               trinuc_context, ref_allele
        FROM {ctx.table_expr}
        WHERE {ctx.where} AND variant_type = 'SNV'
          AND trinuc_context IS NOT NULL AND length(trinuc_context) = 3
    """
    r1_raw = ctx.con.execute(f"""
        SELECT ab.trinuc_context, ab.ref_allele, ab.alt_allele,
               COUNT(*) AS count
        FROM (SELECT * FROM alt_reads {r_reads_filter}) ar
        INNER JOIN ({locus_snv}) ab
        ON  ar.sample_id  = ab.sample_id
        AND ar.chrom      = ab.chrom
        AND ar.pos        = ab.pos
        AND ar.alt_allele = ab.alt_allele
        WHERE ar.is_read1 = TRUE
        GROUP BY ab.trinuc_context, ab.ref_allele, ab.alt_allele
    """).df()
    r2_raw = ctx.con.execute(f"""
        SELECT ab.trinuc_context, ab.ref_allele, ab.alt_allele,
               COUNT(*) AS count
        FROM (SELECT * FROM alt_reads {r_reads_filter}) ar
        INNER JOIN ({locus_snv}) ab
        ON  ar.sample_id  = ab.sample_id
        AND ar.chrom      = ab.chrom
        AND ar.pos        = ab.pos
        AND ar.alt_allele = ab.alt_allele
        WHERE ar.is_read1 = FALSE
        GROUP BY ab.trinuc_context, ab.ref_allele, ab.alt_allele
    """).df()

    for rdf in [r1_raw, r2_raw]:
        if not rdf.empty:
            rdf["sbs_label"] = rdf.apply(
                lambda r: sbs_label(r["trinuc_context"], r["ref_allele"], r["alt_allele"]),
                axis=1,
            )

    r1_s96, n_r1 = to_spec96_strat(r1_raw)
    r2_s96, n_r2 = to_spec96_strat(r2_raw)

    r12_y_max = None
    if r1_s96 is not None and r2_s96 is not None:
        r12_y_max = max(r1_s96["fraction"].max(), r2_s96["fraction"].max())

    ev_r1 = ev_r2 = None
    rc1, rc2 = st.columns(2)
    with rc1:
        if r1_s96 is not None:
            ev_r1 = st.altair_chart(
                strat_sbs96_chart(r1_s96, f"Read 1 (n={n_r1:,})", r12_y_max, sel_name="r1_click"),
                width="stretch", on_select="rerun", key="sbs96_r1",
            )
        else:
            st.info("No R1 SNVs in current selection.")
    with rc2:
        if r2_s96 is not None:
            ev_r2 = st.altair_chart(
                strat_sbs96_chart(r2_s96, f"Read 2 (n={n_r2:,})", r12_y_max, sel_name="r2_click"),
                width="stretch", on_select="rerun", key="sbs96_r2",
            )
        else:
            st.info("No R2 SNVs in current selection.")

    st.caption("Click a bar to drill down to matching loci and open in IGV.")

    _r12_drilldown(ctx, ev_r1, r1_raw, "Read 1", "r1")
    _r12_drilldown(ctx, ev_r2, r2_raw, "Read 2", "r2")


def _r12_drilldown(ctx, event, raw_df, read_label, sel_key_prefix):
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
    sel = ctx.query_records([extra_cond])
    label_str = ", ".join(clicked)
    st.caption(f"**{read_label}** — {len(sel):,} loci · context(s): {label_str}")
    st.dataframe(sel[ctx.table_cols], width="stretch")
    ctx.igv_buttons([extra_cond], sel, key=f"{sel_key_prefix}_sbs_{'_'.join(clicked)}")


# ─────────────────────────────────────────────────────────────────────────────
# Fallback (non-trinuc) and cohort heatmap
# ─────────────────────────────────────────────────────────────────────────────


def _render_simple_fallback(ctx):
    spec = ctx.con.execute(f"""
        SELECT ref_allele || '>' || alt_allele AS substitution, COUNT(*) AS count
        FROM {ctx.table_expr}
        WHERE {ctx.where} AND variant_type = 'SNV'
        GROUP BY substitution
        ORDER BY count DESC
    """).df()

    if spec.empty:
        st.info("No SNVs in current selection.")
        return

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
            sel = ctx.query_records([
                "variant_type = 'SNV'",
                f"ref_allele = '{ref}'",
                f"alt_allele = '{alt_allele}'",
            ])
            st.caption(f"{len(sel):,} records with substitution {sub}")
            st.dataframe(sel[ctx.table_cols], width="stretch")
            ctx.igv_buttons([
                "variant_type = 'SNV'",
                f"ref_allele = '{ref}'",
                f"alt_allele = '{alt_allele}'",
            ], sel, key=f"spectrum_{sub}")


def _render_cohort_heatmap(ctx):
    st.divider()
    st.subheader("SBS96 Heatmap (samples × trinucleotide contexts)")
    if not ctx.has_data("trinuc_context"):
        st.info("Trinucleotide context unavailable — run geac collect with a reference FASTA.")
        return

    has_batch = ctx.has_data("batch")
    id_sql    = "sample_id || ' / ' || batch" if has_batch else "sample_id"
    group_by  = "sample_id, batch" if has_batch else "sample_id"

    raw = ctx.con.execute(f"""
        SELECT {id_sql} AS sample_label,
               trinuc_context, ref_allele, alt_allele, COUNT(*) AS n
        FROM {ctx.table_expr}
        WHERE {ctx.where} AND variant_type = 'SNV'
          AND trinuc_context IS NOT NULL AND length(trinuc_context) = 3
        GROUP BY {group_by}, trinuc_context, ref_allele, alt_allele
    """).df()

    if raw.empty:
        st.info("No SNVs with trinucleotide context in current selection.")
        return

    raw["sbs_label"] = raw.apply(
        lambda row: sbs_label(row["trinuc_context"], row["ref_allele"], row["alt_allele"]),
        axis=1,
    )
    raw = raw.dropna(subset=["sbs_label"])
    agg = raw.groupby(["sample_label", "sbs_label"], as_index=False)["n"].sum()

    totals = agg.groupby("sample_label")["n"].transform("sum")
    agg["fraction"] = agg["n"] / totals

    all_combos = pd.MultiIndex.from_product(
        [agg["sample_label"].unique(), SBS_ORDER],
        names=["sample_label", "sbs_label"],
    )
    full = (
        agg.set_index(["sample_label", "sbs_label"])
        .reindex(all_combos, fill_value=0)
        .reset_index()
    )
    n_samples = agg["sample_label"].nunique()
    n_loci    = int(agg["n"].sum())

    full["mut_type"] = full["sbs_label"].str.extract(r'\[([A-Z]>[A-Z])\]')[0]

    hm_chart = (
        alt.Chart(full)
        .mark_rect()
        .encode(
            alt.X("sbs_label:N", sort=SBS_ORDER, title=None,
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
            height=max(200, 20 * full["sample_label"].nunique()),
            title="Normalised SBS96 profile per sample (fraction of SNVs)",
        )
    )

    hm_label_df = pd.DataFrame([
        {"sbs_label": [l for l in SBS_ORDER if f"[{mt}]" in l][8], "mut_type": mt}
        for mt in SBS_MUT_TYPES
    ])
    hm_label_strip = (
        alt.Chart(hm_label_df)
        .mark_text(align="center", fontSize=11, fontWeight="bold")
        .encode(
            alt.X("sbs_label:N", sort=SBS_ORDER,
                  axis=alt.Axis(labels=False, ticks=False, title=None)),
            alt.Y(value=15),
            alt.Color("mut_type:N", legend=None,
                      scale=alt.Scale(
                          domain=list(SBS_COLORS.keys()),
                          range=list(SBS_COLORS.values()),
                      )),
            alt.Text("mut_type:N"),
        )
        .properties(height=30)
    )
    st.altair_chart(
        alt.vconcat(hm_chart, hm_label_strip, spacing=2)
        .resolve_scale(x="shared"),
        width="stretch",
    )
    st.caption(
        f"{n_samples:,} samples · {n_loci:,} alt loci. "
        "Color = fraction of that sample's SNVs falling in each trinucleotide context. "
        "Contexts ordered by mutation type (C>A, C>G, C>T, T>A, T>C, T>G) then flanking bases."
    )
