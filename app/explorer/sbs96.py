"""SBS96 single-base-substitution helpers.

Used by the Error Spectrum, Cohort, and Duplex/Simplex tabs to render the
canonical 96-channel mutation spectrum (pyrimidine-normalised trinucleotide
context × substitution type), and to look up COSMIC signature etiologies.
"""
from __future__ import annotations

import altair as alt
import pandas as pd
import streamlit as st

COMP = str.maketrans('ACGT', 'TGCA')
SBS_MUT_TYPES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS_COLORS = {
    "C>A": "#1BBDEB", "C>G": "#808080", "C>T": "#E22926",
    "T>A": "#CBCACB", "T>C": "#97D54C", "T>G": "#ECC6C5",
}
SBS_ORDER = [
    f"{b5}[{mt}]{b3}"
    for mt in SBS_MUT_TYPES
    for b5 in "ACGT"
    for b3 in "ACGT"
]


def sbs_label(trinuc_context, ref_allele, alt_allele):
    """Convert raw trinuc/ref/alt into a pyrimidine-normalised SBS96 label."""
    ctx, r, a = trinuc_context, ref_allele, alt_allele
    if not all(b in 'ACGT' for b in (ctx + r + a)):
        return None
    if r in ('A', 'G'):
        ctx = ctx[::-1].translate(COMP)
        r = r.translate(COMP)
        a = a.translate(COMP)
    return f"{ctx[0]}[{r}>{a}]{ctx[2]}"


SBS_ETIOLOGY = {
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
def load_cosmic(p: str) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t", index_col=0)


def to_spec96_strat(raw_df):
    if raw_df.empty:
        return None, 0
    df = raw_df.copy()
    df["sbs_label"] = df.apply(
        lambda r: sbs_label(r["trinuc_context"], r["ref_allele"], r["alt_allele"]),
        axis=1,
    )
    df = df.dropna(subset=["sbs_label"])
    df["mut_type"] = df["sbs_label"].str.extract(r'\[([A-Z]>[A-Z])\]')[0]
    agg = df.groupby(["sbs_label", "mut_type"], as_index=False)["count"].sum()
    full = pd.DataFrame({
        "sbs_label": SBS_ORDER,
        "mut_type":  [lbl[2:5] for lbl in SBS_ORDER],
    })
    s96 = full.merge(agg, on=["sbs_label", "mut_type"], how="left")
    s96["count"] = s96["count"].fillna(0).astype(int)
    total = int(s96["count"].sum())
    s96["fraction"] = s96["count"] / total if total > 0 else 0.0
    return s96, total


def strat_sbs96_chart(spec_df, title, y_max=None, sel_name=None):
    _y_scale = alt.Scale(domain=[0, y_max]) if y_max is not None else alt.Undefined
    _sel = (
        alt.selection_point(name=sel_name, fields=["sbs_label"], on="click")
        if sel_name else None
    )
    panels = []
    for _mt in SBS_MUT_TYPES:
        _s = spec_df[spec_df["mut_type"] == _mt]
        _order = [lbl for lbl in SBS_ORDER if f"[{_mt}]" in lbl]
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
            alt.Chart(_s).mark_bar(color=SBS_COLORS[_mt]).encode(**_enc).properties(
                title=alt.TitleParams(_mt, color=SBS_COLORS[_mt],
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
