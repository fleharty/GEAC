import pandas as pd


def build_unique_pipeline_characterization_df(
    pc_uniq: pd.DataFrame, pipe_a: str, pipe_b: str
) -> pd.DataFrame:
    """Normalize unique-to-pipeline rows into one comparison dataframe."""
    cols = [
        "group",
        "pipeline",
        "sample_id",
        "chrom",
        "pos",
        "alt_allele",
        "variant_type",
        "trinuc_context",
        "ref_allele",
        "vaf",
        "depth",
        "alt_count",
    ]
    if pc_uniq.empty:
        return pd.DataFrame(columns=cols)

    frames = []
    for concordance, pipeline, vaf_col, depth_col, alt_count_col in [
        ("only_a", str(pipe_a), "vaf_a", "depth_a", "alt_count_a"),
        ("only_b", str(pipe_b), "vaf_b", "depth_b", "alt_count_b"),
    ]:
        sub = pc_uniq[pc_uniq["concordance"] == concordance].copy()
        if sub.empty:
            continue
        sub["group"] = f"Only {pipeline}"
        sub["pipeline"] = pipeline
        sub["vaf"] = sub[vaf_col]
        sub["depth"] = sub[depth_col]
        sub["alt_count"] = sub[alt_count_col]
        frames.append(sub[cols])

    if not frames:
        return pd.DataFrame(columns=cols)
    return pd.concat(frames, ignore_index=True)


def summarize_unique_pipeline_groups(
    uniq_cmp: pd.DataFrame, pipe_a: str, pipe_b: str
) -> dict[str, object]:
    """Return per-group summary metrics and a compact descriptive caption."""
    groups = [f"Only {pipe_a}", f"Only {pipe_b}"]
    metrics = {}
    for group in groups:
        sub = uniq_cmp[uniq_cmp["group"] == group]
        metrics[group] = {
            "count": int(len(sub)),
            "median_vaf": float(sub["vaf"].median()) if sub["vaf"].notna().any() else None,
            "median_depth": float(sub["depth"].median()) if sub["depth"].notna().any() else None,
            "median_alt_count": float(sub["alt_count"].median())
            if sub["alt_count"].notna().any()
            else None,
        }

    a_group = groups[0]
    b_group = groups[1]
    a = metrics[a_group]
    b = metrics[b_group]

    if a["count"] == 0 and b["count"] == 0:
        summary = "No unique loci under the current filter selection."
    elif a["count"] == 0:
        summary = f"All unique loci under the current filters are specific to {pipe_b}."
    elif b["count"] == 0:
        summary = f"All unique loci under the current filters are specific to {pipe_a}."
    else:
        parts = []
        if a["median_vaf"] is not None and b["median_vaf"] is not None:
            if a["median_vaf"] < b["median_vaf"]:
                parts.append(f"{pipe_a} unique calls have lower median VAF")
            elif a["median_vaf"] > b["median_vaf"]:
                parts.append(f"{pipe_a} unique calls have higher median VAF")
        if a["median_depth"] is not None and b["median_depth"] is not None:
            if a["median_depth"] < b["median_depth"]:
                parts.append(f"shallower median depth")
            elif a["median_depth"] > b["median_depth"]:
                parts.append(f"deeper median depth")
        if a["median_alt_count"] is not None and b["median_alt_count"] is not None:
            if a["median_alt_count"] < b["median_alt_count"]:
                parts.append(f"and lower median alt support than {pipe_b}")
            elif a["median_alt_count"] > b["median_alt_count"]:
                parts.append(f"and higher median alt support than {pipe_b}")
        summary = " ".join(parts) + "." if parts else "Unique-loci support is broadly similar between pipelines."

    return {"metrics": metrics, "summary": summary}
