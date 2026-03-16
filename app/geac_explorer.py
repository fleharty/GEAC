import marimo

__generated_with = "0.10.0"
app = marimo.App(width="wide", title="GEAC Explorer")


@app.cell
def imports():
    import marimo as mo
    import duckdb
    import altair as alt
    import pandas as pd
    return alt, duckdb, mo, pd


@app.cell
def header(mo):
    mo.md(
        """
        # GEAC Explorer
        **Genomic Evidence Atlas of Cohorts** — inspect alt base metrics from
        per-sample Parquet files or a merged cohort DuckDB database.
        """
    )
    return


@app.cell
def file_input(mo):
    file_path = mo.ui.text(
        placeholder="e.g. sample.parquet  or  cohort.duckdb",
        label="Data file (Parquet or DuckDB)",
        full_width=True,
    )
    return (file_path,)


@app.cell
def load_data(file_path, duckdb, mo):
    mo.stop(
        not file_path.value.strip(),
        mo.callout(mo.md("Enter a Parquet or DuckDB file path above to begin."), kind="info"),
    )

    path = file_path.value.strip()

    if path.endswith(".duckdb"):
        con = duckdb.connect(path, read_only=True)
        table_expr = "alt_bases"
    else:
        con = duckdb.connect()
        table_expr = f"read_parquet('{path}')"

    return con, table_expr


@app.cell
def summary_stats(con, table_expr, mo):
    stats = con.execute(f"""
        SELECT
            COUNT(*)                              AS n_records,
            COUNT(DISTINCT sample_id)             AS n_samples,
            COUNT(DISTINCT chrom)                 AS n_chromosomes,
            COUNT(DISTINCT chrom || ':' || pos)   AS n_positions,
            SUM(alt_count)                        AS total_alt_reads,
            ROUND(AVG(alt_count * 1.0 / total_depth), 4) AS mean_vaf
        FROM {table_expr}
    """).df()

    mo.hstack([
        mo.stat(str(int(stats["n_records"][0])),     label="Alt records"),
        mo.stat(str(int(stats["n_samples"][0])),     label="Samples"),
        mo.stat(str(int(stats["n_chromosomes"][0])), label="Chromosomes"),
        mo.stat(str(int(stats["n_positions"][0])),   label="Positions"),
        mo.stat(str(int(stats["total_alt_reads"][0])), label="Total alt reads"),
        mo.stat(str(stats["mean_vaf"][0]),           label="Mean VAF"),
    ], justify="start")
    return


@app.cell
def filters(con, table_expr, mo):
    chroms = (
        con.execute(f"SELECT DISTINCT chrom FROM {table_expr} ORDER BY chrom")
        .df()["chrom"]
        .tolist()
    )
    samples = (
        con.execute(f"SELECT DISTINCT sample_id FROM {table_expr} ORDER BY sample_id")
        .df()["sample_id"]
        .tolist()
    )

    chrom_filter = mo.ui.dropdown(
        options=["All"] + chroms, value="All", label="Chromosome"
    )
    sample_filter = mo.ui.multiselect(options=samples, label="Samples (blank = all)")
    variant_filter = mo.ui.multiselect(
        options=["SNV", "insertion", "deletion", "MNV"],
        value=["SNV", "insertion", "deletion", "MNV"],
        label="Variant type",
    )
    vaf_filter = mo.ui.range_slider(
        start=0.0, stop=1.0, step=0.01, value=[0.0, 1.0], label="VAF range"
    )
    min_alt_filter = mo.ui.number(start=1, stop=10000, step=1, value=1, label="Min alt count")

    mo.hstack([
        mo.vstack([chrom_filter, sample_filter]),
        mo.vstack([variant_filter, vaf_filter]),
        mo.vstack([min_alt_filter]),
    ], gap=2)
    return chrom_filter, min_alt_filter, sample_filter, variant_filter, vaf_filter


@app.cell
def filtered_data(
    chrom_filter, con, min_alt_filter, mo, sample_filter,
    table_expr, variant_filter, vaf_filter,
):
    conditions = [
        f"alt_count >= {min_alt_filter.value}",
        f"alt_count * 1.0 / total_depth BETWEEN {vaf_filter.value[0]} AND {vaf_filter.value[1]}",
    ]

    if chrom_filter.value != "All":
        conditions.append(f"chrom = '{chrom_filter.value}'")

    if sample_filter.value:
        s = ", ".join(f"'{s}'" for s in sample_filter.value)
        conditions.append(f"sample_id IN ({s})")

    if variant_filter.value:
        t = ", ".join(f"'{t}'" for t in variant_filter.value)
        conditions.append(f"variant_type IN ({t})")

    where = " AND ".join(conditions)

    df = con.execute(f"""
        SELECT *,
               ROUND(alt_count * 1.0 / total_depth, 4) AS vaf,
               ROUND(fwd_alt_count * 1.0 / NULLIF(alt_count, 0), 4) AS fwd_alt_frac
        FROM {table_expr}
        WHERE {where}
        LIMIT 50000
    """).df()

    mo.callout(
        mo.md(f"**{len(df):,}** records shown (capped at 50,000)"),
        kind="info" if len(df) < 50000 else "warn",
    )
    return (df,)


@app.cell
def data_table(df, mo):
    mo.ui.table(
        df[[
            "sample_id", "chrom", "pos", "ref_allele", "alt_allele",
            "variant_type", "vaf", "alt_count", "ref_count", "total_depth",
            "fwd_alt_count", "rev_alt_count", "overlap_alt_agree",
            "overlap_alt_disagree", "variant_called", "variant_filter",
        ]],
        selection=None,
        pagination=True,
    )
    return


@app.cell
def plots(alt, df, mo):
    snvs = df[df["variant_type"] == "SNV"].copy()

    # ── VAF distribution ──────────────────────────────────────────────────────
    vaf_chart = (
        alt.Chart(df)
        .mark_bar(opacity=0.8)
        .encode(
            alt.X("vaf:Q", bin=alt.Bin(maxbins=50), title="VAF"),
            alt.Y("count():Q", title="Count"),
            alt.Color("variant_type:N", title="Variant type"),
            tooltip=["variant_type:N", "count():Q"],
        )
        .properties(title="VAF Distribution", width=380, height=260)
    )

    # ── Error spectrum (SNV substitution types) ───────────────────────────────
    if len(snvs) > 0:
        spectrum_df = (
            snvs.groupby(["ref_allele", "alt_allele"])
            .size()
            .reset_index(name="count")
        )
        spectrum_df["substitution"] = (
            spectrum_df["ref_allele"] + ">" + spectrum_df["alt_allele"]
        )
        spectrum_chart = (
            alt.Chart(spectrum_df)
            .mark_bar()
            .encode(
                alt.X("substitution:N", sort="-y", title="Substitution"),
                alt.Y("count:Q", title="Count"),
                alt.Color("substitution:N", legend=None),
                tooltip=["substitution:N", "count:Q"],
            )
            .properties(title="SNV Error Spectrum", width=380, height=260)
        )
    else:
        spectrum_chart = alt.Chart(snvs).mark_text().encode(
            text=alt.value("No SNVs in current selection")
        ).properties(width=380, height=260)

    # ── Strand bias ───────────────────────────────────────────────────────────
    strand_chart = (
        alt.Chart(df.sample(min(2000, len(df))))
        .mark_point(opacity=0.5, size=30)
        .encode(
            alt.X("fwd_alt_count:Q", title="Forward alt reads"),
            alt.Y("rev_alt_count:Q", title="Reverse alt reads"),
            alt.Color("variant_type:N", title="Variant type"),
            tooltip=["sample_id", "chrom", "pos", "ref_allele", "alt_allele",
                     "fwd_alt_count", "rev_alt_count", "vaf"],
        )
        .properties(title="Strand Bias (up to 2,000 points)", width=380, height=260)
    )

    # ── Overlap agreement ─────────────────────────────────────────────────────
    overlap_df = df[df["overlap_depth"] > 0].copy()
    if len(overlap_df) > 0:
        overlap_df["agree_frac"] = (
            overlap_df["overlap_alt_agree"]
            / (overlap_df["overlap_alt_agree"] + overlap_df["overlap_alt_disagree"]).clip(lower=1)
        )
        overlap_chart = (
            alt.Chart(overlap_df)
            .mark_bar(opacity=0.8)
            .encode(
                alt.X("agree_frac:Q", bin=alt.Bin(maxbins=20), title="Overlap agreement fraction"),
                alt.Y("count():Q", title="Count"),
                tooltip=["count():Q"],
            )
            .properties(title="Overlap Agreement Fraction", width=380, height=260)
        )
    else:
        overlap_chart = alt.Chart(overlap_df).mark_text().encode(
            text=alt.value("No overlapping fragments in selection")
        ).properties(width=380, height=260)

    mo.tabs({
        "VAF distribution":    mo.ui.altair_chart(vaf_chart),
        "Error spectrum":      mo.ui.altair_chart(spectrum_chart),
        "Strand bias":         mo.ui.altair_chart(strand_chart),
        "Overlap agreement":   mo.ui.altair_chart(overlap_chart),
    })
    return


if __name__ == "__main__":
    app.run()
