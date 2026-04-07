import numpy as np
import pandas as pd


def add_read_context_fraction_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Add zero-safe before/after N-fraction metrics for alt-supporting reads."""
    out = df.copy()
    out["frac_n_before_alt"] = np.where(
        out["n_before_alt"] > 0,
        out["n_n_before_alt"] / out["n_before_alt"],
        np.nan,
    )
    out["frac_n_after_alt"] = np.where(
        out["n_after_alt"] > 0,
        out["n_n_after_alt"] / out["n_after_alt"],
        np.nan,
    )
    out["delta_n_fraction"] = out["frac_n_after_alt"] - out["frac_n_before_alt"]
    return out


def compute_locus_n_asymmetry_precomputed(
    con, table_expr: str, where: str
) -> pd.DataFrame:
    """Read per-locus N-asymmetry scores from precomputed columns on alt_bases.

    Available when `geac collect --reads-output` was run with v0.4.7+ and the
    new optional columns (`n_alt_reads_with_n_ctx`, `mean_frac_n_before`,
    `mean_frac_n_after`, `mean_delta_n_frac`, `frac_reads_asymmetric`) are
    present on `alt_bases`.

    Returns one row per (sample_id, chrom, pos, alt_allele, pipeline) so the
    caller can merge against per-pipeline locus rows in multi-pipeline cohorts.
    """
    return con.execute(f"""
        SELECT sample_id,
               chrom,
               pos,
               alt_allele,
               pipeline,
               n_alt_reads_with_n_ctx AS n_alt_reads,
               mean_frac_n_before,
               mean_frac_n_after,
               mean_delta_n_frac,
               frac_reads_asymmetric
        FROM {table_expr}
        WHERE {where}
          AND n_alt_reads_with_n_ctx IS NOT NULL
    """).df()


def compute_locus_n_asymmetry(con, r_join: str) -> pd.DataFrame:
    """Aggregate per-read N-context metrics to per-locus asymmetry scores.

    Parameters
    ----------
    con:
        DuckDB connection.
    r_join:
        The ``_r_join`` FROM fragment already filtered to the current locus
        selection (alt_reads joined to the filtered locus set).

    Returns
    -------
    DataFrame with one row per (sample_id, chrom, pos, alt_allele) containing:

    - ``mean_frac_n_before``  — mean per-read fraction of upstream bases that are N
    - ``mean_frac_n_after``   — mean per-read fraction of downstream bases that are N
    - ``mean_delta_n_frac``   — ``mean_frac_n_after - mean_frac_n_before``
    - ``frac_reads_asymmetric`` — fraction of reads where frac_n_after > 0.5
                                  and frac_n_before < 0.1
    - ``n_alt_reads``         — number of alt reads at the locus
    """
    return con.execute(f"""
        SELECT
            ar.sample_id,
            ar.chrom,
            ar.pos,
            ar.alt_allele,
            COUNT(*)                                                          AS n_alt_reads,
            AVG(CASE WHEN ar.n_before_alt > 0
                     THEN ar.n_n_before_alt * 1.0 / ar.n_before_alt
                     ELSE NULL END)                                           AS mean_frac_n_before,
            AVG(CASE WHEN ar.n_after_alt > 0
                     THEN ar.n_n_after_alt * 1.0 / ar.n_after_alt
                     ELSE NULL END)                                           AS mean_frac_n_after,
            AVG(CASE WHEN ar.n_after_alt > 0
                     THEN ar.n_n_after_alt * 1.0 / ar.n_after_alt
                     ELSE NULL END)
            - AVG(CASE WHEN ar.n_before_alt > 0
                       THEN ar.n_n_before_alt * 1.0 / ar.n_before_alt
                       ELSE NULL END)                                         AS mean_delta_n_frac,
            AVG(CASE
                WHEN ar.n_after_alt  > 0
                 AND ar.n_n_after_alt  * 1.0 / ar.n_after_alt  > 0.5
                 AND (ar.n_before_alt = 0
                      OR ar.n_n_before_alt * 1.0 / ar.n_before_alt < 0.1)
                THEN 1.0 ELSE 0.0 END)                                       AS frac_reads_asymmetric
        FROM {r_join}
        GROUP BY ar.sample_id, ar.chrom, ar.pos, ar.alt_allele
    """).df()
