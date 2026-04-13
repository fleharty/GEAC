from __future__ import annotations

import pandas as pd


def _bait_bias_common_ctes(
    table_expr: str,
    where: str,
    ab_on_tgt: str = "",
    rb_on_tgt: str = "",
    gene_expr: str = "gene",
) -> str:
    return f"""
        WITH hom_alt AS (
            SELECT
                sample_id,
                chrom,
                pos,
                alt_allele,
                variant_type,
                {gene_expr} AS gene,
                total_depth AS observed_depth
            FROM {table_expr}
            WHERE {where}
              AND total_depth > 0
              AND alt_count * 1.0 / total_depth > 0.9
              {ab_on_tgt}
        ),
        sample_baseline AS (
            SELECT
                sample_id,
                MEDIAN(median_target_depth_all) AS sample_baseline_depth
            FROM sample_metrics
            WHERE median_target_depth_all IS NOT NULL
            GROUP BY sample_id
        ),
        noncarrier_rel AS (
            SELECT
                rb.sample_id,
                rb.chrom,
                rb.pos,
                rb.total_depth,
                sb.sample_baseline_depth,
                rb.total_depth * 1.0 / sb.sample_baseline_depth AS relative_depth
            FROM ref_bases rb
            INNER JOIN sample_baseline sb
                    ON rb.sample_id = sb.sample_id
            INNER JOIN (
                SELECT DISTINCT chrom, pos
                FROM hom_alt
            ) hl ON rb.chrom = hl.chrom AND rb.pos = hl.pos
            WHERE rb.total_depth > 0
              AND sb.sample_baseline_depth > 0
              {rb_on_tgt}
        ),
        locus_rel AS (
            SELECT
                chrom,
                pos,
                MEDIAN(relative_depth)                                           AS locus_relative_depth,
                PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY relative_depth)     AS q1_relative_depth,
                PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY relative_depth)     AS q3_relative_depth,
                MIN(relative_depth)                                              AS min_relative_depth,
                MAX(relative_depth)                                              AS max_relative_depth,
                COUNT(*)                                                         AS n_noncarrier_samples
            FROM noncarrier_rel
            GROUP BY chrom, pos
        ),
        carrier_scores AS (
            SELECT
                ha.sample_id,
                ha.chrom,
                ha.pos,
                ha.alt_allele,
                ha.variant_type,
                ha.gene,
                ha.observed_depth,
                sb.sample_baseline_depth,
                lr.locus_relative_depth,
                lr.q1_relative_depth,
                lr.q3_relative_depth,
                lr.min_relative_depth,
                lr.max_relative_depth,
                lr.n_noncarrier_samples,
                sb.sample_baseline_depth * lr.locus_relative_depth               AS expected_depth,
                CASE
                    WHEN ha.observed_depth > 0
                     AND sb.sample_baseline_depth > 0
                     AND lr.locus_relative_depth > 0
                    THEN ha.observed_depth * 1.0
                         / (sb.sample_baseline_depth * lr.locus_relative_depth)
                    ELSE NULL
                END                                                              AS depth_retained
            FROM hom_alt ha
            LEFT JOIN sample_baseline sb
                   ON ha.sample_id = sb.sample_id
            LEFT JOIN locus_rel lr
                   ON ha.chrom = lr.chrom AND ha.pos = lr.pos
        )
    """


def compute_bait_bias_candidates(
    con,
    table_expr: str,
    where: str,
    *,
    ab_on_tgt: str = "",
    rb_on_tgt: str = "",
    gene_expr: str = "gene",
) -> pd.DataFrame:
    """Return per-locus bait-bias candidate rows for the current filtered cohort."""
    return con.execute(
        _bait_bias_common_ctes(
            table_expr,
            where,
            ab_on_tgt=ab_on_tgt,
            rb_on_tgt=rb_on_tgt,
            gene_expr=gene_expr,
        )
        + """
        SELECT
            chrom,
            pos,
            alt_allele,
            variant_type,
            gene,
            COUNT(DISTINCT sample_id)                                            AS n_carrier_samples,
            MAX(n_noncarrier_samples)                                            AS n_noncarrier_samples,
            MEDIAN(observed_depth)                                               AS median_observed_depth,
            MEDIAN(expected_depth)                                               AS median_expected_depth,
            MEDIAN(depth_retained)                                               AS median_depth_retained,
            MEDIAN(sample_baseline_depth)                                        AS median_sample_baseline_depth,
            MAX(locus_relative_depth)                                            AS locus_relative_depth,
            MAX(q1_relative_depth)                                               AS q1_relative_depth,
            MAX(q3_relative_depth)                                               AS q3_relative_depth,
            MAX(min_relative_depth)                                              AS min_relative_depth,
            MAX(max_relative_depth)                                              AS max_relative_depth
        FROM carrier_scores
        GROUP BY chrom, pos, alt_allele, variant_type, gene
        ORDER BY median_depth_retained ASC NULLS LAST, n_noncarrier_samples DESC, chrom, pos
        """
    ).df()


def compute_bait_bias_locus_detail(
    con,
    table_expr: str,
    where: str,
    *,
    chrom: str,
    pos: int,
    alt_allele: str,
    ab_on_tgt: str = "",
    rb_on_tgt: str = "",
    gene_expr: str = "gene",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-carrier details and non-carrier relative-depth rows for one locus."""
    safe_chrom = chrom.replace("'", "''")
    safe_alt = alt_allele.replace("'", "''")
    common = _bait_bias_common_ctes(
        table_expr,
        where,
        ab_on_tgt=ab_on_tgt,
        rb_on_tgt=rb_on_tgt,
        gene_expr=gene_expr,
    )

    carrier_df = con.execute(
        common
        + f"""
        SELECT
            sample_id,
            chrom,
            pos,
            alt_allele,
            variant_type,
            gene,
            observed_depth,
            sample_baseline_depth,
            locus_relative_depth,
            expected_depth,
            depth_retained,
            n_noncarrier_samples,
            q1_relative_depth,
            q3_relative_depth,
            min_relative_depth,
            max_relative_depth
        FROM carrier_scores
        WHERE chrom = '{safe_chrom}'
          AND pos = {int(pos)}
          AND alt_allele = '{safe_alt}'
        ORDER BY depth_retained ASC NULLS LAST, sample_id
        """
    ).df()

    noncarrier_df = con.execute(
        common
        + f"""
        SELECT
            sample_id,
            chrom,
            pos,
            total_depth AS noncarrier_depth,
            sample_baseline_depth,
            relative_depth
        FROM noncarrier_rel
        WHERE chrom = '{safe_chrom}'
          AND pos = {int(pos)}
        ORDER BY relative_depth, sample_id
        """
    ).df()

    return carrier_df, noncarrier_df
