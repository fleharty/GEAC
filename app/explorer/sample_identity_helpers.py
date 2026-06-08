"""Pure-SQL helpers for the Sample Identity / Duplicates tab.

These functions take a live ``duckdb`` connection to a merged cohort DuckDB and
compute pairwise genetic-fingerprint similarity plus per-sample contamination
indicators. They import no Streamlit, so they are unit-testable directly against
a synthetic DuckDB (see ``app/tests/test_sample_identity_helpers.py``).

Method (somalier / NGSCheckMate style): germline het/hom SNVs recorded in
``alt_bases`` form a per-sample genetic fingerprint. Two samples from the same
individual share nearly the same germline-variant set with concordant genotypes;
unrelated samples share only common population SNPs. The primary discriminator is
the Jaccard overlap of each sample's germline-variant set (uses *observed*
non-reference records only, so it is robust to differing coverage and needs no
hom-ref inference); genotype concordance at the intersection confirms identity.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class IdentityParams:
    """Tunable thresholds for marker selection and pair flagging."""

    min_depth: int = 20          # min total_depth for a usable SNV call
    het_lo: float = 0.15         # min VAF to treat an alt record as a germline call
    hom_lo: float = 0.85         # VAF >= this => hom-alt (geno 2), else het (geno 1)
    min_recurrence: int = 2      # marker must appear in >= this many samples
    max_markers: int = 5000      # cap panel size (controls O(S^2 * M) self-join cost)
    use_gnomad: bool = False     # restrict markers to common gnomAD SNPs when available
    gnomad_lo: float = 0.05
    gnomad_hi: float = 0.95
    low_vaf_lo: float = 0.03     # lower bound of the sub-germline "contamination" band
    # Flagging thresholds (applied in classify_flags)
    t_high: float = 0.7          # jaccard >= t_high + high concordance => same individual
    t_low: float = 0.3          # same subject_id but jaccard < t_low => likely swap
    conc_high: float = 0.95
    conc_low: float = 0.90


def _vaf_expr() -> str:
    return "alt_count * 1.0 / NULLIF(total_depth, 0)"


def _marker_ctes(p: IdentityParams) -> str:
    """CTEs producing the germline genotype-at-marker long table ``geno_m``.

    ``geno_m`` has one row per (sample_id, marker_id, geno, vaf) for every
    non-reference germline SNV call that survives marker selection.
    """
    gnomad_clause = (
        f"AND gnomad_af BETWEEN {p.gnomad_lo} AND {p.gnomad_hi}"
        if p.use_gnomad
        else ""
    )
    return f"""
    geno AS (
        SELECT
            sample_id,
            chrom, pos, alt_allele,
            {_vaf_expr()} AS vaf,
            CASE WHEN {_vaf_expr()} >= {p.hom_lo} THEN 2 ELSE 1 END AS geno
        FROM alt_bases
        WHERE variant_type = 'SNV'
          AND total_depth >= {p.min_depth}
          AND {_vaf_expr()} >= {p.het_lo}
          {gnomad_clause}
    ),
    marker_counts AS (
        SELECT chrom, pos, alt_allele, COUNT(DISTINCT sample_id) AS n_samples
        FROM geno
        GROUP BY chrom, pos, alt_allele
    ),
    markers AS (
        SELECT chrom, pos, alt_allele
        FROM marker_counts
        WHERE n_samples >= {p.min_recurrence}
        ORDER BY n_samples DESC, chrom, pos, alt_allele
        LIMIT {p.max_markers}
    ),
    geno_m AS (
        SELECT
            g.sample_id,
            g.chrom || ':' || g.pos || ':' || g.alt_allele AS marker_id,
            g.geno,
            g.vaf
        FROM geno g
        JOIN markers m USING (chrom, pos, alt_allele)
    )
    """


def subject_map(con) -> dict:
    """Return {sample_id -> subject_id (or None)} derived from alt_bases."""
    cols = {r[0] for r in con.execute("DESCRIBE alt_bases").fetchall()}
    if "subject_id" not in cols:
        return {}
    rows = con.execute(
        "SELECT sample_id, ANY_VALUE(subject_id) FROM alt_bases GROUP BY sample_id"
    ).fetchall()
    return {sid: (subj if subj else None) for sid, subj in rows}


def compute_pairwise_identity(con, p: IdentityParams) -> pd.DataFrame:
    """Pairwise germline-fingerprint similarity for all sample pairs sharing >=1 marker.

    Columns: sample_a, sample_b, n_nonref_a, n_nonref_b, n_shared,
    n_shared_concordant, jaccard, concordance.
    """
    sql = f"""
    WITH {_marker_ctes(p)},
    nonref AS (
        SELECT sample_id, COUNT(*) AS n_nonref
        FROM geno_m
        GROUP BY sample_id
    ),
    pairs AS (
        SELECT
            a.sample_id AS sample_a,
            b.sample_id AS sample_b,
            COUNT(*) AS n_shared,
            SUM(CASE WHEN a.geno = b.geno THEN 1 ELSE 0 END) AS n_shared_concordant
        FROM geno_m a
        JOIN geno_m b
          ON a.marker_id = b.marker_id
         AND a.sample_id < b.sample_id
        GROUP BY a.sample_id, b.sample_id
    )
    SELECT
        p.sample_a,
        p.sample_b,
        na.n_nonref AS n_nonref_a,
        nb.n_nonref AS n_nonref_b,
        p.n_shared,
        p.n_shared_concordant,
        p.n_shared * 1.0
            / NULLIF(na.n_nonref + nb.n_nonref - p.n_shared, 0) AS jaccard,
        p.n_shared_concordant * 1.0 / NULLIF(p.n_shared, 0) AS concordance
    FROM pairs p
    JOIN nonref na ON na.sample_id = p.sample_a
    JOIN nonref nb ON nb.sample_id = p.sample_b
    ORDER BY jaccard DESC, concordance DESC
    """
    return con.execute(sql).df()


def compute_all_loci_jaccard(con, samples: list[str], p: IdentityParams) -> pd.DataFrame:
    """All-alt-loci Jaccard (any variant_type, any VAF) for a restricted sample set.

    Used to distinguish technical replicates (which also match at low-VAF/somatic
    loci) from two distinct biological samples of the same individual. Restricted
    to ``samples`` (the candidate set) to keep the second self-join cheap.
    """
    if len(samples) < 2:
        return pd.DataFrame(
            columns=["sample_a", "sample_b", "all_loci_jaccard"]
        )
    in_list = ", ".join("'" + s.replace("'", "''") + "'" for s in samples)
    sql = f"""
    WITH loci AS (
        SELECT DISTINCT
            sample_id,
            chrom || ':' || pos || ':' || alt_allele AS locus_id
        FROM alt_bases
        WHERE total_depth >= {p.min_depth}
          AND sample_id IN ({in_list})
    ),
    counts AS (
        SELECT sample_id, COUNT(*) AS n FROM loci GROUP BY sample_id
    ),
    pairs AS (
        SELECT a.sample_id AS sample_a, b.sample_id AS sample_b,
               COUNT(*) AS n_shared
        FROM loci a
        JOIN loci b ON a.locus_id = b.locus_id AND a.sample_id < b.sample_id
        GROUP BY a.sample_id, b.sample_id
    )
    SELECT
        p.sample_a, p.sample_b,
        p.n_shared * 1.0
            / NULLIF(ca.n + cb.n - p.n_shared, 0) AS all_loci_jaccard
    FROM pairs p
    JOIN counts ca ON ca.sample_id = p.sample_a
    JOIN counts cb ON cb.sample_id = p.sample_b
    """
    return con.execute(sql).df()


def compute_contamination(con, p: IdentityParams) -> pd.DataFrame:
    """Per-sample contamination indicators (heuristic, not a calibrated estimate).

    Columns: sample_id, het_vaf_dispersion, n_het_markers, low_vaf_snv_count.
    het_vaf_dispersion = mean |VAF - 0.5| over germline het markers; drifts up
    when foreign reads skew heterozygous balance. low_vaf_snv_count counts
    sub-germline SNVs (VAF in [low_vaf_lo, het_lo)).
    """
    disp_sql = f"""
    WITH {_marker_ctes(p)}
    SELECT
        sample_id,
        AVG(ABS(vaf - 0.5)) AS het_vaf_dispersion,
        COUNT(*) AS n_het_markers
    FROM geno_m
    WHERE geno = 1
    GROUP BY sample_id
    """
    disp = con.execute(disp_sql).df()

    low_sql = f"""
    SELECT sample_id, COUNT(*) AS low_vaf_snv_count
    FROM alt_bases
    WHERE variant_type = 'SNV'
      AND total_depth >= {p.min_depth}
      AND {_vaf_expr()} >= {p.low_vaf_lo}
      AND {_vaf_expr()} <  {p.het_lo}
    GROUP BY sample_id
    """
    low = con.execute(low_sql).df()

    out = disp.merge(low, on="sample_id", how="outer")
    out["low_vaf_snv_count"] = out["low_vaf_snv_count"].fillna(0).astype(int)
    return out.sort_values("het_vaf_dispersion", ascending=False, na_position="last")


def classify_flags(pairs: pd.DataFrame, subjects: dict, p: IdentityParams) -> pd.DataFrame:
    """Annotate a pairwise dataframe with subject_ids and a flag category.

    Flags (requires a populated subject map; otherwise subject columns are None
    and only UNKNOWN_DUPLICATE/None apply):
      - UNKNOWN_DUPLICATE: high jaccard + high concordance + different/absent subjects
      - POSSIBLE_SWAP:     same subject_id but low jaccard or low concordance
      - EXPECTED_MATCH:    same subject_id and high similarity (sanity check)
    """
    df = pairs.copy()
    df["subject_a"] = df["sample_a"].map(lambda s: subjects.get(s))
    df["subject_b"] = df["sample_b"].map(lambda s: subjects.get(s))

    def _flag(row) -> str:
        sa, sb = row["subject_a"], row["subject_b"]
        jac, conc = row["jaccard"], row["concordance"]
        same_subject = sa is not None and sb is not None and sa == sb
        looks_identical = jac >= p.t_high and conc >= p.conc_high
        if same_subject:
            if jac < p.t_low or conc < p.conc_low:
                return "POSSIBLE_SWAP"
            return "EXPECTED_MATCH"
        # different or missing subject ids
        if looks_identical:
            return "UNKNOWN_DUPLICATE"
        return ""

    df["flag"] = df.apply(_flag, axis=1)
    return df
