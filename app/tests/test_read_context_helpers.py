import math

import duckdb
import pandas as pd

from read_context_helpers import add_read_context_fraction_metrics, compute_locus_n_asymmetry


def test_add_read_context_fraction_metrics_zero_safe():
    df = pd.DataFrame(
        [
            {
                "n_before_alt": 4,
                "n_after_alt": 6,
                "n_n_before_alt": 1,
                "n_n_after_alt": 3,
            },
            {
                "n_before_alt": 0,
                "n_after_alt": 0,
                "n_n_before_alt": 0,
                "n_n_after_alt": 0,
            },
        ]
    )

    out = add_read_context_fraction_metrics(df)

    assert out.loc[0, "frac_n_before_alt"] == 0.25
    assert out.loc[0, "frac_n_after_alt"] == 0.5
    assert out.loc[0, "delta_n_fraction"] == 0.25
    assert math.isnan(out.loc[1, "frac_n_before_alt"])
    assert math.isnan(out.loc[1, "frac_n_after_alt"])
    assert math.isnan(out.loc[1, "delta_n_fraction"])


def test_compute_locus_n_asymmetry():
    """Locus aggregation produces correct per-locus N-asymmetry metrics."""
    con = duckdb.connect()
    # Two reads at the same locus: one highly asymmetric, one not.
    con.execute("""
        CREATE TABLE alt_reads AS
        SELECT * FROM (VALUES
            ('S1', 'chr1', 100, 'A', 10, 1, 90, 45),
            ('S1', 'chr1', 100, 'A', 10, 0,  5,  0)
        ) t(sample_id, chrom, pos, alt_allele,
            n_before_alt, n_n_before_alt, n_after_alt, n_n_after_alt)
    """)
    r_join = "(SELECT * FROM alt_reads) ar INNER JOIN (SELECT DISTINCT 'S1' AS sample_id, 'chr1' AS chrom, 100 AS pos, 'A' AS alt_allele) _filt ON ar.sample_id=_filt.sample_id AND ar.chrom=_filt.chrom AND ar.pos=_filt.pos AND ar.alt_allele=_filt.alt_allele"
    result = compute_locus_n_asymmetry(con, r_join)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["n_alt_reads"] == 2
    # mean_frac_n_before: avg(1/10, 0/10) = 0.05
    assert abs(row["mean_frac_n_before"] - 0.05) < 1e-9
    # mean_frac_n_after: avg(45/90, 0/5) = avg(0.5, 0.0) = 0.25
    assert abs(row["mean_frac_n_after"] - 0.25) < 1e-9
    assert abs(row["mean_delta_n_frac"] - 0.20) < 1e-9
    # frac_reads_asymmetric: read 1 has frac_after=0.5 (not >0.5) and frac_before=0.1 (not <0.1) → 0
    assert row["frac_reads_asymmetric"] == 0.0
