import math

import pandas as pd

from read_context_helpers import add_read_context_fraction_metrics


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
