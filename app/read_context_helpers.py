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
