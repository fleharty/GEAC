"""Tests for NMF-based mutational signature helpers."""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from signature_nmf import compare_signatures_to_cosmic, fit_sbs_nmf


_CONTEXTS = [f"ctx_{i:02d}" for i in range(96)]


def _synthetic_profiles() -> tuple[pd.DataFrame, pd.DataFrame]:
    sig_a = np.zeros(96, dtype=float)
    sig_b = np.zeros(96, dtype=float)
    sig_a[:48] = 1.0
    sig_b[48:] = 1.0
    sig_a /= sig_a.sum()
    sig_b /= sig_b.sum()

    exposures = pd.DataFrame(
        [
            [100.0, 10.0],
            [80.0, 20.0],
            [20.0, 80.0],
            [10.0, 100.0],
        ],
        index=["S1", "S2", "S3", "S4"],
        columns=["SigA", "SigB"],
    )
    counts = exposures.to_numpy() @ np.vstack([sig_a, sig_b]) * 100
    count_df = pd.DataFrame(counts, index=exposures.index, columns=_CONTEXTS)

    cosmic = pd.DataFrame(
        {
            "SBS_A": sig_a,
            "SBS_B": sig_b,
            "SBS_noise": np.full(96, 1 / 96, dtype=float),
        },
        index=_CONTEXTS,
    )
    return count_df, cosmic


def test_fit_sbs_nmf_recovers_two_signature_structure():
    count_df, _ = _synthetic_profiles()

    result = fit_sbs_nmf(count_df, 2, random_state=0, max_iter=2000)

    profiles = result["profiles"]
    exposures = result["exposure_fractions"]

    assert profiles.shape == (2, 96)
    assert exposures.shape == (4, 2)
    assert np.allclose(profiles.sum(axis=1).values, 1.0)
    assert np.allclose(exposures.sum(axis=1).values, 1.0)
    assert result["matrix_cosine"] > 0.999
    assert result["relative_error_pct"] < 0.05


def test_compare_signatures_to_cosmic_finds_expected_best_matches():
    count_df, cosmic = _synthetic_profiles()
    result = fit_sbs_nmf(count_df, 2, random_state=0, max_iter=2000)

    summary, cosine = compare_signatures_to_cosmic(result["profiles"], cosmic)

    assert set(summary["best_cosmic_signature"]) == {"SBS_A", "SBS_B"}
    assert (summary["best_cosine_similarity"] > 0.99).all()
    assert cosine.shape == (2, 3)


def test_fit_sbs_nmf_requires_multiple_samples():
    single = pd.DataFrame([[1.0] * 96], index=["S1"], columns=_CONTEXTS)

    try:
        fit_sbs_nmf(single, 2)
    except ValueError as exc:
        assert "at least two samples" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("Expected ValueError for single-sample NMF input")
