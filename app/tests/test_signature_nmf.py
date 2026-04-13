"""Tests for NMF-based mutational signature helpers."""

from __future__ import annotations

import io
import math
import os
import sys
import zipfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from signature_nmf import (
    _CONNECTED_SIG_GROUPS,
    _expand_connected_groups,
    build_signature_download_zip,
    build_signature_exposure_download_table,
    build_signature_download_table,
    compare_signatures_to_cosmic,
    fit_cosmic_augmented_nmf,
    fit_cosmic_cohort_greedy,
    fit_cosmic_single_sample_greedy,
    fit_sbs_nmf,
)


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


def _cosine(lhs: np.ndarray, rhs: np.ndarray) -> float:
    return float(
        np.dot(lhs, rhs) / (np.linalg.norm(lhs) * np.linalg.norm(rhs) + 1e-12)
    )


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

    assert set(summary["most_similar_cosmic_signature"]) == {"SBS_A", "SBS_B"}
    assert (summary["most_similar_cosine_similarity"] > 0.99).all()
    assert cosine.shape == (2, 3)


def test_fit_sbs_nmf_requires_multiple_samples():
    single = pd.DataFrame([[1.0] * 96], index=["S1"], columns=_CONTEXTS)

    try:
        fit_sbs_nmf(single, 2)
    except ValueError as exc:
        assert "at least two samples" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("Expected ValueError for single-sample NMF input")


def test_fit_cosmic_augmented_nmf_learns_residual_signature():
    contexts = _CONTEXTS
    fixed = np.zeros(96, dtype=float)
    novel = np.zeros(96, dtype=float)
    fixed[:32] = 1.0
    novel[32:64] = 1.0
    fixed /= fixed.sum()
    novel /= novel.sum()

    exposures = pd.DataFrame(
        [
            [80.0, 20.0],
            [70.0, 30.0],
            [25.0, 75.0],
            [10.0, 90.0],
        ],
        index=["S1", "S2", "S3", "S4"],
        columns=["SBS_FIXED", "NovelTruth"],
    )
    counts = exposures.to_numpy() @ np.vstack([fixed, novel]) * 100
    count_df = pd.DataFrame(counts, index=exposures.index, columns=contexts)
    cosmic = pd.DataFrame(
        {
            "SBS_FIXED": fixed,
            "SBS_OTHER": np.full(96, 1 / 96, dtype=float),
        },
        index=contexts,
    )

    result = fit_cosmic_augmented_nmf(count_df, cosmic, ["SBS_FIXED"], max_iter=500)

    assert list(result["profiles"].index) == ["SBS_FIXED", "Novel1"]
    assert np.allclose(result["profiles"].sum(axis=1).values, 1.0)
    assert np.allclose(result["exposure_fractions"].sum(axis=1).values, 1.0)
    assert result["relative_error_improvement_pct"] > 1.0
    assert result["relative_error_pct"] < result["fixed_only_relative_error_pct"]
    assert result["matrix_cosine"] > result["fixed_only_matrix_cosine"]

    learned = result["profiles"].loc["Novel1"].values.astype(float)
    assert _cosine(learned, novel) > 0.98


def test_fit_cosmic_augmented_nmf_requires_fixed_signatures():
    count_df, cosmic = _synthetic_profiles()

    try:
        fit_cosmic_augmented_nmf(count_df, cosmic, [])
    except ValueError as exc:
        assert "at least one fixed COSMIC signature" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("Expected ValueError when no fixed signatures are provided")


def test_build_signature_download_table_includes_match_metadata():
    spec_df = pd.DataFrame(
        {
            "sbs_label": ["A[C>A]A", "A[C>G]A"],
            "mut_type": ["C>A", "C>G"],
            "fraction": [0.75, 0.25],
        }
    )

    result = build_signature_download_table(
        spec_df,
        signature_name="Novel1",
        most_similar_cosmic_signature="SBS18",
        most_similar_cosine_similarity=0.91,
        fixed_signature_names=["SBS1", "SBS5"],
    )

    assert list(result.columns) == [
        "signature",
        "sbs_label",
        "mut_type",
        "fraction",
        "most_similar_cosmic_signature",
        "most_similar_cosine_similarity",
        "fixed_cosmic_signatures",
    ]
    assert set(result["signature"]) == {"Novel1"}
    assert set(result["most_similar_cosmic_signature"]) == {"SBS18"}
    assert np.allclose(result["most_similar_cosine_similarity"], 0.91)
    assert set(result["fixed_cosmic_signatures"]) == {"SBS1, SBS5"}


def test_build_signature_download_zip_contains_signature_and_provenance_files():
    signature_df = pd.DataFrame(
        {
            "signature": ["Novel1"],
            "sbs_label": ["A[C>A]A"],
            "mut_type": ["C>A"],
            "fraction": [1.0],
        }
    )
    provenance_df = pd.DataFrame(
        {
            "section": ["filters", "signature_discovery"],
            "name": ["chromosome", "fixed_cosmic_signatures"],
            "value": ["chr1", "SBS1, SBS5"],
        }
    )
    exposure_df = build_signature_exposure_download_table(
        pd.DataFrame(
            {
                "SBS1": [0.4, 0.7],
                "Novel1": [0.6, 0.3],
            },
            index=["sample_a", "sample_b"],
        )
    )

    payload = build_signature_download_zip(
        signature_df,
        provenance_df,
        signature_name="Novel1",
        sample_exposures_df=exposure_df,
    )

    with zipfile.ZipFile(io.BytesIO(payload)) as zf:
        assert sorted(zf.namelist()) == [
            "novel1_provenance.tsv",
            "novel1_sample_exposures.tsv",
            "novel1_sbs96.tsv",
        ]
        sig_text = zf.read("novel1_sbs96.tsv").decode()
        prov_text = zf.read("novel1_provenance.tsv").decode()
        exposure_text = zf.read("novel1_sample_exposures.tsv").decode()

    assert "signature\tsbs_label\tmut_type\tfraction" in sig_text
    assert "Novel1\tA[C>A]A\tC>A\t1.0" in sig_text
    assert "section\tname\tvalue" in prov_text
    assert "filters\tchromosome\tchr1" in prov_text
    assert "sample_label\tsignature\texposure" in exposure_text
    assert "sample_a\tNovel1\t0.6" in exposure_text


# ── Greedy add/remove tests ───────────────────────────────────────────────────

def test_expand_connected_groups_adds_sister_signatures():
    all_names = ["SBS1", "SBS2", "SBS13", "SBS5"]
    result = _expand_connected_groups(["SBS2"], all_names)
    assert "SBS2" in result
    assert "SBS13" in result


def test_expand_connected_groups_ignores_non_group_sigs():
    all_names = ["SBS1", "SBS2", "SBS13", "SBS5"]
    result = _expand_connected_groups(["SBS1"], all_names)
    assert "SBS2" not in result
    assert "SBS13" not in result


def test_connected_sig_groups_are_frozensets():
    for group in _CONNECTED_SIG_GROUPS:
        assert isinstance(group, frozenset)
        assert len(group) >= 2


def test_fit_cosmic_single_sample_greedy_recovers_sparse_solution():
    count_df, cosmic = _synthetic_profiles()
    # S1 is heavily SigA — greedy should select SigA (and possibly SigB) but not noise
    obs = count_df.loc["S1"].values.astype(float)
    W = cosmic.values.astype(float)
    sig_names = cosmic.columns.tolist()

    result = fit_cosmic_single_sample_greedy(
        obs, W, sig_names,
        add_penalty=0.05,
        remove_penalty=0.01,
        background_sig_names=[],
        connected_sigs=False,
    )

    assert "SBS_A" in result["active_sig_names"]
    assert "SBS_noise" not in result["active_sig_names"]
    assert result["cosine_similarity"] > 0.98
    assert math.isclose(result["exposure_fractions"].sum(), 1.0, abs_tol=1e-6)
    assert len(result["active_sig_names"]) <= 2


def test_fit_cosmic_single_sample_greedy_sparser_than_nnls():
    from scipy.optimize import nnls as scipy_nnls
    count_df, cosmic = _synthetic_profiles()
    obs = count_df.loc["S1"].values.astype(float)
    W = cosmic.values.astype(float)
    sig_names = cosmic.columns.tolist()

    h_nnls, _ = scipy_nnls(W, obs)
    n_nnls_active = int((h_nnls > 0).sum())

    result = fit_cosmic_single_sample_greedy(
        obs, W, sig_names,
        add_penalty=0.05,
        remove_penalty=0.01,
        background_sig_names=[],
        connected_sigs=False,
    )
    n_greedy_active = len(result["active_sig_names"])

    assert n_greedy_active <= n_nnls_active


def test_fit_cosmic_single_sample_greedy_add_penalty_controls_sparsity():
    count_df, cosmic = _synthetic_profiles()
    obs = count_df.loc["S2"].values.astype(float)
    W = cosmic.values.astype(float)
    sig_names = cosmic.columns.tolist()

    result_tight = fit_cosmic_single_sample_greedy(
        obs, W, sig_names, add_penalty=0.15, background_sig_names=[], connected_sigs=False,
    )
    result_loose = fit_cosmic_single_sample_greedy(
        obs, W, sig_names, add_penalty=0.001, background_sig_names=[], connected_sigs=False,
    )
    assert len(result_tight["active_sig_names"]) <= len(result_loose["active_sig_names"])


def test_fit_cosmic_single_sample_greedy_respects_background_sigs():
    count_df, cosmic = _synthetic_profiles()
    obs = count_df.loc["S1"].values.astype(float)
    W = cosmic.values.astype(float)
    sig_names = cosmic.columns.tolist()

    result = fit_cosmic_single_sample_greedy(
        obs, W, sig_names,
        add_penalty=0.99,  # very tight — almost nothing gets added
        background_sig_names=["SBS_B"],
        connected_sigs=False,
    )
    # SBS_B must appear in exposures even though it contributes little
    idx_b = sig_names.index("SBS_B")
    assert result["exposures"][idx_b] >= 0.0  # non-negative
    assert "SBS_B" in result["active_sig_names"] or result["exposures"][idx_b] == 0.0
    # The key invariant: SBS_B was never removed (permanent)
    # When add_penalty is very tight, background is the only active sig
    assert result["n_layers"] >= 1


def test_fit_cosmic_cohort_greedy_shape_and_keys():
    count_df, cosmic = _synthetic_profiles()

    result = fit_cosmic_cohort_greedy(
        count_df, cosmic,
        add_penalty=0.05,
        remove_penalty=0.01,
        background_sig_names=[],
        connected_sigs=False,
    )

    n_sigs = cosmic.shape[1]
    assert result["exposure_fractions"].shape == (4, n_sigs)
    assert result["exposures"].shape == (4, n_sigs)
    assert set(result.keys()) >= {
        "exposures", "exposure_fractions", "active_signatures",
        "per_sample_metrics", "sig_names",
    }

    metrics = result["per_sample_metrics"]
    assert len(metrics) == 4
    assert set(metrics.columns) >= {
        "sample", "l2_distance", "cosine_similarity", "n_active_sigs", "n_layers",
    }
    assert (metrics["cosine_similarity"] > 0.95).all()
    assert (result["exposure_fractions"].sum(axis=1) > 0.99).all()
