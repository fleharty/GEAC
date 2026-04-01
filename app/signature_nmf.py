"""Helpers for de novo mutational-signature discovery with NMF."""

from __future__ import annotations

import numpy as np
import pandas as pd


def _normalize_rows(matrix: np.ndarray) -> np.ndarray:
    row_sums = matrix.sum(axis=1, keepdims=True)
    return np.divide(matrix, row_sums, out=np.zeros_like(matrix), where=row_sums > 0)


def _cosine_similarity_matrix(lhs: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    lhs_norm = np.linalg.norm(lhs, axis=1, keepdims=True)
    rhs_norm = np.linalg.norm(rhs, axis=1, keepdims=True)
    denom = np.clip(lhs_norm, 1e-12, None) * np.clip(rhs_norm.T, 1e-12, None)
    return (lhs @ rhs.T) / denom


def fit_sbs_nmf(
    sample_counts: pd.DataFrame,
    n_components: int,
    *,
    random_state: int = 0,
    max_iter: int = 1000,
) -> dict[str, object]:
    """Factor a sample-by-context count matrix into exposures and signatures.

    ``sample_counts`` must have samples on rows and SBS96 contexts on columns.
    Returns normalized signature profiles and per-sample exposure counts/fractions.
    """
    if sample_counts.empty:
        raise ValueError("No SBS96 counts available for NMF.")
    if sample_counts.shape[0] < 2:
        raise ValueError("NMF signature discovery requires at least two samples.")
    if n_components < 2:
        raise ValueError("NMF signature discovery requires at least two signatures.")
    if n_components > sample_counts.shape[0]:
        raise ValueError("Number of signatures cannot exceed the number of samples.")

    try:
        from sklearn.decomposition import NMF
    except ImportError as exc:  # pragma: no cover - exercised only in missing-dependency envs
        raise RuntimeError(
            "NMF signature discovery requires scikit-learn. "
            "Install it with `pip install scikit-learn`."
        ) from exc

    matrix = sample_counts.to_numpy(dtype=float)
    if np.any(matrix < 0):
        raise ValueError("SBS96 count matrix must be non-negative.")

    model = NMF(
        n_components=n_components,
        init="nndsvda",
        random_state=random_state,
        max_iter=max_iter,
    )
    exposures_raw = model.fit_transform(matrix)
    signatures_raw = model.components_

    signature_scales = signatures_raw.sum(axis=1)
    signatures_norm = _normalize_rows(signatures_raw)
    exposures = exposures_raw * signature_scales[np.newaxis, :]
    exposure_fractions = _normalize_rows(exposures)
    reconstructed = exposures @ signatures_norm

    signature_names = [f"NMF{i + 1}" for i in range(n_components)]
    profile_df = pd.DataFrame(signatures_norm, index=signature_names, columns=sample_counts.columns)
    exposure_df = pd.DataFrame(exposures, index=sample_counts.index, columns=signature_names)
    exposure_frac_df = pd.DataFrame(
        exposure_fractions,
        index=sample_counts.index,
        columns=signature_names,
    )

    flat_obs = matrix.ravel()
    flat_recon = reconstructed.ravel()
    matrix_cosine = float(
        np.dot(flat_obs, flat_recon)
        / (np.linalg.norm(flat_obs) * np.linalg.norm(flat_recon) + 1e-12)
    )
    relative_error_pct = float(
        np.linalg.norm(matrix - reconstructed) / (matrix.sum() + 1e-12) * 100
    )

    return {
        "profiles": profile_df,
        "exposures": exposure_df,
        "exposure_fractions": exposure_frac_df,
        "reconstruction_err": float(model.reconstruction_err_),
        "matrix_cosine": matrix_cosine,
        "relative_error_pct": relative_error_pct,
        "n_iter": int(model.n_iter_),
    }


def compare_signatures_to_cosmic(
    discovered_profiles: pd.DataFrame,
    cosmic_matrix: pd.DataFrame,
    *,
    top_n: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare learned signatures to COSMIC signatures by cosine similarity.

    ``discovered_profiles`` is signature-by-context. ``cosmic_matrix`` is expected to
    be context-by-signature, matching the current Explorer COSMIC loader.
    """
    if discovered_profiles.empty:
        raise ValueError("No discovered signatures available for comparison.")

    cosmic_aligned = cosmic_matrix.reindex(discovered_profiles.columns)
    missing_contexts = cosmic_aligned.isna().any(axis=1)
    if missing_contexts.any():
        missing = missing_contexts.sum()
        raise ValueError(
            f"{missing} context(s) are missing from the COSMIC matrix; cannot compare signatures."
        )

    cosmic_profiles = cosmic_aligned.T
    cosmic_norm = pd.DataFrame(
        _normalize_rows(cosmic_profiles.to_numpy(dtype=float)),
        index=cosmic_profiles.index,
        columns=cosmic_profiles.columns,
    )
    discovered_norm = pd.DataFrame(
        _normalize_rows(discovered_profiles.to_numpy(dtype=float)),
        index=discovered_profiles.index,
        columns=discovered_profiles.columns,
    )

    cosine = pd.DataFrame(
        _cosine_similarity_matrix(
            discovered_norm.to_numpy(dtype=float),
            cosmic_norm.to_numpy(dtype=float),
        ),
        index=discovered_norm.index,
        columns=cosmic_norm.index,
    )

    rows: list[dict[str, object]] = []
    for sig_name in cosine.index:
        ranked = cosine.loc[sig_name].sort_values(ascending=False)
        top_matches = ranked.head(max(1, top_n))
        rows.append(
            {
                "signature": sig_name,
                "best_cosmic_signature": top_matches.index[0],
                "best_cosine_similarity": float(top_matches.iloc[0]),
                "top_matches": ", ".join(
                    f"{name} ({score:.3f})" for name, score in top_matches.items()
                ),
            }
        )

    return pd.DataFrame(rows), cosine
