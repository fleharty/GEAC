"""Helpers for de novo mutational-signature discovery with NMF."""

from __future__ import annotations

import io
import re
import zipfile

import numpy as np
import pandas as pd
from scipy.optimize import nnls


def _normalize_rows(matrix: np.ndarray) -> np.ndarray:
    row_sums = matrix.sum(axis=1, keepdims=True)
    return np.divide(matrix, row_sums, out=np.zeros_like(matrix), where=row_sums > 0)


def _cosine_similarity_matrix(lhs: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    lhs_norm = np.linalg.norm(lhs, axis=1, keepdims=True)
    rhs_norm = np.linalg.norm(rhs, axis=1, keepdims=True)
    denom = np.clip(lhs_norm, 1e-12, None) * np.clip(rhs_norm.T, 1e-12, None)
    return (lhs @ rhs.T) / denom


def _project_to_simplex(values: np.ndarray) -> np.ndarray:
    """Project a vector onto the probability simplex (non-negative, sums to 1)."""
    if values.ndim != 1:
        raise ValueError("Simplex projection expects a 1D vector.")
    if values.size == 0:
        raise ValueError("Cannot project an empty vector onto the simplex.")

    u = np.sort(values)[::-1]
    cssv = np.cumsum(u)
    rho_candidates = u * np.arange(1, len(u) + 1) > (cssv - 1)
    if not np.any(rho_candidates):
        return np.full_like(values, 1.0 / len(values), dtype=float)
    rho = np.nonzero(rho_candidates)[0][-1]
    theta = (cssv[rho] - 1) / (rho + 1)
    return np.maximum(values - theta, 0.0)


def _fit_exposures_nnls(matrix: np.ndarray, profiles: np.ndarray) -> np.ndarray:
    """Fit non-negative per-sample exposures against fixed signature profiles."""
    exposures = np.zeros((matrix.shape[0], profiles.shape[0]), dtype=float)
    design = profiles.T
    for idx, sample in enumerate(matrix):
        coef, _ = nnls(design, sample)
        exposures[idx] = coef
    return exposures


def _fit_metrics(matrix: np.ndarray, reconstruction: np.ndarray) -> tuple[float, float]:
    flat_obs = matrix.ravel()
    flat_recon = reconstruction.ravel()
    matrix_cosine = float(
        np.dot(flat_obs, flat_recon)
        / (np.linalg.norm(flat_obs) * np.linalg.norm(flat_recon) + 1e-12)
    )
    relative_error_pct = float(
        np.linalg.norm(matrix - reconstruction) / (matrix.sum() + 1e-12) * 100
    )
    return matrix_cosine, relative_error_pct


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

    matrix_cosine, relative_error_pct = _fit_metrics(matrix, reconstructed)

    return {
        "profiles": profile_df,
        "exposures": exposure_df,
        "exposure_fractions": exposure_frac_df,
        "reconstruction_err": float(model.reconstruction_err_),
        "matrix_cosine": matrix_cosine,
        "relative_error_pct": relative_error_pct,
        "n_iter": int(model.n_iter_),
    }


def fit_cosmic_augmented_nmf(
    sample_counts: pd.DataFrame,
    cosmic_matrix: pd.DataFrame,
    fixed_signature_names: list[str],
    *,
    max_iter: int = 200,
    tol: float = 1e-7,
) -> dict[str, object]:
    """Fit fixed COSMIC signatures plus one learned residual signature.

    ``sample_counts`` must be sample-by-context counts. ``cosmic_matrix`` must be
    context-by-signature. The returned learned signature is constrained to be
    non-negative and sum to one.
    """
    if sample_counts.empty:
        raise ValueError("No SBS96 counts available for COSMIC-guided NMF.")
    if sample_counts.shape[0] < 2:
        raise ValueError("COSMIC-guided NMF requires at least two samples.")
    if not fixed_signature_names:
        raise ValueError("Select at least one fixed COSMIC signature.")

    missing_signatures = [name for name in fixed_signature_names if name not in cosmic_matrix.columns]
    if missing_signatures:
        raise ValueError(
            "Fixed COSMIC signatures not found in matrix: " + ", ".join(missing_signatures)
        )

    cosmic_aligned = cosmic_matrix.reindex(sample_counts.columns)
    missing_contexts = cosmic_aligned.isna().any(axis=1)
    if missing_contexts.any():
        missing = int(missing_contexts.sum())
        raise ValueError(
            f"{missing} context(s) are missing from the COSMIC matrix; cannot fit COSMIC-guided NMF."
        )

    matrix = sample_counts.to_numpy(dtype=float)
    if np.any(matrix < 0):
        raise ValueError("SBS96 count matrix must be non-negative.")

    fixed_profiles = cosmic_aligned[fixed_signature_names].T.to_numpy(dtype=float)
    fixed_profiles = _normalize_rows(fixed_profiles)

    fixed_exposures = _fit_exposures_nnls(matrix, fixed_profiles)
    fixed_reconstruction = fixed_exposures @ fixed_profiles
    fixed_matrix_cosine, fixed_relative_error_pct = _fit_metrics(matrix, fixed_reconstruction)

    residual = np.maximum(matrix - fixed_reconstruction, 0.0)
    init_profile = residual.sum(axis=0)
    if float(init_profile.sum()) <= 0:
        init_profile = matrix.sum(axis=0)
    if float(init_profile.sum()) <= 0:
        init_profile = np.ones(matrix.shape[1], dtype=float)
    novel_profile = _project_to_simplex(init_profile.astype(float))

    previous_obj = None
    exposures = None
    for iteration in range(1, max_iter + 1):
        profiles = np.vstack([fixed_profiles, novel_profile])
        exposures = _fit_exposures_nnls(matrix, profiles)

        fixed_part = exposures[:, :-1] @ fixed_profiles
        novel_weights = exposures[:, -1]
        weighted_sum_sq = float(np.dot(novel_weights, novel_weights))
        if weighted_sum_sq <= 1e-12:
            residual = np.maximum(matrix - fixed_part, 0.0)
            candidate = residual.sum(axis=0)
        else:
            residual = matrix - fixed_part
            candidate = (novel_weights[:, None] * residual).sum(axis=0) / weighted_sum_sq

        novel_profile = _project_to_simplex(candidate.astype(float))
        reconstruction = exposures @ np.vstack([fixed_profiles, novel_profile])
        current_obj = float(np.linalg.norm(matrix - reconstruction))
        if previous_obj is not None and abs(previous_obj - current_obj) <= tol * max(1.0, previous_obj):
            break
        previous_obj = current_obj

    if exposures is None:
        raise RuntimeError("COSMIC-guided NMF did not run.")

    profiles = np.vstack([fixed_profiles, novel_profile])
    exposures = _fit_exposures_nnls(matrix, profiles)
    reconstruction = exposures @ profiles
    matrix_cosine, relative_error_pct = _fit_metrics(matrix, reconstruction)
    exposure_fractions = _normalize_rows(exposures)

    profile_names = fixed_signature_names + ["Novel1"]
    profile_df = pd.DataFrame(profiles, index=profile_names, columns=sample_counts.columns)
    exposure_df = pd.DataFrame(exposures, index=sample_counts.index, columns=profile_names)
    exposure_frac_df = pd.DataFrame(
        exposure_fractions,
        index=sample_counts.index,
        columns=profile_names,
    )

    return {
        "profiles": profile_df,
        "exposures": exposure_df,
        "exposure_fractions": exposure_frac_df,
        "matrix_cosine": matrix_cosine,
        "relative_error_pct": relative_error_pct,
        "fixed_only_matrix_cosine": fixed_matrix_cosine,
        "fixed_only_relative_error_pct": fixed_relative_error_pct,
        "relative_error_improvement_pct": fixed_relative_error_pct - relative_error_pct,
        "n_iter": iteration,
        "fixed_signature_names": tuple(fixed_signature_names),
        "learned_signature_names": ("Novel1",),
    }


def build_signature_download_table(
    spec_df: pd.DataFrame,
    *,
    signature_name: str,
    most_similar_cosmic_signature: str | None = None,
    most_similar_cosine_similarity: float | None = None,
    fixed_signature_names: list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Build a tidy SBS96 download table for a discovered signature."""
    required = {"sbs_label", "mut_type", "fraction"}
    missing = required.difference(spec_df.columns)
    if missing:
        raise ValueError(
            "Signature download table requires columns: " + ", ".join(sorted(required))
        )

    download_df = spec_df.loc[:, ["sbs_label", "mut_type", "fraction"]].copy()
    download_df.insert(0, "signature", signature_name)
    if most_similar_cosmic_signature is not None:
        download_df["most_similar_cosmic_signature"] = most_similar_cosmic_signature
    if most_similar_cosine_similarity is not None:
        download_df["most_similar_cosine_similarity"] = float(most_similar_cosine_similarity)
    if fixed_signature_names:
        download_df["fixed_cosmic_signatures"] = ", ".join(fixed_signature_names)
    return download_df


def build_signature_exposure_download_table(exposure_fractions: pd.DataFrame) -> pd.DataFrame:
    """Build a tidy per-sample signature exposure table for downloads."""
    if exposure_fractions.empty:
        raise ValueError("Signature exposure download table requires non-empty exposures.")

    return (
        exposure_fractions.reset_index(names="sample_label")
        .melt(id_vars="sample_label", var_name="signature", value_name="exposure")
        .sort_values(["sample_label", "signature"])
        .reset_index(drop=True)
    )


def build_signature_download_zip(
    signature_df: pd.DataFrame,
    provenance_df: pd.DataFrame,
    *,
    signature_name: str,
    sample_exposures_df: pd.DataFrame | None = None,
) -> bytes:
    """Package a discovered signature and its provenance into a zip file."""
    if signature_df.empty:
        raise ValueError("Signature download bundle requires a non-empty signature table.")
    if provenance_df.empty:
        raise ValueError("Signature download bundle requires a non-empty provenance table.")

    slug = re.sub(r"[^A-Za-z0-9]+", "_", signature_name.strip()).strip("_").lower() or "signature"
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(
            f"{slug}_sbs96.tsv",
            signature_df.to_csv(sep="\t", index=False),
        )
        zf.writestr(
            f"{slug}_provenance.tsv",
            provenance_df.to_csv(sep="\t", index=False),
        )
        if sample_exposures_df is not None:
            if sample_exposures_df.empty:
                raise ValueError("Sample exposure table for bundle cannot be empty.")
            zf.writestr(
                f"{slug}_sample_exposures.tsv",
                sample_exposures_df.to_csv(sep="\t", index=False),
            )
    return buffer.getvalue()


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
                "most_similar_cosmic_signature": top_matches.index[0],
                "most_similar_cosine_similarity": float(top_matches.iloc[0]),
                "top_matches": ", ".join(
                    f"{name} ({score:.3f})" for name, score in top_matches.items()
                ),
            }
        )

    return pd.DataFrame(rows), cosine
