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


def _normalized_l2(observed: np.ndarray, reconstructed: np.ndarray) -> float:
    """Normalized L2 distance: ‖obs − recon‖₂ / ‖obs‖₂."""
    denom = np.linalg.norm(observed)
    if denom < 1e-12:
        return 0.0
    return float(np.linalg.norm(observed - reconstructed) / denom)


# Biologically linked signature groups: if any member is selected, all members
# in the same group should be fitted together (mirrors SigProfilerAssignment).
_CONNECTED_SIG_GROUPS: list[frozenset[str]] = [
    frozenset(["SBS2", "SBS13"]),
    frozenset(["SBS7a", "SBS7b", "SBS7c", "SBS7d"]),
    frozenset(["SBS10a", "SBS10b"]),
    frozenset(["SBS17a", "SBS17b"]),
]


def _expand_connected_groups(active_names: list[str], all_names: list[str]) -> list[str]:
    """Add sister signatures for any active member of a connected group.

    Expansion is applied once per call; callers should re-call each outer layer
    so newly activated signatures can unlock further group members.
    """
    all_set = set(all_names)
    expanded = set(active_names)
    for group in _CONNECTED_SIG_GROUPS:
        if group & expanded:
            expanded |= group & all_set
    return [n for n in all_names if n in expanded]


def _greedy_add_phase(
    profiles: np.ndarray,
    sample: np.ndarray,
    present_indices: list[int],
    candidate_indices: list[int],
    add_penalty: float,
) -> tuple[list[int], np.ndarray, float]:
    """Greedy forward selection: add one signature at a time while L2 improves.

    Returns *(new_active_indices, full_exposure_vector, l2_distance)* where the
    exposure vector has length ``profiles.shape[1]`` with zeros at inactive positions.
    """
    n_sigs = profiles.shape[1]
    present = list(present_indices)
    candidates = list(candidate_indices)

    if present:
        h0, _ = nnls(profiles[:, present], sample)
        baseline_l2 = _normalized_l2(sample, profiles[:, present] @ h0)
    else:
        baseline_l2 = _normalized_l2(sample, np.zeros_like(sample))

    improved = True
    while improved and candidates:
        improved = False
        best_l2 = baseline_l2
        best_idx = -1

        for i in candidates:
            trial = present + [i]
            h_trial, _ = nnls(profiles[:, trial], sample)
            l2_trial = _normalized_l2(sample, profiles[:, trial] @ h_trial)
            if l2_trial < best_l2:
                best_l2 = l2_trial
                best_idx = i

        if best_idx >= 0 and (baseline_l2 - best_l2) > add_penalty:
            present.append(best_idx)
            candidates.remove(best_idx)
            baseline_l2 = best_l2
            improved = True

    exposures = np.zeros(n_sigs, dtype=float)
    if present:
        h_final, _ = nnls(profiles[:, present], sample)
        for idx, val in zip(present, h_final):
            exposures[idx] = val

    return present, exposures, _normalized_l2(sample, profiles @ exposures)


def _greedy_remove_phase(
    profiles: np.ndarray,
    sample: np.ndarray,
    present_indices: list[int],
    permanent_indices: list[int],
    remove_penalty: float,
) -> tuple[list[int], np.ndarray, float]:
    """Greedy backward elimination: remove one signature at a time while L2 barely degrades.

    Returns *(surviving_indices, full_exposure_vector, l2_distance)*.
    """
    n_sigs = profiles.shape[1]
    present = list(present_indices)
    permanent = set(permanent_indices)

    if present:
        h0, _ = nnls(profiles[:, present], sample)
        baseline_l2 = _normalized_l2(sample, profiles[:, present] @ h0)
    else:
        return [], np.zeros(n_sigs, dtype=float), _normalized_l2(sample, np.zeros_like(sample))

    improved = True
    while improved:
        improved = False
        removable = [i for i in present if i not in permanent]
        if not removable:
            break

        best_degradation = float("inf")
        best_remove_idx = -1

        for i in removable:
            trial = [j for j in present if j != i]
            if trial:
                h_trial, _ = nnls(profiles[:, trial], sample)
                l2_trial = _normalized_l2(sample, profiles[:, trial] @ h_trial)
            else:
                l2_trial = _normalized_l2(sample, np.zeros_like(sample))
            degradation = l2_trial - baseline_l2
            if degradation < best_degradation:
                best_degradation = degradation
                best_remove_idx = i

        if best_remove_idx >= 0 and best_degradation < remove_penalty:
            present.remove(best_remove_idx)
            if present:
                h_refit, _ = nnls(profiles[:, present], sample)
                baseline_l2 = _normalized_l2(sample, profiles[:, present] @ h_refit)
            else:
                baseline_l2 = _normalized_l2(sample, np.zeros_like(sample))
            improved = True

    exposures = np.zeros(n_sigs, dtype=float)
    if present:
        h_final, _ = nnls(profiles[:, present], sample)
        for idx, val in zip(present, h_final):
            exposures[idx] = val

    return present, exposures, _normalized_l2(sample, profiles @ exposures)


def fit_cosmic_single_sample_greedy(
    sample: np.ndarray,
    cosmic_profiles: np.ndarray,
    sig_names: list[str],
    *,
    add_penalty: float = 0.05,
    remove_penalty: float = 0.01,
    background_sig_names: list[str] | None = None,
    connected_sigs: bool = True,
) -> dict[str, object]:
    """Single-sample greedy add/remove COSMIC signature fitting (SPA-style).

    Mirrors the ``add_remove_signatures`` algorithm in SigProfilerAssignment:
    an outer loop repeatedly proposes each candidate signature, applies a
    greedy add phase followed by a greedy remove phase, and accepts the result
    that most reduces normalized L2 distance.  Iteration stops when no candidate
    improves the current fit.

    Parameters
    ----------
    sample:
        Observed mutation counts, shape ``(n_contexts,)``.
    cosmic_profiles:
        Signature profiles matrix, shape ``(n_contexts, n_sigs)``, columns normalized.
    sig_names:
        Signature names corresponding to columns of *cosmic_profiles*.
    add_penalty:
        Minimum fractional L2 improvement required to accept adding a signature
        (SPA default 0.05).
    remove_penalty:
        Maximum fractional L2 degradation allowed when removing a signature
        (SPA default 0.01).
    background_sig_names:
        Signatures always included and never removed.  ``None`` auto-detects
        SBS1 and SBS5 when present (universal clock-like signatures).
    connected_sigs:
        When ``True``, activating any member of a biologically linked group
        (e.g. SBS2 / SBS13) forces all group members into the active set.
    """
    n_sigs = cosmic_profiles.shape[1]

    # Permanent (background) signatures — always included, never removed
    if background_sig_names is not None:
        perm_names = [n for n in background_sig_names if n in sig_names]
    else:
        perm_names = [n for n in ("SBS1", "SBS5") if n in sig_names]
    perm_indices = [sig_names.index(n) for n in perm_names]
    active_indices = list(perm_indices)

    original_distance = float("inf")
    final_activities: np.ndarray | None = None
    n_layers = 0

    while True:
        n_layers += 1

        if connected_sigs and active_indices:
            expanded = _expand_connected_groups(
                [sig_names[i] for i in active_indices], sig_names
            )
            active_indices = [sig_names.index(n) for n in expanded]

        candidate_indices = [i for i in range(n_sigs) if i not in active_indices]

        layer_best_distance = float("inf")
        layer_best_activities: np.ndarray | None = None

        for i in candidate_indices:
            add_idx, add_exp, add_dist = _greedy_add_phase(
                cosmic_profiles, sample, list(active_indices), [i], add_penalty
            )
            rem_idx, rem_exp, rem_dist = _greedy_remove_phase(
                cosmic_profiles, sample, add_idx, perm_indices, remove_penalty
            )

            # When the two phases produce the same active set prefer the add result;
            # otherwise the remove phase wins (sparser solution same quality).
            if set(np.nonzero(add_exp)[0]) == set(np.nonzero(rem_exp)[0]):
                dist, exp = add_dist, add_exp
            else:
                dist, exp = rem_dist, rem_exp

            if dist < layer_best_distance:
                layer_best_distance = dist
                layer_best_activities = exp

        if layer_best_activities is None or layer_best_distance >= original_distance:
            break

        original_distance = layer_best_distance
        active_indices = list(np.nonzero(layer_best_activities)[0])
        final_activities = layer_best_activities

    # Fallback: if outer loop never improved, use background-only fit
    if final_activities is None:
        final_activities = np.zeros(n_sigs, dtype=float)
        if perm_indices:
            h, _ = nnls(cosmic_profiles[:, perm_indices], sample)
            for idx, val in zip(perm_indices, h):
                final_activities[idx] = val

    total = final_activities.sum()
    exposure_fractions = final_activities / total if total > 0 else final_activities.copy()

    recon = cosmic_profiles @ final_activities
    cos_sim = float(
        np.dot(sample, recon) / (np.linalg.norm(sample) * np.linalg.norm(recon) + 1e-12)
    )

    return {
        "active_sig_names": [sig_names[i] for i in np.nonzero(final_activities)[0]],
        "exposures": final_activities,
        "exposure_fractions": exposure_fractions,
        "l2_distance": _normalized_l2(sample, recon),
        "cosine_similarity": cos_sim,
        "n_layers": n_layers,
    }


def fit_cosmic_cohort_greedy(
    sample_counts: pd.DataFrame,
    cosmic_matrix: pd.DataFrame,
    *,
    add_penalty: float = 0.05,
    remove_penalty: float = 0.01,
    background_sig_names: list[str] | None = None,
    connected_sigs: bool = True,
) -> dict[str, object]:
    """Apply greedy add/remove COSMIC fitting independently to each sample.

    ``sample_counts`` is samples-by-contexts; ``cosmic_matrix`` is contexts-by-signatures.
    Returns per-sample exposure fractions, raw exposures, per-sample metrics, and
    a map of which signatures were selected for each sample.
    """
    cosmic_aligned = cosmic_matrix.reindex(sample_counts.columns)
    W = cosmic_aligned.values.astype(float)
    sig_names = cosmic_aligned.columns.tolist()
    n_sigs = W.shape[1]

    exposure_rows: list[np.ndarray] = []
    metric_rows: list[dict] = []
    active_sigs_map: dict[str, list[str]] = {}

    for sid, row in sample_counts.iterrows():
        obs = row.values.astype(float)
        if obs.sum() == 0:
            exposure_rows.append(np.zeros(n_sigs))
            metric_rows.append({
                "sample": sid, "l2_distance": 0.0,
                "cosine_similarity": 0.0, "n_active_sigs": 0, "n_layers": 0,
            })
            active_sigs_map[str(sid)] = []
            continue

        result = fit_cosmic_single_sample_greedy(
            obs, W, sig_names,
            add_penalty=add_penalty,
            remove_penalty=remove_penalty,
            background_sig_names=background_sig_names,
            connected_sigs=connected_sigs,
        )
        exposure_rows.append(result["exposure_fractions"])
        metric_rows.append({
            "sample": sid,
            "l2_distance": result["l2_distance"],
            "cosine_similarity": result["cosine_similarity"],
            "n_active_sigs": len(result["active_sig_names"]),
            "n_layers": result["n_layers"],
        })
        active_sigs_map[str(sid)] = result["active_sig_names"]

    exposure_frac_df = pd.DataFrame(
        exposure_rows, index=sample_counts.index, columns=sig_names
    )
    exposure_df = exposure_frac_df.multiply(sample_counts.sum(axis=1), axis=0)

    return {
        "exposures": exposure_df,
        "exposure_fractions": exposure_frac_df,
        "active_signatures": active_sigs_map,
        "per_sample_metrics": pd.DataFrame(metric_rows),
        "sig_names": sig_names,
    }


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
