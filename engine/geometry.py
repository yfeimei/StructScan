"""Module 3: Geometry & Metric Engine.

Computes the multi-dimensional structural fingerprint for one aligned variant:

    Global RMSD          sqrt( (1/N) * sum ||r_i - r_i'||^2 )
    Local displacement   per-residue Euclidean shift of the alpha-carbon
    Dihedral delta       per-residue change in backbone phi / psi

Note on the dihedral term: torsion angles are periodic, so the plain difference
|theta_mut - theta_ref| from the design document is wrong at the wrap-around --
it scores 179 deg vs -179 deg as 358 deg apart when they are 2 deg apart. The
implementation below takes the shortest angular separation, giving a value in
[0, 180].
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from .aligner import AlignmentResult
from .parser import StructureData

# Per-residue counting thresholds for the summary fingerprint.
DISPLACEMENT_COUNT_THRESHOLD = 1.0  # angstroms
DIHEDRAL_COUNT_THRESHOLD = 30.0  # degrees


def rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """Root mean square deviation between two paired coordinate sets."""
    coords_a = np.asarray(coords_a, dtype=float)
    coords_b = np.asarray(coords_b, dtype=float)
    if coords_a.shape != coords_b.shape:
        raise ValueError(f"shape mismatch: {coords_a.shape} vs {coords_b.shape}")
    if coords_a.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(np.sum((coords_a - coords_b) ** 2, axis=1))))


def displacement(coords_a: np.ndarray, coords_b: np.ndarray) -> np.ndarray:
    """Per-row Euclidean distance: sqrt(dx^2 + dy^2 + dz^2)."""
    return np.linalg.norm(np.asarray(coords_a, dtype=float) - np.asarray(coords_b, dtype=float), axis=1)


def angular_delta(theta_a: Optional[float], theta_b: Optional[float]) -> Optional[float]:
    """Shortest separation between two angles in degrees, in [0, 180].

    Returns None if either angle is undefined (chain termini, chain breaks).
    """
    if theta_a is None or theta_b is None:
        return None
    delta = abs(float(theta_a) - float(theta_b)) % 360.0
    return delta if delta <= 180.0 else 360.0 - delta


def _nan_max(series: pd.Series) -> float:
    return float(series.max()) if series.notna().any() else float("nan")


def _nan_mean(series: pd.Series) -> float:
    return float(series.mean()) if series.notna().any() else float("nan")


def compute_metrics(
    reference: StructureData,
    variant: StructureData,
    alignment: AlignmentResult,
    displacement_count_threshold: float = DISPLACEMENT_COUNT_THRESHOLD,
    dihedral_count_threshold: float = DIHEDRAL_COUNT_THRESHOLD,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """Return (per-residue metric table, per-variant summary fingerprint)."""
    displacements = displacement(alignment.ref_coords, alignment.var_coords_aligned)

    rows = []
    for index, key in enumerate(alignment.matched_keys):
        ref_res = reference.residues[key]
        var_res = variant.residues[key]
        delta_phi = angular_delta(ref_res.phi, var_res.phi)
        delta_psi = angular_delta(ref_res.psi, var_res.psi)

        finite = [d for d in (delta_phi, delta_psi) if d is not None]
        max_dihedral = max(finite) if finite else None

        rows.append(
            {
                "variant": variant.name,
                "chain": ref_res.chain_id,
                "res_seq": ref_res.res_seq,
                "icode": ref_res.icode.strip(),
                "position": ref_res.label,
                "ref_resname": ref_res.res_name,
                "var_resname": var_res.res_name,
                "is_substitution": ref_res.res_name != var_res.res_name,
                "ca_displacement": float(displacements[index]),
                "delta_phi": delta_phi,
                "delta_psi": delta_psi,
                "max_dihedral_delta": max_dihedral,
                "used_in_fit": bool(alignment.fitted_mask[index]),
            }
        )

    per_residue = pd.DataFrame(rows)

    summary: Dict[str, object] = {
        "variant": variant.name,
        "source_file": variant.path.name,
        "global_rmsd": alignment.rms,
        "global_rmsd_all_matched": alignment.rms_all_matched,
        "n_residues_reference": alignment.n_ref,
        "n_residues_variant": alignment.n_var,
        "n_residues_matched": alignment.n_matched,
        "n_residues_fitted": alignment.n_fitted,
        "coverage": alignment.coverage,
        "max_ca_displacement": _nan_max(per_residue["ca_displacement"]),
        "mean_ca_displacement": _nan_mean(per_residue["ca_displacement"]),
        "n_residues_displaced": int(
            (per_residue["ca_displacement"] > displacement_count_threshold).sum()
        ),
        "max_dihedral_delta": _nan_max(per_residue["max_dihedral_delta"]),
        "mean_dihedral_delta": _nan_mean(per_residue["max_dihedral_delta"]),
        "n_residues_torsioned": int(
            (per_residue["max_dihedral_delta"] > dihedral_count_threshold).sum()
        ),
        "n_substitutions": len(alignment.substitutions),
        "substitutions": alignment.substitution_label,
    }

    return per_residue, summary
