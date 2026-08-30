"""Module 2: Alignment Engine.

Superimposes each variant onto the wild-type reference using Biopython's
Superimposer over matched alpha-carbons, and reports the transform plus the
matched-residue bookkeeping the geometry engine needs.

Residues are matched by position key, not by list index, and only the
intersection of the two structures is used. Point mutants change residue
identity at the mutated site, so matching deliberately ignores residue names --
a substitution is recorded, not treated as a mismatch.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from Bio.PDB import Superimposer

from .parser import ResidueKey, StructureData

MIN_ATOMS_FOR_FIT = 3


class AlignmentError(RuntimeError):
    """Raised when two structures share too little in common to superimpose."""


@dataclass
class AlignmentResult:
    """Outcome of superimposing one variant onto the reference."""

    matched_keys: List[ResidueKey]
    ref_coords: np.ndarray  # (N, 3) reference CA coordinates
    var_coords: np.ndarray  # (N, 3) variant CA coordinates, as read from file
    var_coords_aligned: np.ndarray  # (N, 3) variant CAs after superposition
    rotation: np.ndarray
    translation: np.ndarray
    rms: float  # RMSD over the atoms used in the fit
    rms_all_matched: float  # RMSD over every matched atom under that transform
    n_ref: int
    n_var: int
    fitted_mask: np.ndarray  # bool (N,), False where an atom was trimmed
    substitutions: List[Tuple[ResidueKey, str, str]] = field(default_factory=list)

    @property
    def n_matched(self) -> int:
        return len(self.matched_keys)

    @property
    def n_fitted(self) -> int:
        return int(self.fitted_mask.sum())

    @property
    def coverage(self) -> float:
        """Fraction of the reference's alpha-carbons that had a counterpart."""
        return self.n_matched / self.n_ref if self.n_ref else 0.0

    @property
    def substitution_label(self) -> str:
        return ";".join(
            f"{chain}:{ref_name}{res_seq}{icode.strip()}{var_name}"
            for (chain, res_seq, icode), ref_name, var_name in self.substitutions
        )


def _sorted_key(key: ResidueKey):
    chain_id, res_seq, icode = key
    return (chain_id, res_seq, icode)


def _apply(coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    """Biopython's transform convention: coord . rot + tran."""
    return np.dot(coords, rotation) + translation


def _fit(ref_atoms, var_atoms) -> Tuple[np.ndarray, np.ndarray, float]:
    superimposer = Superimposer()
    superimposer.set_atoms(ref_atoms, var_atoms)
    rotation, translation = superimposer.rotran
    return (
        np.asarray(rotation, dtype=float),
        np.asarray(translation, dtype=float),
        float(superimposer.rms),
    )


def align(
    reference: StructureData,
    variant: StructureData,
    trim_outliers: bool = False,
    trim_cutoff: float = 2.0,
    trim_iterations: int = 3,
    min_coverage: float = 0.5,
) -> AlignmentResult:
    """Superimpose ``variant`` onto ``reference`` over their shared alpha-carbons.

    With ``trim_outliers`` enabled the fit is refined iteratively, discarding
    atoms that deviate by more than ``trim_cutoff`` angstroms so that a single
    mobile loop or disordered terminus cannot drag the whole superposition.
    Trimmed atoms are still measured -- they are excluded from the fit, not from
    the report. Default is off, matching the design document's plain global fit.
    """
    ref_keys = set(reference.ca_keys())
    var_keys = set(variant.ca_keys())
    matched = sorted(ref_keys & var_keys, key=_sorted_key)

    n_ref, n_var = len(ref_keys), len(var_keys)
    if len(matched) < MIN_ATOMS_FOR_FIT:
        raise AlignmentError(
            f"only {len(matched)} shared alpha-carbon(s) between "
            f"{reference.name} ({n_ref}) and {variant.name} ({n_var}); "
            "need at least 3 to superimpose"
        )

    coverage = len(matched) / n_ref if n_ref else 0.0
    if coverage < min_coverage:
        raise AlignmentError(
            f"{variant.name} shares only {coverage:.0%} of the reference's residues "
            f"({len(matched)}/{n_ref}); below the {min_coverage:.0%} minimum. "
            "Different protein, or badly mismatched numbering?"
        )

    ref_atoms = [reference.residues[k].ca_atom for k in matched]
    var_atoms = [variant.residues[k].ca_atom for k in matched]
    ref_coords = np.array([reference.residues[k].ca_coord for k in matched], dtype=float)
    var_coords = np.array([variant.residues[k].ca_coord for k in matched], dtype=float)

    fitted_mask = np.ones(len(matched), dtype=bool)
    rotation, translation, rms = _fit(ref_atoms, var_atoms)

    if trim_outliers:
        for _ in range(max(0, trim_iterations)):
            aligned = _apply(var_coords, rotation, translation)
            deviations = np.linalg.norm(ref_coords - aligned, axis=1)
            candidate = deviations <= trim_cutoff
            if candidate.sum() < MIN_ATOMS_FOR_FIT or np.array_equal(candidate, fitted_mask):
                break
            fitted_mask = candidate
            indices = np.flatnonzero(fitted_mask)
            rotation, translation, rms = _fit(
                [ref_atoms[i] for i in indices], [var_atoms[i] for i in indices]
            )

    var_coords_aligned = _apply(var_coords, rotation, translation)
    rms_all_matched = float(
        np.sqrt(np.mean(np.sum((ref_coords - var_coords_aligned) ** 2, axis=1)))
    )

    substitutions = [
        (k, reference.residues[k].res_name, variant.residues[k].res_name)
        for k in matched
        if reference.residues[k].res_name != variant.residues[k].res_name
    ]

    return AlignmentResult(
        matched_keys=matched,
        ref_coords=ref_coords,
        var_coords=var_coords,
        var_coords_aligned=var_coords_aligned,
        rotation=rotation,
        translation=translation,
        rms=rms,
        rms_all_matched=rms_all_matched,
        n_ref=n_ref,
        n_var=n_var,
        fitted_mask=fitted_mask,
        substitutions=substitutions,
    )
