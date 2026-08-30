"""Module 4: Ranking & Triage Sorter.

Compiles the per-variant fingerprints into a composite impact score, sorts the
batch from highest deviation to lowest, and applies triage tags so benign
background structures can be filtered out.

Two separate mechanisms, deliberately:

  * the composite score is a weighted sum of z-scores, which answers "how does
    this variant rank *against the rest of this batch*" and is the sort key;
  * the triage tag is driven by absolute angstrom / degree thresholds, which
    answers "is this deviation large *in physical terms*".

Ranking alone would tag the top of every batch as interesting even when nothing
in it moved. Thresholds alone would give no ordering within a tag. The report
carries both.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

# Relative contribution of each fingerprint component to the composite score.
DEFAULT_WEIGHTS: Dict[str, float] = {
    "global_rmsd": 1.0,
    "max_ca_displacement": 1.0,
    "max_dihedral_delta": 0.5,
}

# Absolute cut-offs for the triage tag. A variant is tagged at the highest
# level for which it trips *any* criterion.
DEFAULT_THRESHOLDS: Dict[str, float] = {
    "rmsd_hotspot": 1.0,  # angstroms
    "rmsd_moderate": 0.5,
    "displacement_hotspot": 2.0,  # angstroms
    "displacement_moderate": 1.0,
    "dihedral_hotspot": 45.0,  # degrees
    "dihedral_moderate": 20.0,
}

TAG_HOTSPOT = "High Priority Hotspot"
TAG_MODERATE = "Moderate Deviation"
TAG_STABLE = "Stable"

_EPSILON = 1e-9


def _zscore(series: pd.Series) -> pd.Series:
    """Z-score that degrades gracefully on tiny or uniform batches.

    A single variant, or a batch where every value is identical, has no spread
    to normalise against; those cases collapse to 0 rather than NaN/inf so the
    composite stays finite and the absolute-threshold tag carries the meaning.
    """
    values = pd.to_numeric(series, errors="coerce")
    if values.notna().sum() < 2:
        return pd.Series(0.0, index=series.index)
    std = values.std(ddof=0)
    if not np.isfinite(std) or std < _EPSILON:
        return pd.Series(0.0, index=series.index)
    return ((values - values.mean()) / std).fillna(0.0)


def _classify(row: pd.Series, thresholds: Dict[str, float]) -> str:
    rmsd = row.get("global_rmsd", np.nan)
    disp = row.get("max_ca_displacement", np.nan)
    dihedral = row.get("max_dihedral_delta", np.nan)

    def over(value, key) -> bool:
        return bool(np.isfinite(value)) and value > thresholds[key]

    if (
        over(rmsd, "rmsd_hotspot")
        or over(disp, "displacement_hotspot")
        or over(dihedral, "dihedral_hotspot")
    ):
        return TAG_HOTSPOT
    if (
        over(rmsd, "rmsd_moderate")
        or over(disp, "displacement_moderate")
        or over(dihedral, "dihedral_moderate")
    ):
        return TAG_MODERATE
    return TAG_STABLE


def rank_variants(
    summaries: Iterable[Dict[str, object]],
    weights: Optional[Dict[str, float]] = None,
    thresholds: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    """Score, tag and sort the batch. Highest deviation first."""
    weights = dict(weights or DEFAULT_WEIGHTS)
    thresholds = {**DEFAULT_THRESHOLDS, **(thresholds or {})}

    frame = pd.DataFrame(list(summaries))
    if frame.empty:
        return frame

    total_weight = sum(abs(w) for w in weights.values()) or 1.0
    composite = pd.Series(0.0, index=frame.index)
    for column, weight in weights.items():
        if column not in frame.columns:
            continue
        z = _zscore(frame[column])
        frame[f"z_{column}"] = z
        composite = composite + weight * z

    frame["composite_score"] = composite / total_weight
    frame["triage_tag"] = frame.apply(lambda row: _classify(row, thresholds), axis=1)

    frame = frame.sort_values(
        ["composite_score", "global_rmsd"], ascending=False, kind="mergesort"
    ).reset_index(drop=True)
    frame.insert(0, "rank", np.arange(1, len(frame) + 1))
    return frame


def top_residues_per_variant(
    per_residue: pd.DataFrame, top_n: int = 5
) -> Dict[str, str]:
    """Compact 'A:45(2.31A)' summary of each variant's worst-deviating residues."""
    if per_residue.empty:
        return {}

    labels: Dict[str, str] = {}
    for variant, group in per_residue.groupby("variant", sort=False):
        worst = group.nlargest(top_n, "ca_displacement")
        labels[str(variant)] = "; ".join(
            f"{row.position}({row.ca_displacement:.2f}A)" for row in worst.itertuples()
        )
    return labels


def aggregate_hotspots(
    per_residue: pd.DataFrame,
    displacement_threshold: float = 1.0,
    dihedral_threshold: float = 30.0,
) -> pd.DataFrame:
    """Batch-level view: which positions move across many variants?

    A residue that shifts in one variant is that variant's problem. A residue
    that shifts across a large fraction of the batch is a structural hotspot,
    which is what the pipeline is ultimately looking for.
    """
    if per_residue.empty:
        return pd.DataFrame()

    frame = per_residue.copy()
    frame["_displaced"] = frame["ca_displacement"] > displacement_threshold
    frame["_torsioned"] = frame["max_dihedral_delta"] > dihedral_threshold

    grouped = frame.groupby(["chain", "res_seq", "icode", "position"], sort=False)
    hotspots = grouped.agg(
        ref_resname=("ref_resname", "first"),
        n_variants=("variant", "nunique"),
        n_variants_displaced=("_displaced", "sum"),
        n_variants_torsioned=("_torsioned", "sum"),
        max_ca_displacement=("ca_displacement", "max"),
        mean_ca_displacement=("ca_displacement", "mean"),
        max_dihedral_delta=("max_dihedral_delta", "max"),
        n_substituted=("is_substitution", "sum"),
    ).reset_index()

    hotspots["fraction_displaced"] = (
        hotspots["n_variants_displaced"] / hotspots["n_variants"].replace(0, np.nan)
    ).fillna(0.0)

    return hotspots.sort_values(
        ["fraction_displaced", "max_ca_displacement"], ascending=False, kind="mergesort"
    ).reset_index(drop=True)
