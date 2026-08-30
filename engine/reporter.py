"""Module 5: Report Generator.

Terminal progress logging during the batch run, and CSV export of the finished
tables for direct lab review.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from .ranker import TAG_HOTSPOT, TAG_MODERATE, TAG_STABLE

RANKED_CSV = "ranked_hotspots.csv"
PER_RESIDUE_CSV = "per_residue_detail.csv"
HOTSPOT_CSV = "residue_hotspot_frequency.csv"
FAILURES_CSV = "failed_structures.csv"

# Column order for the headline report, so the file opens legibly in Excel.
RANKED_COLUMNS = [
    "rank",
    "variant",
    "triage_tag",
    "composite_score",
    "global_rmsd",
    "max_ca_displacement",
    "mean_ca_displacement",
    "max_dihedral_delta",
    "n_residues_displaced",
    "n_residues_torsioned",
    "top_deviating_residues",
    "n_substitutions",
    "substitutions",
    "n_residues_matched",
    "coverage",
    "source_file",
]


def _out(message: str = "") -> None:
    print(message, file=sys.stdout, flush=True)


def log_header(reference_name: str, reference_desc: str, n_variants: int) -> None:
    _out("=" * 78)
    _out("StructScan - batch structural variant triage")
    _out("=" * 78)
    _out(f"  Reference : {reference_name}  [{reference_desc}]")
    _out(f"  Variants  : {n_variants} file(s)")
    _out("-" * 78)


def log_variant(index: int, total: int, name: str, summary: Dict[str, object]) -> None:
    _out(
        f"  [{index:>3}/{total}] {name:<28} "
        f"RMSD {summary['global_rmsd']:6.3f} A   "
        f"max shift {summary['max_ca_displacement']:6.3f} A   "
        f"max dtheta {summary['max_dihedral_delta']:6.1f} deg   "
        f"matched {summary['n_residues_matched']}"
    )


def log_skip(index: int, total: int, name: str, reason: str) -> None:
    _out(f"  [{index:>3}/{total}] {name:<28} SKIPPED - {reason}")


def log_summary(ranked: pd.DataFrame, hotspots: pd.DataFrame, top_n: int = 10) -> None:
    _out("-" * 78)
    if ranked.empty:
        _out("  No variants were successfully analysed.")
        return

    counts = ranked["triage_tag"].value_counts()
    _out(
        f"  Triage: {counts.get(TAG_HOTSPOT, 0)} hotspot, "
        f"{counts.get(TAG_MODERATE, 0)} moderate, "
        f"{counts.get(TAG_STABLE, 0)} stable  "
        f"(of {len(ranked)} analysed)"
    )
    _out("")
    _out(f"  Top {min(top_n, len(ranked))} by composite impact score:")
    _out(f"    {'#':>3}  {'variant':<28} {'score':>7}  {'RMSD':>6}  {'tag'}")
    for row in ranked.head(top_n).itertuples():
        _out(
            f"    {row.rank:>3}  {str(row.variant):<28} "
            f"{row.composite_score:>7.3f}  {row.global_rmsd:>6.3f}  {row.triage_tag}"
        )

    if not hotspots.empty:
        recurrent = hotspots[hotspots["n_variants_displaced"] > 0].head(top_n)
        if not recurrent.empty:
            _out("")
            _out(f"  Most frequently perturbed positions across the batch:")
            for row in recurrent.itertuples():
                _out(
                    f"    {row.position:<10} {row.ref_resname:<4} "
                    f"displaced in {int(row.n_variants_displaced):>3}/{int(row.n_variants)} variants   "
                    f"max {row.max_ca_displacement:.2f} A"
                )


def write_reports(
    output_dir,
    ranked: pd.DataFrame,
    per_residue: pd.DataFrame,
    hotspots: pd.DataFrame,
    failures: Optional[List[Dict[str, str]]] = None,
) -> Dict[str, Path]:
    """Write the CSV tables. Returns a map of label -> path written."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    written: Dict[str, Path] = {}

    if not ranked.empty:
        ordered = [c for c in RANKED_COLUMNS if c in ranked.columns]
        remaining = [c for c in ranked.columns if c not in ordered]
        path = output_dir / RANKED_CSV
        ranked[ordered + remaining].to_csv(path, index=False, float_format="%.4f")
        written["ranked"] = path

    if not per_residue.empty:
        path = output_dir / PER_RESIDUE_CSV
        per_residue.to_csv(path, index=False, float_format="%.4f")
        written["per_residue"] = path

    if not hotspots.empty:
        path = output_dir / HOTSPOT_CSV
        hotspots.to_csv(path, index=False, float_format="%.4f")
        written["hotspots"] = path

    if failures:
        path = output_dir / FAILURES_CSV
        pd.DataFrame(failures).to_csv(path, index=False)
        written["failures"] = path

    _out("")
    for label, path in written.items():
        _out(f"  wrote {label:<12} -> {path}")
    _out("=" * 78)
    return written
