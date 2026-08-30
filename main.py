"""StructScan - execution entry point for the batch analysis.

    python main.py
    python main.py --reference data/wildtype_reference.pdb --variants data/variants
    python main.py --trim-outliers --top-residues 8

Runs every structure in the variants directory against the wild-type reference
and writes a prioritised triage report to the output directory.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from engine import __version__
from engine.aligner import AlignmentError, align
from engine.geometry import compute_metrics
from engine.parser import StructureParseError, load_structure
from engine.ranker import (
    DEFAULT_THRESHOLDS,
    DEFAULT_WEIGHTS,
    aggregate_hotspots,
    rank_variants,
    top_residues_per_variant,
)
from engine.reporter import log_header, log_skip, log_summary, log_variant, write_reports

STRUCTURE_SUFFIXES = (".pdb", ".ent", ".cif", ".mmcif")

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_REFERENCE = REPO_ROOT / "data" / "wildtype_reference.pdb"
DEFAULT_VARIANTS = REPO_ROOT / "data" / "variants"
DEFAULT_OUTPUT = REPO_ROOT / "output"


class BatchResult(Tuple):
    pass


def discover_variants(
    variants_dir, reference_path: Optional[Path] = None
) -> List[Path]:
    """Every structure file in the directory, excluding the reference itself."""
    variants_dir = Path(variants_dir)
    if not variants_dir.is_dir():
        raise NotADirectoryError(f"Variants directory not found: {variants_dir}")

    reference_resolved = Path(reference_path).resolve() if reference_path else None
    found = [
        path
        for path in sorted(variants_dir.iterdir())
        if path.is_file()
        and path.suffix.lower() in STRUCTURE_SUFFIXES
        and (reference_resolved is None or path.resolve() != reference_resolved)
    ]
    return found


def run_batch(
    reference_path,
    variants_dir,
    top_residues: int = 5,
    trim_outliers: bool = False,
    trim_cutoff: float = 2.0,
    min_coverage: float = 0.5,
    displacement_threshold: float = 1.0,
    dihedral_threshold: float = 30.0,
    weights: Optional[Dict[str, float]] = None,
    thresholds: Optional[Dict[str, float]] = None,
    chains: Optional[Sequence[str]] = None,
    quiet: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[Dict[str, str]]]:
    """Run the full pipeline. Returns (ranked, per_residue, hotspots, failures).

    One unreadable or mismatched file does not abort the batch -- it is recorded
    in ``failures`` and the run continues.
    """
    reference_path = Path(reference_path)
    chain_list = list(chains) if chains else None

    reference = load_structure(reference_path, chains=chain_list)
    variant_paths = discover_variants(variants_dir, reference_path)

    if not quiet:
        log_header(reference.name, reference.describe(), len(variant_paths))

    per_residue_frames: List[pd.DataFrame] = []
    summaries: List[Dict[str, object]] = []
    failures: List[Dict[str, str]] = []

    total = len(variant_paths)
    for index, path in enumerate(variant_paths, start=1):
        try:
            variant = load_structure(path, chains=chain_list)
            alignment = align(
                reference,
                variant,
                trim_outliers=trim_outliers,
                trim_cutoff=trim_cutoff,
                min_coverage=min_coverage,
            )
            frame, summary = compute_metrics(
                reference,
                variant,
                alignment,
                displacement_count_threshold=displacement_threshold,
                dihedral_count_threshold=dihedral_threshold,
            )
        except (StructureParseError, AlignmentError) as exc:
            failures.append({"file": path.name, "reason": str(exc)})
            if not quiet:
                log_skip(index, total, path.stem, str(exc))
            continue
        except Exception as exc:  # noqa: BLE001 - never let one file kill the batch
            failures.append({"file": path.name, "reason": f"unexpected error: {exc}"})
            if not quiet:
                log_skip(index, total, path.stem, f"unexpected error: {exc}")
            continue

        per_residue_frames.append(frame)
        summaries.append(summary)
        if not quiet:
            log_variant(index, total, variant.name, summary)

    per_residue = (
        pd.concat(per_residue_frames, ignore_index=True)
        if per_residue_frames
        else pd.DataFrame()
    )
    ranked = rank_variants(summaries, weights=weights, thresholds=thresholds)
    hotspots = aggregate_hotspots(
        per_residue,
        displacement_threshold=displacement_threshold,
        dihedral_threshold=dihedral_threshold,
    )

    if not ranked.empty:
        top_labels = top_residues_per_variant(per_residue, top_n=top_residues)
        ranked["top_deviating_residues"] = ranked["variant"].map(top_labels).fillna("")

    return ranked, per_residue, hotspots, failures


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="structscan",
        description=(
            "Batch-screen macromolecular variant structures against a wild-type "
            "reference and rank them by structural deviation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE,
                        help="wild-type reference structure (.pdb/.cif)")
    parser.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS,
                        help="directory of variant structures")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                        help="directory for the CSV reports")
    parser.add_argument("--chains", nargs="*", default=None,
                        help="restrict the analysis to these chain IDs")
    parser.add_argument("--top-residues", type=int, default=5,
                        help="how many worst-deviating residues to name per variant")

    fit = parser.add_argument_group("superposition")
    fit.add_argument("--trim-outliers", action="store_true",
                     help="iteratively refit, excluding high-deviation atoms, so one "
                          "mobile loop cannot distort the whole superposition")
    fit.add_argument("--trim-cutoff", type=float, default=2.0,
                     help="angstrom cutoff for --trim-outliers")
    fit.add_argument("--min-coverage", type=float, default=0.5,
                     help="minimum fraction of reference residues a variant must match")

    metrics = parser.add_argument_group("metric thresholds")
    metrics.add_argument("--displacement-threshold", type=float, default=1.0,
                         help="angstroms; counts a residue as displaced")
    metrics.add_argument("--dihedral-threshold", type=float, default=30.0,
                         help="degrees; counts a residue as torsioned")

    score = parser.add_argument_group("composite score weights")
    score.add_argument("--weight-rmsd", type=float, default=DEFAULT_WEIGHTS["global_rmsd"])
    score.add_argument("--weight-displacement", type=float,
                       default=DEFAULT_WEIGHTS["max_ca_displacement"])
    score.add_argument("--weight-dihedral", type=float,
                       default=DEFAULT_WEIGHTS["max_dihedral_delta"])

    triage = parser.add_argument_group("triage tag cutoffs")
    triage.add_argument("--rmsd-hotspot", type=float, default=DEFAULT_THRESHOLDS["rmsd_hotspot"])
    triage.add_argument("--displacement-hotspot", type=float,
                        default=DEFAULT_THRESHOLDS["displacement_hotspot"])
    triage.add_argument("--dihedral-hotspot", type=float,
                        default=DEFAULT_THRESHOLDS["dihedral_hotspot"])

    parser.add_argument("--quiet", action="store_true", help="suppress the progress log")
    parser.add_argument("--version", action="version", version=f"StructScan {__version__}")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    weights = {
        "global_rmsd": args.weight_rmsd,
        "max_ca_displacement": args.weight_displacement,
        "max_dihedral_delta": args.weight_dihedral,
    }
    thresholds = {
        "rmsd_hotspot": args.rmsd_hotspot,
        "displacement_hotspot": args.displacement_hotspot,
        "dihedral_hotspot": args.dihedral_hotspot,
    }

    try:
        ranked, per_residue, hotspots, failures = run_batch(
            reference_path=args.reference,
            variants_dir=args.variants,
            top_residues=args.top_residues,
            trim_outliers=args.trim_outliers,
            trim_cutoff=args.trim_cutoff,
            min_coverage=args.min_coverage,
            displacement_threshold=args.displacement_threshold,
            dihedral_threshold=args.dihedral_threshold,
            weights=weights,
            thresholds=thresholds,
            chains=args.chains,
            quiet=args.quiet,
        )
    except (StructureParseError, NotADirectoryError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        print(
            "hint: populate data/ first -- 'python tools/make_synthetic.py' for offline "
            "test structures, or 'python tools/fetch_data.py' to download the T4 "
            "lysozyme mutant series from the RCSB PDB.",
            file=sys.stderr,
        )
        return 2

    if not args.quiet:
        log_summary(ranked, hotspots)

    write_reports(args.output, ranked, per_residue, hotspots, failures)

    if ranked.empty:
        print("error: no variants could be analysed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
