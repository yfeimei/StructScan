"""Correctness tests for the StructScan pipeline.

Every assertion here is anchored to a perturbation that tools/make_synthetic.py
planted deliberately, so "the tool found something" can be distinguished from
"the tool found the right thing".

    python -m pytest tests -q          # with pytest installed
    python tests/test_pipeline.py      # standalone, no pytest needed
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from engine.geometry import angular_delta, rmsd  # noqa: E402
from engine.ranker import TAG_HOTSPOT, TAG_STABLE  # noqa: E402
from main import run_batch  # noqa: E402
from tools.make_synthetic import generate  # noqa: E402

LENGTH = 60

_FIXTURES: dict = {}


def fixtures() -> dict:
    """Build the synthetic battery once and run the pipeline over it."""
    if _FIXTURES:
        return _FIXTURES

    tmpdir = Path(tempfile.mkdtemp(prefix="structscan_test_"))
    manifest = generate(tmpdir, length=LENGTH)
    ranked, per_residue, hotspots, failures = run_batch(
        tmpdir / "wildtype_reference.pdb", tmpdir / "variants", quiet=True
    )
    _FIXTURES.update(
        tmpdir=tmpdir,
        manifest=manifest,
        ranked=ranked,
        per_residue=per_residue,
        hotspots=hotspots,
        failures=failures,
    )
    return _FIXTURES


def _row(variant: str):
    ranked = fixtures()["ranked"]
    matches = ranked[ranked["variant"] == variant]
    assert not matches.empty, f"{variant} missing from the ranked report"
    return matches.iloc[0]


def _residues(variant: str):
    per_residue = fixtures()["per_residue"]
    return per_residue[per_residue["variant"] == variant]


# --------------------------------------------------------------------------
# Unit-level maths
# --------------------------------------------------------------------------

def test_angular_delta_wraps_around_the_periodic_boundary():
    # The plain |a - b| in the design document returns 358 here. It is 2.
    assert angular_delta(179.0, -179.0) == 2.0
    assert angular_delta(-179.0, 179.0) == 2.0
    assert angular_delta(350.0, 10.0) == 20.0
    assert angular_delta(0.0, 0.0) == 0.0
    assert angular_delta(0.0, 180.0) == 180.0
    # Never exceeds a half turn.
    for a in np.linspace(-360, 360, 97):
        for b in np.linspace(-360, 360, 97):
            assert 0.0 <= angular_delta(float(a), float(b)) <= 180.0


def test_angular_delta_is_none_when_an_angle_is_undefined():
    assert angular_delta(None, 12.0) is None
    assert angular_delta(12.0, None) is None


def test_rmsd_of_identical_coordinates_is_zero():
    coords = np.random.default_rng(0).normal(size=(25, 3))
    assert rmsd(coords, coords) == 0.0


def test_rmsd_matches_the_closed_form():
    a = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    b = np.array([[3.0, 4.0, 0.0], [0.0, 0.0, 0.0]])
    assert abs(rmsd(a, b) - np.sqrt(25.0 / 2)) < 1e-12


# --------------------------------------------------------------------------
# Pipeline behaviour against planted ground truth
# --------------------------------------------------------------------------

def test_every_synthetic_variant_is_analysed():
    data = fixtures()
    assert not data["failures"], f"unexpected failures: {data['failures']}"
    assert len(data["ranked"]) == len(data["manifest"]["variants"])


def test_rigid_body_motion_is_removed_by_superposition():
    # v06 is the reference rotated and translated. Any residual RMSD beyond
    # file precision is a bug in the superposition, not a property of the
    # structure. PDB stores coordinates to 3 decimal places, so a round trip
    # through the format costs up to ~0.0005 A per axis (~0.001 A in 3D);
    # the tolerance below is that floor, not a fudge factor.
    row = _row("v06_rigid_body")
    assert row["global_rmsd"] < 2e-3, row["global_rmsd"]
    assert row["max_ca_displacement"] < 5e-3, row["max_ca_displacement"]
    assert row["triage_tag"] == TAG_STABLE


def test_coordinate_noise_alone_is_tagged_stable():
    row = _row("v01_stable")
    assert row["triage_tag"] == TAG_STABLE
    assert row["global_rmsd"] < 0.2


def test_planted_loop_shift_is_localised_to_the_right_residues():
    row = _row("v02_loop_shift_20_26")
    assert row["triage_tag"] == TAG_HOTSPOT

    residues = _residues("v02_loop_shift_20_26")
    top = residues.nlargest(7, "ca_displacement")["res_seq"].tolist()
    assert set(top) <= set(range(20, 27)), f"expected residues 20-26, got {sorted(top)}"

    # Residues far from the perturbation should barely move.
    quiet = residues[(residues["res_seq"] < 18) | (residues["res_seq"] > 28)]
    assert quiet["ca_displacement"].max() < 1.5


def test_larger_planted_shift_is_also_localised():
    row = _row("v03_core_shift_40_48")
    assert row["triage_tag"] == TAG_HOTSPOT
    residues = _residues("v03_core_shift_40_48")
    top = residues.nlargest(9, "ca_displacement")["res_seq"].tolist()
    assert set(top) <= set(range(40, 49)), f"expected residues 40-48, got {sorted(top)}"


def test_bigger_perturbation_outranks_smaller_one():
    # 5 A over 9 residues must score above 2.5 A over 7.
    assert _row("v03_core_shift_40_48")["rank"] < _row("v02_loop_shift_20_26")["rank"]
    assert _row("v02_loop_shift_20_26")["rank"] < _row("v01_stable")["rank"]


def test_torsion_kink_is_detected_at_the_planted_residue():
    residues = _residues("v05_dihedral_kink_30")
    worst = residues.nlargest(1, "max_dihedral_delta").iloc[0]
    assert worst["res_seq"] == 30, f"kink reported at {worst['res_seq']}, planted at 30"
    # phi was moved +60 and psi -60; the shortest-separation metric recovers 60.
    assert abs(worst["max_dihedral_delta"] - 60.0) < 1.0
    assert _row("v05_dihedral_kink_30")["triage_tag"] == TAG_HOTSPOT


def test_missing_residues_do_not_break_alignment():
    # v07 has residues 10-14 deleted; matching is by position key, so the
    # remaining residues must still line up with their true counterparts.
    row = _row("v07_gapped")
    assert row["n_residues_matched"] == LENGTH - 5
    assert row["triage_tag"] == TAG_STABLE
    assert row["global_rmsd"] < 0.2

    residues = _residues("v07_gapped")
    assert not set(range(10, 15)) & set(residues["res_seq"])


def test_substitutions_are_recorded_not_treated_as_mismatches():
    row = _row("v08_substitutions")
    assert row["n_residues_matched"] == LENGTH
    assert row["n_substitutions"] == 3
    for position in (15, 33, 50):
        assert f"{position}TRP" in row["substitutions"]


def test_outlier_trimming_isolates_a_mobile_terminus():
    data = fixtures()
    plain = _row("v04_terminal_flex")["global_rmsd"]

    trimmed_ranked, _, _, _ = run_batch(
        data["tmpdir"] / "wildtype_reference.pdb",
        data["tmpdir"] / "variants",
        trim_outliers=True,
        quiet=True,
    )
    trimmed = trimmed_ranked[trimmed_ranked["variant"] == "v04_terminal_flex"].iloc[0]

    # Excluding the mobile tail from the fit must tighten the core superposition.
    assert trimmed["global_rmsd"] < plain
    # ...without hiding the displacement itself, which is still measured.
    assert trimmed["max_ca_displacement"] > 3.0


def test_hotspot_aggregation_surfaces_recurrently_perturbed_positions():
    hotspots = fixtures()["hotspots"]
    assert not hotspots.empty
    top = hotspots.nlargest(10, "max_ca_displacement")["res_seq"].tolist()
    planted = set(range(20, 27)) | set(range(40, 49)) | set(range(LENGTH - 5, LENGTH + 1))
    assert set(top) & planted, f"no planted positions in the top hotspots: {sorted(top)}"


def test_report_columns_are_present_and_finite():
    ranked = fixtures()["ranked"]
    for column in ("rank", "variant", "triage_tag", "composite_score",
                   "global_rmsd", "max_ca_displacement", "top_deviating_residues"):
        assert column in ranked.columns, f"missing column {column}"
    assert np.isfinite(ranked["composite_score"]).all()
    assert np.isfinite(ranked["global_rmsd"]).all()


# --------------------------------------------------------------------------

def _run_standalone() -> int:
    tests = [
        (name, obj)
        for name, obj in sorted(globals().items())
        if name.startswith("test_") and callable(obj)
    ]
    failed = []
    for name, test in tests:
        try:
            test()
        except AssertionError as exc:
            failed.append((name, exc))
            print(f"FAIL  {name}\n        {exc}")
        except Exception as exc:  # noqa: BLE001
            failed.append((name, exc))
            print(f"ERROR {name}\n        {type(exc).__name__}: {exc}")
        else:
            print(f"ok    {name}")

    tmpdir = _FIXTURES.get("tmpdir")
    if tmpdir:
        shutil.rmtree(tmpdir, ignore_errors=True)

    print(f"\n{len(tests) - len(failed)}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
