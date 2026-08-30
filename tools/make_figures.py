"""Render the two static figures used in the README and the project page.

    python tools/make_figures.py

Figure 1 (validation)  - per-residue signal for three synthetic variants, with the
                         perturbation window read from data/ground_truth.json shaded
                         behind the measured curve. Shows the pipeline recovers a
                         planted perturbation at the position it was planted.
Figure 2 (T4 ranking)  - the real-crystal batch ordered by composite impact score,
                         coloured by triage tag.

Each figure is written twice, for light and dark page surfaces, so the README can
serve the right one with <picture>. Inputs are the CSVs that main.py already
writes; this script does no analysis of its own.

Colours are not chosen by eye. The triage ramp is a single-hue ordinal ramp
(validated: monotone lightness, adjacent dL >= 0.06, light end clears the
surface). A red/amber/green severity ramp was measured first and rejected -
green #0ca30c against red #d03b3b is dE 4.1 under simulated deuteranopia, which
is indistinguishable for roughly 1 in 12 male readers.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

REPO_ROOT = Path(__file__).resolve().parent.parent
GROUND_TRUTH = REPO_ROOT / "data" / "ground_truth.json"
SYNTH_DETAIL = REPO_ROOT / "output" / "per_residue_detail.csv"
T4_RANKED = REPO_ROOT / "output" / "t4_lysozyme" / "ranked_hotspots.csv"
FIGURE_DIR = REPO_ROOT / "docs" / "figures"

# Palette slots, per mode. Both columns are selected for their own surface -
# the dark column is the same hue re-stepped, not an automatic inversion.
THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "ink": "#0b0b0b",
        "ink_secondary": "#52514e",
        "muted": "#898781",
        "grid": "#e1e0d9",
        "axis": "#c3c2b7",
        "series": "#2a78d6",
        "band": "#e8e7e1",
        # Stable -> Moderate -> High Priority, light to dark.
        "triage": ("#86b6ef", "#2a78d6", "#0d366b"),
    },
    "dark": {
        "surface": "#1a1a19",
        "ink": "#ffffff",
        "ink_secondary": "#c3c2b7",
        "muted": "#898781",
        "grid": "#2c2c2a",
        "axis": "#383835",
        "series": "#3987e5",
        "band": "#282826",
        "triage": ("#184f95", "#3987e5", "#9ec5f4"),
    },
}

TRIAGE_ORDER = ["Stable", "Moderate Deviation", "High Priority Hotspot"]

# Which synthetic variants to show, and which measured column carries their signal.
# A rigid window shift lands in alpha-carbon displacement; a torsion kink lands in
# the dihedral delta, so each panel plots the metric its perturbation actually moves.
PANELS = [
    ("v02_loop_shift_20_26", "ca_displacement", "CA displacement", "A"),
    ("v03_core_shift_40_48", "ca_displacement", "CA displacement", "A"),
    ("v05_dihedral_kink_30", "max_dihedral_delta", "backbone torsion change", "deg"),
]


def style(theme: dict) -> None:
    """Recessive chrome, system sans, no top/right spines."""
    plt.rcParams.update({
        "figure.facecolor": theme["surface"],
        "axes.facecolor": theme["surface"],
        "savefig.facecolor": theme["surface"],
        "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"],
        "text.color": theme["ink"],
        "axes.labelcolor": theme["ink_secondary"],
        "axes.edgecolor": theme["axis"],
        "xtick.color": theme["muted"],
        "ytick.color": theme["muted"],
        "xtick.labelcolor": theme["ink_secondary"],
        "ytick.labelcolor": theme["ink_secondary"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": theme["grid"],
        "grid.linewidth": 0.8,
        "font.size": 9,
        "axes.titlesize": 10,
        "figure.dpi": 200,
    })


def planted_window(truth: dict, variant: str) -> tuple[list[int], str]:
    """The residues the generator perturbed, plus a human label for the band."""
    spec = truth["variants"][variant]
    residues = spec.get("hotspot_residues", [])
    if spec.get("kind") == "torsion":
        return residues, f"planted {spec['delta_degrees']:.0f} deg kink at {spec['kink_resseq']}"
    lo, hi = spec["window"]
    return residues, f"planted {spec['magnitude']:.1f} A shift, residues {lo}-{hi}"


def figure_validation(mode: str) -> Path:
    truth = json.loads(GROUND_TRUTH.read_text())
    detail = pd.read_csv(SYNTH_DETAIL)
    theme = THEMES[mode]
    style(theme)

    fig, axes = plt.subplots(3, 1, figsize=(7.6, 7.4), sharex=True)

    for ax, (variant, column, metric, unit) in zip(axes, PANELS):
        rows = detail[detail["variant"] == variant].sort_values("res_seq")
        residues, band_label = planted_window(truth, variant)
        values = rows[column].fillna(0.0)

        # Ground-truth window behind the measured curve, so overlap is the claim.
        # The band label sits at the floor of the axes; the peak label sits at the
        # peak. Keeping them on opposite edges is what stops them colliding.
        if residues:
            ax.axvspan(min(residues) - 0.5, max(residues) + 0.5,
                       color=theme["band"], zorder=0, linewidth=0)
            # A one-residue band is narrower than its own label, and the spike would
            # run straight through centred text - anchor it beside the band instead.
            narrow = max(residues) - min(residues) < 4
            ax.annotate(band_label,
                        xy=(max(residues) + 1.5 if narrow
                            else (min(residues) + max(residues)) / 2, 0.04),
                        xycoords=("data", "axes fraction"),
                        ha="left" if narrow else "center", va="bottom",
                        fontsize=8, color=theme["ink_secondary"])

        ax.plot(rows["res_seq"], values, color=theme["series"], linewidth=2, zorder=3)

        # One direct label on the peak, not a number on every point. Flip the label
        # inward when the peak is in the right half, or it runs off the figure.
        peak = values.idxmax()
        px, py = rows.loc[peak, "res_seq"], values.loc[peak]
        ax.plot([px], [py], "o", markersize=8, color=theme["series"],
                markeredgecolor=theme["surface"], markeredgewidth=2, zorder=4)
        midpoint = (rows["res_seq"].min() + rows["res_seq"].max()) / 2
        flip = px > midpoint
        ax.annotate(f"measured peak {py:.2f} {unit} at residue {int(px)}",
                    xy=(px, py), xytext=(-10 if flip else 10, 6),
                    textcoords="offset points",
                    ha="right" if flip else "left", va="center",
                    fontsize=8, color=theme["ink"])

        ax.set_title(f"{variant}   -   {metric}", loc="left", color=theme["ink"], pad=6)
        ax.set_ylabel(f"{metric} ({unit})")
        # Headroom so the peak label clears the panel title.
        ax.set_ylim(bottom=0, top=float(values.max()) * 1.18)
        ax.set_xlim(rows["res_seq"].min(), rows["res_seq"].max())
        ax.set_axisbelow(True)
        ax.grid(axis="x", visible=False)

    axes[-1].set_xlabel("residue number")

    fig.suptitle("StructScan recovers planted perturbations at the planted position",
                 x=0.055, y=0.995, ha="left", va="top", fontsize=13, color=theme["ink"])
    fig.text(0.055, 0.963,
             "Shaded band is ground truth from the generator; the line is what the pipeline measured, "
             "blind to it.",
             ha="left", va="top", fontsize=9, color=theme["ink_secondary"])
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    out = FIGURE_DIR / f"validation-{mode}.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def figure_t4_ranking(mode: str) -> Path:
    ranked = pd.read_csv(T4_RANKED).sort_values("composite_score")
    theme = THEMES[mode]
    style(theme)

    colours = dict(zip(TRIAGE_ORDER, theme["triage"]))
    bar_colours = [colours.get(tag, theme["series"]) for tag in ranked["triage_tag"]]

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    positions = range(len(ranked))
    # height < 1 leaves a surface gap between adjacent bars
    ax.barh(list(positions), ranked["composite_score"], height=0.72,
            color=bar_colours, zorder=3)

    ax.set_yticks(list(positions))
    ax.set_yticklabels(ranked["variant"])
    ax.axvline(0, color=theme["axis"], linewidth=1, zorder=2)

    # RMSD is the physical quantity a reader can interpret; the bar is the
    # batch-relative score. Label both rather than making one stand for the other.
    span = ranked["composite_score"].max() - ranked["composite_score"].min()
    for y, (score, rmsd) in enumerate(zip(ranked["composite_score"], ranked["global_rmsd"])):
        offset = 0.012 * span
        ax.annotate(f"{rmsd:.2f} A",
                    xy=(score + (offset if score >= 0 else -offset), y),
                    ha="left" if score >= 0 else "right", va="center",
                    fontsize=8, color=theme["ink_secondary"])

    ax.set_xlabel("composite impact score  (z-units across the batch; 0 = batch mean)")
    ax.set_axisbelow(True)
    ax.grid(axis="y", visible=False)
    ax.margins(x=0.14)

    present = [tag for tag in TRIAGE_ORDER if tag in set(ranked["triage_tag"])]
    ax.legend(
        handles=[Line2D([], [], marker="s", linestyle="none", markersize=9,
                        color=colours[tag], label=tag) for tag in present],
        title="triage tag (absolute cutoffs)", loc="lower right", frameon=False,
        fontsize=8, title_fontsize=8, labelcolor=theme["ink_secondary"],
    )
    ax.get_legend().get_title().set_color(theme["ink_secondary"])

    # va="top" so the two-line subtitle grows downward and cannot climb into the title.
    fig.suptitle("T4 lysozyme mutants ranked by structural deviation",
                 x=0.055, y=0.985, ha="left", va="top", fontsize=13, color=theme["ink"])
    fig.text(0.055, 0.935,
             f"{len(ranked)} crystal structures against wild-type 2LZM. Bar length is the batch-relative "
             "score; colour is the\nabsolute triage tag; the label is global RMSD. Independent crystals - "
             "packing and resolution contribute.",
             ha="left", va="top", fontsize=8.5, color=theme["ink_secondary"])
    fig.tight_layout(rect=(0, 0, 1, 0.855))

    out = FIGURE_DIR / f"t4-ranking-{mode}.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    missing = [p for p in (GROUND_TRUTH, SYNTH_DETAIL) if not p.exists()]
    if missing:
        print("missing input(s): " + ", ".join(str(p.relative_to(REPO_ROOT)) for p in missing))
        print("run: python tools/make_synthetic.py && python main.py")
        return 1

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    written = []
    for mode in ("light", "dark"):
        written.append(figure_validation(mode))
        if T4_RANKED.exists():
            written.append(figure_t4_ranking(mode))

    if not T4_RANKED.exists():
        print(f"skipped the T4 figure - {T4_RANKED.relative_to(REPO_ROOT)} not found.")
        print("run: python tools/fetch_data.py && python main.py --reference "
              "data/t4_lysozyme/wildtype_reference.pdb --variants data/t4_lysozyme/variants "
              "--output output/t4_lysozyme --trim-outliers")

    for path in written:
        print(f"  wrote {path.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
