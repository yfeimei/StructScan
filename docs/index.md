---
title: StructScan
description: Screening macromolecular variant structures for structural hotspots.
---

<p align="center">
  <a href="https://github.com/yfeimei/StructScan/actions/workflows/tests.yml"><img alt="tests" src="https://img.shields.io/github/actions/workflow/status/yfeimei/StructScan/tests.yml?branch=main&label=tests&style=flat-square"></a>
  <img alt="16 tests" src="https://img.shields.io/badge/tests-16%20passing-1a7f37?style=flat-square">
  <img alt="validation" src="https://img.shields.io/badge/ground%20truth-8%2F8%20recovered-1a7f37?style=flat-square">
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-3776ab?style=flat-square">
  <img alt="licence" src="https://img.shields.io/badge/licence-MIT-blue?style=flat-square">
</p>

<p align="center">
  <b><a href="https://github.com/yfeimei/StructScan">Source on GitHub</a></b> ·
  <b><a href="https://github.com/yfeimei/StructScan/blob/main/GETTING_STARTED.md">Getting started</a></b> ·
  <b><a href="https://github.com/yfeimei/StructScan/blob/main/PROJECT_GUIDE.md">Technical guide</a></b>
</p>

---

## What it is

A protein's shape determines what it does. When a single amino acid is
substituted, the backbone may barely flinch or may distort badly — and which of
those happened is the interesting question.

Labs now generate variant structures in bulk, whether by crystallography or by
prediction. Opening each one in PyMOL or ChimeraX to look for differences works
for five structures and collapses at five hundred. **StructScan automates that
first pass: point it at a directory, and it reports which variants deviate from
the reference, by how much, and at which residues.**

### What it does

1. **Parses** every structure in the directory and extracts backbone coordinates.
2. **Matches residues** between reference and variant by `(chain, residue number,
   insertion code)` — never by list position.
3. **Superimposes** them and computes three geometric metrics per variant:
   global RMSD, per-residue displacement, and backbone torsion change.
4. **Ranks** the batch, assigns each variant a triage tag, and writes four CSV
   reports.

It is a **triage tool, not a functional-impact predictor.** It measures geometry.
Whether a given deviation matters biologically is a separate question this
project does not answer.

---

## Results

A ranking over real structures can only show the pipeline produced *a* result —
not the *right* one. So it is scored against structures whose distortion is known
exactly. `tools/make_synthetic.py` builds an idealised 60-residue backbone, then
derives variants by planting a specific, measured distortion at a specific
residue — a rigid 2.5 Å shift of residues 20–26, a 60° torsion kink at residue
30, and others. The ground truth goes into a file **the pipeline never reads.**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/validation-dark.png">
  <img alt="Three panels of per-residue signal. In each panel the measured curve spikes inside the shaded band marking the residues that were deliberately perturbed." src="figures/validation-light.png">
</picture>

The shaded band is what was planted. The line is what the pipeline measured,
blind to it. The signal lands inside the window every time.

**Eight out of eight variants are triaged correctly.** The four built to be
hotspots are flagged; the four built to be benign are called Stable — including
the two designed to *look* alarming, a structure missing five residues and a
structure with three mutated side chains. The torsion kink is recovered at
**60.02° against a planted 60°**, and the variant that was merely rotated and
translated through space — identical shape, completely different coordinates —
returns **0.001 Å**, because superposition is meant to remove motion that isn't a
shape change.

The whole batch runs from two commands and finishes in about a second, and
**16 tests pass** alongside it. Full ranking, per-residue detail and the raw run
log are in
[the repository](https://github.com/yfeimei/StructScan#results).

### The two results worth reading closely are the imperfect ones

The rigid shifts measure **low** — 2.19 Å for a planted 2.5 Å, 4.03 Å for a
planted 5.0 Å. Least-squares superposition distributes part of any local shift
across the whole fit, so a real displacement is systematically under-reported.
That is a property of superposition-based metrics rather than a bug, and it is
why [the limitations](#what-this-does-not-do) below are not boilerplate.

Ranks 2–4 also show **negative composite scores while still tagged High Priority
Hotspot**, which looks like a contradiction and isn't. The composite score is
*relative* to the batch, and `v05` is extreme enough to pull the mean above them.
The triage tag is *absolute*, so it doesn't care. Both ship because neither
answers the other's question.

---

## On real crystal structures

T4 lysozyme is the standard model system here: the Matthews lab deposited
hundreds of single-point mutants against a common wild-type (`2LZM`). A
sequence-identity search at 90% currently returns **499 entries** — a genuine
mutant batch rather than a contrived one.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/t4-ranking-dark.png">
  <img alt="Horizontal bar chart of 12 T4 lysozyme crystal structures ordered by composite impact score, coloured by triage tag, each labelled with its global RMSD." src="figures/t4-ranking-light.png">
</picture>

Bar length is the batch-relative score; colour is the absolute triage tag. That
the colour blocks do not split cleanly at the zero line is the point — `3LZM` and
`5LZM` fall below the batch mean and still clear the absolute hotspot cutoffs.
The top-ranked entry, `149L`, shows 6.04 Å of displacement concentrated around
residues 36–41.

**A caveat that belongs beside that figure:** these are independent crystals.
Crystal packing, resolution, bound ligands, and crystallisation conditions all
contribute to the measured deviation alongside the mutation itself. A
high-ranking entry is not automatically a mutation effect.

---

## Install and run it locally

**Python 3.10 or newer.** Everything below runs offline.

### 1. Download

```bash
git clone https://github.com/yfeimei/StructScan.git
cd StructScan
```

No git? Download the ZIP from the GitHub page, unpack it, and `cd` into the
folder.

### 2. Set up an environment

Keeps StructScan's dependencies out of your system Python:

```bash
python -m venv .venv
source .venv/bin/activate        # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 3. Run it

```bash
python tools/make_synthetic.py   # generate structures with known perturbations
python main.py                   # run the batch, write output/
python tests/test_pipeline.py    # 16 tests, no pytest required
```

That reproduces the terminal output and the ranking shown above. Results land in
`output/ranked_hotspots.csv`, sorted worst-first — open it in Excel, Numbers,
Google Sheets, or pandas.

### 4. Run it on your own structures

```bash
python main.py --reference /path/to/wildtype.pdb \
               --variants  /path/to/variants/ \
               --output    output/my_run
```

Or download the real T4 lysozyme crystal series shown above:

```bash
python tools/fetch_data.py --limit 60
python main.py --reference data/t4_lysozyme/wildtype_reference.pdb \
               --variants  data/t4_lysozyme/variants \
               --output    output/t4_lysozyme --trim-outliers
```

> **[Getting Started](https://github.com/yfeimei/StructScan/blob/main/GETTING_STARTED.md)**
> is the complete walkthrough — per-shell setup instructions, what every report
> column means, how to trace a flagged variant to its residues, the full option
> table, and troubleshooting.

Structure coordinate data from the RCSB PDB is public domain (CC0).

---

## What the reports contain

| File | What it answers |
|---|---|
| `ranked_hotspots.csv` | *Which variants should I look at first?* One row per variant, sorted by composite score, with triage tag |
| `per_residue_detail.csv` | *Where inside this variant did it move?* Per-residue displacement and torsion delta |
| `residue_hotspot_frequency.csv` | *Which positions move across the whole batch?* Per position: how many variants perturb it |
| `failed_structures.csv` | *What couldn't be processed, and why?* Written only on failure |

One unreadable or mismatched file does not abort the batch — it is logged and the
run continues.

---

## The metrics

| Metric | What it captures |
|---|---|
| Global RMSD over matched α-carbons | Overall deviation after superposition |
| Per-residue displacement | Where the warping is localised |
| Backbone dihedral delta (φ/ψ) | Local torsion change, invisible to displacement |

Torsion angles are periodic, which is a genuine trap. A naive
`|θ_variant − θ_reference|` scores 179° against −179° as **358° apart when they
are 2° apart**. The pipeline takes the shortest angular separation:

```
dθ = min( |a − b| mod 360,  360 − (|a − b| mod 360) )
```

A test pins this behaviour so it cannot regress.

Residues are matched by `(chain, residue number, insertion code)` rather than by
list position, because real PDB entries have missing residues and insertion codes
that positional indexing silently misaligns. Only the intersection is compared,
and coverage is reported.

### Two scores, deliberately independent

- **Composite score** — a weighted sum of z-scores across the batch. Answers
  *how does this variant rank against its peers?* This is the sort key.
- **Triage tag** — absolute Å and degree cutoffs. Answers *is this deviation
  large in physical terms?*

Neither alone is sufficient. Ranking alone flags the top of every batch as
interesting even when nothing moved. Thresholds alone give no ordering within a
tag.

---

## What this does not do

- **Superposition-based RMSD is dominated by rigid-body motion and flexible
  loops.** A buried point mutation can be functionally devastating with
  near-zero RMSD; a floppy terminus can produce large RMSD and mean nothing.
  Superposition-free metrics (lDDT, distance-difference matrices) are more
  robust for this purpose.
- **The composite score is an unvalidated heuristic.** The weights are chosen,
  not fitted, and nothing here has been benchmarked against known-pathogenic or
  known-benign variants. The ranking orders structures by geometric deviation;
  it does not predict functional impact.
- **Backbone only.** α-carbon displacement and φ/ψ say nothing about side-chain
  repacking, which is where much of a point mutation's local effect lands.
- **A single extreme outlier compresses the z-scores** of everything else, since
  the composite is normalised against batch spread. The absolute triage tags are
  unaffected.
- **On predicted structures** (AlphaFold, ESMFold), backbones of point mutants
  are often nearly identical regardless of the mutation's real effect, so a
  geometric screen over predicted mutants may largely track prediction
  confidence rather than biology. Check pLDDT before trusting a ranking.

---

*Built with Biopython, NumPy, pandas, and Matplotlib. MIT licensed.*
