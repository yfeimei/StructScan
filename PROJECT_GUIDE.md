# StructScan — Project Guide

Complete reference for the design, data and limitations of StructScan. For a
short overview see [`README.md`](README.md); to install and run it, see
[`GETTING_STARTED.md`](GETTING_STARTED.md). This document is the long form of
*why* it works the way it does.

**Contents**

1. [What StructScan is](#1-what-structscan-is)
2. [Design](#2-design)
3. [Input data](#3-input-data)
4. [Output data](#4-output-data)
5. [Concerns and limitations](#5-concerns-and-limitations)
6. [How to run](#6-how-to-run)
7. [Deploying to GitHub](#7-deploying-to-github)

---

## 1. What StructScan is

A headless batch pipeline that screens a directory of protein variant structures
against a wild-type reference, computes geometric deviation metrics for each, and
emits a ranked triage report.

**The bottleneck it addresses.** When a lab models or predicts dozens to hundreds
of protein variants, inspecting each one in PyMOL or ChimeraX does not scale.
StructScan replaces the first pass of that manual review with an automated screen
that says which variants moved, by how much, and where.

**What it is not.** It is not a functional-impact predictor. It measures geometry.
Whether a geometric deviation implies a functional consequence is a separate
question that this tool does not answer — see [Concerns](#5-concerns-and-limitations).

---

## 2. Design

### 2.1 Data flow

```
[ Input Folder: Wild-type + Batch Variant PDB Files ]
                           |
                           v
 +------------------------------------------------------+
 |             StructScan Analysis Engine                |
 +------------------------------------------------------+
 | 1. Structure Parser & Coordinate Extraction           |
 | 2. Structural Alignment & Pairwise Mapping            |
 | 3. Multi-Metric Geometry Engine (RMSD, Shift, dtheta) |
 | 4. Prioritization Scoring & Anomaly Ranking           |
 +----------------------------+-------------------------+
                              |
                              v
 +------------------------------------------------------+
 |        Generated Output & Prioritization Report       |
 |        (Sorted CSV Sheets & Terminal Progress Log)    |
 +------------------------------------------------------+
```

The reference is parsed once. Each variant is then processed independently:
parse → align → measure → collect. Ranking is a batch-level operation that runs
after every variant has been measured, because the composite score is normalised
against the spread of the batch.

### 2.2 Repository layout

```
StructScan/
├── main.py                          Entry point and CLI
├── requirements.txt                 biopython, pandas, numpy
├── README.md
├── PROJECT_GUIDE.md                 This document
│
├── engine/
│   ├── __init__.py
│   ├── parser.py                    Module 1 — coordinate extraction
│   ├── aligner.py                   Module 2 — superposition, residue mapping
│   ├── geometry.py                  Module 3 — RMSD, displacement, torsion
│   ├── ranker.py                    Module 4 — scoring, triage, aggregation
│   └── reporter.py                  Module 5 — terminal log, CSV export
│
├── tools/
│   ├── make_synthetic.py            Offline fixtures with known ground truth
│   └── fetch_data.py                RCSB PDB downloader
│
├── tests/
│   └── test_pipeline.py             16 tests
│
├── data/                            Input structures
└── output/                          Generated reports
```

### 2.3 Module reference

**`engine/parser.py` — Structure Parser**

Reads `.pdb` / `.cif` via Biopython and normalises each file into a flat table of
residues keyed by `(chain_id, res_seq, insertion_code)`. Extracts backbone N, CA
and C coordinates, residue names, and backbone phi/psi angles (computed with
`PPBuilder`, which correctly returns `None` at chain termini and across chain
breaks rather than inventing an angle).

Waters and ligands are dropped. Only the first model of a multi-model file is
read. For disordered atoms, Biopython's default highest-occupancy altloc is used.

Key types: `ResidueRecord`, `StructureData`. Raises `StructureParseError`.

**`engine/aligner.py` — Alignment Engine**

Matches residues between reference and variant by position key, takes the
intersection, and superimposes the variant onto the reference over matched
alpha-carbons using Biopython's `Superimposer`.

Optional iterative outlier trimming (`--trim-outliers`) refits while excluding
atoms deviating beyond a cutoff, so a single mobile loop cannot distort the whole
superposition. Trimmed atoms are excluded from the *fit*, not from the *report*.

Returns `AlignmentResult` carrying the transform, matched keys, the fitted mask,
coverage, and the list of residue substitutions. Raises `AlignmentError` when
fewer than 3 alpha-carbons match or coverage falls below `--min-coverage`.

**`engine/geometry.py` — Geometry & Metric Engine**

Computes the per-residue metric table and the per-variant summary fingerprint.
Pure functions `rmsd()`, `displacement()` and `angular_delta()` are separately
testable.

**`engine/ranker.py` — Ranking & Triage Sorter**

Builds the composite score, assigns triage tags, sorts the batch, and aggregates
per-residue results into the batch-level hotspot table.

**`engine/reporter.py` — Report Generator**

Terminal progress logging and CSV export. Holds the column ordering for the
headline report so the file opens legibly in Excel.

### 2.4 The metrics

**Global RMSD** over matched alpha-carbons after superposition:

```
RMSD = sqrt( (1/N) * sum_i || r_i - r_i' ||^2 )
```

**Local displacement** — per-residue Euclidean shift, to localise warping:

```
d = sqrt( (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 )
```

**Dihedral delta** — change in backbone phi/psi. Torsion angles are periodic, so
this is the *shortest* angular separation, bounded to `[0, 180]`:

```
raw    = |theta_mut - theta_ref| mod 360
dtheta = min(raw, 360 - raw)
```

A plain `|theta_mut - theta_ref|` scores 179° against −179° as 358° apart when
they are 2° apart. This is pinned by a test.

### 2.5 Scoring and triage

Two independent mechanisms, both present in the report:

| | Composite score | Triage tag |
|---|---|---|
| Question | How does this variant rank *against its peers*? | Is the deviation large *in physical terms*? |
| Basis | Weighted sum of z-scores across the batch | Absolute Å / degree cutoffs |
| Role | Sort key | Filter |

Ranking alone would flag the top of every batch as interesting even when nothing
moved. Thresholds alone would give no ordering within a tag. Both are needed.

Default weights — RMSD 1.0, max displacement 1.0, max dihedral 0.5. Default
hotspot cutoffs — RMSD > 1.0 Å, displacement > 2.0 Å, dihedral > 45°. A variant is
tagged at the highest level for which it trips *any* criterion.

### 2.6 Design decisions worth knowing

**Residues are matched by key, not by list index.** Real PDB entries have missing
residues, insertion codes and chain-specific numbering. Positional indexing
silently misaligns them and produces confidently wrong numbers. Only the
intersection is compared, and coverage is reported so partial matches are visible.

**Residue identity is ignored during matching.** A point mutant differs in residue
name at the mutated site by definition. Substitutions are recorded in the report,
not treated as a mismatch.

**Alpha-carbons only for the superposition.** Side-chain atoms differ between
wild-type and mutant residues, so they cannot be paired.

**One bad file does not abort the batch.** Parse and alignment failures are logged
to `failed_structures.csv` with a reason, and the run continues.

**Trimming is off by default.** The default behaviour is the plain global fit
described in the original design. Trimming is available but opt-in, so results are
reproducible against the specification unless explicitly changed.

---

## 3. Input data

### 3.1 Format requirements

- `.pdb`, `.ent`, `.cif` or `.mmcif`
- One reference file, plus a directory of variant files
- Variants must share chain IDs and residue numbering with the reference —
  matching is by `(chain, resseq, icode)`
- At least 3 shared alpha-carbons, and by default at least 50% coverage of the
  reference (`--min-coverage`)

Missing residues, insertion codes, waters, ligands, multiple chains and
alternate conformations are all handled. Renumbered chains are not — if a variant
uses different residue numbering than the reference, matching will fail or
mismatch, and the coverage figure in the report is the signal to check.

### 3.2 Dataset 1 — synthetic (default)

```
data/
├── wildtype_reference.pdb      14 KB   idealised 60-residue N-CA-C backbone
├── ground_truth.json          2.1 KB   what was planted in each variant
└── variants/                  113 KB   8 files
```

| File | Planted perturbation | Expected tag |
|---|---|---|
| `v01_stable.pdb` | coordinate noise, sigma 0.03 Å | Stable |
| `v02_loop_shift_20_26.pdb` | 2.5 Å rigid shift, residues 20–26 | Hotspot |
| `v03_core_shift_40_48.pdb` | 5.0 Å rigid shift, residues 40–48 | Hotspot |
| `v04_terminal_flex.pdb` | 4.0 Å shift, last 6 residues | Hotspot |
| `v05_dihedral_kink_30.pdb` | 60° torsion kink at residue 30 | Hotspot |
| `v06_rigid_body.pdb` | rotated + translated, otherwise identical | Stable |
| `v07_gapped.pdb` | residues 10–14 deleted | Stable |
| `v08_substitutions.pdb` | 3 residues renamed to TRP | Stable |

**Provenance.** Generated locally by `tools/make_synthetic.py`. Not experimental
structures — an idealised backbone built from ideal bond geometry (N-CA 1.458 Å,
CA-C 1.525 Å, C-N 1.329 Å) and specified phi/psi, using NeRF placement.
Deterministic under seed `20240517`, so regeneration is byte-identical. Rebuilds
in 0.2 seconds, offline.

**Why it exists.** Real structures can confirm the pipeline produced *a* ranking.
Only planted perturbations confirm it produced the *right* one. `tests/` asserts
recovery of each perturbation at the correct residue; `data/ground_truth.json`
records what was planted.

Note on `v05`: the torsion change is at residue 30, but the largest *displacement*
appears near residue 60. That is correct — a backbone kink propagates downstream
like a lever arm. The 60° shows up in the dihedral column at residue 30.

### 3.3 Dataset 2 — real crystal structures

```
data/t4_lysozyme/
├── wildtype_reference.pdb     146 KB   PDB entry 2LZM, T4 lysozyme wild-type
└── variants/                  2.5 MB   149L 150L 1DYC 1L17 1LYD 1T6H
                                        3LZM 4LZM 4S0W 5LZM 6LZM 7LZM
```

**Provenance.** Downloaded unmodified from the RCSB PDB by `tools/fetch_data.py`.
T4 lysozyme was chosen because the Matthews lab solved an unusually large series
of single-point mutants against a common wild-type — a sequence-identity search at
90% currently returns **499 entries**, which is a genuine mutant batch rather than
a contrived one. PDB coordinate data is public domain (CC0).

The 12 files present are a capped sample (`--limit 12`). To pull more:

```bash
python tools/fetch_data.py --limit 200      # ~35 MB, a few minutes
```

Already-downloaded files are skipped, so re-running resumes rather than refetching.

**Other usable series:** `--reference 1STN` (staphylococcal nuclease),
`--reference 1HHP` (HIV-1 protease drug-resistance mutants).

### 3.4 Using your own data

Nothing in the pipeline is specific to these datasets:

```bash
python main.py --reference /path/to/wt.pdb --variants /path/to/variants/
```

If the variants come from structure prediction, read the note on predicted
structures in [Concerns](#51-scientific).

---

## 4. Output data

Written to `output/` by default (`--output`).

### 4.1 `ranked_hotspots.csv` — the headline report

One row per variant, sorted by composite score, highest deviation first. 24 columns.

| Column | Meaning |
|---|---|
| `rank` | 1 = most deviant in this batch |
| `variant` | Variant name (filename stem) |
| `triage_tag` | `High Priority Hotspot` / `Moderate Deviation` / `Stable` |
| `composite_score` | Weighted sum of z-scores; the sort key |
| `global_rmsd` | Å, over atoms used in the fit |
| `max_ca_displacement` | Å, largest single-residue shift |
| `mean_ca_displacement` | Å, averaged over matched residues |
| `max_dihedral_delta` | degrees, largest backbone torsion change |
| `n_residues_displaced` | count above `--displacement-threshold` |
| `n_residues_torsioned` | count above `--dihedral-threshold` |
| `top_deviating_residues` | e.g. `A:40(4.03A); A:43(4.01A); ...` |
| `n_substitutions` | residues whose identity differs from reference |
| `substitutions` | e.g. `A:ALA15TRP;A:GLY33TRP` |
| `n_residues_matched` | residues compared |
| `coverage` | matched / reference residues; low values mean check the input |
| `source_file` | original filename |
| `global_rmsd_all_matched` | RMSD over every matched atom, including trimmed ones |
| `n_residues_reference` / `n_residues_variant` | alpha-carbon counts |
| `n_residues_fitted` | atoms used in the superposition after trimming |
| `mean_dihedral_delta` | degrees |
| `z_global_rmsd`, `z_max_ca_displacement`, `z_max_dihedral_delta` | the z-scores feeding the composite |

**Reading it.** Sort order answers *what to look at first*. `triage_tag` answers
*is this worth looking at at all*. If `coverage` is well below 1.0, treat the row
with suspicion — the structures may not correspond as well as assumed.

### 4.2 `per_residue_detail.csv` — the localisation table

One row per (variant, residue). 13 columns: `variant`, `chain`, `res_seq`,
`icode`, `position`, `ref_resname`, `var_resname`, `is_substitution`,
`ca_displacement`, `delta_phi`, `delta_psi`, `max_dihedral_delta`, `used_in_fit`.

This is where you find *where* in the structure a flagged variant moved.
`used_in_fit` is `False` for atoms excluded by `--trim-outliers` — still measured,
just not part of the superposition.

### 4.3 `residue_hotspot_frequency.csv` — the batch-level view

One row per position, aggregated across all variants. 13 columns, notably
`n_variants`, `n_variants_displaced`, `fraction_displaced`, `max_ca_displacement`,
`mean_ca_displacement`, `n_substituted`.

A residue that moves in one variant is that variant's problem. A residue that
moves across a large fraction of the batch is a structural hotspot — which is
what the pipeline is ultimately looking for. Sorted by `fraction_displaced`.

Note that `n_variants` differs per position: a residue absent from some structures
(truncated termini, disordered loops) is only counted against variants where it
exists, so `3/6` and `3/12` can both appear in one report.

### 4.4 `failed_structures.csv`

Written only when something failed. Columns `file`, `reason`. Absence of this file
means every input was analysed.

### 4.5 Terminal output

A per-file progress line during the run, then a triage summary, the top variants
by composite score, and the most frequently perturbed positions. Suppress with
`--quiet`.

---

## 5. Concerns and limitations

This section is deliberately blunt. Every item here is something a reviewer could
raise; better to have the answer written down than discovered mid-conversation.

### 5.1 Scientific

**Superposition RMSD is a weak proxy for variant impact.** It is dominated by
rigid-body motion and flexible loops or termini. A buried point mutation can be
functionally devastating with near-zero RMSD; a floppy tail can produce large RMSD
and mean nothing. Superposition-free metrics — lDDT, distance-difference matrices,
TM-score — are more robust for this purpose. This is the single largest
limitation of the approach.

**The composite score is an unvalidated heuristic.** Weights are chosen, not
fitted. Thresholds are chosen, not derived from a null distribution. Nothing has
been benchmarked against known-pathogenic or known-benign variants. The ranking
orders structures by geometric deviation; it does not predict functional impact,
and no claim about its accuracy can be supported from what is in this repository.

**Backbone only.** Alpha-carbon displacement and phi/psi say nothing about
side-chain repacking, which is where much of a point mutation's local effect
actually lands. Burial and solvent accessibility — strong predictors of mutational
tolerance — are not considered at all.

**Crystal-packing confound in real data.** The T4 lysozyme structures are
independent crystals. Packing contacts, resolution, bound ligands and
crystallisation conditions all contribute to measured deviation alongside the
mutation. A high-ranking entry is not automatically a mutation effect.

**Predicted structures behave differently.** On AlphaFold/ESMFold-type models,
backbones of point mutants are frequently near-identical regardless of the
mutation's real effect. A geometric screen over predicted mutants may largely
track prediction confidence (pLDDT) and disorder rather than biology. Check
confidence scores before trusting any such ranking.

**Thresholds need retuning per dataset.** The defaults are calibrated for the
synthetic set. On the real T4 batch they produced 6 hotspot / 6 moderate / 0
stable — the 45° dihedral cutoff is tripped by ordinary crystal-to-crystal
backbone variation. Raise `--dihedral-hotspot` or set `--weight-dihedral 0` for
real crystal data.

**A single extreme outlier squashes the batch.** Because the composite is
normalised against batch spread, one very large deviation pushes every other
z-score negative. In the synthetic run `v05` scored 2.61 and everything else went
negative. The absolute triage tags are unaffected, which is why both mechanisms
are reported.

**Single-reference assumption.** Everything is measured against one wild-type. If
the reference itself is a poor or unrepresentative structure, every number
inherits that.

### 5.2 Implementation

- **First model only.** NMR ensembles and other multi-model files are read at
  model 0. Conformational spread within an ensemble is invisible.
- **Highest-occupancy altloc only.** Alternate conformations are not compared.
- **Heteroatoms are dropped.** Ligand-induced or cofactor-dependent differences
  are invisible by construction.
- **Chain IDs must correspond.** Renamed or reordered chains break matching.
  Use `--chains` to restrict the comparison when they do not.
- **PDB precision floor.** PDB stores coordinates to 3 decimal places, so a round
  trip through the format costs up to ~0.001 Å. Deviations at that scale are
  format noise, not signal — the rigid-body test asserts against this floor.
- **`icode` round-trips as a float.** Empty insertion codes become `NaN` when the
  CSV is re-read by pandas. Cosmetic; the `position` column carries the full
  label.
- **Single-threaded.** 12 real structures take 1.3 s, so this is not currently a
  constraint. At tens of thousands of files it would be.
- **All per-residue frames are held in memory** before concatenation. Fine for
  hundreds of small proteins; would need streaming for very large batches.

### 5.3 Data and licensing

PDB coordinate data is released public domain (CC0) and is safe to redistribute.
Synthetic data is generated by this repository and carries its licence. Any
unpublished laboratory data requires the PI's permission before it is committed
to a public repository — see [Deploying to GitHub](#7-deploying-to-github).

### 5.4 Authorship and claims

This code was written with AI assistance. That is unremarkable and increasingly
standard, but it places a real obligation on whoever presents the work: be able to
explain every design decision in it, particularly the choice of metrics, the
residue-matching scheme, and the validation strategy. The tests and this document
exist partly to make that possible.

Be conservative about biological claims. In the T4 run, positions 19–23 recur
across variants and sit near a known hinge region, and residue 164 is the
C-terminus. The correct statement is "these positions recur in this batch"; the
incorrect one is "these are functional hotspots". The tool cannot distinguish
those, and saying so is a stronger position than overclaiming.

---

## 6. How to run

Operating instructions live in one place: **[`GETTING_STARTED.md`](GETTING_STARTED.md)**.
It covers installation and virtual environments, running the pipeline, reading
each output file, running against your own structures, the complete flag table,
and troubleshooting.

The short version, for readers who already have the dependencies installed:

```bash
python tools/make_synthetic.py       # 0.2 s — generates data/
python main.py                       # 1.1 s — writes output/
python tests/test_pipeline.py        # 16 tests, ~1.2 s, no pytest required
```

Exit codes: `0` success, `1` no variants analysable, `2` bad input paths or
network failure.

---

## 7. Deploying to GitHub

### 7.1 What to commit

**Commit:** all source (`engine/`, `tools/`, `tests/`, `main.py`),
`requirements.txt`, both markdown documents, a licence, `.gitignore`, and the CI
workflow.

**Commit the synthetic data?** Optional — it is 130 KB and regenerates
deterministically in 0.2 s. Committing it means the repo runs immediately after
clone with no generation step, which is friendlier for a reader who just wants to
look. Either choice is defensible; the `.gitignore` below excludes it, so delete
those two lines if you would rather include it.

**Do not commit:** `output/` (regenerated on every run), `data/t4_lysozyme/`
(2.5 MB of files that `fetch_data.py` re-downloads on demand), `__pycache__/`,
virtual environments, and — most importantly — **any unpublished laboratory data**
without the PI's explicit permission. Once pushed to a public repository, content
can be cached and indexed even if deleted afterwards.

### 7.2 `.gitignore`

```gitignore
# Generated reports
output/

# Downloaded structures — regenerate with tools/fetch_data.py
data/t4_lysozyme/
data/*/variants/*.pdb
data/*/wildtype_reference.pdb

# Synthetic fixtures — regenerate with tools/make_synthetic.py
# (delete these three lines to ship the demo data with the repo)
data/variants/
data/wildtype_reference.pdb
data/ground_truth.json

# Python
__pycache__/
*.py[cod]
.pytest_cache/
.venv/
venv/
*.egg-info/

# Editors / OS
.vscode/
.idea/
.DS_Store
Thumbs.db
*.swp
```

### 7.3 Licence

MIT is the conventional choice for a small research utility and imposes no
constraints on a reader. Save as `LICENSE`:

```
MIT License

Copyright (c) 2026 <YOUR NAME>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### 7.4 Continuous integration

A green "tests passing" badge is cheap and does real work — it shows a reader the
project is verified without them running anything. Save as
`.github/workflows/tests.yml`:

```yaml
name: tests

on:
  push:
    branches: [main]
  pull_request:

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.10", "3.12"]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - run: pip install -r requirements.txt
      - run: python tests/test_pipeline.py
```

The tests generate their own fixtures in a temp directory and need no network, so
CI passes on a bare checkout regardless of whether data is committed.

Add the badge to the top of `README.md`:

```markdown
![tests](https://github.com/<USER>/structscan/actions/workflows/tests.yml/badge.svg)
```

### 7.5 Publishing

```bash
cd StructScan
git init -b main
git add .
git status                  # confirm no data/ or output/ before committing
git commit -m "StructScan: batch structural variant triage pipeline"
```

Then either with the GitHub CLI:

```bash
gh repo create structscan --public --source=. --push \
  --description "Automated batch pipeline for screening protein variant structures and ranking structural hotspots"
```

or create an empty repository on github.com and:

```bash
git remote add origin https://github.com/<USER>/structscan.git
git push -u origin main
```

### 7.6 Repository polish

Small things that disproportionately affect how the repo reads:

- **Description and topics.** Set the description; add topics
  `bioinformatics`, `structural-biology`, `protein-structure`, `pdb`, `biopython`.
- **Paste real output into the README.** A reader who will not clone anything
  still sees it work. The ground-truth comparison table is the strongest single
  artefact — planted perturbation beside recovered perturbation.
- **Commit history.** A handful of meaningful commits reads better than one
  `initial commit` containing everything. If the work genuinely happened in
  stages, let the history show that.
- **Do not inflate the claims** in the description. "Screens variant structures
  and ranks them by geometric deviation" is accurate and sufficient. "Predicts
  mutation impact" is not supported by anything in the repository.
