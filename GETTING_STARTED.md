# Getting Started with StructScan

Everything needed to go from a fresh download to reading a finished report. No
prior knowledge of the codebase assumed.

If you only want to know *what* StructScan is and *why* it works the way it does,
read [`README.md`](README.md) or the [project page](docs/index.md) instead. For
the design rationale, metric derivations and scientific limitations, see
[`PROJECT_GUIDE.md`](PROJECT_GUIDE.md).

**Contents**

1. [Before you start](#1-before-you-start)
2. [Setup](#2-setup)
3. [Run it](#3-run-it)
4. [Read the report](#4-read-the-report)
5. [Run it on your own structures](#5-run-it-on-your-own-structures)
6. [All options](#6-all-options)
7. [Troubleshooting](#7-troubleshooting)

---

## 1. Before you start

**Python 3.10 or newer.** Check:

```bash
python --version
```

If that prints 3.9 or lower, or `command not found`, install a current Python
from [python.org](https://www.python.org/downloads/). On Windows, tick *"Add
Python to PATH"* in the installer. On macOS and Linux `python` may need to be
`python3` — if so, use `python3` everywhere below.

**Roughly 5 minutes and 200 MB** of disk for the dependencies. Everything in
sections 2–4 runs offline; only section 5's crystal-structure download needs a
network connection.

---

## 2. Setup

### 2.1 Get the code

```bash
git clone https://github.com/yfeimei/StructScan.git
cd StructScan
```

No git? Download the ZIP from the GitHub page, unpack it, and `cd` into the
unpacked folder.

### 2.2 Create a virtual environment

This keeps StructScan's dependencies out of your system Python. Skipping it is a
common cause of `externally-managed-environment` errors on macOS and Linux.

```bash
python -m venv .venv
```

Then activate it:

| Shell | Command |
|---|---|
| macOS / Linux (bash, zsh) | `source .venv/bin/activate` |
| Windows PowerShell | `.venv\Scripts\Activate.ps1` |
| Windows cmd.exe | `.venv\Scripts\activate.bat` |
| Windows Git Bash | `source .venv/Scripts/activate` |

Your prompt should now be prefixed with `(.venv)`. That prefix is how you know
the environment is active — you need to re-activate it in every new terminal.

> **Windows PowerShell:** if activation fails with *"running scripts is disabled
> on this system"*, run
> `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` once, then try again.

### 2.3 Install the dependencies

```bash
pip install -r requirements.txt
```

Installs Biopython, pandas and NumPy for the pipeline, plus Matplotlib for the
figure script. Takes a minute or two.

### 2.4 Confirm it works

```bash
python main.py --version
```

Prints the version. If instead you get `ModuleNotFoundError`, the virtual
environment is not active — go back to 2.2.

---

## 3. Run it

Two commands. The first generates test structures; the second analyses them.

```bash
python tools/make_synthetic.py    # ~0.2 s — writes data/
python main.py                    # ~1 s   — writes output/
```

`make_synthetic.py` builds an idealised 60-residue protein backbone and derives
eight variants from it, each with a **known, deliberately planted distortion** —
a 2.5 Å rigid shift of residues 20–26, a 60° torsion kink at residue 30, and so
on. Because the answer is known in advance, you can check the pipeline against
it. (The data ships with the repo, so this step is only strictly necessary if
you have deleted or modified `data/`.)

### What you should see

```
==============================================================================
StructScan - batch structural variant triage
==============================================================================
  Reference : wildtype_reference  [60 residues (60 with CA) across chain(s) A]
  Variants  : 8 file(s)
------------------------------------------------------------------------------
  [  1/8] v01_stable                   RMSD  0.051 A   max shift  0.101 A   max dtheta    8.7 deg   matched 60
  [  2/8] v02_loop_shift_20_26         RMSD  0.796 A   max shift  2.189 A   max dtheta    0.0 deg   matched 60
  [  3/8] v03_core_shift_40_48         RMSD  1.722 A   max shift  4.027 A   max dtheta    0.0 deg   matched 60
  [  4/8] v04_terminal_flex            RMSD  0.931 A   max shift  2.611 A   max dtheta    0.0 deg   matched 60
  [  5/8] v05_dihedral_kink_30         RMSD 13.719 A   max shift 27.197 A   max dtheta   60.0 deg   matched 60
  [  6/8] v06_rigid_body               RMSD  0.001 A   max shift  0.001 A   max dtheta    0.1 deg   matched 60
  [  7/8] v07_gapped                   RMSD  0.050 A   max shift  0.111 A   max dtheta   11.2 deg   matched 55
  [  8/8] v08_substitutions            RMSD  0.054 A   max shift  0.104 A   max dtheta    9.8 deg   matched 60
------------------------------------------------------------------------------
  Triage: 4 hotspot, 0 moderate, 4 stable  (of 8 analysed)

  Top 8 by composite impact score:
      #  variant                        score    RMSD  tag
      1  v05_dihedral_kink_30           2.607  13.719  High Priority Hotspot
      2  v03_core_shift_40_48          -0.182   1.722  High Priority Hotspot
      3  v04_terminal_flex             -0.319   0.931  High Priority Hotspot
      4  v02_loop_shift_20_26          -0.351   0.796  High Priority Hotspot
      5  v07_gapped                    -0.397   0.050  Stable
      6  v08_substitutions             -0.411   0.054  Stable
      7  v01_stable                    -0.423   0.051  Stable
      8  v06_rigid_body                -0.523   0.001  Stable

  Most frequently perturbed positions across the batch:
    A:60       ASP  displaced in   3/8 variants   max 27.20 A
    ...

  wrote ranked       -> .../output/ranked_hotspots.csv
  wrote per_residue  -> .../output/per_residue_detail.csv
  wrote hotspots     -> .../output/residue_hotspot_frequency.csv
==============================================================================
```

Three sanity checks you can make on that output immediately:

- **`v06_rigid_body` reports RMSD 0.001 Å.** That variant is the reference
  rotated and translated in space — same shape, different coordinates.
  Superposition is supposed to remove motion that isn't a shape change, and it
  does.
- **`v05_dihedral_kink_30` reports `max dtheta 60.0 deg`** against a planted
  60°. The torsion is recovered almost exactly.
- **`v07_gapped` matched 55 of 60 residues**, because five were deleted from it.
  The pipeline compares the intersection instead of failing.

### Verify against the tests

```bash
python tests/test_pipeline.py     # standalone, no pytest needed
python -m pytest tests -q         # if you have pytest
```

16 tests, ~1.2 s. They assert that each planted perturbation is recovered at the
correct residue, along with rigid-body invariance, the periodic dihedral
boundary, gap tolerance and substitution handling.

### Redraw the figures (optional)

```bash
python tools/make_figures.py      # regenerates docs/figures/ from output/
```

---

## 4. Read the report

The run writes four CSV files to `output/`. Open them in Excel, Numbers, Google
Sheets, or `pandas` — they are plain comma-separated text with a header row.

| File | What it answers | Written |
|---|---|---|
| **`ranked_hotspots.csv`** | *Which variants should I look at first?* | always |
| `per_residue_detail.csv` | *Where inside this variant did it move?* | always |
| `residue_hotspot_frequency.csv` | *Which positions move across the whole batch?* | always |
| `failed_structures.csv` | *What couldn't be processed, and why?* | only on failure |

**Start with `ranked_hotspots.csv`.** One row per variant, already sorted worst
first.

### The columns that matter

| Column | Read it as |
|---|---|
| `rank` | 1 = most deviant **in this batch** |
| `variant` | The variant's filename, minus the extension |
| `triage_tag` | `High Priority Hotspot` / `Moderate Deviation` / `Stable` |
| `composite_score` | The sort key — weighted sum of z-scores |
| `global_rmsd` | Å. Overall backbone deviation after superposition |
| `max_ca_displacement` | Å. The single worst-moving residue |
| `max_dihedral_delta` | Degrees. The single worst torsion change |
| `top_deviating_residues` | e.g. `A:40(4.03A); A:43(4.01A)` — go look at these |
| `substitutions` | e.g. `A:ALA15TRP` — which residues actually changed identity |
| `coverage` | Matched ÷ reference residues. **Below ~0.9, distrust the row** |

There are 24 columns in total; the rest are intermediate values (per-metric
z-scores, pre-trim RMSD, residue counts). [`PROJECT_GUIDE.md`
§4](PROJECT_GUIDE.md#4-output-data) documents every one.

### Why there are two independent scores

This trips people up, so it is worth being explicit:

- **`composite_score` is relative.** It asks *how does this variant compare to
  its peers in this batch?* It is a z-score, so it is negative for anything
  below the batch mean — including rows that are genuinely bad in absolute
  terms.
- **`triage_tag` is absolute.** It asks *is this deviation physically large?*,
  using fixed Å and degree cutoffs that don't care about the rest of the batch.

In the run above, ranks 2–4 have **negative composite scores and still read
`High Priority Hotspot`** — they sit below the batch mean only because
`v05_dihedral_kink_30` is so extreme it drags the mean up. Both columns are
correct; they answer different questions. Neither alone is sufficient: ranking
alone flags the top of every batch even when nothing moved, and thresholds alone
give no ordering within a tag.

### Following a hotspot to its location

Once `ranked_hotspots.csv` tells you *which* variant, `per_residue_detail.csv`
tells you *where*. Filter it to that variant and sort by `ca_displacement`:

```bash
python -c "import pandas as pd; d=pd.read_csv('output/per_residue_detail.csv'); print(d[d.variant=='v03_core_shift_40_48'].nlargest(10,'ca_displacement').to_string(index=False))"
```

For `v03_core_shift_40_48` the top rows are residues 40–48 — exactly the window
the generator perturbed.

Useful columns there: `ca_displacement` (Å), `delta_phi` / `delta_psi`
(degrees), `is_substitution`, and `used_in_fit` (`False` for atoms excluded by
`--trim-outliers` — still measured, just not part of the superposition).

### The batch-level view

`residue_hotspot_frequency.csv` aggregates across every variant, sorted by
`fraction_displaced`. A residue that moves in one variant is that variant's
problem; a residue that moves across a large fraction of the batch is a
structural hotspot, which is what the pipeline is ultimately looking for.

Note `n_variants` differs per position — a residue absent from some structures
is only counted against variants where it exists, so `3/6` and `3/12` can both
appear in one report.

### One honest caveat about the numbers

The rigid shifts read **low**: 2.19 Å measured for a planted 2.5 Å, 4.03 Å for a
planted 5.0 Å. Least-squares superposition distributes part of any local shift
across the whole fit, so a real displacement is systematically under-reported.
That is a property of superposition-based metrics rather than a bug — but it
means the measured magnitude is a lower bound, not an exact figure. See
[`PROJECT_GUIDE.md` §5](PROJECT_GUIDE.md#5-concerns-and-limitations) for the
full list of what this tool cannot tell you.

---

## 5. Run it on your own structures

### Your own files

```bash
python main.py --reference /path/to/wildtype.pdb \
               --variants  /path/to/variants/ \
               --output    output/my_run
```

Requirements for input:

- Format `.pdb`, `.ent`, `.cif` or `.mmcif`
- One reference file, plus a **directory** of variant files
- Variants must share chain IDs and residue numbering with the reference —
  matching is by `(chain, residue number, insertion code)`, never by position
- At least 3 shared α-carbons and, by default, 50% coverage (`--min-coverage`)

Missing residues, insertion codes, waters, ligands, multiple chains and
alternate conformations are all handled. **Renumbered chains are not** — if a
variant numbers its residues differently from the reference, matching will fail
or silently mismatch, and the `coverage` column is your signal to check.

### Real crystal structures from the PDB

Downloads the T4 lysozyme mutant series — the Matthews lab deposited hundreds of
single-point mutants against a common wild-type (`2LZM`):

```bash
python tools/fetch_data.py --limit 60

python main.py --reference data/t4_lysozyme/wildtype_reference.pdb \
               --variants  data/t4_lysozyme/variants \
               --output    output/t4_lysozyme \
               --trim-outliers --dihedral-hotspot 90
```

`--trim-outliers` is strongly recommended on real crystal data — without it,
flexible termini dominate the superposition and swamp everything else.

Other usable series: `--reference 1STN` (staphylococcal nuclease), `--reference
1HHP` (HIV-1 protease drug-resistance mutants).

**Caveat:** these are independent crystals. Packing, resolution, bound ligands
and crystallisation conditions all contribute to the measured deviation
alongside the mutation. A high-ranking entry is not automatically a mutation
effect.

---

## 6. All options

### `main.py`

| Flag | Default | Purpose |
|---|---|---|
| `--reference` | `data/wildtype_reference.pdb` | Wild-type reference structure |
| `--variants` | `data/variants` | Directory of variant structures |
| `--output` | `output` | Directory for the CSV reports |
| `--chains` | all | Restrict analysis to these chain IDs |
| `--top-residues` | 5 | Worst-deviating residues named per variant |
| `--trim-outliers` | off | Iterative refit excluding high-deviation atoms |
| `--trim-cutoff` | 2.0 | Å cutoff for trimming |
| `--min-coverage` | 0.5 | Minimum fraction of reference residues matched |
| `--displacement-threshold` | 1.0 | Å; counts a residue as displaced |
| `--dihedral-threshold` | 30.0 | Degrees; counts a residue as torsioned |
| `--weight-rmsd` | 1.0 | Composite score weight |
| `--weight-displacement` | 1.0 | Composite score weight |
| `--weight-dihedral` | 0.5 | Composite score weight |
| `--rmsd-hotspot` | 1.0 | Å; triage tag cutoff |
| `--displacement-hotspot` | 2.0 | Å; triage tag cutoff |
| `--dihedral-hotspot` | 45.0 | Degrees; triage tag cutoff |
| `--quiet` | off | Suppress the progress log |
| `--version` | — | Print version and exit |

Exit codes: `0` success, `1` no variants analysable, `2` bad input paths or
network failure.

### Helper scripts

| Script | Options |
|---|---|
| `tools/make_synthetic.py` | `--outdir`, `--length`, `--seed` |
| `tools/fetch_data.py` | `--reference` (PDB ID, default `2LZM`), `--entity`, `--outdir`, `--limit`, `--identity` |
| `tools/make_figures.py` | redraws `docs/figures/` from `output/` |
| `tools/preview_page.py` | `--port`; previews `docs/index.md` locally (needs `pip install grip`) |

---

## 7. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `ModuleNotFoundError: No module named 'Bio'` | Virtual environment not active, or dependencies not installed. Re-run §2.2 and §2.3 |
| `python: command not found` | Try `python3` instead, or install Python — see §1 |
| `running scripts is disabled on this system` (PowerShell) | `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`, then re-activate |
| `externally-managed-environment` from pip | You skipped the virtual environment — see §2.2 |
| `error: No such coordinate file` | Run `python tools/make_synthetic.py` first |
| `error: no variants could be analysed` | Check `output/failed_structures.csv` for the per-file reason |
| Everything tagged `Stable` | Thresholds too high for this batch; lower `--rmsd-hotspot` and friends |
| Nothing tagged `Stable` | Thresholds too low — common on real crystal data; raise `--dihedral-hotspot` |
| Variants skipped with a coverage error | Residue numbering differs from the reference, or it is a different protein |
| Huge RMSD across the whole batch | A flexible terminus is dominating the fit; add `--trim-outliers` |
| `could not reach the RCSB PDB` | Network or firewall. The synthetic path in §3 works entirely offline |
| Figures look stale | `tools/make_figures.py` reads `output/` — re-run `main.py` first |
