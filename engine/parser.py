"""Module 1: Structure Parser & Coordinate Extraction.

Reads a structural coordinate file and normalises it into a flat, position-keyed
table of residues so that a wild-type reference and any number of variants can be
compared side by side.

Residues are keyed by (chain_id, res_seq, insertion_code) rather than by their
index in the file. Real PDB entries have missing residues, insertion codes and
chain-specific numbering, so positional indexing silently misaligns them.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PPBuilder
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Backbone atoms required for the geometry engine.
BACKBONE_ATOMS: Tuple[str, ...] = ("N", "CA", "C")

# (chain_id, residue sequence number, insertion code)
ResidueKey = Tuple[str, int, str]


class StructureParseError(RuntimeError):
    """Raised when a coordinate file cannot be turned into usable residue data."""


@dataclass
class ResidueRecord:
    """A single residue, normalised for cross-structure comparison."""

    key: ResidueKey
    chain_id: str
    res_seq: int
    icode: str
    res_name: str
    ca_coord: Optional[np.ndarray]
    ca_atom: object  # Bio.PDB.Atom, retained so Superimposer can consume it
    backbone: Dict[str, np.ndarray]
    phi: Optional[float]  # degrees, None at chain starts / breaks
    psi: Optional[float]  # degrees, None at chain ends / breaks

    @property
    def has_ca(self) -> bool:
        return self.ca_coord is not None

    @property
    def label(self) -> str:
        return f"{self.chain_id}:{self.res_seq}{self.icode.strip()}"


@dataclass
class StructureData:
    """All residues extracted from one coordinate file."""

    name: str
    path: Path
    residues: Dict[ResidueKey, ResidueRecord]
    order: List[ResidueKey]
    chains: List[str]

    def __len__(self) -> int:
        return len(self.order)

    def ca_keys(self) -> List[ResidueKey]:
        """Keys of residues that carry an alpha-carbon, in file order."""
        return [k for k in self.order if self.residues[k].has_ca]

    def describe(self) -> str:
        n_ca = len(self.ca_keys())
        chains = ",".join(self.chains) if self.chains else "-"
        return f"{len(self.order)} residues ({n_ca} with CA) across chain(s) {chains}"


def _residue_key(chain_id: str, residue) -> ResidueKey:
    _, res_seq, icode = residue.get_id()
    return (chain_id, int(res_seq), icode)


def _pick_parser(path: Path):
    if path.suffix.lower() in (".cif", ".mmcif"):
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)


def _collect_dihedrals(model) -> Dict[ResidueKey, Tuple[Optional[float], Optional[float]]]:
    """Backbone phi/psi in degrees, keyed by residue.

    PPBuilder splits a chain into separate peptides at chain breaks, so residues
    flanking a gap correctly receive None instead of an angle computed across the
    break.
    """
    builder = PPBuilder()
    angles: Dict[ResidueKey, Tuple[Optional[float], Optional[float]]] = {}
    for chain in model:
        for peptide in builder.build_peptides(chain):
            for residue, (phi, psi) in zip(peptide, peptide.get_phi_psi_list()):
                key = _residue_key(chain.get_id(), residue)
                angles[key] = (
                    None if phi is None else float(np.degrees(phi)),
                    None if psi is None else float(np.degrees(psi)),
                )
    return angles


def load_structure(
    path,
    name: Optional[str] = None,
    model_index: int = 0,
    chains: Optional[List[str]] = None,
) -> StructureData:
    """Parse a .pdb / .cif file into a StructureData record.

    Waters, ligands and other heteroatoms are dropped: only standard polymer
    residues participate in the comparison. For NMR ensembles and other
    multi-model files, only one model is read (the first, by default).
    """
    path = Path(path)
    if not path.is_file():
        raise StructureParseError(f"No such coordinate file: {path}")

    parser = _pick_parser(path)
    try:
        with warnings.catch_warnings():
            # Real PDB entries routinely trip these (discontinuous chains,
            # duplicate atoms); they are informational, not failures.
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = parser.get_structure(name or path.stem, str(path))
    except Exception as exc:  # noqa: BLE001 - surface any parser failure uniformly
        raise StructureParseError(f"Could not parse {path.name}: {exc}") from exc

    models = list(structure)
    if not models:
        raise StructureParseError(f"{path.name} contains no models")
    if model_index >= len(models):
        raise StructureParseError(
            f"{path.name} has {len(models)} model(s); model_index={model_index} is out of range"
        )
    model = models[model_index]

    dihedrals = _collect_dihedrals(model)

    residues: Dict[ResidueKey, ResidueRecord] = {}
    order: List[ResidueKey] = []
    seen_chains: List[str] = []

    for chain in model:
        chain_id = chain.get_id()
        if chains and chain_id not in chains:
            continue
        if chain_id not in seen_chains:
            seen_chains.append(chain_id)

        for residue in chain:
            # residue.id == (hetflag, resseq, icode); a non-blank hetflag means
            # water ('W') or a ligand ('H_XXX').
            if residue.get_id()[0] != " ":
                continue

            key = _residue_key(chain_id, residue)
            if key in residues:
                continue  # defensive: duplicate keys in malformed files

            backbone: Dict[str, np.ndarray] = {}
            for atom_name in BACKBONE_ATOMS:
                if atom_name in residue:
                    # For a disordered atom this returns the highest-occupancy
                    # altloc, which is Biopython's default selection.
                    backbone[atom_name] = np.asarray(
                        residue[atom_name].get_coord(), dtype=float
                    )

            ca_atom = residue["CA"] if "CA" in residue else None
            phi, psi = dihedrals.get(key, (None, None))

            record = ResidueRecord(
                key=key,
                chain_id=chain_id,
                res_seq=key[1],
                icode=key[2],
                res_name=residue.get_resname().strip(),
                ca_coord=backbone.get("CA"),
                ca_atom=ca_atom,
                backbone=backbone,
                phi=phi,
                psi=psi,
            )
            residues[key] = record
            order.append(key)

    if not residues:
        raise StructureParseError(
            f"{path.name} yielded no standard polymer residues "
            "(is it a ligand-only or nucleic-acid entry?)"
        )

    return StructureData(
        name=name or path.stem,
        path=path,
        residues=residues,
        order=order,
        chains=seen_chains,
    )
