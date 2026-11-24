"""Redocking: does the setup reproduce the crystal pose?

Before trusting a single docking score, take the ligand that came with the
structure, dock it back into its own site, and measure how far the top pose is
from the crystal one. Under 2.0 A is the convention for "the setup works";
above it, fix the setup before screening anything.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

from .box import Box
from .dock import DockingResult, dock
from .errors import LigandError
from .ligand import ligand_from_pdb_block, prepare_ligand
from .results import Pose, poses_to_rdkit
from .structure import Residue

#: The usual success threshold in the docking literature.
RMSD_SUCCESS = 2.0


@dataclass
class RedockResult:
    """A redocking check: the poses, their RMSD to the crystal pose, and a verdict."""

    ligand: str
    docking: DockingResult
    rmsd_per_mode: list[float]
    best_rmsd: float
    top_pose_rmsd: float
    top_affinity: float
    threshold: float = RMSD_SUCCESS
    symmetry_corrected: bool = True
    notes: list[str] = field(default_factory=list)

    @property
    def success(self) -> bool:
        """True when Vina's top-ranked pose reproduces the crystal pose."""
        return self.top_pose_rmsd <= self.threshold

    @property
    def found_in_top_poses(self) -> bool:
        """True when any returned pose is close to the crystal pose."""
        return self.best_rmsd <= self.threshold

    def to_dataframe(self):
        import pandas as pd

        frame = self.docking.to_dataframe().copy()
        frame["rmsd_to_crystal"] = self.rmsd_per_mode
        return frame

    def summary(self) -> str:
        verdict = "PASS" if self.success else "CHECK THE SETUP"
        lines = [
            f"Redocking {self.ligand}: {verdict}",
            f"  top pose affinity : {self.top_affinity:.2f} kcal/mol",
            f"  top pose RMSD     : {self.top_pose_rmsd:.2f} A "
            f"(threshold {self.threshold:.1f} A)",
            f"  best pose RMSD    : {self.best_rmsd:.2f} A "
            f"(mode {self.rmsd_per_mode.index(self.best_rmsd) + 1})",
            f"  symmetry-corrected: {'yes' if self.symmetry_corrected else 'no'}",
        ]
        if not self.success and self.found_in_top_poses:
            lines.append(
                "  The crystal pose is among the returned modes but is not ranked first. "
                "The box and the preparation are fine; the scoring function is not."
            )
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


def pose_rmsd(probe_mol, reference_mol, conf_id: int = 0, symmetry: bool = True) -> float:
    """Symmetry-corrected, in-place RMSD between a pose and the crystal ligand.

    No alignment is done: the two molecules are already in the same frame, and
    aligning them would hide exactly the error we want to measure. Symmetry
    correction matters for rings and for equivalent substituents, where a
    chemically identical pose can otherwise score 1-2 A worse.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    probe = Chem.RemoveHs(Chem.Mol(probe_mol))
    reference = Chem.RemoveHs(Chem.Mol(reference_mol))
    if probe.GetNumAtoms() != reference.GetNumAtoms():
        raise LigandError(
            "The docked pose and the crystal ligand have different heavy-atom counts "
            f"({probe.GetNumAtoms()} vs {reference.GetNumAtoms()}). The template SMILES "
            "probably does not match the crystal ligand."
        )
    if symmetry:
        return float(rdMolAlign.CalcRMS(probe, reference, prbId=conf_id, refId=0))
    import numpy as np

    a = probe.GetConformer(conf_id).GetPositions()
    b = reference.GetConformer(0).GetPositions()
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def rmsd_to_crystal(
    out_pdbqt: str | os.PathLike[str],
    crystal_mol,
    symmetry: bool = True,
) -> list[float]:
    """RMSD of every docked pose in a Vina output file to the crystal ligand."""
    poses = poses_to_rdkit(out_pdbqt)
    values = []
    for conf in poses.GetConformers():
        values.append(pose_rmsd(poses, crystal_mol, conf_id=conf.GetId(), symmetry=symmetry))
    return values


def redock(
    receptor_pdbqt: str | os.PathLike[str],
    crystal_ligand: Residue,
    template_smiles: str,
    box: Box,
    out_dir: str | os.PathLike[str] = "redock",
    name: str | None = None,
    threshold: float = RMSD_SUCCESS,
    symmetry: bool = True,
    seed: int = 42,
    **dock_kwargs,
) -> RedockResult:
    """Dock the co-crystal ligand back into its own structure and score the result.

    ``template_smiles`` supplies the bond orders the PDB file does not have.
    Look it up once from the PDB chemical component (for example AQ4 in 1M17 is
    erlotinib) and keep it next to the structure.
    """
    ligand_name = name or crystal_ligand.res_name
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    crystal_mol = ligand_from_pdb_block(
        crystal_ligand.to_pdb(), template_smiles=template_smiles, name=f"{ligand_name}_crystal"
    )
    prepared = prepare_ligand(
        smiles=template_smiles, name=f"{ligand_name}_redock", out_dir=out_dir, seed=seed
    )
    docking = dock(
        receptor_pdbqt,
        prepared.pdbqt,
        box,
        out_pdbqt=out_dir / f"{ligand_name}_redock_out.pdbqt",
        name=ligand_name,
        seed=seed,
        **dock_kwargs,
    )
    values = rmsd_to_crystal(docking.out_pdbqt, crystal_mol, symmetry=symmetry)
    notes = []
    if len(values) != len(docking.poses):  # pragma: no cover - defensive
        notes.append("Pose count and RMSD count differ; check the output file.")
    return RedockResult(
        ligand=ligand_name,
        docking=docking,
        rmsd_per_mode=values,
        best_rmsd=min(values),
        top_pose_rmsd=values[0],
        top_affinity=docking.poses[0].affinity,
        threshold=threshold,
        symmetry_corrected=symmetry,
        notes=notes,
    )
