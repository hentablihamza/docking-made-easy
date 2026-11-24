"""Parse Vina output into a tidy table, and get the poses back as SDF or PDB.

Vina's own output is a PDBQT file with ``REMARK VINA RESULT`` lines and a
progress table on stdout. Neither is convenient. Both are parsed here into one
row per pose.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Sequence

from .errors import LigandError, ParseError

VINA_RESULT = re.compile(
    r"^REMARK\s+VINA\s+RESULT:\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)", re.IGNORECASE
)
STDOUT_ROW = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")
NAME_REMARK = re.compile(r"^REMARK\s+Name\s*=\s*(.+)$", re.IGNORECASE)


@dataclass
class Pose:
    """One docked pose."""

    ligand: str
    mode: int
    affinity: float          # kcal/mol
    rmsd_lb: float           # A, to the best mode
    rmsd_ub: float           # A, to the best mode
    pdbqt_block: str = ""

    def as_row(self) -> dict:
        row = asdict(self)
        row.pop("pdbqt_block")
        return row


def parse_vina_pdbqt(
    path: str | os.PathLike[str], ligand_name: str | None = None
) -> list[Pose]:
    """Parse a Vina output PDBQT (``--out``) into poses, keeping the coordinates."""
    p = Path(path)
    if not p.is_file():
        raise ParseError(f"No such Vina output file: {p}")
    text = p.read_text()
    if "VINA RESULT" not in text.upper():
        raise ParseError(
            f"{p} has no 'REMARK VINA RESULT' lines. Is it a Vina output file, "
            "or the input ligand PDBQT?"
        )
    name = ligand_name or p.stem.replace("_out", "")
    poses: list[Pose] = []
    current: list[str] = []
    scores: tuple[float, float, float] | None = None
    embedded_name: str | None = None

    def flush() -> None:
        nonlocal current, scores
        if scores is not None and current:
            poses.append(
                Pose(
                    ligand=embedded_name or name,
                    mode=len(poses) + 1,
                    affinity=scores[0],
                    rmsd_lb=scores[1],
                    rmsd_ub=scores[2],
                    pdbqt_block="\n".join(current) + "\n",
                )
            )
        current = []
        scores = None

    for line in text.splitlines():
        if line.startswith("MODEL"):
            flush()
            current = [line]
            continue
        match = VINA_RESULT.match(line)
        if match:
            scores = (float(match.group(1)), float(match.group(2)), float(match.group(3)))
        name_match = NAME_REMARK.match(line)
        if name_match and embedded_name is None:
            embedded_name = name_match.group(1).strip()
        current.append(line)
        if line.startswith("ENDMDL"):
            flush()
    flush()
    if not poses:
        raise ParseError(f"No poses could be parsed from {p}.")
    for i, pose in enumerate(poses, start=1):
        pose.mode = i
        if embedded_name and ligand_name is None:
            pose.ligand = embedded_name
    return poses


def parse_vina_log(text: str, ligand_name: str = "ligand") -> list[Pose]:
    """Parse the score table Vina prints to stdout (no coordinates)."""
    poses = []
    for line in text.splitlines():
        match = STDOUT_ROW.match(line)
        if match:
            poses.append(
                Pose(
                    ligand=ligand_name,
                    mode=int(match.group(1)),
                    affinity=float(match.group(2)),
                    rmsd_lb=float(match.group(3)),
                    rmsd_ub=float(match.group(4)),
                )
            )
    if not poses:
        raise ParseError("No score table found in the Vina output.")
    return poses


def poses_to_dataframe(poses: Iterable[Pose]):
    """Tidy table: one row per pose."""
    import pandas as pd

    rows = [p.as_row() for p in poses]
    if not rows:
        raise ParseError("No poses to tabulate.")
    frame = pd.DataFrame(rows)
    return frame[["ligand", "mode", "affinity", "rmsd_lb", "rmsd_ub"]]


def results_table(results: Iterable["object"]):
    """Combine several docking runs into one table, ranked by best affinity."""
    import pandas as pd

    frames = []
    for result in results:
        poses = getattr(result, "poses", result)
        frames.append(poses_to_dataframe(poses))
    if not frames:
        raise ParseError("No docking results to tabulate.")
    table = pd.concat(frames, ignore_index=True)
    best = table.groupby("ligand")["affinity"].transform("min")
    table = table.assign(_best=best).sort_values(
        ["_best", "ligand", "mode"]
    ).drop(columns="_best")
    return table.reset_index(drop=True)


def summary_table(results: Iterable["object"]):
    """One row per ligand: the best affinity and how the modes spread out."""
    import pandas as pd

    table = results_table(results)
    grouped = table.groupby("ligand", sort=False)
    out = grouped.agg(
        best_affinity=("affinity", "min"),
        mean_affinity=("affinity", "mean"),
        worst_affinity=("affinity", "max"),
        n_modes=("mode", "count"),
    ).reset_index()
    top = table[table["mode"] == 1][["ligand"]]
    out = out.merge(top.drop_duplicates(), on="ligand", how="left")
    return out.sort_values("best_affinity").reset_index(drop=True)


def save_table(frame, path: str | os.PathLike[str]) -> Path:
    """Write a table as ``.csv``, ``.tsv`` or ``.parquet``."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    suffix = out.suffix.lower()
    if suffix == ".csv":
        frame.to_csv(out, index=False)
    elif suffix in (".tsv", ".tab"):
        frame.to_csv(out, index=False, sep="\t")
    elif suffix == ".parquet":
        frame.to_parquet(out, index=False)
    else:
        raise ParseError(f"Unknown table format '{suffix}'. Use .csv, .tsv or .parquet.")
    return out


# -- getting poses back into chemistry formats -------------------------------


def poses_to_rdkit(path: str | os.PathLike[str], only_best: bool = False):
    """Read a Vina output PDBQT back into RDKit molecules (one per pose).

    Meeko restores the bond orders and the hydrogens from the ``REMARK SMILES``
    lines it wrote during ligand preparation, so the result is a real molecule
    rather than a bag of atoms.
    """
    try:
        from meeko import PDBQTMolecule, RDKitMolCreate
    except ImportError as exc:
        raise LigandError(
            "Converting poses back to SDF needs Meeko (`pip install meeko`)."
        ) from exc
    from rdkit import Chem

    pdbqt_mol = PDBQTMolecule.from_file(str(path), skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    mols = [m for m in mols if m is not None]
    if not mols:
        raise ParseError(
            f"Meeko could not rebuild a molecule from {path}. The ligand PDBQT was "
            "probably not written by Meeko, so the SMILES remark is missing."
        )
    mol = mols[0]
    if only_best and mol.GetNumConformers() > 1:
        keep = Chem.Mol(mol)
        keep.RemoveAllConformers()
        keep.AddConformer(mol.GetConformer(0), assignId=True)
        return keep
    return mol


def poses_to_sdf(
    path: str | os.PathLike[str],
    sdf_out: str | os.PathLike[str],
    only_best: bool = False,
    affinities: Sequence[float] | None = None,
    name: str | None = None,
) -> Path:
    """Write the docked poses to SDF, with the affinity as an SD property."""
    from rdkit import Chem

    mol = poses_to_rdkit(path, only_best=only_best)
    out = Path(sdf_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(out))
    for i, conf in enumerate(mol.GetConformers()):
        pose = Chem.Mol(mol)
        pose.RemoveAllConformers()
        pose.AddConformer(conf, assignId=True)
        pose.SetProp("_Name", f"{name or Path(path).stem}_mode{i + 1}")
        pose.SetProp("mode", str(i + 1))
        if affinities is not None and i < len(affinities):
            pose.SetProp("vina_affinity_kcal_per_mol", f"{affinities[i]:.3f}")
        writer.write(pose)
        if only_best:
            break
    writer.close()
    return out


def pose_to_pdb(
    path: str | os.PathLike[str], pdb_out: str | os.PathLike[str], mode: int = 1
) -> Path:
    """Write one pose as a PDB file (handy for PyMOL and for complex building)."""
    from rdkit import Chem

    mol = poses_to_rdkit(path)
    if mode < 1 or mode > mol.GetNumConformers():
        raise ParseError(
            f"Mode {mode} does not exist; the file has {mol.GetNumConformers()} pose(s)."
        )
    out = Path(pdb_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    Chem.MolToPDBFile(mol, str(out), confId=mol.GetConformer(mode - 1).GetId())
    return out


def write_complex(
    receptor_pdb: str | os.PathLike[str],
    pose_pdb: str | os.PathLike[str],
    out_path: str | os.PathLike[str],
) -> Path:
    """Concatenate a receptor and a docked pose into one PDB complex."""
    receptor_lines = [
        ln for ln in Path(receptor_pdb).read_text().splitlines()
        if ln.startswith(("ATOM", "HETATM"))
    ]
    pose_lines = [
        "HETATM" + ln[6:] for ln in Path(pose_pdb).read_text().splitlines()
        if ln.startswith(("ATOM", "HETATM"))
    ]
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(receptor_lines + ["TER"] + pose_lines + ["END"]) + "\n")
    return out
