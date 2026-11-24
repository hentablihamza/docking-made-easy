"""Turn a crystal structure into a receptor PDBQT that Vina can read.

The steps are the ones everyone does by hand, in the order that matters:

1. keep one model and one alternate location,
2. drop the ligand, the waters and the buffer molecules,
3. add the missing heavy atoms and the hydrogens at a chosen pH (PDBFixer),
4. assign AutoDock atom types and charges (Meeko), with Open Babel as a fallback.
"""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

from .errors import MissingToolError, StructureError
from .structure import Atom, Structure, write_pdb
from .tools import require, which


@dataclass
class ReceptorPrepResult:
    """Everything the preparation produced, so a run can be reported."""

    pdbqt: Path
    protonated_pdb: Path
    raw_pdb: Path
    method: str
    ph: float
    n_atoms: int
    n_hydrogens: int
    removed: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"Receptor PDBQT: {self.pdbqt}",
            f"  method          : {self.method}",
            f"  protonation pH  : {self.ph}",
            f"  atoms           : {self.n_atoms} ({self.n_hydrogens} hydrogens)",
        ]
        for key, value in self.removed.items():
            lines.append(f"  removed {key:<8}: {value}")
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


def add_hydrogens(
    pdb_in: str | os.PathLike[str],
    pdb_out: str | os.PathLike[str],
    ph: float = 7.4,
    add_missing_residues: bool = False,
    keep_oxt: bool = False,
) -> tuple[Path, list[str]]:
    """Add missing heavy atoms and hydrogens with PDBFixer.

    ``keep_oxt=False`` removes the terminal OXT oxygens that PDBFixer adds.
    Meeko's residue templates treat a ragged chain end as blunt and refuse a
    terminal residue that carries an extra oxygen, which is the most common
    reason receptor preparation fails on a crystal structure with chain breaks.
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError as exc:
        raise MissingToolError(
            "pdbfixer",
            "Install it with `conda install -c conda-forge pdbfixer openmm`, "
            "or protonate the receptor yourself and pass --protonated.",
        ) from exc

    notes: list[str] = []
    fixer = PDBFixer(filename=str(pdb_in))
    fixer.findMissingResidues()
    if not add_missing_residues:
        if fixer.missingResidues:
            notes.append(
                f"{len(fixer.missingResidues)} chain gap(s) left unmodelled "
                "(use --add-missing-residues to let PDBFixer build them, but loops "
                "modelled from nothing are a guess)."
            )
        fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        names = ", ".join(sorted({r[0].name for r in fixer.nonstandardResidues}))
        notes.append(f"Non-standard residues replaced by their standard form: {names}")
        fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    n_missing = sum(len(v) for v in fixer.missingAtoms.values())
    if n_missing:
        notes.append(f"{n_missing} missing heavy atom(s) rebuilt in incomplete side chains.")
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    out = Path(pdb_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    if not keep_oxt:
        lines = out.read_text().splitlines()
        kept = [ln for ln in lines if not (ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() in ("OXT", "HXT"))]
        removed = len(lines) - len(kept)
        if removed:
            out.write_text("\n".join(kept) + "\n")
            notes.append(
                f"{removed} terminal OXT/HXT atom(s) removed for Meeko compatibility "
                "(pass --keep-oxt to keep them)."
            )
    return out, notes


def _run(cmd: Sequence[str], what: str) -> subprocess.CompletedProcess:
    result = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if result.returncode != 0:
        tail = (result.stderr or result.stdout or "").strip().splitlines()[-6:]
        raise StructureError(f"{what} failed:\n  " + "\n  ".join(tail))
    return result


def pdbqt_with_meeko(
    protonated_pdb: str | os.PathLike[str],
    output_basename: str | os.PathLike[str],
    extra_args: Sequence[str] = (),
) -> Path:
    """Assign AutoDock types and charges with ``mk_prepare_receptor.py``."""
    exe = require("mk_prepare_receptor.py")
    base = Path(output_basename)
    base.parent.mkdir(parents=True, exist_ok=True)
    _run(
        [exe, "--read_pdb", str(protonated_pdb), "-o", str(base), "-p", *extra_args],
        "Meeko receptor preparation",
    )
    pdbqt = base.with_suffix(".pdbqt")
    if not pdbqt.is_file():
        raise StructureError(f"Meeko did not write {pdbqt}.")
    return pdbqt


def pdbqt_with_openbabel(
    protonated_pdb: str | os.PathLike[str],
    pdbqt_out: str | os.PathLike[str],
) -> Path:
    """Fallback route: Open Babel writes a rigid receptor PDBQT (``-xr``)."""
    exe = require("obabel")
    out = Path(pdbqt_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    _run(
        [exe, str(protonated_pdb), "-O", str(out), "-xr", "--partialcharge", "gasteiger"],
        "Open Babel receptor conversion",
    )
    if not out.is_file() or out.stat().st_size == 0:
        raise StructureError(f"Open Babel did not write {out}.")
    return out


def prepare_receptor(
    source: str | os.PathLike[str],
    out_dir: str | os.PathLike[str] = "receptor",
    name: str = "receptor",
    keep_waters: bool = False,
    keep_ions: bool = False,
    keep_residues: Sequence[str] = (),
    chains: Sequence[str] = (),
    ph: float = 7.4,
    method: str = "auto",
    add_missing_residues: bool = False,
    keep_oxt: bool = False,
    protonated: str | os.PathLike[str] | None = None,
    structure: Structure | None = None,
) -> ReceptorPrepResult:
    """Prepare a receptor PDBQT from a PDB ID or a PDB file.

    ``method`` is ``"meeko"``, ``"openbabel"`` or ``"auto"`` (Meeko first, Open
    Babel if Meeko is missing or refuses the structure).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    if structure is None:
        source_str = str(source)
        if len(source_str) == 4 and not Path(source_str).exists():
            structure = Structure.fetch(source_str)
        else:
            structure = Structure.from_file(source_str)

    before = structure.summary()
    receptor_atoms: list[Atom] = structure.receptor_atoms(
        keep_residues=keep_residues, keep_waters=keep_waters, keep_ions=keep_ions, chains=chains
    )
    raw_pdb = write_pdb(
        receptor_atoms,
        out / f"{name}_raw.pdb",
        remarks=[
            "Receptor prepared with docking-made-easy",
            f"Source: {structure.pdb_id or structure.path}",
        ],
    )
    removed = {
        "ligands": [lig["key"] for lig in before["ligands"]] if not keep_residues else [],
        "waters": 0 if keep_waters else before["n_waters"],
        "ions": [] if keep_ions else before["ions"],
    }

    notes: list[str] = []
    if protonated is not None:
        protonated_pdb = Path(protonated)
        if not protonated_pdb.is_file():
            raise StructureError(f"No such protonated receptor: {protonated_pdb}")
        notes.append(f"Used the supplied protonated receptor {protonated_pdb}.")
    else:
        protonated_pdb, fixer_notes = add_hydrogens(
            raw_pdb,
            out / f"{name}_h.pdb",
            ph=ph,
            add_missing_residues=add_missing_residues,
            keep_oxt=keep_oxt,
        )
        notes.extend(fixer_notes)

    used = method
    if method in ("auto", "meeko") and which("mk_prepare_receptor.py"):
        try:
            pdbqt = pdbqt_with_meeko(protonated_pdb, out / name)
            used = "meeko"
        except (StructureError, MissingToolError) as exc:
            if method == "meeko":
                raise
            notes.append(f"Meeko refused the structure ({exc}); fell back to Open Babel.")
            pdbqt = pdbqt_with_openbabel(protonated_pdb, out / f"{name}.pdbqt")
            used = "openbabel"
    elif method == "meeko":
        raise MissingToolError(
            "mk_prepare_receptor.py", "Install Meeko (`pip install meeko`) or use --method openbabel."
        )
    else:
        pdbqt = pdbqt_with_openbabel(protonated_pdb, out / f"{name}.pdbqt")
        used = "openbabel"

    text = pdbqt.read_text().splitlines()
    atoms = [ln for ln in text if ln.startswith(("ATOM", "HETATM"))]
    hydrogens = [ln for ln in atoms if ln[12:16].strip().startswith("H")]
    return ReceptorPrepResult(
        pdbqt=pdbqt,
        protonated_pdb=Path(protonated_pdb),
        raw_pdb=raw_pdb,
        method=used,
        ph=ph,
        n_atoms=len(atoms),
        n_hydrogens=len(hydrogens),
        removed=removed,
        notes=notes,
    )
