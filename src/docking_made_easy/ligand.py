"""Ligand preparation: SMILES, SDF or a crystal pose to a Vina-ready PDBQT.

RDKit builds the 3D structure, Meeko assigns the AutoDock types, the Gasteiger
charges and the rotatable bonds.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from .errors import LigandError, MissingToolError


def _rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:  # pragma: no cover - RDKit is a hard dependency
        raise MissingToolError("rdkit", "Install it with `pip install rdkit`.") from exc
    return Chem, AllChem


def _meeko():
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except ImportError as exc:
        raise MissingToolError("meeko", "Install it with `pip install meeko`.") from exc
    return MoleculePreparation, PDBQTWriterLegacy


@dataclass
class PreparedLigand:
    name: str
    smiles: str
    pdbqt: Path
    n_heavy_atoms: int
    n_rotatable_bonds: int
    sdf: Path | None = None


def mol_from_smiles(smiles: str, name: str = "ligand"):
    """Parse a SMILES string, with a clear message when it is wrong."""
    Chem, _ = _rdkit()
    if not smiles or not smiles.strip():
        raise LigandError("Empty SMILES string.")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise LigandError(
            f"RDKit could not parse the SMILES {smiles!r} for {name!r}. "
            "Check the brackets, ring closures and aromatic atoms."
        )
    if mol.GetNumAtoms() == 0:
        raise LigandError(f"The SMILES {smiles!r} contains no atoms.")
    mol.SetProp("_Name", name)
    return mol


def embed_3d(mol, seed: int = 42, n_tries: int = 10, optimize: bool = True):
    """Add hydrogens and generate one reasonable 3D conformer (ETKDGv3 + MMFF)."""
    Chem, AllChem = _rdkit()
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useSmallRingTorsions = True
    status = AllChem.EmbedMolecule(mol, params)
    tries = 0
    while status != 0 and tries < n_tries:
        tries += 1
        params.randomSeed = seed + tries
        params.useRandomCoords = True
        status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        raise LigandError(
            f"RDKit could not generate a 3D conformer after {n_tries + 1} attempts. "
            "Very strained macrocycles sometimes need `--no-optimize` and a hand-built SDF."
        )
    if optimize:
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception:  # pragma: no cover - force-field edge cases
            pass
    return mol


def protonate(mol, ph: float = 7.4):
    """Apply a small set of pH 7.4 rules: acids deprotonated, strong bases protonated.

    This is deliberately simple and explicit. For anything subtle, prepare the
    ligand in a dedicated tool and pass the SDF with the right protonation.
    """
    Chem, _ = _rdkit()
    if ph is None:
        return mol
    rules = [
        ("[CX3](=O)[OX2H1]", "[CX3](=O)[O-]"),          # carboxylic acid
        ("[SX4](=O)(=O)[OX2H1]", "[SX4](=O)(=O)[O-]"),  # sulfonic acid
        ("[PX4](=O)([OX2H1])", "[PX4](=O)([O-])"),      # phosphate/phosphonate
        ("[NX3;H2,H1;!$(N[!C]);!$(N=*);!$(NC=[O,N,S])]", None),  # aliphatic amine
    ]
    editable = Chem.RWMol(mol)
    for smarts, _ in rules[:3]:
        patt = Chem.MolFromSmarts(smarts)
        for match in editable.GetSubstructMatches(patt):
            for idx in match:
                atom = editable.GetAtomWithIdx(idx)
                if atom.GetSymbol() == "O" and atom.GetTotalNumHs() == 1:
                    atom.SetFormalCharge(-1)
                    atom.SetNoImplicit(True)
                    atom.SetNumExplicitHs(0)
    amine = Chem.MolFromSmarts(rules[3][0])
    for match in editable.GetSubstructMatches(amine):
        atom = editable.GetAtomWithIdx(match[0])
        if atom.GetFormalCharge() == 0 and not atom.GetIsAromatic():
            atom.SetFormalCharge(1)
            atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
            atom.SetNoImplicit(True)
    out = editable.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:  # pragma: no cover - fall back to the neutral form
        return mol
    return out


def mol_to_pdbqt(mol, path: str | os.PathLike[str], rigid: bool = False) -> Path:
    """Write a 3D molecule to PDBQT with Meeko."""
    MoleculePreparation, PDBQTWriterLegacy = _meeko()
    if mol.GetNumConformers() == 0:
        raise LigandError("The molecule has no 3D conformer; call embed_3d() first.")
    prep = MoleculePreparation(rigid_macrocycles=rigid)
    try:
        setups = prep.prepare(mol)
    except Exception as exc:
        raise LigandError(f"Meeko could not prepare the ligand: {exc}") from exc
    if not setups:
        raise LigandError("Meeko returned no molecule setup.")
    text, ok, error = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise LigandError(f"Meeko could not write the PDBQT: {error}")
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    return out


def prepare_ligand(
    smiles: str | None = None,
    mol=None,
    name: str = "ligand",
    out_dir: str | os.PathLike[str] = "ligands",
    seed: int = 42,
    ph: float | None = 7.4,
    optimize: bool = True,
    write_sdf: bool = True,
) -> PreparedLigand:
    """SMILES (or an RDKit molecule) to a docking-ready PDBQT."""
    Chem, AllChem = _rdkit()
    from rdkit.Chem import Descriptors, rdMolDescriptors

    if mol is None:
        if smiles is None:
            raise LigandError("Give either a SMILES string or an RDKit molecule.")
        mol = mol_from_smiles(smiles, name=name)
    if ph is not None:
        mol = protonate(mol, ph=ph)
    smiles_out = Chem.MolToSmiles(mol)
    mol3d = embed_3d(mol, seed=seed, optimize=optimize)
    mol3d.SetProp("_Name", name)

    out_dir = Path(out_dir)
    pdbqt = mol_to_pdbqt(mol3d, out_dir / f"{name}.pdbqt")
    sdf_path = None
    if write_sdf:
        sdf_path = out_dir / f"{name}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol3d)
        writer.close()

    heavy = Chem.RemoveHs(mol3d)
    return PreparedLigand(
        name=name,
        smiles=smiles_out,
        pdbqt=pdbqt,
        n_heavy_atoms=heavy.GetNumAtoms(),
        n_rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(heavy),
        sdf=sdf_path,
    )


def read_ligands(path: str | os.PathLike[str]) -> Iterator[tuple[str, object | None, str | None]]:
    """Read ligands from ``.smi``/``.csv``/``.txt`` (SMILES) or ``.sdf``.

    Yields ``(name, mol_or_None, smiles_or_None)``.
    """
    Chem, _ = _rdkit()
    p = Path(path)
    if not p.is_file():
        raise LigandError(f"No such ligand file: {p}")
    suffix = p.suffix.lower()
    if suffix in (".sdf", ".mol", ".sd"):
        supplier = Chem.SDMolSupplier(str(p), removeHs=False)
        found = False
        for i, mol in enumerate(supplier):
            if mol is None:
                continue
            found = True
            name = mol.GetProp("_Name") if mol.HasProp("_Name") and mol.GetProp("_Name") else f"{p.stem}_{i + 1}"
            yield name, mol, None
        if not found:
            raise LigandError(f"{p} contains no readable molecules.")
        return

    rows = [ln.strip() for ln in p.read_text().splitlines() if ln.strip()]
    if not rows:
        raise LigandError(f"{p} is empty.")
    header = rows[0].lower()
    start = 1 if ("smiles" in header and ("," in header or "\t" in header or " " in header)) else 0
    count = 0
    for i, row in enumerate(rows[start:], start=1):
        parts = [c for c in row.replace("\t", ",").replace(" ", ",").split(",") if c]
        if not parts:
            continue
        smiles = parts[0]
        name = parts[1] if len(parts) > 1 else f"{p.stem}_{i}"
        count += 1
        yield name, None, smiles
    if count == 0:
        raise LigandError(f"No SMILES found in {p}.")


def prepare_ligands(
    path: str | os.PathLike[str],
    out_dir: str | os.PathLike[str] = "ligands",
    seed: int = 42,
    ph: float | None = 7.4,
    skip_failures: bool = True,
) -> tuple[list[PreparedLigand], list[tuple[str, str]]]:
    """Prepare every ligand in a file. Returns (prepared, failures)."""
    prepared: list[PreparedLigand] = []
    failures: list[tuple[str, str]] = []
    seen: set[str] = set()
    for name, mol, smiles in read_ligands(path):
        safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in name)[:60] or "ligand"
        unique = safe
        n = 1
        while unique in seen:
            n += 1
            unique = f"{safe}_{n}"
        seen.add(unique)
        try:
            prepared.append(
                prepare_ligand(
                    smiles=smiles, mol=mol, name=unique, out_dir=out_dir, seed=seed, ph=ph
                )
            )
        except Exception as exc:
            if not skip_failures:
                raise
            failures.append((unique, str(exc)))
    if not prepared:
        raise LigandError(
            f"No ligand from {path} could be prepared. First error: "
            + (failures[0][1] if failures else "unknown")
        )
    return prepared, failures


def ligand_from_pdb_block(block: str, template_smiles: str | None = None, name: str = "crystal"):
    """Rebuild a crystal ligand from PDB text, using a SMILES template for bond orders.

    A PDB file has no bond orders, so without a template the molecule is only
    good for coordinates (which is all the box needs, but not enough for RMSD).
    """
    Chem, AllChem = _rdkit()
    mol = Chem.MolFromPDBBlock(block, removeHs=True, sanitize=False)
    if mol is None:
        raise LigandError("RDKit could not read the crystal ligand PDB block.")
    if template_smiles:
        template = mol_from_smiles(template_smiles, name=name)
        if template.GetNumAtoms() != mol.GetNumAtoms():
            raise LigandError(
                "The template SMILES does not match the crystal ligand: "
                f"{mol.GetNumAtoms()} heavy atoms in the PDB block vs "
                f"{template.GetNumAtoms()} in the template."
            )
        try:
            mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
        except Exception as exc:
            raise LigandError(
                "The template SMILES does not match the crystal ligand "
                f"({mol.GetNumAtoms()} heavy atoms vs {template.GetNumAtoms()}): {exc}"
            ) from exc
    mol.SetProp("_Name", name)
    return mol
