"""Read a PDB file and split it into receptor, ligands, waters and ions.

Everything here is plain Python plus NumPy, so it works without AutoDock Vina,
Open Babel or Meeko installed. That keeps the box calculation and the structure
inspection usable (and testable) on any machine.
"""

from __future__ import annotations

import os
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from .errors import StructureError

#: Residue names that are solvent, buffer or cryoprotectant rather than a ligand.
SOLVENT_RESIDUES = frozenset(
    {
        "HOH", "WAT", "DOD", "H2O", "SOL", "TIP3", "TIP",
        "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "DMS", "MPD",
        "TRS", "EPE", "MES", "ACT", "ACY", "FMT", "SO4", "PO4", "IPA", "BME",
        "NO3", "CIT", "TLA", "MLI", "BTB", "IMD", "AZI", "CO3", "SCN", "BU3",
    }
)

#: Monoatomic ions. They are not docking ligands, but they can matter in the site.
ION_RESIDUES = frozenset(
    {
        "NA", "K", "LI", "CS", "RB", "MG", "CA", "ZN", "MN", "FE", "FE2", "CO",
        "NI", "CU", "CU1", "CD", "HG", "CL", "BR", "IOD", "F", "SR", "BA", "AU",
        "PT", "AG", "PB", "TL", "YB", "SM", "GD", "EU", "W", "MO", "V",
    }
)

#: The 20 standard residues plus the usual alternatives, and nucleotides.
STANDARD_RESIDUES = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "MSE", "SEC", "PYL", "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN",
        "A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DI", "DU",
        "ACE", "NME", "NMA", "UNK",
    }
)

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


@dataclass(frozen=True)
class Atom:
    """One ATOM/HETATM record."""

    record: str
    serial: int
    name: str
    alt_loc: str
    res_name: str
    chain: str
    res_seq: int
    i_code: str
    x: float
    y: float
    z: float
    occupancy: float
    b_factor: float
    element: str
    line: str

    @property
    def coord(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)

    @property
    def is_hydrogen(self) -> bool:
        return self.element.upper() == "H" or self.element.upper() == "D"


@dataclass(frozen=True)
class Residue:
    """A residue (or a HETATM group such as a bound ligand)."""

    res_name: str
    chain: str
    res_seq: int
    i_code: str
    atoms: tuple[Atom, ...] = field(repr=False)

    @property
    def key(self) -> str:
        """Short selector, e.g. ``AQ4:A:1``."""
        return f"{self.res_name}:{self.chain}:{self.res_seq}{self.i_code}".strip()

    @property
    def n_heavy_atoms(self) -> int:
        return sum(1 for a in self.atoms if not a.is_hydrogen)

    @property
    def coords(self) -> np.ndarray:
        return np.array([[a.x, a.y, a.z] for a in self.atoms], dtype=float)

    @property
    def heavy_coords(self) -> np.ndarray:
        heavy = [a for a in self.atoms if not a.is_hydrogen]
        if not heavy:
            raise StructureError(f"{self.key} has no heavy atoms.")
        return np.array([[a.x, a.y, a.z] for a in heavy], dtype=float)

    def to_pdb(self) -> str:
        return "".join(a.line if a.line.endswith("\n") else a.line + "\n" for a in self.atoms)


def _f(line: str, start: int, end: int, default: float = 0.0) -> float:
    chunk = line[start:end].strip()
    if not chunk:
        return default
    try:
        return float(chunk)
    except ValueError:
        return default


def _i(line: str, start: int, end: int, default: int = 0) -> int:
    chunk = line[start:end].strip()
    if not chunk:
        return default
    try:
        return int(chunk)
    except ValueError:
        return default


def _element_of(line: str, name: str) -> str:
    element = line[76:78].strip()
    if element:
        return element.capitalize()
    # Fall back to the atom name, as old PDB files often leave columns 77-78 empty.
    stripped = name.strip().lstrip("0123456789")
    return stripped[:1].capitalize() if stripped else "C"


def parse_pdb_atoms(text: str, keep_altloc: str = "A", model: int = 1) -> list[Atom]:
    """Parse ATOM/HETATM records from PDB text.

    Only the first model is kept by default, and alternate locations other than
    the blank one and ``keep_altloc`` are dropped (a very common silent source
    of duplicated atoms in docking setups).
    """
    atoms: list[Atom] = []
    current_model = 1
    for line in text.splitlines():
        if line.startswith("MODEL"):
            current_model = _i(line, 10, 14, default=current_model)
            continue
        if line.startswith("ENDMDL"):
            if current_model >= model:
                break
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if current_model != model:
            continue
        alt_loc = line[16:17].strip()
        if alt_loc and alt_loc.upper() != keep_altloc.upper():
            continue
        name = line[12:16].strip()
        atoms.append(
            Atom(
                record=line[:6].strip(),
                serial=_i(line, 6, 11),
                name=name,
                alt_loc=alt_loc,
                res_name=line[17:20].strip(),
                chain=(line[21:22].strip() or "A"),
                res_seq=_i(line, 22, 26),
                i_code=line[26:27].strip(),
                x=_f(line, 30, 38),
                y=_f(line, 38, 46),
                z=_f(line, 46, 54),
                occupancy=_f(line, 54, 60, 1.0),
                b_factor=_f(line, 60, 66),
                element=_element_of(line, name),
                line=line.rstrip("\n"),
            )
        )
    if not atoms:
        raise StructureError(
            "No ATOM/HETATM records found. Is this really a PDB file? "
            "mmCIF input is not supported; download the .pdb format instead."
        )
    return atoms


def group_residues(atoms: Iterable[Atom]) -> list[Residue]:
    """Group atoms into residues, keeping the order they appear in the file."""
    order: list[tuple[str, str, int, str]] = []
    buckets: dict[tuple[str, str, int, str], list[Atom]] = {}
    for atom in atoms:
        key = (atom.res_name, atom.chain, atom.res_seq, atom.i_code)
        if key not in buckets:
            buckets[key] = []
            order.append(key)
        buckets[key].append(atom)
    return [
        Residue(res_name=k[0], chain=k[1], res_seq=k[2], i_code=k[3], atoms=tuple(buckets[k]))
        for k in order
    ]


@dataclass
class Structure:
    """A parsed PDB structure, split into its parts."""

    atoms: list[Atom]
    path: Path | None = None
    pdb_id: str | None = None

    @classmethod
    def from_file(cls, path: str | os.PathLike[str], **kwargs) -> "Structure":
        p = Path(path)
        if not p.is_file():
            raise StructureError(f"No such file: {p}")
        text = p.read_text(errors="replace")
        return cls(atoms=parse_pdb_atoms(text, **kwargs), path=p)

    @classmethod
    def from_text(cls, text: str, **kwargs) -> "Structure":
        return cls(atoms=parse_pdb_atoms(text, **kwargs))

    @classmethod
    def fetch(
        cls,
        pdb_id: str,
        cache_dir: str | os.PathLike[str] | None = None,
        timeout: float = 60.0,
        **kwargs,
    ) -> "Structure":
        """Download a structure from the RCSB PDB (cached on disk)."""
        path = fetch_pdb(pdb_id, cache_dir=cache_dir, timeout=timeout)
        structure = cls.from_file(path, **kwargs)
        structure.pdb_id = pdb_id.upper()
        return structure

    # -- parts ---------------------------------------------------------------

    @property
    def residues(self) -> list[Residue]:
        return group_residues(self.atoms)

    @property
    def polymer_residues(self) -> list[Residue]:
        return [r for r in self.residues if r.res_name in STANDARD_RESIDUES]

    @property
    def waters(self) -> list[Residue]:
        return [r for r in self.residues if r.res_name in ("HOH", "WAT", "DOD", "H2O")]

    @property
    def ions(self) -> list[Residue]:
        return [r for r in self.residues if r.res_name in ION_RESIDUES]

    def hetero_ligands(self, min_heavy_atoms: int = 6) -> list[Residue]:
        """HETATM groups that look like a real small-molecule ligand.

        Solvent, buffers, cryoprotectants and ions are filtered out, and groups
        smaller than ``min_heavy_atoms`` are ignored.
        """
        out = []
        for res in self.residues:
            if res.res_name in STANDARD_RESIDUES or res.res_name in SOLVENT_RESIDUES:
                continue
            if res.res_name in ION_RESIDUES:
                continue
            if not all(a.record == "HETATM" for a in res.atoms):
                continue
            if res.n_heavy_atoms < min_heavy_atoms:
                continue
            out.append(res)
        return sorted(out, key=lambda r: -r.n_heavy_atoms)

    def select_ligand(self, selector: str | None = None, min_heavy_atoms: int = 6) -> Residue:
        """Pick the co-crystal ligand.

        ``selector`` may be a residue name (``AQ4``), ``NAME:CHAIN`` or
        ``NAME:CHAIN:RESSEQ``. Without a selector the largest hetero group wins,
        and an error lists the candidates when the choice is ambiguous.
        """
        candidates = self.hetero_ligands(min_heavy_atoms=min_heavy_atoms)
        if selector:
            parts = [p for p in selector.split(":")]
            name = parts[0].strip().upper()
            chain = parts[1].strip() if len(parts) > 1 and parts[1].strip() else None
            seq = int(parts[2]) if len(parts) > 2 and parts[2].strip() else None
            pool = [r for r in self.residues if r.res_name.upper() == name]
            if chain is not None:
                pool = [r for r in pool if r.chain == chain]
            if seq is not None:
                pool = [r for r in pool if r.res_seq == seq]
            if not pool:
                known = ", ".join(sorted({r.key for r in candidates})) or "none found"
                raise StructureError(
                    f"No residue matches '{selector}'. Ligand-like residues here: {known}"
                )
            if len(pool) > 1:
                options = ", ".join(r.key for r in pool)
                raise StructureError(
                    f"'{selector}' matches {len(pool)} residues ({options}). "
                    "Add the chain and residue number, e.g. AQ4:A:1."
                )
            return pool[0]
        if not candidates:
            raise StructureError(
                "No co-crystal ligand found (only protein, solvent and ions). "
                "Give the box centre manually with --center, or lower --min-heavy-atoms."
            )
        if len(candidates) > 1 and candidates[0].n_heavy_atoms == candidates[1].n_heavy_atoms:
            options = ", ".join(f"{r.key} ({r.n_heavy_atoms} heavy atoms)" for r in candidates)
            raise StructureError(
                f"Several ligands are equally large: {options}. Choose one with --ligand-res."
            )
        return candidates[0]

    def receptor_atoms(
        self,
        keep_residues: Sequence[str] = (),
        keep_waters: bool = False,
        keep_ions: bool = False,
        chains: Sequence[str] = (),
    ) -> list[Atom]:
        """Everything that should go into the receptor file."""
        keep = {r.upper() for r in keep_residues}
        chain_set = {c.upper() for c in chains}
        out: list[Atom] = []
        for atom in self.atoms:
            if chain_set and atom.chain.upper() not in chain_set:
                continue
            name = atom.res_name.upper()
            if name in STANDARD_RESIDUES:
                out.append(atom)
                continue
            if name in keep:
                out.append(atom)
                continue
            if name in ("HOH", "WAT", "DOD", "H2O"):
                if keep_waters:
                    out.append(atom)
                continue
            if name in ION_RESIDUES:
                if keep_ions:
                    out.append(atom)
                continue
            # Any other HETATM (ligand, buffer, sugar) is dropped.
        if not out:
            raise StructureError(
                "The receptor selection is empty. Check the --chains option: "
                f"chains present are {sorted({a.chain for a in self.atoms})}."
            )
        return out

    def summary(self) -> dict:
        ligands = self.hetero_ligands()
        return {
            "pdb_id": self.pdb_id,
            "n_atoms": len(self.atoms),
            "chains": sorted({a.chain for a in self.atoms}),
            "n_polymer_residues": len(self.polymer_residues),
            "n_waters": len(self.waters),
            "ions": sorted({r.res_name for r in self.ions}),
            "ligands": [
                {"key": r.key, "res_name": r.res_name, "n_heavy_atoms": r.n_heavy_atoms}
                for r in ligands
            ],
        }


def write_pdb(atoms: Sequence[Atom], path: str | os.PathLike[str], remarks: Sequence[str] = ()) -> Path:
    """Write atoms back out as a PDB file."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"REMARK   {r}" for r in remarks]
    lines += [a.line for a in atoms]
    lines.append("END")
    out.write_text("\n".join(lines) + "\n")
    return out


def fetch_pdb(
    pdb_id: str,
    cache_dir: str | os.PathLike[str] | None = None,
    timeout: float = 60.0,
) -> Path:
    """Download ``pdb_id`` from the RCSB PDB, reusing the cached copy if present."""
    code = pdb_id.strip().upper()
    if len(code) != 4 or not code[0].isdigit():
        raise StructureError(f"'{pdb_id}' is not a 4-character PDB ID (for example 1M17).")
    directory = Path(cache_dir) if cache_dir else Path.home() / ".cache" / "docking-made-easy"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{code}.pdb"
    if path.is_file() and path.stat().st_size > 0:
        return path
    url = RCSB_URL.format(pdb_id=code)
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:  # noqa: S310
            data = response.read()
    except Exception as exc:  # pragma: no cover - network failure path
        raise StructureError(
            f"Could not download {code} from the RCSB PDB ({exc}). "
            "Check the ID and your network, or pass a local .pdb file instead."
        ) from exc
    if not data:
        raise StructureError(f"The RCSB PDB returned an empty file for {code}.")
    path.write_bytes(data)
    return path


def split_structure(
    structure: Structure,
    ligand_selector: str | None = None,
    keep_waters: bool = False,
    keep_ions: bool = False,
    chains: Sequence[str] = (),
    min_heavy_atoms: int = 6,
) -> tuple[list[Atom], Residue]:
    """Split a structure into (receptor atoms, co-crystal ligand residue)."""
    ligand = structure.select_ligand(ligand_selector, min_heavy_atoms=min_heavy_atoms)
    receptor = structure.receptor_atoms(
        keep_waters=keep_waters, keep_ions=keep_ions, chains=chains
    )
    return receptor, ligand
