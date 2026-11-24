"""The docking box: where Vina is allowed to look.

The box is the single setting that decides whether a docking run means
anything. Too small and the right pose cannot fit; too large and Vina spreads
its search thin and returns nonsense with a good-looking score.

When a structure already contains a co-crystal ligand, that ligand *is* the
answer: it marks the pocket experimentally. This module turns it into a box.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from .errors import StructureError
from .structure import Atom, Residue

#: AutoDock Vina's default grid spacing (A). Sizes are rounded to a multiple of it.
GRID_SPACING = 0.375

#: Padding added around the ligand, in angstrom. 4 A is enough room for a
#: ligand of similar size to rotate and translate without hitting the wall.
DEFAULT_PADDING = 4.0

#: Vina's search slows down sharply past this; warn rather than fail.
LARGE_BOX_EDGE = 30.0


def _as_array(points) -> np.ndarray:
    array = np.asarray(points, dtype=float)
    if array.ndim != 2 or array.shape[1] != 3:
        raise StructureError("Coordinates must be an (N, 3) array.")
    if array.shape[0] == 0:
        raise StructureError("Cannot build a box from an empty set of atoms.")
    return array


@dataclass(frozen=True)
class Box:
    """A Vina search box: a centre and the edge lengths, in angstrom."""

    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    source: str = "manual"

    # -- convenience ---------------------------------------------------------

    @property
    def center(self) -> tuple[float, float, float]:
        return (self.center_x, self.center_y, self.center_z)

    @property
    def size(self) -> tuple[float, float, float]:
        return (self.size_x, self.size_y, self.size_z)

    @property
    def volume(self) -> float:
        return self.size_x * self.size_y * self.size_z

    @property
    def min_corner(self) -> tuple[float, float, float]:
        return tuple(c - s / 2 for c, s in zip(self.center, self.size))  # type: ignore[return-value]

    @property
    def max_corner(self) -> tuple[float, float, float]:
        return tuple(c + s / 2 for c, s in zip(self.center, self.size))  # type: ignore[return-value]

    def contains(self, points, margin: float = 0.0) -> np.ndarray:
        """Boolean mask of which points fall inside the box."""
        array = _as_array(points)
        low = np.array(self.min_corner) - margin
        high = np.array(self.max_corner) + margin
        return np.all((array >= low) & (array <= high), axis=1)

    def fraction_inside(self, points, margin: float = 0.0) -> float:
        mask = self.contains(points, margin=margin)
        return float(mask.sum()) / float(len(mask))

    def warnings(self) -> list[str]:
        """Human-readable warnings about a box that is likely to misbehave."""
        notes = []
        for axis, value in zip("xyz", self.size):
            if value < 10.0:
                notes.append(
                    f"size_{axis} is {value:.1f} A, which is small: a ligand needs room "
                    "to rotate, so give at least ~10 A per axis."
                )
            if value > LARGE_BOX_EDGE:
                notes.append(
                    f"size_{axis} is {value:.1f} A. Boxes above {LARGE_BOX_EDGE:.0f} A per side "
                    "slow Vina down and weaken the search; raise --exhaustiveness or shrink the box."
                )
        return notes

    # -- serialisation -------------------------------------------------------

    def to_dict(self) -> dict:
        return asdict(self)

    def to_config(self) -> str:
        """Vina ``--config`` text."""
        return (
            f"center_x = {self.center_x:.3f}\n"
            f"center_y = {self.center_y:.3f}\n"
            f"center_z = {self.center_z:.3f}\n"
            f"size_x = {self.size_x:.3f}\n"
            f"size_y = {self.size_y:.3f}\n"
            f"size_z = {self.size_z:.3f}\n"
        )

    def write_config(self, path: str | os.PathLike[str]) -> Path:
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(self.to_config())
        return out

    def write_json(self, path: str | os.PathLike[str]) -> Path:
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(self.to_dict(), indent=2) + "\n")
        return out

    @classmethod
    def from_config(cls, path: str | os.PathLike[str]) -> "Box":
        values: dict[str, float] = {}
        for line in Path(path).read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            key, _, raw = line.partition("=")
            key = key.strip()
            if key in {
                "center_x", "center_y", "center_z", "size_x", "size_y", "size_z"
            }:
                values[key] = float(raw.strip())
        missing = {
            "center_x", "center_y", "center_z", "size_x", "size_y", "size_z"
        } - values.keys()
        if missing:
            raise StructureError(
                f"{path} is missing {', '.join(sorted(missing))}."
            )
        return cls(**values, source=f"config:{Path(path).name}")

    def to_pymol(self, name: str = "docking_box") -> str:
        """A PyMOL script that draws the box, so you can check it by eye."""
        cx, cy, cz = self.center
        sx, sy, sz = self.size
        return (
            "from pymol.cgo import *\n"
            "from pymol import cmd\n"
            f"minX, minY, minZ = {cx - sx / 2:.3f}, {cy - sy / 2:.3f}, {cz - sz / 2:.3f}\n"
            f"maxX, maxY, maxZ = {cx + sx / 2:.3f}, {cy + sy / 2:.3f}, {cz + sz / 2:.3f}\n"
            "box = [LINEWIDTH, 2.0, BEGIN, LINES, COLOR, 1.0, 1.0, 0.0]\n"
            "corners = [(minX, minY, minZ), (maxX, minY, minZ), (maxX, maxY, minZ), (minX, maxY, minZ),\n"
            "           (minX, minY, maxZ), (maxX, minY, maxZ), (maxX, maxY, maxZ), (minX, maxY, maxZ)]\n"
            "edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]\n"
            "for i, j in edges:\n"
            "    box += [VERTEX, *corners[i], VERTEX, *corners[j]]\n"
            "box.append(END)\n"
            f"cmd.load_cgo(box, '{name}')\n"
        )

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"center ({self.center_x:.3f}, {self.center_y:.3f}, {self.center_z:.3f}) "
            f"size ({self.size_x:.1f} x {self.size_y:.1f} x {self.size_z:.1f}) A"
        )


def _round_to_grid(value: float, spacing: float = GRID_SPACING) -> float:
    if spacing <= 0:
        return value
    return math.ceil(value / spacing) * spacing


def box_from_points(
    points,
    padding: float = DEFAULT_PADDING,
    mode: str = "extent",
    cubic: bool = False,
    min_size: float = 12.0,
    max_size: float | None = None,
    snap_to_grid: bool = True,
    source: str = "points",
) -> Box:
    """Build a box around a set of coordinates.

    ``mode="extent"`` (the default) centres the box on the midpoint of the
    bounding box and sizes it from the ligand's own extent plus ``padding`` on
    every side. ``mode="centroid"`` centres it on the mean atom position
    instead, which differs when a ligand is lopsided, and sizes it so that every
    atom still sits ``padding`` inside the wall.
    """
    array = _as_array(points)
    low = array.min(axis=0)
    high = array.max(axis=0)
    extent = high - low

    if mode == "extent":
        center = (low + high) / 2.0
        size = extent + 2.0 * padding
    elif mode == "centroid":
        center = array.mean(axis=0)
        reach = np.maximum(np.abs(high - center), np.abs(center - low))
        size = 2.0 * (reach + padding)
    else:
        raise StructureError("mode must be 'extent' or 'centroid'.")

    size = np.maximum(size, float(min_size))
    if cubic:
        size = np.repeat(size.max(), 3)
    if max_size is not None:
        size = np.minimum(size, float(max_size))
    if snap_to_grid:
        size = np.array([_round_to_grid(float(s)) for s in size])

    return Box(
        center_x=float(center[0]),
        center_y=float(center[1]),
        center_z=float(center[2]),
        size_x=float(size[0]),
        size_y=float(size[1]),
        size_z=float(size[2]),
        source=source,
    )


def box_from_ligand(
    ligand: Residue | Iterable[Atom],
    padding: float = DEFAULT_PADDING,
    mode: str = "extent",
    heavy_atoms_only: bool = True,
    **kwargs,
) -> Box:
    """Build the box from a co-crystal ligand. This is the recommended route."""
    if isinstance(ligand, Residue):
        coords = ligand.heavy_coords if heavy_atoms_only else ligand.coords
        source = f"ligand:{ligand.key}"
    else:
        atoms = [a for a in ligand]
        if heavy_atoms_only:
            atoms = [a for a in atoms if not a.is_hydrogen] or atoms
        if not atoms:
            raise StructureError("Cannot build a box from an empty ligand.")
        coords = np.array([[a.x, a.y, a.z] for a in atoms], dtype=float)
        source = "ligand:atoms"
    kwargs.setdefault("source", source)
    return box_from_ligand_coords(coords, padding=padding, mode=mode, **kwargs)


def box_from_ligand_coords(coords, padding: float = DEFAULT_PADDING, **kwargs) -> Box:
    return box_from_points(coords, padding=padding, **kwargs)


def box_from_residues(
    residues: Sequence[Residue],
    padding: float = DEFAULT_PADDING,
    **kwargs,
) -> Box:
    """Build a box around a set of pocket residues (no ligand needed)."""
    if not residues:
        raise StructureError("No residues selected for the box.")
    coords = np.vstack([r.heavy_coords for r in residues])
    kwargs.setdefault("source", "residues:" + ",".join(r.key for r in residues))
    return box_from_points(coords, padding=padding, **kwargs)


def box_from_center(
    center: Sequence[float],
    size: float | Sequence[float] = 22.5,
    snap_to_grid: bool = False,
) -> Box:
    """A manual box: an explicit centre and size."""
    if len(center) != 3:
        raise StructureError("center must be three numbers: x y z.")
    sizes = [float(size)] * 3 if isinstance(size, (int, float)) else [float(s) for s in size]
    if len(sizes) != 3:
        raise StructureError("size must be one number or three.")
    if snap_to_grid:
        sizes = [_round_to_grid(s) for s in sizes]
    if any(s <= 0 for s in sizes):
        raise StructureError("Box sizes must be positive.")
    return Box(
        center_x=float(center[0]),
        center_y=float(center[1]),
        center_z=float(center[2]),
        size_x=sizes[0],
        size_y=sizes[1],
        size_z=sizes[2],
        source="manual",
    )


def pocket_residues(
    receptor_atoms: Iterable[Atom],
    ligand: Residue,
    cutoff: float = 5.0,
) -> list[tuple[str, float]]:
    """Residues with any heavy atom within ``cutoff`` A of the ligand.

    Useful as a sanity check: the box should enclose these residues' side
    chains, and they are what a reviewer will ask about.
    """
    ligand_coords = ligand.heavy_coords
    groups: dict[tuple[str, str, int, str], list[Atom]] = {}
    for atom in receptor_atoms:
        if atom.is_hydrogen:
            continue
        groups.setdefault((atom.res_name, atom.chain, atom.res_seq, atom.i_code), []).append(atom)
    out = []
    for key, atoms in groups.items():
        coords = np.array([[a.x, a.y, a.z] for a in atoms], dtype=float)
        distance = float(
            np.sqrt(((coords[:, None, :] - ligand_coords[None, :, :]) ** 2).sum(-1)).min()
        )
        if distance <= cutoff:
            out.append((f"{key[0]}{key[2]}{key[3]}".strip() + f":{key[1]}", distance))
    return sorted(out, key=lambda item: item[1])


def check_box(box: Box, ligand: Residue, receptor_atoms: Iterable[Atom] | None = None) -> dict:
    """Report whether the box actually covers the ligand, and how much is protein.

    A box full of solvent wastes the search; a box that clips the ligand cannot
    reproduce the crystal pose.
    """
    coords = ligand.heavy_coords
    report = {
        "ligand": ligand.key,
        "ligand_atoms_inside": float(box.fraction_inside(coords)),
        "ligand_extent": [float(v) for v in (coords.max(axis=0) - coords.min(axis=0))],
        "box_volume": box.volume,
        "warnings": box.warnings(),
    }
    if report["ligand_atoms_inside"] < 1.0:
        report["warnings"].append(
            "The box does not contain every ligand atom. Increase --padding."
        )
    if receptor_atoms is not None:
        heavy = np.array(
            [[a.x, a.y, a.z] for a in receptor_atoms if not a.is_hydrogen], dtype=float
        )
        if len(heavy):
            inside = int(box.contains(heavy).sum())
            report["receptor_heavy_atoms_in_box"] = inside
            report["receptor_atom_density"] = inside / box.volume * 1000.0
            if inside < 50:
                report["warnings"].append(
                    f"Only {inside} receptor heavy atoms are inside the box. "
                    "Is the box centred on the pocket, or out in the solvent?"
                )
    return report
