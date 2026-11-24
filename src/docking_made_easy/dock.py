"""Run AutoDock Vina: one ligand, or a whole library in parallel."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Sequence

from .box import Box
from .errors import VinaError
from .results import Pose, parse_vina_pdbqt
from .tools import require, which

#: Vina's own defaults. 8 is fine for a first pass; 16-32 for a careful run.
DEFAULT_EXHAUSTIVENESS = 8
DEFAULT_NUM_MODES = 9
DEFAULT_SEED = 42


@dataclass
class DockingResult:
    """The outcome of one Vina run."""

    ligand: str
    poses: list[Pose]
    out_pdbqt: Path
    log: str = ""
    seconds: float = 0.0
    command: list[str] = field(default_factory=list)
    box: Box | None = None

    @property
    def best_affinity(self) -> float:
        return min(p.affinity for p in self.poses)

    @property
    def best_pose(self) -> Pose:
        return min(self.poses, key=lambda p: p.affinity)

    def to_dataframe(self):
        from .results import poses_to_dataframe

        return poses_to_dataframe(self.poses)


def vina_executable() -> str:
    """Path to the Vina binary (``DME_VINA`` overrides the search)."""
    return require("vina", env_var="DME_VINA")


def vina_available() -> bool:
    return which("vina", env_var="DME_VINA") is not None


def build_command(
    receptor: str | os.PathLike[str],
    ligand: str | os.PathLike[str],
    box: Box,
    out_pdbqt: str | os.PathLike[str],
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    num_modes: int = DEFAULT_NUM_MODES,
    seed: int = DEFAULT_SEED,
    cpu: int | None = None,
    energy_range: float = 3.0,
    scoring: str = "vina",
    extra_args: Sequence[str] = (),
    executable: str = "vina",
) -> list[str]:
    """Assemble the Vina command line (kept separate so it is easy to test)."""
    if exhaustiveness < 1:
        raise VinaError("exhaustiveness must be at least 1.")
    if num_modes < 1:
        raise VinaError("num_modes must be at least 1.")
    cmd = [
        str(executable),
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--out", str(out_pdbqt),
        "--center_x", f"{box.center_x:.3f}",
        "--center_y", f"{box.center_y:.3f}",
        "--center_z", f"{box.center_z:.3f}",
        "--size_x", f"{box.size_x:.3f}",
        "--size_y", f"{box.size_y:.3f}",
        "--size_z", f"{box.size_z:.3f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--energy_range", str(energy_range),
        "--seed", str(seed),
    ]
    if scoring and scoring != "vina":
        cmd += ["--scoring", scoring]
    if cpu:
        cmd += ["--cpu", str(cpu)]
    cmd += [str(a) for a in extra_args]
    return cmd


def dock(
    receptor: str | os.PathLike[str],
    ligand: str | os.PathLike[str],
    box: Box,
    out_pdbqt: str | os.PathLike[str] | None = None,
    name: str | None = None,
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    num_modes: int = DEFAULT_NUM_MODES,
    seed: int = DEFAULT_SEED,
    cpu: int | None = None,
    energy_range: float = 3.0,
    scoring: str = "vina",
    extra_args: Sequence[str] = (),
    timeout: float | None = 3600.0,
    overwrite: bool = True,
) -> DockingResult:
    """Dock one ligand into one receptor and parse the poses."""
    receptor_path = Path(receptor)
    ligand_path = Path(ligand)
    for path, label in ((receptor_path, "receptor"), (ligand_path, "ligand")):
        if not path.is_file():
            raise VinaError(f"No such {label} PDBQT: {path}")
        if path.stat().st_size == 0:
            raise VinaError(f"The {label} PDBQT {path} is empty.")

    ligand_name = name or ligand_path.stem
    out = Path(out_pdbqt) if out_pdbqt else ligand_path.with_name(f"{ligand_name}_out.pdbqt")
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.is_file() and not overwrite:
        poses = parse_vina_pdbqt(out, ligand_name=ligand_name)
        return DockingResult(ligand_name, poses, out, log="(cached)", box=box)

    cmd = build_command(
        receptor_path, ligand_path, box, out,
        exhaustiveness=exhaustiveness, num_modes=num_modes, seed=seed, cpu=cpu,
        energy_range=energy_range, scoring=scoring, extra_args=extra_args,
        executable=vina_executable(),
    )
    start = time.perf_counter()
    try:
        process = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise VinaError(
            f"Vina did not finish within {timeout:.0f} s for {ligand_name}. "
            "Lower --exhaustiveness, shrink the box, or raise --timeout."
        ) from exc
    seconds = time.perf_counter() - start

    if process.returncode != 0:
        tail = (process.stderr or process.stdout or "").strip().splitlines()[-8:]
        raise VinaError(
            f"Vina failed for {ligand_name} (exit code {process.returncode}):\n  "
            + "\n  ".join(tail)
        )
    if not out.is_file():
        raise VinaError(
            f"Vina reported success but wrote no output for {ligand_name}. "
            "Usually the ligand lies entirely outside the box."
        )
    poses = parse_vina_pdbqt(out, ligand_name=ligand_name)
    return DockingResult(
        ligand=ligand_name, poses=poses, out_pdbqt=out,
        log=process.stdout, seconds=seconds, command=cmd, box=box,
    )


def _dock_one(args: tuple) -> tuple[str, DockingResult | None, str | None]:
    name, receptor, ligand, box, out, kwargs = args
    try:
        result = dock(receptor, ligand, box, out_pdbqt=out, name=name, **kwargs)
        return name, result, None
    except Exception as exc:  # returned rather than raised: one bad ligand must not stop a screen
        return name, None, str(exc)


def dock_many(
    receptor: str | os.PathLike[str],
    ligands: Sequence[tuple[str, str | os.PathLike[str]]],
    box: Box,
    out_dir: str | os.PathLike[str] = "poses",
    n_workers: int = 1,
    progress: Callable[[str, int, int], None] | None = None,
    **kwargs,
) -> tuple[list[DockingResult], list[tuple[str, str]]]:
    """Dock a library. Returns (results, failures).

    Vina already uses every core for one ligand, so ``n_workers`` splits the
    cores between ligands instead: with 8 cores, ``n_workers=4`` runs 4 ligands
    on 2 cores each, which is faster overall for a library of small ligands.
    """
    if not ligands:
        raise VinaError("No ligands to dock.")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_workers = max(1, int(n_workers))
    if n_workers > 1 and "cpu" not in kwargs:
        cores = os.cpu_count() or 1
        kwargs["cpu"] = max(1, cores // n_workers)

    jobs = [
        (name, str(receptor), str(path), box, str(out_dir / f"{name}_out.pdbqt"), kwargs)
        for name, path in ligands
    ]
    results: list[DockingResult] = []
    failures: list[tuple[str, str]] = []
    done = 0
    total = len(jobs)

    if n_workers == 1:
        for job in jobs:
            name, result, error = _dock_one(job)
            done += 1
            if result is not None:
                results.append(result)
            else:
                failures.append((name, error or "unknown error"))
            if progress:
                progress(name, done, total)
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(_dock_one, job): job[0] for job in jobs}
            for future in as_completed(futures):
                name, result, error = future.result()
                done += 1
                if result is not None:
                    results.append(result)
                else:
                    failures.append((name, error or "unknown error"))
                if progress:
                    progress(name, done, total)

    results.sort(key=lambda r: r.best_affinity)
    return results, failures


def clear_cache(directory: str | os.PathLike[str]) -> int:
    """Delete cached pose files from a directory. Returns how many were removed."""
    path = Path(directory)
    if not path.is_dir():
        return 0
    removed = 0
    for file in path.glob("*_out.pdbqt"):
        file.unlink()
        removed += 1
    return removed
