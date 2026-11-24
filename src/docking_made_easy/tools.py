"""Find the external programs docking needs, and say what to do when they are missing."""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass

from .errors import MissingToolError

INSTALL_HINTS = {
    "vina": (
        "Install AutoDock Vina with `conda install -c conda-forge vina` "
        "(or `pip install vina`), or download the binary from "
        "https://github.com/ccsb-scripps/AutoDock-Vina/releases and put it on your PATH. "
        "Set DME_VINA to the executable if it has another name."
    ),
    "mk_prepare_receptor.py": (
        "Install Meeko with `pip install meeko` (version 0.6 or newer provides "
        "mk_prepare_receptor.py)."
    ),
    "pdbfixer": (
        "Install PDBFixer with `conda install -c conda-forge pdbfixer`. "
        "It adds missing heavy atoms and hydrogens before receptor typing."
    ),
    "obabel": (
        "Install Open Babel with `conda install -c conda-forge openbabel`. "
        "It is the fallback route for receptor PDBQT files."
    ),
    "meeko": "Install Meeko with `pip install meeko`.",
    "rdkit": "Install RDKit with `conda install -c conda-forge rdkit` or `pip install rdkit`.",
}


@dataclass(frozen=True)
class ToolStatus:
    name: str
    available: bool
    path: str | None = None
    version: str | None = None
    hint: str = ""


def which(program: str, env_var: str | None = None) -> str | None:
    """Locate ``program``, honouring an override environment variable."""
    if env_var:
        override = os.environ.get(env_var)
        if override:
            return override if shutil.which(override) or os.path.isfile(override) else None
    return shutil.which(program)


def require(program: str, env_var: str | None = None) -> str:
    path = which(program, env_var=env_var)
    if not path:
        raise MissingToolError(program, INSTALL_HINTS.get(program, ""))
    return path


def have_python_module(name: str) -> bool:
    import importlib.util

    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):  # pragma: no cover - odd import states
        return False


def _version(cmd: list[str]) -> str | None:
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    except Exception:  # pragma: no cover - environment dependent
        return None
    text = (result.stdout or result.stderr or "").strip().splitlines()
    return text[0].strip() if text else None


def check_tools() -> list[ToolStatus]:
    """What is installed right now. Drives the ``dme doctor`` command."""
    statuses: list[ToolStatus] = []

    vina = which("vina", env_var="DME_VINA")
    statuses.append(
        ToolStatus(
            "AutoDock Vina (vina)",
            vina is not None,
            vina,
            _version([vina, "--version"]) if vina else None,
            INSTALL_HINTS["vina"],
        )
    )

    for module, label in (("rdkit", "RDKit"), ("meeko", "Meeko"), ("pdbfixer", "PDBFixer")):
        ok = have_python_module(module)
        version = None
        if ok:
            try:
                version = __import__(module).__version__
            except Exception:
                version = "installed"
        statuses.append(
            ToolStatus(label, ok, None, version, INSTALL_HINTS.get(module, ""))
        )

    receptor_cli = which("mk_prepare_receptor.py")
    statuses.append(
        ToolStatus(
            "Meeko receptor prep (mk_prepare_receptor.py)",
            receptor_cli is not None,
            receptor_cli,
            None,
            INSTALL_HINTS["mk_prepare_receptor.py"],
        )
    )

    obabel = which("obabel")
    statuses.append(
        ToolStatus(
            "Open Babel (obabel, fallback)",
            obabel is not None,
            obabel,
            _version([obabel, "-V"]) if obabel else None,
            INSTALL_HINTS["obabel"],
        )
    )
    return statuses
