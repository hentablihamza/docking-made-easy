"""Exception types with actionable messages."""

from __future__ import annotations


class DockingError(Exception):
    """Base class for every docking-made-easy error."""


class MissingToolError(DockingError):
    """A required external program or Python package is not installed."""

    def __init__(self, tool: str, hint: str = "") -> None:
        message = f"{tool} is not available."
        if hint:
            message += f" {hint}"
        super().__init__(message)
        self.tool = tool


class StructureError(DockingError):
    """A PDB/SDF structure could not be read, or a selection is empty."""


class LigandError(DockingError):
    """A ligand could not be built from SMILES/SDF or converted to PDBQT."""


class VinaError(DockingError):
    """AutoDock Vina failed to run or returned no poses."""


class ParseError(DockingError):
    """An output file could not be parsed."""
