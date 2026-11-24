"""docking-made-easy: AutoDock Vina docking without the setup pain.

Protein preparation, a docking box taken from the co-crystal ligand, ligand
preparation from SMILES, batch docking and tidy result tables.
"""

from .version import __version__
from .errors import (
    DockingError,
    LigandError,
    MissingToolError,
    ParseError,
    StructureError,
    VinaError,
)
from .structure import (
    Atom,
    Residue,
    Structure,
    fetch_pdb,
    split_structure,
    write_pdb,
)
from .box import (
    Box,
    box_from_center,
    box_from_ligand,
    box_from_points,
    box_from_residues,
    check_box,
    pocket_residues,
)
from .receptor import ReceptorPrepResult, add_hydrogens, prepare_receptor
from .ligand import PreparedLigand, prepare_ligand, prepare_ligands, read_ligands
from .dock import DockingResult, dock, dock_many, vina_available
from .results import (
    Pose,
    parse_vina_log,
    parse_vina_pdbqt,
    pose_to_pdb,
    poses_to_dataframe,
    poses_to_sdf,
    results_table,
    save_table,
    summary_table,
    write_complex,
)
from .validate import RedockResult, pose_rmsd, redock, rmsd_to_crystal
from .tools import check_tools

__all__ = [
    "__version__",
    "Atom", "Residue", "Structure", "fetch_pdb", "split_structure", "write_pdb",
    "Box", "box_from_center", "box_from_ligand", "box_from_points",
    "box_from_residues", "check_box", "pocket_residues",
    "ReceptorPrepResult", "add_hydrogens", "prepare_receptor",
    "PreparedLigand", "prepare_ligand", "prepare_ligands", "read_ligands",
    "DockingResult", "dock", "dock_many", "vina_available",
    "Pose", "parse_vina_pdbqt", "parse_vina_log", "poses_to_dataframe",
    "results_table", "summary_table", "save_table", "poses_to_sdf",
    "pose_to_pdb", "write_complex",
    "RedockResult", "redock", "rmsd_to_crystal", "pose_rmsd",
    "check_tools",
    "DockingError", "StructureError", "LigandError", "VinaError",
    "ParseError", "MissingToolError",
]
