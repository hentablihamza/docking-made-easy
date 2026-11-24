"""Redocking checks. The RMSD maths runs everywhere; the docking itself needs Vina."""

import pytest

from docking_made_easy import Structure, pose_rmsd, rmsd_to_crystal
from docking_made_easy.errors import LigandError
from docking_made_easy.ligand import ligand_from_pdb_block
from docking_made_easy.results import poses_to_rdkit

ERLOTINIB = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"


@pytest.fixture
def crystal(mini_pdb):
    ligand = Structure.from_file(mini_pdb).select_ligand("AQ4")
    return ligand_from_pdb_block(ligand.to_pdb(), template_smiles=ERLOTINIB)


def test_a_pose_has_zero_rmsd_to_itself(crystal):
    assert pose_rmsd(crystal, crystal) == pytest.approx(0.0, abs=1e-6)


def test_rmsd_detects_a_translation(crystal):
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    moved = Chem.Mol(crystal)
    conf = moved.GetConformer()
    for i in range(moved.GetNumAtoms()):
        position = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(position.x + 2.0, position.y, position.z))
    assert pose_rmsd(moved, crystal) == pytest.approx(2.0, abs=1e-6)


def test_rmsd_needs_matching_molecules(crystal):
    from docking_made_easy.ligand import embed_3d, mol_from_smiles

    other = embed_3d(mol_from_smiles("CCO"))
    with pytest.raises(LigandError, match="different heavy-atom counts"):
        pose_rmsd(other, crystal)


def test_poses_are_rebuilt_from_the_vina_output(vina_out):
    mol = poses_to_rdkit(vina_out)
    from rdkit import Chem

    assert Chem.RemoveHs(mol).GetNumAtoms() == 29
    assert mol.GetNumConformers() == 3


def test_redocked_poses_are_scored_against_the_crystal_pose(vina_out, crystal):
    """The committed fixture is the real 1M17 redocking run.

    Vina's top-ranked pose sits 1.37 A from the crystal pose, under the 2.0 A
    threshold, so the setup reproduces the experiment.
    """
    values = rmsd_to_crystal(vina_out, crystal)
    assert len(values) == 3
    assert values[0] == pytest.approx(1.371, abs=0.05)
    assert values[0] < 2.0
    assert all(v > values[0] for v in values[1:])


def test_symmetry_correction_never_increases_the_rmsd(vina_out, crystal):
    with_symmetry = rmsd_to_crystal(vina_out, crystal, symmetry=True)
    without = rmsd_to_crystal(vina_out, crystal, symmetry=False)
    assert all(a <= b + 1e-6 for a, b in zip(with_symmetry, without))
