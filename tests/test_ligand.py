import pytest

from docking_made_easy import Structure, prepare_ligand, prepare_ligands, read_ligands
from docking_made_easy.errors import LigandError
from docking_made_easy.ligand import embed_3d, ligand_from_pdb_block, mol_from_smiles, protonate

ERLOTINIB = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
GEFITINIB = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"


def test_a_good_smiles_parses():
    assert mol_from_smiles(ERLOTINIB).GetNumAtoms() == 29


@pytest.mark.parametrize("bad", ["", "   ", "not-a-molecule", "c1ccccc", "C(C"])
def test_a_bad_smiles_says_so(bad):
    with pytest.raises(LigandError):
        mol_from_smiles(bad)


def test_embedding_is_reproducible_for_a_given_seed():
    a = embed_3d(mol_from_smiles(ERLOTINIB), seed=7).GetConformer().GetPositions()
    b = embed_3d(mol_from_smiles(ERLOTINIB), seed=7).GetConformer().GetPositions()
    assert (a == b).all()


def test_protonation_charges_an_acid_and_a_base():
    from rdkit import Chem

    acid = protonate(mol_from_smiles("CC(=O)O"))
    assert Chem.GetFormalCharge(acid) == -1
    base = protonate(mol_from_smiles("CCCN"))
    assert Chem.GetFormalCharge(base) == 1
    neutral = protonate(mol_from_smiles("c1ccccc1"))
    assert Chem.GetFormalCharge(neutral) == 0


def test_prepare_ligand_writes_a_pdbqt_and_an_sdf(tmp_path):
    prepared = prepare_ligand(ERLOTINIB, name="erlotinib", out_dir=tmp_path)
    text = prepared.pdbqt.read_text()
    assert text.startswith("REMARK SMILES")
    assert "ROOT" in text and "TORSDOF" in text
    assert prepared.n_heavy_atoms == 29
    assert prepared.n_rotatable_bonds == 10
    assert prepared.sdf is not None and prepared.sdf.is_file()


def test_atom_count_in_the_pdbqt_matches_the_molecule(tmp_path):
    prepared = prepare_ligand(ERLOTINIB, name="erlotinib", out_dir=tmp_path)
    atoms = [ln for ln in prepared.pdbqt.read_text().splitlines() if ln.startswith("ATOM")]
    # Heavy atoms plus the polar hydrogens AutoDock keeps.
    assert len(atoms) >= 29


def test_reading_a_smi_file(tmp_path):
    path = tmp_path / "ligands.smi"
    path.write_text(f"{ERLOTINIB} erlotinib\n{GEFITINIB} gefitinib\n")
    items = list(read_ligands(path))
    assert [name for name, _, _ in items] == ["erlotinib", "gefitinib"]


def test_reading_a_csv_with_a_header(tmp_path):
    path = tmp_path / "ligands.csv"
    path.write_text(f"smiles,name\n{ERLOTINIB},erlotinib\n{GEFITINIB},gefitinib\n")
    assert [name for name, _, _ in read_ligands(path)] == ["erlotinib", "gefitinib"]


def test_reading_an_sdf(tmp_path):
    prepared = prepare_ligand(ERLOTINIB, name="erlotinib", out_dir=tmp_path)
    assert [name for name, _, _ in read_ligands(prepared.sdf)] == ["erlotinib"]


def test_an_empty_or_missing_ligand_file_is_an_error(tmp_path):
    empty = tmp_path / "empty.smi"
    empty.write_text("")
    with pytest.raises(LigandError, match="empty"):
        list(read_ligands(empty))
    with pytest.raises(LigandError, match="No such ligand file"):
        list(read_ligands(tmp_path / "nope.smi"))


def test_a_bad_ligand_is_skipped_not_fatal(tmp_path):
    path = tmp_path / "ligands.smi"
    path.write_text(f"{ERLOTINIB} erlotinib\nnot-a-molecule broken\n")
    prepared, failures = prepare_ligands(path, out_dir=tmp_path / "out")
    assert [p.name for p in prepared] == ["erlotinib"]
    assert failures[0][0] == "broken"


def test_all_ligands_failing_is_an_error(tmp_path):
    path = tmp_path / "ligands.smi"
    path.write_text("not-a-molecule broken\n")
    with pytest.raises(LigandError, match="No ligand"):
        prepare_ligands(path, out_dir=tmp_path / "out")


def test_duplicate_names_do_not_overwrite_each_other(tmp_path):
    path = tmp_path / "ligands.smi"
    path.write_text(f"{ERLOTINIB} lig\n{GEFITINIB} lig\n")
    prepared, _ = prepare_ligands(path, out_dir=tmp_path / "out")
    assert [p.name for p in prepared] == ["lig", "lig_2"]


def test_crystal_ligand_gets_its_bond_orders_from_the_template(mini_pdb):
    from rdkit import Chem

    ligand = Structure.from_file(mini_pdb).select_ligand("AQ4")
    mol = ligand_from_pdb_block(ligand.to_pdb(), template_smiles=ERLOTINIB)
    assert mol.GetNumAtoms() == 29
    assert any(b.GetBondType() == Chem.BondType.TRIPLE for b in mol.GetBonds())
    assert mol.GetConformer().GetPositions().shape == (29, 3)


def test_a_wrong_template_is_reported(mini_pdb):
    ligand = Structure.from_file(mini_pdb).select_ligand("AQ4")
    with pytest.raises(LigandError, match="does not match"):
        ligand_from_pdb_block(ligand.to_pdb(), template_smiles="CCO")
