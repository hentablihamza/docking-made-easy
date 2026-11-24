import pytest

from docking_made_easy import Structure, split_structure, write_pdb
from docking_made_easy.errors import StructureError
from docking_made_easy.structure import group_residues, parse_pdb_atoms

PDB_LINE = (
    "ATOM      1  N   ALA A  12      11.104  13.207  10.000  1.00 20.00           N"
)


def test_parse_basic_fields():
    atom = parse_pdb_atoms(PDB_LINE)[0]
    assert (atom.name, atom.res_name, atom.chain, atom.res_seq) == ("N", "ALA", "A", 12)
    assert atom.coord.tolist() == [11.104, 13.207, 10.0]
    assert atom.element == "N"
    assert not atom.is_hydrogen


def test_element_falls_back_to_the_atom_name():
    line = "ATOM      2  CB  ALA A  12      11.000  13.000  10.000  1.00 20.00"
    assert parse_pdb_atoms(line)[0].element == "C"


def test_alternate_locations_are_filtered():
    text = (
        "ATOM      1  N  AALA A  12      11.000  13.000  10.000  0.60 20.00           N\n"
        "ATOM      2  N  BALA A  12      11.500  13.500  10.500  0.40 20.00           N\n"
    )
    atoms = parse_pdb_atoms(text)
    assert len(atoms) == 1
    assert atoms[0].alt_loc == "A"


def test_only_the_first_model_is_kept():
    text = (
        "MODEL        1\n" + PDB_LINE + "\nENDMDL\n"
        "MODEL        2\n" + PDB_LINE + "\nENDMDL\n"
    )
    assert len(parse_pdb_atoms(text)) == 1


def test_empty_input_is_an_error():
    with pytest.raises(StructureError, match="No ATOM/HETATM"):
        parse_pdb_atoms("HEADER something\nEND\n")


def test_fixture_is_split_into_its_parts(mini_pdb):
    structure = Structure.from_file(mini_pdb)
    summary = structure.summary()
    assert summary["chains"] == ["A"]
    assert summary["n_waters"] == 2
    assert summary["ions"] == ["ZN"]
    assert [lig["res_name"] for lig in summary["ligands"]] == ["AQ4"]
    assert summary["ligands"][0]["n_heavy_atoms"] == 29


def test_select_ligand_picks_the_largest_hetero_group(mini_pdb):
    ligand = Structure.from_file(mini_pdb).select_ligand()
    assert ligand.res_name == "AQ4"
    assert ligand.key == "AQ4:A:999"
    assert ligand.n_heavy_atoms == 29


@pytest.mark.parametrize("selector", ["AQ4", "aq4:A", "AQ4:A:999"])
def test_select_ligand_accepts_several_selector_forms(mini_pdb, selector):
    assert Structure.from_file(mini_pdb).select_ligand(selector).res_name == "AQ4"


def test_unknown_selector_lists_the_candidates(mini_pdb):
    with pytest.raises(StructureError, match="AQ4:A:999"):
        Structure.from_file(mini_pdb).select_ligand("XYZ")


def test_receptor_excludes_ligand_water_and_ions_by_default(mini_pdb):
    structure = Structure.from_file(mini_pdb)
    receptor, ligand = split_structure(structure)
    names = {a.res_name for a in receptor}
    assert "AQ4" not in names and "HOH" not in names and "ZN" not in names
    assert ligand.res_name == "AQ4"
    kept = structure.receptor_atoms(keep_waters=True, keep_ions=True)
    assert len(kept) == len(receptor) + 2 + 1


def test_empty_chain_selection_is_an_error(mini_pdb):
    with pytest.raises(StructureError, match="selection is empty"):
        Structure.from_file(mini_pdb).receptor_atoms(chains=["Z"])


def test_write_pdb_round_trips(mini_pdb, tmp_path):
    structure = Structure.from_file(mini_pdb)
    out = write_pdb(structure.receptor_atoms(), tmp_path / "rec.pdb", remarks=["test"])
    text = out.read_text()
    assert text.startswith("REMARK   test")
    assert text.rstrip().endswith("END")
    assert len(Structure.from_file(out).atoms) == len(structure.receptor_atoms())


def test_group_residues_keeps_file_order(mini_pdb):
    residues = group_residues(Structure.from_file(mini_pdb).atoms)
    assert residues[-1].res_name == "ZN"
    assert len({r.key for r in residues}) == len(residues)


def test_fetch_rejects_a_bad_pdb_id():
    with pytest.raises(StructureError, match="not a 4-character PDB ID"):
        Structure.fetch("erlotinib")
