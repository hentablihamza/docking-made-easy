import json

import pytest

from docking_made_easy.cli import build_parser, main

ERLOTINIB = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"


def test_help_lists_every_command(capsys):
    with pytest.raises(SystemExit):
        main(["--help"])
    text = capsys.readouterr().out
    for command in ("info", "box", "prep-receptor", "prep-ligand", "dock", "parse", "redock", "doctor"):
        assert command in text


def test_version_flag():
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0


def test_info_on_the_fixture(mini_pdb, capsys):
    assert main(["info", str(mini_pdb)]) == 0
    out = capsys.readouterr().out
    assert "AQ4:A:999" in out
    assert "29 heavy atoms" in out


def test_info_as_json(mini_pdb, capsys):
    assert main(["info", str(mini_pdb), "--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["ligands"][0]["res_name"] == "AQ4"


def test_box_prints_a_vina_config(mini_pdb, capsys):
    assert main(["box", str(mini_pdb), "--ligand-res", "AQ4"]) == 0
    out = capsys.readouterr().out
    assert "center_x = 21.697" in out
    assert "ligand atoms inside: 100%" in out


def test_box_writes_config_json_and_pymol(mini_pdb, tmp_path, capsys):
    config = tmp_path / "box.txt"
    assert main([
        "box", str(mini_pdb), "-o", str(config),
        "--json-out", str(tmp_path / "box.json"),
        "--pymol", str(tmp_path / "box.pml"),
        "--pocket",
    ]) == 0
    assert "size_x" in config.read_text()
    assert json.loads((tmp_path / "box.json").read_text())["size_x"] > 0
    assert "load_cgo" in (tmp_path / "box.pml").read_text()
    assert "MET769" in capsys.readouterr().out


def test_box_with_a_manual_centre(mini_pdb, capsys):
    assert main(["box", str(mini_pdb), "--center", "1", "2", "3", "--size", "20"]) == 0
    out = capsys.readouterr().out
    assert "center_x = 1.000" in out and "size_z = 20.000" in out


def test_box_with_an_unknown_ligand_fails_cleanly(mini_pdb, capsys):
    assert main(["box", str(mini_pdb), "--ligand-res", "XYZ"]) == 1
    assert "error:" in capsys.readouterr().err


def test_prep_ligand_from_smiles(tmp_path, capsys):
    assert main(["prep-ligand", ERLOTINIB, "--name", "erlotinib", "-o", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "heavy atoms     : 29" in out
    assert (tmp_path / "erlotinib.pdbqt").is_file()


def test_prep_ligand_rejects_a_bad_smiles(capsys, tmp_path):
    assert main(["prep-ligand", "not-a-molecule", "-o", str(tmp_path)]) == 1
    assert "could not parse" in capsys.readouterr().err


def test_parse_prints_the_table(vina_out, tmp_path, capsys):
    csv = tmp_path / "poses.csv"
    assert main(["parse", str(vina_out), "--csv", str(csv), "--name", "erlotinib"]) == 0
    assert "-7.179" in capsys.readouterr().out
    assert "erlotinib" in csv.read_text()


def test_dock_without_a_box_is_an_error(capsys):
    assert main(["dock", "-r", "rec.pdbqt", "--smiles", "CCO"]) == 1
    assert "error:" in capsys.readouterr().err


def test_doctor_runs(capsys):
    assert main(["doctor"]) == 0
    assert "RDKit" in capsys.readouterr().out


def test_parser_defaults():
    args = build_parser().parse_args(["box", "1M17"])
    assert args.padding == 4.0
    assert args.box_mode == "extent"
