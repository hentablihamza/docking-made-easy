import pytest

from docking_made_easy import (
    parse_vina_log,
    parse_vina_pdbqt,
    poses_to_dataframe,
    results_table,
    save_table,
    summary_table,
)
from docking_made_easy.errors import ParseError
from docking_made_easy.results import Pose, write_complex

LOG = """
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -7.179      0.000      0.000
   2       -7.089      2.822      7.841
"""


def test_parse_a_real_vina_output_file(vina_out):
    poses = parse_vina_pdbqt(vina_out)
    assert [p.mode for p in poses] == [1, 2, 3]
    assert poses[0].affinity == pytest.approx(-7.179)
    assert poses[1].rmsd_lb == pytest.approx(2.822)
    assert poses[2].rmsd_ub == pytest.approx(8.554)
    assert all(p.pdbqt_block.startswith("MODEL") for p in poses)
    assert all("ENDMDL" in p.pdbqt_block for p in poses)


def test_affinities_are_ranked_best_first(vina_out):
    affinities = [p.affinity for p in parse_vina_pdbqt(vina_out)]
    assert affinities == sorted(affinities)


def test_ligand_name_can_be_overridden(vina_out):
    assert parse_vina_pdbqt(vina_out, ligand_name="erlotinib")[0].ligand == "erlotinib"


def test_parsing_a_file_without_results_is_an_error(tmp_path):
    path = tmp_path / "ligand.pdbqt"
    path.write_text("ROOT\nATOM      1  C   UNL     1       0.000   0.000   0.000\nENDROOT\n")
    with pytest.raises(ParseError, match="no 'REMARK VINA RESULT'"):
        parse_vina_pdbqt(path)


def test_missing_file_is_an_error(tmp_path):
    with pytest.raises(ParseError, match="No such Vina output"):
        parse_vina_pdbqt(tmp_path / "nope.pdbqt")


def test_parse_the_stdout_table():
    poses = parse_vina_log(LOG, ligand_name="erlotinib")
    assert len(poses) == 2
    assert poses[0].affinity == pytest.approx(-7.179)
    assert poses[0].ligand == "erlotinib"
    with pytest.raises(ParseError):
        parse_vina_log("nothing useful here")


def test_dataframe_columns_and_order(vina_out):
    frame = poses_to_dataframe(parse_vina_pdbqt(vina_out))
    assert list(frame.columns) == ["ligand", "mode", "affinity", "rmsd_lb", "rmsd_ub"]
    assert len(frame) == 3


class _Result:
    def __init__(self, poses):
        self.poses = poses


def _poses(name, affinities):
    return [
        Pose(ligand=name, mode=i + 1, affinity=a, rmsd_lb=0.0, rmsd_ub=0.0)
        for i, a in enumerate(affinities)
    ]


def test_results_table_sorts_the_best_ligand_first():
    weak = _Result(_poses("weak", [-5.0, -4.5]))
    strong = _Result(_poses("strong", [-9.1, -8.0]))
    table = results_table([weak, strong])
    assert table.iloc[0]["ligand"] == "strong"
    assert list(table["mode"]) == [1, 2, 1, 2]


def test_summary_table_has_one_row_per_ligand():
    summary = summary_table([_Result(_poses("a", [-9.1, -8.0])), _Result(_poses("b", [-5.0]))])
    assert list(summary["ligand"]) == ["a", "b"]
    assert summary.iloc[0]["best_affinity"] == pytest.approx(-9.1)
    assert summary.iloc[0]["n_modes"] == 2


def test_save_table_formats(tmp_path, vina_out):
    frame = poses_to_dataframe(parse_vina_pdbqt(vina_out))
    csv = save_table(frame, tmp_path / "poses.csv")
    tsv = save_table(frame, tmp_path / "poses.tsv")
    assert "affinity" in csv.read_text().splitlines()[0]
    assert "\t" in tsv.read_text().splitlines()[0]
    with pytest.raises(ParseError, match="Unknown table format"):
        save_table(frame, tmp_path / "poses.xyz")


def test_write_complex_merges_receptor_and_pose(tmp_path):
    receptor = tmp_path / "rec.pdb"
    receptor.write_text(
        "ATOM      1  N   ALA A  12      11.000  13.000  10.000  1.00 20.00           N\nEND\n"
    )
    pose = tmp_path / "pose.pdb"
    pose.write_text(
        "ATOM      1  C   UNL A   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n"
    )
    text = write_complex(receptor, pose, tmp_path / "complex.pdb").read_text()
    assert text.count("ATOM") == 1
    assert text.count("HETATM") == 1
    assert "TER" in text
