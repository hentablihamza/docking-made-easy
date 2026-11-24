import math

import numpy as np
import pytest

from docking_made_easy import (
    Box,
    Structure,
    box_from_center,
    box_from_ligand,
    box_from_points,
    box_from_residues,
    check_box,
    pocket_residues,
)
from docking_made_easy.box import GRID_SPACING
from docking_made_easy.errors import StructureError

CUBE = [
    [-5.0, -5.0, -5.0],
    [5.0, 5.0, 5.0],
    [0.0, 0.0, 0.0],
]


def test_box_from_points_centres_on_the_bounding_box():
    box = box_from_points(CUBE, padding=4.0, snap_to_grid=False)
    assert box.center == (0.0, 0.0, 0.0)
    assert box.size == (18.0, 18.0, 18.0)   # 10 A extent + 2 x 4 A padding


def test_padding_is_applied_on_both_sides():
    points = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
    box = box_from_points(points, padding=3.0, snap_to_grid=False, min_size=0.0)
    assert box.size_x == pytest.approx(16.0)
    assert box.center_x == pytest.approx(5.0)


def test_centroid_mode_differs_for_a_lopsided_ligand():
    points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.0, 0.0], [20.0, 0.0, 0.0]]
    extent = box_from_points(points, padding=2.0, mode="extent", snap_to_grid=False, min_size=0)
    centroid = box_from_points(points, padding=2.0, mode="centroid", snap_to_grid=False, min_size=0)
    assert extent.center_x == pytest.approx(10.0)
    assert centroid.center_x == pytest.approx(5.625)
    # Either way, every atom still has the requested clearance.
    assert centroid.fraction_inside(points) == 1.0
    assert extent.fraction_inside(points) == 1.0


def test_sizes_snap_up_to_the_vina_grid():
    box = box_from_points(CUBE, padding=4.1, snap_to_grid=True)
    for size in box.size:
        assert math.isclose(size / GRID_SPACING, round(size / GRID_SPACING), abs_tol=1e-9)
        assert size >= 18.2


def test_minimum_size_protects_a_tiny_ligand():
    box = box_from_points([[0, 0, 0], [1, 0, 0]], padding=1.0, min_size=12.0, snap_to_grid=False)
    assert box.size == (12.0, 12.0, 12.0)


def test_cubic_option_uses_the_longest_edge():
    box = box_from_points([[0, 0, 0], [20, 2, 2]], padding=2.0, cubic=True, snap_to_grid=False)
    assert box.size_x == box.size_y == box.size_z == pytest.approx(24.0)


def test_max_size_caps_the_box():
    box = box_from_points([[0, 0, 0], [60, 0, 0]], padding=4.0, max_size=30.0, snap_to_grid=False)
    assert max(box.size) == 30.0


def test_box_from_points_rejects_empty_input():
    with pytest.raises(StructureError):
        box_from_points(np.zeros((0, 3)))


def test_box_mode_must_be_known():
    with pytest.raises(StructureError, match="extent"):
        box_from_points(CUBE, mode="middle")


def test_manual_box_accepts_one_or_three_sizes():
    assert box_from_center([1, 2, 3], 20).size == (20.0, 20.0, 20.0)
    assert box_from_center([1, 2, 3], [20, 22, 24]).size == (20.0, 22.0, 24.0)
    with pytest.raises(StructureError):
        box_from_center([1, 2], 20)
    with pytest.raises(StructureError):
        box_from_center([1, 2, 3], -5)


def test_contains_and_fraction_inside():
    box = box_from_center([0, 0, 0], 10)
    points = [[0, 0, 0], [4, 4, 4], [20, 0, 0]]
    assert box.contains(points).tolist() == [True, True, False]
    assert box.fraction_inside(points) == pytest.approx(2 / 3)


def test_corners_match_the_centre_and_size():
    box = box_from_center([1, 2, 3], [10, 20, 30])
    assert box.min_corner == (-4.0, -8.0, -12.0)
    assert box.max_corner == (6.0, 12.0, 18.0)
    assert box.volume == pytest.approx(6000.0)


def test_warnings_flag_boxes_that_are_too_small_or_too_large():
    assert any("small" in w for w in box_from_center([0, 0, 0], 8).warnings())
    assert any("slow Vina down" in w for w in box_from_center([0, 0, 0], 40).warnings())
    assert box_from_center([0, 0, 0], 22.5).warnings() == []


def test_config_round_trip(tmp_path):
    box = box_from_center([1.2345, -2.5, 3.0], [20, 21, 22])
    path = box.write_config(tmp_path / "box.txt")
    text = path.read_text()
    assert "center_x = 1.234" in text or "center_x = 1.235" in text
    again = Box.from_config(path)
    assert again.center == pytest.approx(box.center, abs=1e-3)
    assert again.size == pytest.approx(box.size, abs=1e-3)
    assert again.source.startswith("config:")


def test_config_with_a_missing_key_is_an_error(tmp_path):
    path = tmp_path / "bad.txt"
    path.write_text("center_x = 1\ncenter_y = 2\n")
    with pytest.raises(StructureError, match="missing"):
        Box.from_config(path)


def test_json_and_pymol_output(tmp_path):
    box = box_from_center([0, 0, 0], 20)
    assert "center_x" in box.write_json(tmp_path / "box.json").read_text()
    script = box.to_pymol("mybox")
    assert "load_cgo" in script and "mybox" in script


# -- the real EGFR pocket ----------------------------------------------------


def test_box_from_the_1m17_co_crystal_ligand(mini_pdb):
    structure = Structure.from_file(mini_pdb)
    ligand = structure.select_ligand("AQ4")
    box = box_from_ligand(ligand, padding=4.0)

    # Erlotinib's centroid in 1M17, to three decimals.
    assert box.center == pytest.approx((21.697, 0.303, 52.093), abs=1e-3)
    assert box.size_x == pytest.approx(24.0, abs=0.4)
    assert box.fraction_inside(ligand.heavy_coords) == 1.0
    assert box.source == "ligand:AQ4:A:999"


def test_larger_padding_gives_a_larger_box(mini_pdb):
    ligand = Structure.from_file(mini_pdb).select_ligand("AQ4")
    small = box_from_ligand(ligand, padding=2.0)
    large = box_from_ligand(ligand, padding=8.0)
    assert large.volume > small.volume
    assert large.center == pytest.approx(small.center, abs=1e-6)


def test_check_box_reports_coverage_and_density(mini_pdb):
    structure = Structure.from_file(mini_pdb)
    ligand = structure.select_ligand("AQ4")
    report = check_box(box_from_ligand(ligand), ligand, structure.receptor_atoms())
    assert report["ligand_atoms_inside"] == 1.0
    assert report["receptor_heavy_atoms_in_box"] > 0


def test_check_box_warns_when_the_box_misses_the_ligand(mini_pdb):
    structure = Structure.from_file(mini_pdb)
    ligand = structure.select_ligand("AQ4")
    report = check_box(box_from_center([200, 200, 200], 20), ligand, structure.receptor_atoms())
    assert report["ligand_atoms_inside"] == 0.0
    assert any("does not contain every ligand atom" in w for w in report["warnings"])
    assert any("out in the solvent" in w for w in report["warnings"])


def test_pocket_residues_find_the_egfr_hinge(mini_pdb):
    structure = Structure.from_file(mini_pdb)
    ligand = structure.select_ligand("AQ4")
    residues = pocket_residues(structure.receptor_atoms(), ligand, cutoff=5.0)
    names = [name for name, _ in residues]
    # Met769 is the hinge residue erlotinib hydrogen-bonds to.
    assert "MET769:A" in names
    assert "THR766:A" in names and "LEU694:A" in names
    assert residues == sorted(residues, key=lambda item: item[1])


def test_box_from_pocket_residues_needs_no_ligand(mini_pdb):
    structure = Structure.from_file(mini_pdb)
    residues = [r for r in structure.polymer_residues if r.res_seq in (769, 766, 694, 820)]
    box = box_from_residues(residues, padding=4.0)
    ligand = structure.select_ligand("AQ4")
    # The pocket-residue box still covers most of the crystal ligand.
    assert box.fraction_inside(ligand.heavy_coords) > 0.5
    with pytest.raises(StructureError):
        box_from_residues([])
