import pytest

from docking_made_easy import box_from_center, dock, dock_many
from docking_made_easy.dock import build_command, clear_cache, vina_available
from docking_made_easy.errors import MissingToolError, VinaError
from docking_made_easy.tools import check_tools, require, which

BOX = box_from_center([21.697, 0.303, 52.093], [24.0, 13.0, 15.5])


def test_command_contains_the_box_and_the_settings():
    cmd = build_command("rec.pdbqt", "lig.pdbqt", BOX, "out.pdbqt",
                        exhaustiveness=16, num_modes=5, seed=7, cpu=4)
    joined = " ".join(cmd)
    assert "--center_x 21.697" in joined
    assert "--size_y 13.000" in joined
    assert "--exhaustiveness 16" in joined
    assert "--num_modes 5" in joined
    assert "--seed 7" in joined
    assert "--cpu 4" in joined


def test_command_omits_cpu_when_not_set():
    assert "--cpu" not in " ".join(build_command("r", "l", BOX, "o"))


def test_extra_arguments_are_passed_through():
    cmd = build_command("r", "l", BOX, "o", extra_args=["--min_rmsd", "1.5"])
    assert cmd[-2:] == ["--min_rmsd", "1.5"]


@pytest.mark.parametrize("kwargs", [{"exhaustiveness": 0}, {"num_modes": 0}])
def test_nonsense_settings_are_rejected(kwargs):
    with pytest.raises(VinaError):
        build_command("r", "l", BOX, "o", **kwargs)


def test_missing_input_files_give_a_clear_error(tmp_path):
    receptor = tmp_path / "rec.pdbqt"
    receptor.write_text("ATOM\n")
    with pytest.raises(VinaError, match="No such ligand PDBQT"):
        dock(receptor, tmp_path / "missing.pdbqt", BOX)
    with pytest.raises(VinaError, match="No such receptor PDBQT"):
        dock(tmp_path / "missing.pdbqt", receptor, BOX)


def test_empty_input_file_is_caught(tmp_path):
    receptor = tmp_path / "rec.pdbqt"
    receptor.write_text("ATOM\n")
    empty = tmp_path / "lig.pdbqt"
    empty.write_text("")
    with pytest.raises(VinaError, match="is empty"):
        dock(receptor, empty, BOX)


def test_docking_an_empty_library_is_an_error():
    with pytest.raises(VinaError, match="No ligands"):
        dock_many("rec.pdbqt", [], BOX)


def test_missing_tool_message_names_the_fix():
    with pytest.raises(MissingToolError) as excinfo:
        require("a-program-that-does-not-exist")
    assert "not available" in str(excinfo.value)


def test_which_honours_the_override_variable(monkeypatch, tmp_path):
    fake = tmp_path / "vina"
    fake.write_text("#!/bin/sh\n")
    fake.chmod(0o755)
    monkeypatch.setenv("DME_VINA", str(fake))
    assert which("vina", env_var="DME_VINA") == str(fake)


def test_doctor_reports_every_tool():
    names = [status.name for status in check_tools()]
    assert any("Vina" in n for n in names)
    assert any("RDKit" in n for n in names)
    assert all(isinstance(s.available, bool) for s in check_tools())


def test_vina_available_is_a_boolean():
    assert isinstance(vina_available(), bool)


def test_clear_cache_removes_pose_files(tmp_path):
    (tmp_path / "a_out.pdbqt").write_text("x")
    (tmp_path / "keep.txt").write_text("x")
    assert clear_cache(tmp_path) == 1
    assert (tmp_path / "keep.txt").exists()
    assert clear_cache(tmp_path / "nope") == 0
