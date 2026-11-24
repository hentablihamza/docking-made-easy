from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures() -> Path:
    return FIXTURES


@pytest.fixture
def mini_pdb(fixtures) -> Path:
    return fixtures / "mini_complex.pdb"


@pytest.fixture
def vina_out(fixtures) -> Path:
    return fixtures / "vina_out.pdbqt"
