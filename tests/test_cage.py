import os

import pytest
from MDAnalysis import Universe

from AtomPacker import Cage, Cavity, Cluster

DATADIR = os.path.abspath(os.path.dirname(__file__))


# Define a fixture for different cage file formats
@pytest.fixture
def mmcif():
    return os.path.join(DATADIR, "data", "ZOCXOH.pdb")


@pytest.fixture
def mol2():
    return os.path.join(DATADIR, "data", "ZOCXOH.pdb")


@pytest.fixture
def pdb():
    return os.path.join(DATADIR, "data", "ZOCXOH.pdb")


@pytest.fixture
def xyz():
    return os.path.join(DATADIR, "data", "ZOCXOH.pdb")


def test_mmcif(mmcif):
    cage = Cage()
    cage.load(mmcif)
    assert isinstance(cage.universe, Universe)


def test_mol2(mol2):
    cage = Cage()
    cage.load(mol2)
    assert isinstance(cage.universe, Universe)


def test_pdb(pdb):
    cage = Cage()
    cage.load(pdb)
    assert isinstance(cage.universe, Universe)


def test_xyz(xyz):
    cage = Cage()
    cage.load(xyz)
    assert isinstance(cage.universe, Universe)


def test_atomic(pdb):
    cage = Cage()
    cage.load(pdb)
    # [resnum, chain, resname, element, x, y, z, radius]
    assert cage.atomic.shape == (228, 8)


def test_coordinates(pdb):
    cage = Cage()
    cage.load(pdb)
    # [x, y, z]
    assert cage.coordinates.shape == (228, 3)


def test_universe_for_unloaded_cage():
    cage = Cage()
    assert cage.universe is None


def test_atomic_for_unloaded_cage():
    cage = Cage()
    assert cage.atomic is None


def test_coordinates_for_unloaded_cage():
    cage = Cage()
    assert cage.coordinates is None


if __name__ == "__main__":
    pytest.main()
