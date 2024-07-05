import os

import numpy
import pytest
import warnings

from AtomPacker import Cage, Cavity

DATADIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")


@pytest.fixture
def pdb():
    return os.path.join(DATADIR, "ZOCXOH.pdb")


def test_cavity_volume(pdb):
    cage = Cage()
    cage.load(pdb)

    # Create dummy cavity
    grid = numpy.full((3, 3, 3), 0)

    # Create dummy vertices
    daxis = 3 * 0.6
    vertices = numpy.array(
        [
            cage.centroid,
            cage.centroid + [daxis, 0, 0],
            cage.centroid + [0, daxis, 0],
            cage.centroid + [0, 0, daxis],
        ]
    )
    # Create grid
    grid[0, 0, 0] = 2

    # Dummy Cavity
    cage.cavity = Cavity(
        grid=grid,
        step=0.6,
        probe_in=1.4,
        probe_out=10.0,
        removal_distance=1.0,
        volume_cutoff=5.0,
        vertices=vertices,
        surface="SES",
    )

    # Check volume
    assert cage.cavity.volume == round(cage.cavity._step**3, 2)


def test_cavity_coordinates(pdb):
    cage = Cage()
    cage.load(pdb)

    # Create dummy cavity
    grid = numpy.full((3, 3, 3), 0)

    # Create dummy vertices
    daxis = 3 * 0.6
    vertices = numpy.array(
        [
            cage.centroid,
            cage.centroid + [daxis, 0, 0],
            cage.centroid + [0, daxis, 0],
            cage.centroid + [0, 0, daxis],
        ]
    )
    # Create grid
    grid[0, 0, 0] = 2

    # Dummy Cavity
    cage.cavity = Cavity(
        grid=grid,
        step=0.6,
        probe_in=1.4,
        probe_out=10.0,
        removal_distance=1.0,
        volume_cutoff=5.0,
        vertices=vertices,
        surface="SES",
    )

    # Check coordinates
    assert cage.cavity.coordinates.shape == (1, 3)
    assert (cage.cavity.coordinates == cage.centroid).all()


def test_cavity_universe(pdb):
    cage = Cage()
    cage.load(pdb)

    # Create dummy cavity
    grid = numpy.full((3, 3, 3), 0)

    # Create dummy vertices
    daxis = 3 * 0.6
    vertices = numpy.array(
        [
            cage.centroid,
            cage.centroid + [daxis, 0, 0],
            cage.centroid + [0, daxis, 0],
            cage.centroid + [0, 0, daxis],
        ]
    )
    # Create grid
    grid[0, 0, 0] = 2

    # Dummy Cavity
    cage.cavity = Cavity(
        grid=grid,
        step=0.6,
        probe_in=1.4,
        probe_out=10.0,
        removal_distance=1.0,
        volume_cutoff=5.0,
        vertices=vertices,
        surface="SES",
    )

    # Check universe
    assert len(cage.cavity.universe.atoms) == 1
    assert (cage.cavity.universe.atoms.positions[0] == cage.centroid).all()


def test_more_than_one_cavity(pdb):
    cage = Cage()
    cage.load(pdb)

    # Create dummy cavity
    grid = numpy.full((3, 3, 3), 0)

    # Create dummy vertices
    daxis = 3 * 0.6
    vertices = numpy.array(
        [
            cage.centroid,
            cage.centroid + [daxis, 0, 0],
            cage.centroid + [0, daxis, 0],
            cage.centroid + [0, 0, daxis],
        ]
    )
    # Create grid
    grid[0, 0, 0] = 2
    grid[1, 1, 1] = 3

    # Dummy Cavity
    cage.cavity = Cavity(
        grid=grid,
        step=0.6,
        probe_in=1.4,
        probe_out=10.0,
        removal_distance=1.0,
        volume_cutoff=5.0,
        vertices=vertices,
        surface="SES",
    )

    with pytest.warns(UserWarning):
        cage.cavity.volume


def test_select_cavity(pdb):
    cage = Cage()
    cage.load(pdb)

    # Create dummy cavity
    grid = numpy.full((3, 3, 3), 0)

    # Create dummy vertices
    daxis = 3 * 0.6
    vertices = numpy.array(
        [
            cage.centroid,
            cage.centroid + [daxis, 0, 0],
            cage.centroid + [0, daxis, 0],
            cage.centroid + [0, 0, daxis],
        ]
    )
    # Create grid
    grid[0, 0, 0] = 2
    grid[1, 1, 1] = 3

    # Dummy Cavity
    cage.cavity = Cavity(
        grid=grid,
        step=0.6,
        probe_in=1.4,
        probe_out=10.0,
        removal_distance=1.0,
        volume_cutoff=5.0,
        vertices=vertices,
        surface="SES",
    )

    # Select cavity
    cage.cavity.select_cavity([2])

    # Check volume
    assert cage.cavity.volume == round(cage.cavity._step**3, 2)

    # Check coordinates
    assert cage.cavity.coordinates.shape == (1, 3)
    assert (cage.cavity.coordinates == cage.centroid).all()


def test_select_cavity_invalid_indices(pdb):
    cage = Cage()
    cage.load(pdb)

    # Create dummy cavity
    grid = numpy.full((3, 3, 3), 0)

    # Create dummy vertices
    daxis = 3 * 0.6
    vertices = numpy.array(
        [
            cage.centroid,
            cage.centroid + [daxis, 0, 0],
            cage.centroid + [0, daxis, 0],
            cage.centroid + [0, 0, daxis],
        ]
    )
    # Create grid
    grid[0, 0, 0] = 2
    grid[1, 1, 1] = 3

    # Dummy Cavity
    cage.cavity = Cavity(
        grid=grid,
        step=0.6,
        probe_in=1.4,
        probe_out=10.0,
        removal_distance=1.0,
        volume_cutoff=5.0,
        vertices=vertices,
        surface="SES",
    )

    # Select cavity
    with pytest.raises(TypeError):
        cage.cavity.select_cavity(2)


if __name__ == "__main__":
    pytest.main()
