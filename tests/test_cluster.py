import os

import numpy
import pytest

from AtomPacker import Cage, Cavity, Cluster
from ase.cluster import FaceCenteredCubic

DATADIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")


@pytest.fixture
def pdb():
    return os.path.join(DATADIR, "ZOCXOH.pdb")


@pytest.fixture
def grid():
    return numpy.full((20, 20, 20), 2)


@pytest.fixture
def vertices():
    return numpy.array(
        [
            [-10 * 0.6, -10 * 0.6, -10 * 0.6],
            [10 * 0.6, 0.0, 0.0],
            [0.0, 10 * 0.6, 0.0],
            [0.0, 0.0, 10 * 0.6],
        ]
    )


@pytest.fixture
def cavity(grid, vertices):
    return Cavity(
        grid=grid,
        step=0.6,
        probe_in=1.4,
        probe_out=10.0,
        removal_distance=1.0,
        volume_cutoff=5.0,
        vertices=vertices,
        surface="SES",
    )


@pytest.fixture  # 13 atoms
def cluster():
    _cluster = FaceCenteredCubic(
        symbols="Au",
        surfaces=[(1, 0, 0), (0, 1, 0), (0, 0, 1)],
        layers=[1, 1, 1],
        latticeconstant=4.08,
    )
    # Add radius to ase.cluster
    distances = _cluster.get_all_distances()
    radii = (distances[distances > 0] / 2).min()
    _cluster.info.update({"radii": radii})
    # Add lattice_constants to ase.cluster
    _cluster.info.update({"lattice_constants": 4.08})
    return _cluster


def test_cluster_distances(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity, log=None)

    # Check distances for multi-atom cluster
    # Get distances from positions
    positions = cage.cluster._cluster.positions
    distances = numpy.linalg.norm(
        positions[:, None, :] - positions[None, :, :], axis=-1
    )
    assert (distances == cage.cluster._cluster.get_all_distances()).all()
    assert distances.shape == (len(cage.cluster._cluster), len(cage.cluster._cluster))
    assert (distances >= 0).all()


def test_cluster_radii(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity, log=None)

    # Calculate radii
    distances = cage.cluster._cluster.get_all_distances()
    radii = distances[distances > 0].min() / 2

    # Check radii for multi-atom cluster
    assert cage.cluster.radii == radii


def test_cluster_coordinates(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity, log=None)

    # Check coordinates
    assert (cage.cluster.coordinates - cage.cluster._cluster.positions < 0.001).all()


def test_cluster_lattice_constants(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity, log=None)
    print(cage.cluster.lattice_constants)

    # Check lattice constants
    assert cage.cluster.lattice_constants == 4.08


def test_cluster_number_of_atoms(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity, log=None)

    # Check maximum number of atoms
    assert cage.cluster.number_of_atoms == 13


def test_cluster_maximum_number_of_atoms(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity, log=None)

    # Get maximum number of atoms
    maximum_number_of_atoms = numpy.ceil(
        numpy.prod(cage.cavity.grid.shape)
        * 0.6**3
        / (4 / 3 * numpy.pi * (cage.cluster.radii**3)),
    ).astype(int)

    # Check maximum number of atoms
    assert cage.cluster.maximum_number_of_atoms == maximum_number_of_atoms
