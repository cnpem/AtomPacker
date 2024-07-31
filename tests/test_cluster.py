import os

import numpy
import pytest

from AtomPacker import Cage, Cavity, Cluster
from ase.cluster import FaceCenteredCubic
from ase.data import atomic_numbers, covalent_radii

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


@pytest.fixture  # 1 atom
def atom():
    return FaceCenteredCubic(
        symbols="Au",
        surfaces=[(1, 0, 0), (0, 1, 0), (0, 0, 1)],
        layers=[0, 0, 0],
        latticeconstant=4.08,
    )


@pytest.fixture  # 13 atoms
def cluster():
    return FaceCenteredCubic(
        symbols="Au",
        surfaces=[(1, 0, 0), (0, 1, 0), (0, 0, 1)],
        layers=[1, 1, 1],
        latticeconstant=4.08,
    )


def test_one_atom_cluster_distances(pdb, cavity, atom):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(atom, cavity)

    # Check distances for one-atom cluster
    assert cage.cluster._get_distances() == numpy.array([0])


def test_multi_atom_cluster_distances(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity)

    # Check distances for multi-atom cluster
    assert (
        cage.cluster._cluster.get_all_distances() == cage.cluster._get_distances()
    ).all()


def test_one_atom_cluster_radii(pdb, cavity, atom):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(atom, cavity)

    # Check radii for one-atom cluster
    assert cage.cluster._get_radii() == [covalent_radii[atomic_numbers["Au"]]]


def test_multi_atom_cluster_radii(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity)

    # Calculate radii
    distances = cage.cluster._cluster.get_all_distances()
    radii = distances[distances > 0].min() / 2

    # Check radii for multi-atom cluster
    assert cage.cluster._get_radii() == radii


def test_one_atom_cluster_coordinates(pdb, cavity, atom):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(atom, cavity)

    # Check coordinates
    assert (cage.cluster.coordinates == numpy.array([0, 0, 0])).all()


def test_multi_atom_cluster_coordinates(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity)

    # Check coordinates
    assert (cage.cluster.coordinates - cage.cluster._cluster.positions < 0.001).all()


def test_cluster_lattice_constants(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity)

    # Check lattice constants
    assert (cage.cluster.lattice_constants == 4.08).all()


def test_one_atom_cluster_number_of_atoms(pdb, cavity, atom):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(atom, cavity)

    # Check number of atoms
    assert cage.cluster.number_of_atoms == 1


def test_multi_atom_cluster_number_of_atoms(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity)

    # Check maximum number of atoms
    assert cage.cluster.number_of_atoms == 13


def test_one_atom_cluster_maximum_number_of_atoms(pdb, cavity, atom):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(atom, cavity)

    # Get maximum number of atoms
    maximum_number_of_atoms = numpy.ceil(
        numpy.prod(cage.cavity.grid.shape)
        * 0.6**3
        / (4 / 3 * numpy.pi * (cage.cluster._get_radii() ** 3)),
    ).astype(int)

    # Check maximum number of atoms
    assert cage.cluster.maximum_number_of_atoms == maximum_number_of_atoms


def test_multi_atom_cluster_maximum_number_of_atoms(pdb, cavity, cluster):
    cage = Cage()
    cage.load(pdb)

    # Dummy Cavity
    cage.cavity = cavity

    # Dummy Cluster
    cage.cluster = Cluster(cluster, cavity)

    # Get maximum number of atoms
    maximum_number_of_atoms = numpy.ceil(
        numpy.prod(cage.cavity.grid.shape)
        * 0.6**3
        / (4 / 3 * numpy.pi * (cage.cluster._get_radii() ** 3)),
    ).astype(int)

    # Check maximum number of atoms
    assert cage.cluster.maximum_number_of_atoms == maximum_number_of_atoms
