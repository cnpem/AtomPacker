import os

import ase
import numpy
import pytest
from ase.cluster.cubic import FaceCenteredCubic
from ase.data import atomic_numbers, covalent_radii
from MDAnalysis import Universe

from AtomPacker import Cage, Cavity, Cluster
from AtomPacker.core.io.vdw.file import read_vdw

DATADIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")


# Define a fixture for different cage file formats
@pytest.fixture
def mmcif():
    return os.path.join(DATADIR, "ZOCXOH.cif")


@pytest.fixture
def mol2():
    return os.path.join(DATADIR, "ZOCXOH.mol2")


@pytest.fixture
def pdb():
    return os.path.join(DATADIR, "ZOCXOH.pdb")


@pytest.fixture
def xyz():
    return os.path.join(DATADIR, "ZOCXOH.xyz")


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


def test_inexistent_file():
    cage = Cage()
    with pytest.raises(FileNotFoundError):
        cage.load("inexistent_file.pdb")


def test_custom_vdw(pdb):
    cage = Cage()
    cage.load(pdb, vdw={"H": 0.91, "C": 1.66, "Si": 2.1, "O": 1.69})
    assert isinstance(cage.universe, Universe)


def test_custom_vdw_with_missing_atom(pdb):
    cage = Cage()
    with pytest.raises(ValueError):
        cage.load(pdb, vdw={"H": 0.91, "C": 1.66, "Si": 2.1})


def test_atomic(pdb):
    cage = Cage()
    cage.load(pdb)
    with open(pdb, "r") as f:
        for line in f.readlines():
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                # [resnum, chain, resname, element, x, y, z, radius]
                resnum = int(line[6:11])
                chain = line[21].strip()
                resname = line[17:20].strip()
                atom = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                radius = read_vdw()["GEN"][atom]
                break

    assert (
        (cage.atomic[0, 0] == resnum)
        & (cage.atomic[0, 1] == chain)
        & (cage.atomic[0, 2] == resname)
        & (cage.atomic[0, 3] == atom)
        & (round(cage.atomic[0, 4], 3) == x)
        & (round(cage.atomic[0, 5], 3) == y)
        & (round(cage.atomic[0, 6], 3) == z)
        & (cage.atomic[0, 7] == radius)
    )


def test_atomic_shape(pdb):
    cage = Cage()
    cage.load(pdb)
    # [resnum, chain, resname, element, x, y, z, radius]
    assert cage.atomic.shape == (228, 8)


def test_coordinates_shape(pdb):
    cage = Cage()
    cage.load(pdb)
    # [x, y, z]
    assert cage.coordinates.shape == (228, 3)


def test_centroid_shape(pdb):
    cage = Cage()
    cage.load(pdb)
    assert cage.centroid.shape == (3,)


def test_universe_for_unloaded_cage():
    cage = Cage()
    assert cage.universe is None


def test_atomic_for_unloaded_cage():
    cage = Cage()
    assert cage.atomic is None


def test_coordinates_for_unloaded_cage():
    cage = Cage()
    assert cage.coordinates is None


def test_centroid_for_unloaded_cage():
    cage = Cage()
    assert cage.centroid is None


def test_cavity_for_unloaded_cage():
    cage = Cage()
    assert cage.cavity is None


def test_cluster_for_unloaded_cage():
    cage = Cage()
    assert cage.cluster is None


# TODO: This test is based on an older box size calculation method.
#       It should be updated to reflect the current implementation.
# def test_get_cluster_layers(pdb):
#     # Load cage structure
#     cage = Cage()
#     cage.load(pdb)

#     # Start min max coordinates
#     xmin, ymin, zmin = float("inf"), float("inf"), float("inf")
#     xmax, ymax, zmax = float("-inf"), float("-inf"), float("-inf")

#     # Read pdb
#     with open(pdb, "r") as f:
#         for line in f.readlines():
#             if line.startswith("ATOM  ") or line.startswith("HETATM"):
#                 x = float(line[30:38])
#                 y = float(line[38:46])
#                 z = float(line[46:54])
#                 xmin = min(x, xmin)
#                 ymin = min(y, ymin)
#                 zmin = min(z, zmin)
#                 xmax = max(x, xmax)
#                 ymax = max(y, ymax)
#                 zmax = max(z, zmax)

#     # Get xrange, yrange, zrange
#     xrange = abs(xmax - xmin)
#     yrange = abs(ymax - ymin)
#     zrange = abs(zmax - zmin)

#     # Get van der Waals radii of atom
#     atomic_number = atomic_numbers["Au"]
#     atom_size = covalent_radii[atomic_number] * 2
#     layers = numpy.ceil([xrange, yrange, zrange] / atom_size).astype(int)

#     # Assert
#     assert (cage._get_cluster_layers("Au", factor=0.0) == layers).all()


def test_detect_cavity(pdb):
    cage = Cage()
    cage.load(pdb)

    # Detect cavity
    cage.detect_cavity(
        step=0.6, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0
    )

    assert isinstance(cage.cavity, Cavity)
    assert isinstance(cage.cavity.grid, numpy.ndarray)
    assert isinstance(cage.cavity.coordinates, numpy.ndarray)


def test_filter_clashing_atoms(pdb):
    cage = Cage()
    cage.load(pdb)

    # Detect cavity
    cage.detect_cavity(
        step=0.6, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0
    )

    # Cluster with atom clashing with cage
    cluster = FaceCenteredCubic(
        symbols="Au",
        surfaces=[(1, 0, 0), (0, 1, 0), (0, 0, 1)],
        layers=[0, 0, 0],
        latticeconstant=4.08,
    )
    # Force clash with first atom of the cage
    cluster.positions += cage.coordinates[0]

    # Add radii to cluster info (Au)
    # (a: 4.0874, radius: 1.441932148 Å)
    cluster.info.update({"radii": 1.441932148})

    # Filter clashing atoms
    cluster = cage._filter_clashing_atoms(cluster, clashing_tolerance=0.0)

    assert len(cluster) == 0


def test_filter_not_clashing_atoms(pdb):
    cage = Cage()
    cage.load(pdb)

    # Detect cavity
    cage.detect_cavity(
        step=0.6, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0
    )

    # Cluster with atom clashing with cage
    cluster = FaceCenteredCubic(
        symbols="Au",
        surfaces=[(1, 0, 0), (0, 1, 0), (0, 0, 1)],
        layers=[0, 0, 0],
        latticeconstant=4.08,
    )
    # Place on center of the cage
    cluster.positions += cage.centroid

    # Add radii to cluster info (Au)
    # (a: 4.0874, radius: 1.441932148 Å)
    cluster.info.update({"radii": 1.441932148})

    # Filter clashing atoms
    cluster = cage._filter_clashing_atoms(cluster, clashing_tolerance=0.0)

    assert len(cluster) == 1


def test_filter_outside_cavity(pdb):
    cage = Cage()
    cage.load(pdb)

    # Detect cavity
    cage.detect_cavity(
        step=0.6, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0
    )

    # Cluster outside the cavity
    cluster = FaceCenteredCubic(
        symbols="Au",
        surfaces=[(1, 0, 0), (0, 1, 0), (0, 0, 1)],
        layers=[0, 0, 0],
        latticeconstant=4.08,
    )
    # Place outside cavity
    cluster.positions += cage.cavity.coordinates.min(axis=0) - 1.0

    # Filter outside the cavity
    cluster = cage._filter_outside_cavity(cluster)

    assert len(cluster) == 0


def test_not_filter_inside_cavity(pdb):
    cage = Cage()
    cage.load(pdb)

    # Detect cavity
    cage.detect_cavity(
        step=0.6, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0
    )

    # Cluster outside the cavity
    cluster = FaceCenteredCubic(
        symbols="Au",
        surfaces=[(1, 0, 0), (0, 1, 0), (0, 0, 1)],
        layers=[0, 0, 0],
        latticeconstant=4.08,
    )
    # Place on center of the cage (inside the cavity)
    cluster.positions += cage.centroid

    # Filter outside the cavity
    cluster = cage._filter_outside_cavity(cluster)

    assert len(cluster) == 1


def test_build_cluster(pdb):
    cage = Cage()
    cage.load(pdb)

    # Detect cavity
    cage.detect_cavity(
        step=0.6, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0
    )

    # Build cluster
    _cluster, log = cage._build_cluster(
        atom_type="Au",
        lattice_type="fcc",
        lattice_constants=4.08,
        center=cage.centroid,
        clashing_tolerance=0.0,
        angles=[0.0],
        translations=[0.0],
        optsave=False,
        optdir=None,
    )
    assert isinstance(_cluster, ase.cluster.Cluster)

    # Convert ase.cluster.Cluster to AtomPacker.Cluster
    cage.cluster = Cluster(_cluster, cavity=cage.cavity, log=log)
    assert isinstance(cage.cluster, Cluster)

    # Check number of atoms in AtomPacker.Cluster is the same as in
    # ase.cluster.Cluster
    assert cage.cluster.number_of_atoms == len(_cluster)
    # Check number of atoms in AtomPacker.Cluster is 17 (known value)
    assert cage.cluster.number_of_atoms == 17


def test_pack(pdb):
    cage = Cage()
    cage.load(pdb)

    # Detect cavity
    cage.detect_cavity(
        step=0.6, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0
    )

    # Pack cluster
    cage.pack(
        "Au",
        "fcc",
        4.08,
        clashing_tolerance=0.0,
        angles=[0.0],
        translations=[0.0],
        optsave=False,
        optdir=None,
    )

    assert isinstance(cage.cluster, Cluster)
    # Check number of atoms in AtomPacker.Cluster is 17 (known value)
    assert cage.cluster.number_of_atoms == 17


def test_pack_with_invalid_lattice_type(pdb):
    cage = Cage()
    cage.load(pdb)

    with pytest.raises(ValueError):
        cage.pack(
            "Au",
            "invalid",
            4.08,
            clashing_tolerance=0.0,
            angles=[0.0],
            translations=[0.0],
            optsave=False,
            optdir=None,
        )


def test_pack_with_negative_clashing_tolerance(pdb):
    cage = Cage()
    cage.load(pdb)

    with pytest.raises(ValueError):
        cage.pack(
            "Au",
            "fcc",
            4.08,
            clashing_tolerance=-1.0,
            angles=[0.0],
            translations=[0.0],
            optsave=False,
            optdir=None,
        )


def test_pack_for_unloaded_cage():
    cage = Cage()
    with pytest.raises(ValueError):
        cage.pack(
            "Au",
            "fcc",
            4.08,
            clashing_tolerance=0.0,
            angles=[0.0],
            translations=[0.0],
            optsave=False,
            optdir=None,
        )


def test_pack_without_cavity(pdb):
    cage = Cage()
    cage.load(pdb)

    with pytest.raises(ValueError):
        cage.pack(
            "Au",
            "fcc",
            4.08,
            clashing_tolerance=0.0,
            angles=[0.0],
            translations=[0.0],
            optsave=False,
            optdir=None,
        )


if __name__ == "__main__":
    pytest.main()
