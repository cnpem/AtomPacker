import os

import pytest

from AtomPacker import Cage

DATADIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")


@pytest.fixture
def pdb():
    return os.path.join(DATADIR, "ZOCXOH.pdb")


def test_pipeline_without_optimization(pdb):
    cage = Cage()
    cage.load(pdb)
    cage.detect_cavity(
        step=0.6, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0
    )
    cage.pack(
        atom_type="Au",
        lattice_type="fcc",
        a=4.08,
        b=None,
        c=None,
        clashing_tolerance=0.0,
        angles=[0.0],
        translations=[0.0],
        optsave=False,
        optdir=None,
    )

    assert cage.cluster is not None
    assert cage.cluster.number_of_atoms > 0


def test_pipeline_with_optimization(pdb):
    cage = Cage()
    cage.load(pdb)
    cage.detect_cavity(
        step=0.6, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0
    )
    cage.pack(
        atom_type="Au",
        lattice_type="fcc",
        a=4.08,
        b=None,
        c=None,
        clashing_tolerance=0.0,
        angles=[-45, 0, 45],
        translations=[0.0, 0.2],
        optsave=False,
        optdir=None,
    )

    assert cage.cluster is not None
    assert cage.cluster.number_of_atoms > 0


if __name__ == "__main__":
    pytest.main()
