import os
import pytest
from AtomPacker import *

HERE = os.path.abspath(os.path.dirname(__file__))
REPLICATES = 10


# Define a fixture for different np_atom and np_atom_radius combinations
@pytest.fixture(
    params=[
        {"np_atom_file": os.path.join(HERE, "data", "Au.pdb"), "np_atom_radius": 1.36},
        {"np_atom_file": os.path.join(HERE, "data", "Pd.pdb"), "np_atom_radius": 1.39},
    ]
)
def np_atom_params(request):
    return request.param


# Define a test function to test AtomPacker with different np_atom and np_atom_radius
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_pipeline_1(np_atom_params):
    # Define a basename for the temporary directory
    basename = np_atom_params["np_atom_file"].split("/")[-1].replace(".pdb", "")

    # Extract parameters from the fixture
    np_atom_file = np_atom_params["np_atom_file"]
    np_atom_radius = np_atom_params["np_atom_radius"]

    # Load supramolecular cage into PackmolStructure object
    smc = PackmolStructure(
        os.path.join(HERE, "data", "C1.pdb"),
        number=1,
        instructions=["center", "fixed 0. 0. 0. 0. 0. 0."],
    )

    # Load nanoparticle atoms into PackmolStructure object (specified by np_atom_file)
    np_atom = PackmolStructure(
        np_atom_file,
        number=40,
        instructions=["inside sphere 0. 0. 0. 7."],
    )

    # Create a CavityDetector with detection parameters appropriate for the supramolecular cage
    cd = CavityDetector(
        step=0.25,
        probe_in=1.4,
        probe_out=10.0,
        removal_distance=1.0,
        volume_cutoff=5.0,
        vdw=None,
    )

    # Create the AtomPacker object with specified np_atom_radius
    ap = AtomPacker(
        smc,
        np_atom,
        np_atom_radius=np_atom_radius,
        cavity_detector=cd,
        basedir=os.path.join(HERE, "pipeline", basename),
    )

    # Run Packing algorithm
    ap.packing(replicates=REPLICATES)

    # Perform assertions to check if files were created for the specified np_atom and np_atom_radius
    assert os.path.exists(os.path.join(HERE, "pipeline", basename, "PackedAtoms.csv"))
    assert os.path.exists(os.path.join(HERE, "pipeline", basename, "cavity.pdb"))
    assert os.path.exists(os.path.join(HERE, "pipeline", basename, "packmol.inp"))
    assert os.path.exists(os.path.join(HERE, "pipeline", basename, "packmol.stdout"))
    for replicate in range(REPLICATES):
        assert os.path.exists(
            os.path.join(HERE, "pipeline", basename, f"packed{replicate}.pdb")
        )


# Run the pytest framework
if __name__ == "__main__":
    pytest.main([__file__])
