# %%
import os
import pandas

from AtomPacker import Cage

# Create directories to store the results
if not os.path.exists("results"):
    os.makedirs("results")

RUNCONFIG = [
    # Co: FCC a: 3.5445 Å  # Ok
    ["PCC-2a.mol2", "Co", "fcc", 3.5445, None, None, (0.6, 4, 12, 1.8, 200)],
    # Co: HCP a: 2.5071 Å, c: 4.0694 Å  # Ok
    ["PCC-2a.mol2", "Co", "hcp", 2.5071, None, 4.0694, (0.6, 4, 12, 1.8, 200)],
    # Pt: FCC a: 3.9233 Å  # Ok
    ["CIAC-121.mol2", "Pt", "fcc", 3.9233, None, None, (0.6, 3.3, 12, 1.8, 200)],
    # Pd: FCC a: 3.8901 Å  # Ok
    ["MTC1.mol2", "Pd", "fcc", 3.8901, None, None, (0.6, 1.4, 12, 2.4, 200)],
    # Pt: FCC a: 3.9233 Å  # Ok
    ["Cu-MOC.mol2", "Pt", "fcc", 3.9233, None, None, (0.6, 1.4, 12, 3.6, 200)],
    # Ru: HCP a: 2.7053 Å, c: 4.2814 Å  # Ok
    ["Cu-MOC.mol2", "Ru", "hcp", 2.7053, None, 4.2814, (0.6, 1.4, 12, 3.6, 200)],
    # Pd: FCC a: 3.8901 Å  # Ok
    ["PIC-T.xyz", "Pd", "fcc", 3.8901, None, None, (0.6, 1.4, 12, 1.2, 200)],
    # Ag: FCC a: 4.0861 Å  # Ok | probe_out can be 50
    ["OB4r.xyz", "Ag", "fcc", 4.0861, None, None, (0.6, 1.4, 12, 2.4, 200)],
    # Pd: FCC a: 3.8901 Å   # Ok
    ["Phos-cage.mol2", "Pd", "fcc", 3.8901, None, None, (0.6, 1.4, 12, 1.8, 5)],
    # Au: FCC a: 4.0784 Å  # Ok
    ["PIC-1.xyz", "Au", "fcc", 4.0784, None, None, (0.6, 1.4, 12, 1.2, 200)],
    # Au: FCC a: 4.0784 Å  # Ok
    ["PIC-2.xyz", "Au", "fcc", 4.0784, None, None, (0.6, 4, 18, 1.8, 200)],
    # Au: FCC a: 4.0784 Å (new cage)
    ["PIC-3.xyz", "Au", "fcc", 4.0784, None, None, (0.6, 4, 30, 1.8, 200)],
]

# TESTMOLECULES
TESTMOLECULES = [
    # Pd: FCC a: 3.8901 Å
    ["Phos-cage.mol2", "Pd", "fcc", 3.8901, None, None, (0.6, 1.4, 12, 1.8, 5)],
    # Pd: FCC a: 3.8901 Å
    [
        "Phos-cage-decamer.mol2",
        "Pd",
        "fcc",
        3.8901,
        None,
        None,
        (0.6, 2.3, 12.0, 2.4, 200.0),
    ],
    # Pd: FCC a: 3.8901 Å
    ["MTC1.mol2", "Pd", "fcc", 3.8901, None, None, (0.6, 1.4, 12, 2.4, 200)],
    # Pt: FCC a: 3.9233 Å
    ["Cu-MOC.mol2", "Pt", "fcc", 3.9233, None, None, (0.6, 1.4, 12, 3.6, 200)],
    # Ru: HCP a: 2.7053 Å, c: 4.2814 Å
    ["Cu-MOC.mol2", "Ru", "hcp", 2.7053, None, 4.2814, (0.6, 1.4, 12, 3.6, 200)],
]


# %%
print(TESTMOLECULES[2])
smc_file, atom_type, lattice_type, a, b, c, (s, pi, po, rd, vc) = TESTMOLECULES[2]
print(smc_file, atom_type, lattice_type, a, b, c, (s, pi, po, rd, vc))

# Create a DataFrame to store the results
summary = pandas.DataFrame()

# Get cage name and
cage_name = os.path.basename(smc_file).split(".")[0]

# Create basedir to store the results
if not os.path.exists(os.path.join("results", cage_name)):
    os.makedirs(os.path.join("results", cage_name))

# Print running message
print(f"Running {cage_name} with {atom_type} for {lattice_type} lattice ...")

# Step 0: Create a Cage object
cage = Cage()

# Step 1: Load cage
cage.load(smc_file)

# Step 2: Detect cavities
cage.detect_cavity(
    step=s,
    probe_in=pi,
    probe_out=po,
    removal_distance=rd,
    volume_cutoff=vc,
)

# Step 3: Pack atoms (No optimization)
cage.pack(
    atom_type=atom_type,
    lattice_type=lattice_type,
    a=a,
    b=b,
    c=c,
    clashing_tolerance=0.0,
    optimize=False,
)

# %%

# Step 3.1: Get cluster paramters
from ase.cluster.cubic import FaceCenteredCubic

surfaces = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
layers = cage._get_cluster_layers("Pd", 0.5)

# Step 3.2 Build raw cluster
cluster = FaceCenteredCubic(
    symbols=atom_type,
    surfaces=surfaces,
    layers=layers,
    latticeconstant=a,
)

# Step 3.3: Move cluster to centroid of cage
cluster.positions += cage.centroid

# Step 3.4: Save raw cluster
cluster.write("raw-cluster.pdb")

# Step 3.5: Filter atoms outside the cavity
f1 = cage._filter_outside_cavity(cluster)
f1.write("filtered-cluster-1.pdb")

# %%

# Step 3.6: Filter atoms clashing with cage
f2 = cage._filter_clashing_atoms(cluster)  # -> Not working properly
f2.write("filtered-cluster-2.pdb")
# %%
import numpy

cluster_distances = numpy.linalg.norm(
    cluster.positions[:, numpy.newaxis, :] - cluster.positions[numpy.newaxis, :, :],
    axis=-1,
)
radius = (cluster_distances[cluster_distances > 0] / 2).min()
cluster_radii = numpy.full(len(cluster), radius)
limits = cluster_radii[:, numpy.newaxis] + cage.universe.atoms.radii

# %%

# Save experiment conditions
exp = f"{atom_type}({lattice_type}-{cage.cluster.lattice_constants})@{cage_name}".strip(
    " "
)

# Step 4: Save cavity and packed atoms
cage.cavity.save(f"results/{cage_name}/cavity.pdb")
cage.cluster.save(os.path.join("results", cage_name, f"{exp}.pdb"))

# Step 5: Save summary
summary = pandas.concat([summary, cage.cluster.summary], axis=1)

# Change column names to reflect the experiment conditions
summary.columns = [exp]

# Save summary to file
summary.to_csv("results/summary.csv")

display(summary)
# %%
