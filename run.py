# %%
from AtomPacker import Cage

import numpy as np
from sklearn.decomposition import PCA
from ase.cluster import FaceCenteredCubic


def obb(points):
    """
    Calculate the Oriented Bounding Box (OBB) for a set of 3D points.

    Parameters:
        points (numpy.ndarray): Nx3 array of 3D coordinates for the points.

    Returns:
        center (numpy.ndarray): The center of the OBB.
        axes (numpy.ndarray): 3x3 array where each row is an axis of the OBB.
        extents (numpy.ndarray): Half-lengths (extents) along each OBB axis.
    """
    # Step 1: Compute the PCA to get the principal axes
    pca = PCA(n_components=3)
    pca.fit(points)

    # Principal axes (directions of the OBB)
    axes = pca.components_

    # Step 2: Project the points onto each principal axis to find the extents
    # Transform points to the PCA space
    transformed_points = pca.transform(points)

    # Find min and max projections along each axis
    min_proj = np.min(transformed_points, axis=0)
    max_proj = np.max(transformed_points, axis=0)

    # Step 3: Calculate the center and extents of the OBB
    center_pca = (min_proj + max_proj) / 2.0  # Center in the PCA space
    extents = (max_proj - min_proj) / 2.0  # Half-lengths (extents) along each axis

    # Convert the center from PCA space back to the original space
    center = pca.inverse_transform(center_pca)

    return center, axes, extents


cage = Cage()
cage.load("MTC1.mol2")
cage.detect_cavity(
    step=0.6, probe_in=1.4, probe_out=12, removal_distance=2.4, volume_cutoff=200
)

# Calculate OBB
center, axes, extents = obb(cage.coordinates)

# Cluster parameters
symbols = "Pd"
surfaces = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
layers = cage._get_cluster_layers(symbols, factor=0.2)
lattice_constant = 3.8901
cluster = FaceCenteredCubic(
    symbols="Pd",
    surfaces=surfaces,
    layers=layers,
    latticeconstant=lattice_constant,
)
cluster.positions += cage.centroid
cluster.write("raw_cluster.pdb")


def align_cluster_to_obb(cluster, obb_center, obb_axes):
    """
    Aligns a cluster to an Oriented Bounding Box (OBB) defined by a center and axes.

    Parameters:
        cluster (ase.Atoms): The cluster to align.
        obb_center (array-like): 3D coordinates of the OBB center.
        obb_axes (numpy.ndarray): 3x3 array where each row is a unit vector for the box orientation.

    Returns:
        aligned_cluster (ase.Atoms): The aligned and translated atomic cluster.
    """
    # Step 1: Translate the cluster to the origin (centered around its own center of mass)
    com = cluster.get_center_of_mass()
    cluster.translate(-com)
    
    # Step 2: Rotate the cluster to align with the OBB's orientation
    # We assume the default cluster orientation is along [1, 0, 0], [0, 1, 0], [0, 0, 1]
    rotation_matrix = obb_axes  # OBB axes act as the rotation matrix
    
    # Apply rotation matrix to each atom's position
    rotated_positions = np.dot(cluster.get_positions(), rotation_matrix)
    cluster.set_positions(rotated_positions)
    
    # Step 3: Translate the rotated cluster to the OBB center
    cluster.translate(obb_center)

    return cluster

cluster = align_cluster_to_obb(cluster, center, axes)

cluster.write("obb_cluster.pdb")

# Step 3: Translate the cluster to the OBB center


# %% FINAL PROCEDURE
# Cage
cage = Cage()
cage.load("MTC1.mol2")
cage.detect_cavity(
    step=0.6, probe_in=1.4, probe_out=12, removal_distance=2.4, volume_cutoff=200
)

# OBB
pca = PCA(n_components=3)
pca.fit(cage.coordinates)
obb_axes = pca.components_
obb_extents = pca.transform(cage.coordinates).ptp(axis=0)
obb_center = (cage.coordinates.min(axis=0) + cage.coordinates.max(axis=0)) / 2.0

# Atom type
symbols = atom_type = "Pd"

# Layers
from ase.data import atomic_numbers, covalent_radii
obb_extents += 0.2 * obb_extents
atomic_number = atomic_numbers[atom_type]
atom_size = covalent_radii[atomic_number] * 2
layers = numpy.ceil(obb_extents / atom_size).astype(int)

# Surface 
surfaces = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

# Lattice constant
lattice_constant = 3.8901

# Cluster
cluster = FaceCenteredCubic(
    symbols="Pd",
    surfaces=surfaces,
    layers=layers,
    latticeconstant=lattice_constant,
)
cluster.write("raw_cluster.pdb")

# Apply rotation matrix to each atom's position
rotated_positions = numpy.dot(cluster.positions, obb_axes)
cluster.set_positions(rotated_positions)

# Translate to cage centroid
cluster.translate(cage.centroid)
cluster.write("cluster.pdb")
