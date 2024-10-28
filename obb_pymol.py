# %%
from AtomPacker import Cage

import numpy as np
from pymol import cmd
from pymol.cgo import LINEWIDTH, BEGIN, LINES, COLOR, VERTEX, END

import numpy as np
from sklearn.decomposition import PCA


def calculate_obb(points):
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


def create_obb_box(center, axes, extents, name="obb_box"):
    """
    Create an oriented bounding box in PyMOL.

    Parameters:
        center (list or np.ndarray): The center of the box (3 elements).
        axes (np.ndarray): 3x3 array where each row is a unit vector for the box orientation.
        extents (list or np.ndarray): Half-lengths (extents) along each box axis (3 elements).
        name (str): The name for the PyMOL object.
    """
    # Calculate the 8 corners of the box
    corners = []
    for i in [-1, 1]:
        for j in [-1, 1]:
            for k in [-1, 1]:
                corner = (
                    center
                    + i * extents[0] * axes[0]
                    + j * extents[1] * axes[1]
                    + k * extents[2] * axes[2]
                )
                corners.append(corner)

    # Define edges of the box by specifying pairs of corner indices
    edges = [
        (0, 1),
        (1, 3),
        (3, 2),
        (2, 0),  # Bottom face
        (4, 5),
        (5, 7),
        (7, 6),
        (6, 4),  # Top face
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),  # Vertical edges
    ]

    # Create the CGO object to draw the box
    box_cgo = [
        LINEWIDTH,
        2.0,
        BEGIN,
        LINES,
        COLOR,
        1.0,
        1.0,
        0.0,  # Yellow color for the box lines
    ]

    for start, end in edges:
        box_cgo.extend([VERTEX, *corners[start]])
        box_cgo.extend([VERTEX, *corners[end]])

    box_cgo.append(END)

    # Load the CGO object into PyMOL
    cmd.load_cgo(box_cgo, name)


# %%

cage = Cage()
cage.load("MTC1.mol2")
cage.detect_cavity(
    step=0.6, probe_in=1.4, probe_out=12, removal_distance=2.4, volume_cutoff=200
)

center, axes, extents = calculate_obb(cage.coordinates)

print("OBB Center:", center)
print("OBB Axes:\n", axes)
print("OBB Extents:", extents)


create_obb_box(center, axes, extents, name="obb_box")
