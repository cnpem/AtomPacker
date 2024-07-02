# This source code is part of the AtomPacker package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
A subpackage of *AtomPacker* that contains the `get_coordinates` function.

The `get_coordinates` function is used to convert a grid of points to a numpy
.ndarray of coordinates.
"""

__all__ = ["get_coordinates"]

import numpy

from pyKVFinder.grid import _get_sincos


def get_coordinates(
    grid: numpy.ndarray, step: float, vertices: numpy.ndarray
) -> numpy.ndarray:
    """
    Convert a grid of points to a numpy.ndarray of xyz coordinates.

    Parameters
    ----------
    grid : numpy.ndarray
        The grid points.
    step : float
        The step size of the grid.
    vertices : numpy.ndarray
        The vertices (origin, X-axis, Y-axis, Z-axis) of the grid.

    Returns
    -------
    coordinates : numpy.ndarray
        The xyz coordinates of the grid points.
    """
    # Get the indexes of the grid points
    indexes = numpy.argwhere(grid > 1)

    # Get origin from vertices
    P1, _, _, _ = vertices

    # Calculate sin and cos for each axis
    sincos = _get_sincos(vertices)

    # Convert indexes to xyz coordinates
    xaux, yaux, zaux = (indexes * step).T
    x = (
        (xaux * sincos[3])
        + (yaux * sincos[0] * sincos[2])
        - (zaux * sincos[1] * sincos[2])
        + P1[0]
    )
    y = (yaux * sincos[1]) + (zaux * sincos[0]) + P1[1]
    z = (
        (xaux * sincos[2])
        - (yaux * sincos[0] * sincos[3])
        + (zaux * sincos[1] * sincos[3])
        + P1[2]
    )

    # Merge xyz coordinates into a numpy.ndarray
    coordinates = numpy.array([x, y, z]).T

    return coordinates
