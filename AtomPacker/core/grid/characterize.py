# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
A subpackage of *AtomPacker* that contains the `get_depth` function.

The `get_depth` function is used to convert a grid of points to a numpy
.ndarray of coordinates.
"""

__all__ = ["get_depths"]

import numpy

from pyKVFinder import depth


def get_depths(grid: numpy.ndarray, step: float) -> numpy.ndarray:
    """
    Get the depth of a grid of cavity points.

    Parameters
    ----------
    grid : numpy.ndarray
        The grid of cavity points.
    step : float
        The step size of the grid.

    Returns
    -------
    depths : numpy.ndarray
        The depths of the cavity points.
    """
    # Get depth of grid points
    depths, _, _ = depth(grid, step)

    # Get the indexes of the grid points
    indexes = numpy.argwhere(grid > 1)

    return depths[indexes[:, 0], indexes[:, 1], indexes[:, 2]]
