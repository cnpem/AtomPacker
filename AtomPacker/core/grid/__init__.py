# This source code is part of the AtomPacker package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
A subpackage of *AtomPacker* that contains the `get_coordinates` and
`get_depth` functions.

The `get_coordinates` function is used to convert a grid of points to a numpy
.ndarray of coordinates. The `get_depth` function is used to calculate the depth
of a grid of points.
"""

from .coordinates import get_coordinates
from .characterize import get_depths
