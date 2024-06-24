# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for looking up the van der Waals radii of atoms.
"""

__all__ = ["_lookup_radii"]

import numpy
from pyKVFinder import read_vdw


def _lookup_radii(names: numpy.ndarray) -> numpy.ndarray:
    """
    Look up the van der Waals radii for a list of atom names.

    Parameters
    ----------
    vdw : Dict[str, float]
        A dictionary containing the van der Waals radii for each atom type.
    names : numpy.ndarray
        An array containing the names of the atoms.

    Returns
    -------
    radii : numpy.ndarray
        An array containing the van der Waals radii for each atom.
    """
    # Read van der Waals radii
    vdw = read_vdw()["GEN"]

    return numpy.vectorize(lambda x: vdw[x.upper()])(names)
