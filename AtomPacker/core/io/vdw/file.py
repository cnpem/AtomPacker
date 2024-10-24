# This source code is part of the AtomPacker package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for looking up the van der Waals radii of atoms.
"""

__all__ = ["_lookup_radii"]

from typing import Dict, Optional

import numpy
from pyKVFinder import read_vdw


def _lookup_radii(
    names: numpy.ndarray, vdw: Optional[Dict[str, float]]
) -> numpy.ndarray:
    """
    Look up the van der Waals radii for a list of atom names.

    Parameters
    ----------
    names : numpy.ndarray
        An array containing the names of the atoms.
    vdw : Dict[str, float]
        A dictionary containing the van der Waals radii for each atom type.
        If None, the van der Waals radii are looked up from the `pyKVFinder`
        package.

    Returns
    -------
    radii : numpy.ndarray
        An array containing the van der Waals radii for each atom.

    Raises
    ------
    ValueError
        If not all elements are provided in the `vdw` dictionary.
    """
    # Read van der Waals radii
    if vdw is None:
        vdw = read_vdw()["GEN"]
        return numpy.vectorize(lambda x: vdw[x.upper()])(names)
    elif isinstance(vdw, dict):
        # Check if all elements are provided in vdw
        if not all([name in vdw for name in names]):
            raise ValueError(
                f"Missing atoms in `vdw` dictionary: {set(names) - set(vdw.keys())}"
            )
        return numpy.vectorize(lambda x: vdw[x])(names)
