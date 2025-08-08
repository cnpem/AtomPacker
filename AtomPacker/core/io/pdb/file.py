# This source code is part of the AtomPacker package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading PDB file format.
"""

__all__ = ["load_pdb"]

from typing import Dict, Optional

from MDAnalysis import Universe

from ..vdw import _lookup_radii


def load_pdb(filename: str, vdw: Optional[Dict[str, float]] = None) -> Universe:
    """
    Load a PDB file into the :class:`MDAnalysis.Univese` object.

    Parameters
    ----------
    filename : str
        The filename of the structure file.
    vdw : Dict[str, float]
        A dictionary containing the van der Waals radii for each atom type,
        by default None. If None, the radii will be looked up from the
        `pyKVFinder` package.

    Returns
    -------
    universe: MDAnalysis.Univese
        The MDAnalysis universe object containing the structure.
    """
    # Read PDB file in MDAnalysis.Universe
    universe = Universe(filename)

    # Get van der Waals radii
    try:
        radii = _lookup_radii(universe.atoms.names, vdw)
    except KeyError:
        radii = _lookup_radii(universe.atoms.elements, vdw)

    # Add radii to topology
    universe.add_TopologyAttr("radii", radii)

    return universe
