# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading XYZ file format.
"""

__all__ = ["load_xyz"]

from MDAnalysis import Universe

from ..vdw import _lookup_radii


def load_xyz(filename: str) -> Universe:
    """
    Load a XYZ file into the :class:`MDAnalysis.Univese` object.

    Parameters
    ----------
    filename : str
        The filename of the structure file.

    Returns
    -------
    universe: MDAnalysis.Univese
        The MDAnalysis universe object containing the structure.
    """
    # Read XYZ file in MDAnalysis.Universe
    universe = Universe(filename)

    # Add radii to topology
    universe.add_TopologyAttr(
        "radii",
        _lookup_radii(universe.atoms.names),
    )

    return universe
