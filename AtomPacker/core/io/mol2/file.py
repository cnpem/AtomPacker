# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading MOL2 file format.
"""

__all__ = ["load_mol2"]


from MDAnalysis import Universe

from ..vdw import _lookup_radii


def load_mol2(filename: str) -> Universe:
    """
    Load a MOL2 file into the :class:`MDAnalysis.Univese` object.

    Parameters
    ----------
    filename : str
        The filename of the structure file.

    Returns
    -------
    universe: MDAnalysis.Univese
        The MDAnalysis universe object containing the structure.
    """
    # Read PDB file in MDAnalysis.Universe
    universe = Universe(filename)

    # Add radii to topology
    universe.add_TopologyAttr(
        "radii",
        _lookup_radii(universe.atoms.elements),
    )

    return universe
