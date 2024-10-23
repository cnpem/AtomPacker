# This source code is part of the AtomPacker package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading XYZ file format.
"""

__all__ = ["load_xyz"]

from typing import Dict, Optional

from MDAnalysis import Universe

from ..vdw import _lookup_radii


def load_xyz(filename: str, vdw: Optional[Dict[str, float]] = None) -> Universe:
    """
    Load a XYZ file into the :class:`MDAnalysis.Univese` object.

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
    # Read XYZ file in MDAnalysis.Universe
    universe = Universe(filename)

    # if chainIDs not provided, set them to 'X'
    if "chainIDs" not in universe.atoms._SETATTR_WHITELIST:
        universe.add_TopologyAttr(
            "chainIDs", ["X"] * universe.atoms.n_atoms
        )  # chainIDs

    # if resnames not provided, set them to 'X'
    if "resnames" not in universe.atoms._SETATTR_WHITELIST:
        universe.add_TopologyAttr("resnames", ["UNK"])  # resnames

    # Add radii to topology
    universe.add_TopologyAttr(
        "radii",
        _lookup_radii(universe.atoms.names, vdw),
    )

    return universe
