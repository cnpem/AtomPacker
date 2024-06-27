# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading MOL2 file format.
"""

__all__ = ["load_mol2"]

import warnings
from string import digits

from MDAnalysis import Universe

from ..vdw import _lookup_radii

# Suppress MDAnalysis warnings
warnings.filterwarnings("ignore", category=UserWarning)


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

    # if chainIDs not provided, set them to 'X'
    if "chainIDs" not in universe.atoms._SETATTR_WHITELIST:
        universe.add_TopologyAttr(
            "chainIDs", ["X"] * universe.atoms.n_atoms
        )  # chainIDs

    # Check if all elements are provided
    if (universe.atoms.elements == "").any():
        # Remove digits from element names
        elements = [
            atom.translate(str.maketrans("", "", digits))
            # Convert any hydrogens to H
            .replace("HA", "H").replace("HB", "H").replace("HC", "H")
            for atom in universe.atoms.names
        ]
        # Add elements to topology
        universe.add_TopologyAttr("elements", elements)

    # Add radii to topology
    universe.add_TopologyAttr(
        "radii",
        _lookup_radii(universe.atoms.elements),
    )

    return universe
