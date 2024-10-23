# This source code is part of the AtomPacker package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading MOL2 file format.
"""

__all__ = ["load_mol2"]

import warnings
from string import digits
from typing import Dict, Optional

from MDAnalysis import Universe

from ..vdw import _lookup_radii

# Suppress MDAnalysis warnings
warnings.filterwarnings("ignore", category=UserWarning)


def load_mol2(filename: str, vdw: Optional[Dict[str, float]] = None) -> Universe:
    """
    Load a MOL2 file into the :class:`MDAnalysis.Univese` object.

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
        _lookup_radii(universe.atoms.elements, vdw),
    )

    return universe
