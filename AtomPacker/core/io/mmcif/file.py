# This source code is part of the AtomPacker package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading mmCIF file format.
"""

__all__ = ["load_mmcif"]


from typing import Dict, Optional

import numpy
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from MDAnalysis import Universe

from ..vdw import _lookup_radii


def load_mmcif(filename: str, vdw: Optional[Dict[str, float]] = None) -> Universe:
    """
    Load a PDBx/mmCIF file into the :class:`MDAnalysis.Univese` object.

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
    # Read mmCIF file into a dictionary
    mmcif = MMCIF2Dict(filename)

    # Number of atoms and residues
    n_atoms = len(mmcif["_atom_site.id"])

    # Create an empty Universe
    universe = Universe.empty(n_atoms=n_atoms, trajectory=True)

    # Add topology attributes
    universe.add_TopologyAttr("names", mmcif["_atom_site.label_atom_id"])  # atom name
    universe.atoms.positions = numpy.c_[
        numpy.asarray(mmcif["_atom_site.Cartn_x"], dtype=float),
        numpy.asarray(mmcif["_atom_site.Cartn_y"], dtype=float),
        numpy.asarray(mmcif["_atom_site.Cartn_z"], dtype=float),
    ]  # xyz coordinates
    universe.add_TopologyAttr("elements", mmcif["_atom_site.type_symbol"])  # atom type

    # Add radii to topology
    universe.add_TopologyAttr(
        "radii",
        _lookup_radii(universe.atoms.names, vdw),
    )

    return universe
