# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
The *AtomPacker* is a Python package for packing nanoclusters into \
supramolecular cages. It provides a set of tools for detecting cavities in \
molecular structures, packing nanoclusters into these cavities, and exporting \
the resulting structures.
"""

__version__ = "0.2.0"
__name__ = "AtomPacker"
license = "GNU GPL-3.0 License"

from .core import (get_coordinates, get_depths, load_mmcif, load_mol2,
                   load_pdb, load_xyz)
from .structure import (Cage, Cavity, Cluster, get_lattice_constants,
                        lattice_constants)
