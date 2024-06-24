# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
A subpackage of *AtomPacker* that contains the `Cluster` class.

The `Cluster` class is used to represent a cluster of atoms, following a
lattice, inside a supramolecular cages. It provides methods for calculating
cluster properties, and exporting cavity data.
"""

import os
from typing import Any, Dict, Tuple

import ase
import numpy
import pandas
from MDAnalysis import Universe
from plotly.express import scatter_3d

from .Cavity import Cavity


class Cluster:
    """
    A class representing a cluster of atoms inside a supramolecular cage.

    The `Cluster` class provides methods for calculating cluster properties and
    exporting cluster data.
    """

    def __init__(
        self,
        cluster: ase.cluster.Cluster,
        cavity: Cavity,
    ):
        """
        Create a new `Cluster` object.

        Parameters
        ----------
        cluster : ase.cluster.Cluster
            The cluster of atoms from `ase` package.
        cavity : Cavity
            The cavity object representing the cavity in the supramolecular
            cage.
        """
        self._cavity = cavity
        self._cluster = cluster

        # Cluster information
        self.atom_type = cluster.get_chemical_symbols()[0]
        self.lattice_type = cluster.symmetry
        self.lattice_constants = self._get_lattice_constants()

        # Atomic information: MDAnalysis Universe
        self.universe = self._get_universe()

    def __repr__(self) -> str:
        """
        Return a string representation of the :class:`AtomPacker.structure
        .Cluster` object.

        Returns
        -------
        str
            A string representation of the object.
        """
        return f"<AtomPacker.structure.Cluster at {hex(id(self))}>"

    def diameter(self, method: str = "maximum") -> float:
        """
        Calculate the diameter of the cluster.

        Parameters
        ----------
        method : str, optional
            The method to use for calculating the diameter of the cluster.
            Supported methods are 'maximum' (default), 'shape', and 'volume'.

            * 'maximum': calculates the maximum distance between any two \
            atoms in the cluster, considering their radius.

            * 'shape': uses `ase.cluster.Cluster.get_diameter('shape')` \
            to return the diameter of a sphere with the same volume as \
            the atoms.

            * 'volume': uses `ase.cluster.Cluster.get_diameter('volume')` \
            to return the averaged diameter calculated from the \
            directions given by the defined surfaces.

        Returns
        -------
        float
            The diameter of the cluster.
        """
        # Calculate distances between atoms
        if method == "maximum":
            diameter = self._get_distances().max()
        elif method == "shape":
            diameter = self._cluster.get_diameter(method="shape")
        elif method == "volume":
            diameter = self._cluster.get_diameter(method="volume")

        return diameter

    def preview(
        self, renderer: str = "browser", opacity: float = 1.0, **kwargs: Dict[str, Any]
    ) -> None:
        """
        Preview the cavity in a 3D viewer.

        Parameters
        ----------
        renderer : str, optional
            The renderer to use for the 3D viewer. Supported renderers are
            'browser' (default), 'notebook' and 'png'.
        opacity : float, optional
            The opacity of the atoms in the 3D viewer, by default 1.0. The
            opacity value ranges from 0.0 (completely transparent) to 1.0
            (completely opaque).
        **kwargs : Dict[str, Any]
            Additional keyword arguments to pass to the scatter_3d function of
            the plotly.express package.
        """
        if self._cluster is not None:
            x, y, z = self.coordinates.T
            radii = self.universe.atoms.radii
            fig = scatter_3d(
                x=x,
                y=y,
                z=z,
                size=(radii * 2),
                color=numpy.full(x.shape[0], "Cluster"),
                color_discrete_map={
                    "Cage": "rgba(99, 110, 255, 1)",  # blue
                    "Cavity": "rgba(102, 102, 102, 0.25)",  # light gray
                    "Cluster": "rgba(239, 85, 59, 1)",  # red
                },
                labels={"color": "Structures"},
                opacity=opacity,
                **kwargs,
            )
            # Update layout
            factor = 0.5
            fig.update_layout(
                scene=dict(
                    xaxis=dict(
                        showticklabels=False,
                        range=[
                            (x.min() - factor * x.ptp()),
                            (x.max() + factor * x.ptp()),
                        ],
                    ),
                    yaxis=dict(
                        showticklabels=False,
                        range=[
                            (y.min() - factor * y.ptp()),
                            (y.max() + factor * y.ptp()),
                        ],
                    ),
                    zaxis=dict(
                        showticklabels=False,
                        range=[
                            (z.min() - factor * z.ptp()),
                            (z.max() + factor * z.ptp()),
                        ],
                    ),
                )
            )
            # Update marker
            fig.update_traces(
                marker=dict(
                    sizeref=0.001,
                )
            )
            fig.show(renderer)

    def save(self, filename: str = "cluster.pdb") -> None:
        """
        Export the cluster data to a file.

        This method exports the cluster data represented by the `Cluster`
        object to a file in a specified format. The supported formats include
        PDB, and XYZ.

        Parameters
        ----------
        filename : str, optional
            The name of the file to export the cavity data to, by default
            "cluster.pdb". The file format is determined by the suffix.
            Supported formats are: .pdb, .xyz.
        """
        # We only need the suffix here
        _, suffix = os.path.splitext(filename)

        # Check if filename is a supported format
        if suffix not in [".pdb", ".xyz"]:
            raise ValueError(
                f"Unsupported file format: {suffix}. Supported formats are: \
.pdb, .xyz."
            )

        # Save the cavity to a file
        self.universe.atoms.write(filename)

    @property
    def summary(self) -> pandas.DataFrame:
        """
        Print a summary of the cluster properties.
        """
        return pandas.DataFrame.from_dict(
            {
                "Atom Type": self.atom_type,
                "Atom Radius": self._get_radii(),
                "Cavity Volume (Å³)": self._cavity.volume,
                "Diameter (maximum)": self.diameter(method="maximum"),
                "Diameter (shape)": self.diameter(method="shape"),
                "Diameter (volume)": self.diameter(method="volume"),
                "Lattice Constants": self.lattice_constants,
                "Lattice Type": self.lattice_type,
                "Maximum Number of Atoms": self.maximum_number_of_atoms,
                "Number of Atoms": self.number_of_atoms,
            },
            orient="index",
            columns=[self._cluster.get_chemical_formula()],
        )

    def _get_distances(self) -> numpy.ndarray:
        """
        Get the distances between atoms in the cluster.

        Returns
        -------
        numpy.ndarray
            A 2D array of distances between atoms in the
            cluster.
        """
        # Calculate distances between atoms
        distances = numpy.linalg.norm(
            self._cluster.positions[:, numpy.newaxis, :]
            - self._cluster.positions[numpy.newaxis, :, :],
            axis=-1,
        )

        return distances

    def _get_lattice_constants(
        self,
    ) -> Tuple[float, float, float] | Tuple[float, float] | float:
        """
        Get the lattice constants of the cluster.

        Returns
        -------
        Tuple[float, float, float] | Tuple[float, float] | float
            The lattice constants of the cluster.
        """
        if self.lattice_type == "hcp":
            return (
                self._cluster.lattice_basis[0, 0],
                self._cluster.lattice_basis[2, 2],
            )
        elif self.lattice_type == "fcc":
            return self._cluster.lattice_basis[0, 0]
        elif self.lattice_type == "bcc":
            return self._cluster.lattice_basis[0, 0]
        elif self.lattice_type == "sc":
            return self._cluster.lattice_basis[0, 0]

    def _get_radii(self) -> float:
        """
        Get the radius of the atoms in the cluster.

        Returns
        -------
        float
            The radius of the atoms in the cluster.
        """
        # Calculate distances between atoms
        distances = self._get_distances()

        # Calculate radii, ignoring self-distances
        radii = (distances[distances > 0] / 2).min()

        return radii

    def _get_universe(self) -> Universe:
        """
        Get the `MDAnalysis.Universe` object of the _cluster.

        This method creates an `MDAnalysis.Universe` object from the
        `ase.cluster.Cluster` object.

        Returns
        -------
        Universe
            The `MDAnalysis.Universe` object of the cluster.
        """
        # Create an empty Universe
        universe = Universe.empty(n_atoms=len(self._cluster), trajectory=True)

        # Add coordinates
        universe.atoms.positions = self._cluster.positions

        # Add CRYST1 record
        universe.dimensions = [1.0, 1.0, 1.0, 90.0, 90.0, 90.0]

        # Add topology attributes
        universe.add_TopologyAttr(
            "record_types", ["ATOM"] * universe.atoms.n_atoms
        )  # record type
        universe.add_TopologyAttr(
            "names", [self.atom_type] * universe.atoms.n_atoms
        )  # atom name
        universe.add_TopologyAttr(
            "elements", [self.atom_type] * universe.atoms.n_atoms
        )  # atom type
        universe.add_TopologyAttr(
            "radii",
            [self._get_radii()] * universe.atoms.n_atoms,  # atom radius
        )
        universe.add_TopologyAttr("resids", [1])  # resids
        universe.add_TopologyAttr("resnums", [1])  # resnums
        universe.add_TopologyAttr("resnames", ["NC"])  # resnames
        universe.add_TopologyAttr(
            "chainIDs", ["X"] * universe.atoms.n_atoms
        )  # chainIDs
        universe.add_TopologyAttr("segids", [" "])  # segids
        universe.add_TopologyAttr("icodes", [" "])  # icodes
        universe.add_TopologyAttr("altLocs", [" "] * universe.atoms.n_atoms)  # altLocs
        universe.add_TopologyAttr(
            "occupancies", [1.0] * universe.atoms.n_atoms
        )  # occupancy
        universe.add_TopologyAttr(
            "tempfactors", [0.0] * universe.atoms.n_atoms
        )  # temperature factor
        universe.add_TopologyAttr(
            "formalcharges", [0] * universe.atoms.n_atoms
        )  # formalcharge

        return universe

    @property
    def coordinates(self):
        """
        Return the coordinates of the atoms in the cluster.
        """
        return self.universe.atoms.positions

    @property
    def maximum_number_of_atoms(self) -> int:
        """
        Calculate the maximum number of atoms in the cluster.

        The maximum number of atoms is calculated based on the volume of the
        cavity divided by the volume of the atom in the cluster. The formula
        is given by:

        .. math::
            N_{max} = \\left \\lceil \\frac{V_{cav}}{V_a} \\right \\rceil
                    = \\left \\lceil \\frac{V_{cav}}{\\frac{4}{3}\\pi R_a^3} \\
                        \\right \\rceil

        where :math:`V_{cav}` is the volume of the cavity, :math:`V_a` is the
        volume of the atom in the cluster, :math:`R_a` is the radius of the
        atom in the cluster, and :math:`\\lceil \\cdot \\rceil` is the ceiling
        function.
        """
        return numpy.ceil(
            self._cavity.volume
            / ((4 / 3) * numpy.pi * self.universe.atoms.radii[0] ** 3)
        ).astype(int)

    @property
    def number_of_atoms(self):
        """
        Get the number of atoms in the cluster.
        """
        return len(self._cluster)
