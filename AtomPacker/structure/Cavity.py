# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
A subpackage of *AtomPacker* that contains the `Cavity` class.

The `Cavity` class is used to represent a cavity in a molecular structure. It
provides methods for detecting cavities, calculating cavity properties, and
exporting cavity data.
"""

import os
import warnings
from typing import Any, Dict, List

import numpy
from MDAnalysis import Universe
from plotly.express import scatter_3d

from ..core.grid import get_coordinates, get_depths


class Cavity:
    """
    A class representing a cavity in a molecular structure.

    The `Cavity` class provides methods for detecting cavities, calculating
    cavity properties, and exporting cavity data.
    """

    def __init__(
        self,
        grid: numpy.ndarray,
        step: float,
        probe_in: float,
        probe_out: float,
        removal_distance: float,
        volume_cutoff: float,
        vertices: numpy.ndarray,
        surface: str = "SES",
    ):
        """
        Create a new `AtomPacker.Cage` object.

        Parameters
        ----------
        grid : numpy.ndarray
            The grid points of the cavity.
        step : float
            The step size of the grid.
        probe_in : float
            The radius of the Probe In.
        probe_out : float
            The radius of the Probe Out.
        removal_distance : float
            The removal distance.
        volume_cutoff : float
            The volume cutoff.
        vertices : numpy.ndarray
            The vertices (origin, X-axis, Y-axis, Z-axis) of the grid.
        surface : str, optional
            Surface representation. Keywords options are SES (Solvent Excluded
            Surface) or SAS (Solvent Accessible Surface), by default SES.
        """
        self.grid = grid
        self._step = step
        self._probe_in = probe_in
        self._probe_out = probe_out
        self._removal_distance = removal_distance
        self._volume_cutoff = volume_cutoff
        self._vertices = vertices
        self._surface = surface
        self.universe = self._get_universe()

    def __repr__(self) -> str:
        """
        Return a string representation of the :class:`AtomPacker.structure
        .Cavity` object.

        Returns
        -------
        str
            A string representation of the object.
        """
        return f"<AtomPacker.structure.Cavity at {hex(id(self))}>"

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
        if self.grid is not None:
            x, y, z = self.coordinates.T
            radii = self.universe.atoms.radii
            labels = self.universe.atoms.tempfactors
            fig = scatter_3d(
                x=x,
                y=y,
                z=z,
                size=(radii * 2),
                color=labels,
                color_continuous_scale="rainbow",
                opacity=opacity,
                labels={"color": "Depth"},
                **kwargs,
            )
            # Update layout
            factor = 0.3
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

    def select_cavity(self, indices: List[int]) -> None:
        """
        Select a cavity from the list of detected cavities.

        This method selects a cavity from the list of detected cavities based
        on the indices of the cavity in the list.

        Parameters
        ----------
        indices : List[int]
            The indices of the cavities to select.
            Cavity points in the 3D grid (grid[nx][ny][nz]). Cavities
            array has integer labels in each position, that are:

            * -1: solvent points;

            * 0: cage points;

            * 1: empty space points;

            * >=2: cavity points.
        """
        # Copy cavities
        selected = numpy.copy(self.grid)

        # When outside selection, change cavities tags to 1
        for cav in range(2, self.grid.max() + 1):
            if cav not in indices:
                selected[self.grid == cav] = 1

        # Update grid
        self.grid = selected
        self.universe = self._get_universe()

    def save(self, filename: str = "cavity.pdb") -> None:
        """
        Export the cavity data to a file.

        This method exports the cavity data represented by the `Cavity` object
        to a file in a specified format. The supported formats include PDB, and
        XYZ.

        Parameters
        ----------
        filename : str
            The name of the file to export the cavity data to. The filename of
            the structure file. The file format is determined by the suffix.
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

    def _get_universe(self) -> Universe:
        """
        Get the `MDAnalysis.Universe` object of the cavity.

        This method creates an `MDAnalysis.Universe` object from the cavity grid
        data.

        Returns
        -------
        Universe
            The `MDAnalysis.Universe` object of the cavity.
        """
        # Create an empty Universe
        universe = Universe.empty(n_atoms=(self.grid > 1).sum(), trajectory=True)

        # Add coordinates
        universe.atoms.positions = get_coordinates(
            self.grid, self._step, self._vertices
        )

        # Add CRYST1 record
        universe.dimensions = [1.0, 1.0, 1.0, 90.0, 90.0, 90.0]

        # Add topology attributes
        universe.add_TopologyAttr(
            "record_types", ["ATOM"] * universe.atoms.n_atoms
        )  # record type
        universe.add_TopologyAttr(
            "names", ["DUMMY"] * universe.atoms.n_atoms
        )  # atom name
        universe.add_TopologyAttr(
            "elements", ["DUMMY"] * universe.atoms.n_atoms
        )  # atom type
        universe.add_TopologyAttr(
            "radii",
            [self._step / 2] * universe.atoms.n_atoms,  # atom radius
        )
        universe.add_TopologyAttr("resids", [1])  # resids
        universe.add_TopologyAttr("resnums", [1])  # resnums
        universe.add_TopologyAttr("resnames", ["CAV"])  # resnames
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
            "tempfactors", get_depths(self.grid, self._step)
        )  # depth
        universe.add_TopologyAttr(
            "formalcharges", [0] * universe.atoms.n_atoms
        )  # formalcharge

        return universe

    @property
    def coordinates(self) -> numpy.ndarray:
        """
        Get the xyz coordinates of the cavity.
        """
        return self.universe.atoms.positions

    @property
    def volume(self) -> numpy.float64:
        """
        Calculate the volume of the cavity.

        This method calculates the volume of the cavity represented by the
        `Cavity` object. It uses a grid-based algorithm to estimate the volume
        of the cavity based on the number of grid points inside the cavity.
        """
        if self.grid.max() > 2:
            warnings.warn(
                "Cavity has more than one cavity. Returning total volume.\
Users can select cavities using select_cavity method."
            )
        return ((self.grid > 1).sum() * (self._step**3)).round(2)
