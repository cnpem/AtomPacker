# This source code is part of the AtomPacker package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
A subpackage of *AtomPacker* that contains the `Openings` class.

The `Openings` class is used to represent openings in a cavity structure.
It provides methods for detecting openings, calculating their properties, and
exporting opening data.
"""

import os
import warnings

import numpy
from MDAnalysis import Universe
from plotly.express import scatter_3d
from pyKVFinder import openings, depth

from ..core.grid import get_coordinates, get_depths


class Openings:
    """
    A class representing openings in a cavity structure.

    The `Openings` class provides methods for detecting openings, calculating
    their properties, and exporting openings data.
    """

    def __init__(
        self,
        cavities: numpy.ndarray,
        step: float,
        vertices: numpy.ndarray,
        openings_cutoff: float = 1,
        verbose: bool = False,
    ):
        """
        Create a new `AtomPacker.Openings` object.

        Parameters
        ----------
        cavities : numpy.ndarray
            The cavity points.
        step : float
            The step size of the grid.
        vertices : numpy.ndarray
            The vertices (origin, X-axis, Y-axis, Z-axis) of the grid.
        openings_cutoff : float, optional
            The cutoff value for detecting openings (default is 1). The minimum
            number of points in an opening to be considered valid.
        verbose : bool, optional
            If True, print detailed information during processing (default is False).
        """
        self.nopenings: int
        self.grid: numpy.ndarray
        self.areas: dict[str, float]

        self.nopenings, self.grid, self.areas = self._detect(
            cavities, step, openings_cutoff=openings_cutoff, verbose=verbose
        )
        self.diameters: dict[str, float] = self._get_diameter()

        self._step: float = step
        self._vertices: numpy.ndarray = vertices

        self.universe: Universe = self._get_universe()

    def __repr__(self) -> str:
        return f"<AtomPacker.structure.Openings at {hex(id(self))}>"

    @property
    def coordinates(self) -> numpy.ndarray:
        """
        Get the xyz coordinates of the cavity.
        """
        return self.universe.atoms.positions

    def preview(
        self,
        renderer: str = "browser",
        opacity: float = 1.0,
        **kwargs: dict[str, object],
    ) -> None:
        """
        Preview the openings in a 3D viewer.

        Parameters
        ----------
        renderer : str, optional
            The renderer to use for the 3D viewer. Supported renderers are
            'browser' (default), 'notebook' and 'png'.
        opacity : float, optional
            The opacity of the atoms in the 3D viewer, by default 1.0. The
            opacity value ranges from 0.0 (completely transparent) to 1.0
            (completely opaque).
        **kwargs : dict[str, object], optional
            Additional keyword arguments to pass to the scatter_3d function of
            the plotly.express package.
        """
        if self.grid is not None:
            x, y, z = self.coordinates.T
            radii = self.universe.atoms.radii
            fig = scatter_3d(
                x=x,
                y=y,
                z=z,
                size=(radii * 2),
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
                            (x.min() - factor * numpy.ptp(x)),
                            (x.max() + factor * numpy.ptp(x)),
                        ],
                    ),
                    yaxis=dict(
                        showticklabels=False,
                        range=[
                            (y.min() - factor * numpy.ptp(y)),
                            (y.max() + factor * numpy.ptp(y)),
                        ],
                    ),
                    zaxis=dict(
                        showticklabels=False,
                        range=[
                            (z.min() - factor * numpy.ptp(z)),
                            (z.max() + factor * numpy.ptp(z)),
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

    def save(self, filename: str = "openings.pdb") -> None:
        """
        Save the openings to a PDB file.

        Parameters
        ----------
        filename : str, optional
            The name of the output PDB file (default is "openings.pdb").
        """
        # We only need the suffix here
        _, ext = os.path.splitext(filename)

        # Check if filename is a supported format
        if ext not in [".pdb", ".xyz"]:
            raise ValueError(
                f"Unsupported file format: {ext}. Supported formats are: \
.pdb, .xyz."
            )

        # Save the openings to a file
        self.universe.atoms.write(filename)

    def _detect(
        self,
        cavities: numpy.ndarray,
        step: float,
        openings_cutoff: float = 1,
        verbose: bool = False,
    ) -> tuple[int, numpy.ndarray, dict[str, float]]:
        """
        Detect openings in the cavity grid.

        Parameters
        ----------
        cavities : numpy.ndarray
            The cavity points.
        step : float
            The step size of the grid.
        openings_cutoff : float, optional
            The cutoff value for detecting openings (default is 1). The minimum
            number of points in an opening to be considered valid.
        verbose : bool, optional
            If True, print detailed information during processing (default is False).

        Returns
        -------
        nopenings : int
            The number of openings detected.
        ogrid : numpy.ndarray
            The grid of openings.
        oarea : dict[str, float]
            The area of each opening.
        """
        if cavities.max() > 2:
            warnings.warn("Cavity has more than one cavity.")

        # Calculate depth of cavity points
        depths, _, _ = depth(cavities, step, verbose=verbose)

        # Calculate openings and area of openings
        nopenings, grid, aopenings = openings(
            cavities,
            depths,
            step=step,
            openings_cutoff=openings_cutoff,
            verbose=verbose,
        )

        # Flatten and sort opening areas
        areas = {
            key: opening
            for cavity in aopenings.values()
            for key, opening in cavity.items()
        }
        areas = dict(sorted(areas.items()))

        return nopenings, grid, areas

    def _get_diameter(self) -> dict[str, float]:
        """
        Get the diameter of each opening.

        The diameter is calculated from the area using the formula:

            diameter = 2 * sqrt(area / π)

        Returns
        -------
        diameters : dict[str, float]
            The diameter of each opening.
        """
        diameters = {}
        for key, area in self.areas.items():
            diameters[key] = float(2 * numpy.sqrt(area / numpy.pi))

        return diameters

    def _get_universe(self) -> Universe:
        """
        Get the `MDAnalysis.Universe` object of the openings.

        This method creates an `MDAnalysis.Universe` object from the cavity
        grid data.

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

        # Assign chainIDs based on opening indices
        indexes = (
            self.grid[
                tuple(
                    ((universe.atoms.positions - self._vertices[0]) / self._step)
                    .round()
                    .astype(int)
                    .T
                )
            ]
            - 1
        ).astype(int)
        chain_ids = [chr(65 + idx) for idx in indexes]

        # Add CRYST1 record
        universe.dimensions = [1.0, 1.0, 1.0, 90.0, 90.0, 90.0]

        # Add topology attributes
        n = universe.atoms.n_atoms
        universe.add_TopologyAttr("record_types", ["ATOM"] * n)  # record type
        universe.add_TopologyAttr("names", ["DUMMY"] * n)  # atom name
        universe.add_TopologyAttr("elements", ["DUMMY"] * n)  # atom type
        universe.add_TopologyAttr(
            "radii",
            [self._step / 2] * n,  # atom radius
        )
        universe.add_TopologyAttr("resids", [1])  # resids
        universe.add_TopologyAttr("resnums", [1])  # resnums
        universe.add_TopologyAttr("resnames", ["OPE"])  # resnames
        universe.add_TopologyAttr("chainIDs", chain_ids)  # chainIDs
        universe.add_TopologyAttr("segids", [" "])  # segids
        universe.add_TopologyAttr("icodes", [" "])  # icodes
        universe.add_TopologyAttr("altLocs", [" "] * n)  # altLocs
        universe.add_TopologyAttr("occupancies", [1.0] * n)  # occupancy
        universe.add_TopologyAttr(
            "tempfactors", get_depths(self.grid, self._step)
        )  # depth
        universe.add_TopologyAttr("formalcharges", [0] * n)  # formalcharge

        return universe
