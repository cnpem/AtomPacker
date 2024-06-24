# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
A subpackage of *AtomPacker* that contains the `Cage` class.

Macromolecular structure files (PDB, PDBx/mmCIF, XYZ, etc.) can be used to
create a :class:`AtomPacker.structure.Cage` object, which is a container for
the atoms of the structure.
"""

import os
import warnings
from typing import Any, Dict, Optional, Tuple

import ase
import numpy
from ase.cluster.cubic import BodyCenteredCubic, FaceCenteredCubic, SimpleCubic
from ase.cluster.hexagonal import HexagonalClosedPacked
from ase.data import atomic_numbers, covalent_radii
from plotly.express import scatter_3d
from pyKVFinder import detect, get_vertices

from ..core import load_mmcif, load_mol2, load_pdb, load_xyz
from .Cavity import Cavity
from .Cluster import Cluster
from .data import get_lattice_constants


class Cage:
    """
    A container for the atoms of a macromolecular structure.

    The :class:`AtomPacker.Cage` class is used to store the atoms of a
    macromolecular structure.
    """

    def __init__(self):
        """
        Create a new :class:`AtomPacker.structure.Cage` object.
        """
        self.universe = None
        self.cavity = None
        self.cluster = None

    def __repr__(self) -> str:
        """
        Return a string representation of the :class:`AtomPacker.structure
        .Cage` object.

        Returns
        -------
        str
            A string representation of the object.
        """
        return f"<AtomPacker.structure.Cage at {hex(id(self))}>"

    def detect_cavity(
        self,
        step: float = 0.6,
        probe_in: float = 1.4,
        probe_out: float = 10.0,
        removal_distance: float = 2.4,
        volume_cutoff: float = 100.0,
        surface: str = "SES",
        nthreads: Optional[int] = None,
        verbose: bool = False,
        **kwargs: Dict[str, Any],
    ) -> Cavity:
        """
        Detect the cavity in the macromolecular structure.

        Parameters
        ----------
        step : float, optional
            The grid spacing for the cavity detection algorithm.
        probe_in : float, optional
            The radius of the inner probe.
        probe_out : float, optional
            The radius of the outer probe.
        removal_distance : float, optional
            The distance to remove atoms from the cavity.
        volume_cutoff : float, optional
            The volume cutoff for the cavity.
        surface : str, optional
            Surface representation. Keywords options are SES (Solvent Excluded
            Surface) or SAS (Solvent Accessible Surface), by default SES.
        nthreads : int, optional
            Number of threads, by default None. If None, the number of threads
            is `os.cpu_count() - 1`.
        verbose : bool, optional
            Print extra information to standard output, by default False.
        kwargs : Dict[str, Any]
            Additional keyword arguments to pass to the cavity detection
            algorithm (pyKVFinder.detect) of pyKVFinder package.

        Returns
        -------
        Cavity
            The cavity object. The cavity object contains the grid points of
            the cavity, the step size, the Probe In radius, the Probe Out
            radius, the removal distance, the volume cutoff, and the vertices
            (origin, X-axis, Y-axis, Z-axis) of the grid.
            Cavity points in the 3D grid (grid[nx][ny][nz]). Cavities
            array has integer labels in each position, that are:

            * -1: solvent points;

            * 0: cage points;

            * 1: empty space points;

            * >=2: cavity points.

            The empty space points are regions that do not meet the chosen
            volume cutoff to be considered a cavity.

        Raises
        ------
        ValueError
            If the cage is not loaded.
        """
        if self.universe is None:
            raise ValueError("No cage loaded. Please run load() first.")

        # Get the vertices (origin, X-axis, Y-axis, Z-axis) of the cavity grid
        vertices = get_vertices(self.atomic, probe_out=probe_out, step=step)

        # Detect the cavity in the cage structure
        ncav, cavities = detect(
            atomic=self.atomic,
            vertices=vertices,
            step=step,
            probe_in=probe_in,
            probe_out=probe_out,
            removal_distance=removal_distance,
            volume_cutoff=volume_cutoff,
            surface=surface,
            nthreads=nthreads,
            verbose=verbose,
            **kwargs,
        )

        # Check the number of cavities detected
        if ncav < 1:
            warnings.warn(
                "No cavity detected. Please run detect_cavity() \
first."
            )

        if ncav > 1:
            warnings.warn(
                f"{ncav} cavities detected. Keeping all cavities to futher \
analysis."
            )

        # Assign the cavity object
        self.cavity = Cavity(
            grid=cavities,
            step=step,
            probe_in=probe_in,
            probe_out=probe_out,
            removal_distance=removal_distance,
            volume_cutoff=volume_cutoff,
            vertices=vertices,
            surface=surface,
        )

        # As cavity was successfully detected, clean the cluster object
        if self.cluster is not None:
            self.cluster = None

    def load(self, filename: str) -> None:
        """
        Load a supramolecular cage structure file into the :class:`MDAnalysis
        .Univese` object.

        Parameters
        ----------
        filename : str
            The filename of the structure file. The file format is determined
            by the suffix.  Supported formats are: .cif, .mol2, .pdb, .xyz.

        Returns
        -------
        universe: MDAnalysis.Univese
            The MDAnalysis universe object containing the structure.

        Raises
        ------
        ValueError
            If the file format is not supported.
        """
        # We only need the suffix here
        _, suffix = os.path.splitext(filename)

        # Match the suffix to the appropriate file format
        match suffix:
            case ".cif":
                self.universe = load_mmcif(filename)
            case ".mol2":
                self.universe = load_mol2(filename)
            case ".pdb":
                self.universe = load_pdb(filename)
            case ".xyz":
                self.universe = load_xyz(filename)
            case _:
                raise ValueError(
                    f"Unsupported file format: {suffix}. Supported formats \
are: .cif, .pdb, .xyz, .mol2."
                )

    def pack(
        self,
        atom_type: str,
        lattice_type: str,
        a: Optional[float] = None,
        b: Optional[float] = None,
        c: Optional[float] = None,
    ) -> Cluster:
        """
        Pack the cluster of atoms into the cage structure.

        Parameters
        ----------
        atom_type : str
            The type of atom in the cluster.
        lattice_type : str
            The type of lattice in the cluster. The available lattice types are
            'bcc', 'fcc', 'hcp', 'graphite', 'sc', 'decahedron', 'icosahedron',
            and 'octahedron'. The lattices are based on the `ase.cluster`
            module.
            The `ase.cluster` developers state that the module works properly
            for the three cubic crystal structures: FaceCenteredCubic ('fcc'),
            BodyCenteredCubic ('bcc'), and SimpleCubic ('sc'). Other structures
            like HexagonalClosedPacked ('hcp'), Graphite ('graphite'),
            Decahedron ('decahedron'), Icosahedron ('icosahedron'), and
            Octahedron ('octahedron') are implemented, but currently do not
            work correctly.
        a : float, optional
            The lattice constant `a`. If not specified, the experimental value
            from `ase.data` is used.
        b : float, optional
            The lattice constant `b`. If not specified, the experimental value
            from `ase.data` is used.
        c : float, optional
            The lattice constant `c`. If not specified, the experimental value
            from `ase.data` is used.

        Raises
        ------
        ValueError
            If the cage is not loaded.
        ValueError
            If the cavity is not detected.
        """
        if self.universe is None:
            raise ValueError("No cage loaded. Please run load() first.")

        if self.cavity is None:
            raise ValueError(
                "No cavity detected. Please run detect_cavity() \
first."
            )

        # Get lattice constants
        lattice_constants = get_lattice_constants(
            atom_type, lattice_type, a=a, b=b, c=c
        )

        # Make cluster
        _cluster = self._build_cluster(
            atom_type,
            lattice_type,
            lattice_constants,
            center=self.centroid,
        )

        # Create `AtomPacker.structure.Cluster` object
        self.cluster = Cluster(cluster=_cluster, cavity=self.cavity)

    def preview(
        self,
        show_cavity: bool = False,
        show_cluster: bool = False,
        renderer: str = "browser",
        **kwargs: Dict[str, Any],
    ) -> None:
        """
        Preview the cage system (cage, cavity, and cluster) in a 3D viewer.

        Parameters
        ----------
        show_cavity : bool, optional
            Show the cavity in the 3D viewer, by default False.
        show_cluster : bool, optional
            Show the cluster in the 3D viewer, by default False.
        renderer : str, optional
            The renderer to use for the 3D viewer. Supported renderers are
            'browser' (default), 'notebook' and 'png'.
        **kwargs : Dict[str, Any]
            Additional keyword arguments to pass to the scatter_3d function of
            the plotly.express package.

        Raises
        ------
        ValueError
            If the cage is not loaded.

        Warns
        -----
        UserWarning
            If the cavity is not detected, a warning is issued.
        UserWarning
            If the cluster is not packed, a warning is issued.
        """
        # Check if cage is loaded
        if self.universe is None:
            raise ValueError("No cage loaded. Please run load() first.")
        else:
            coordinates = self.coordinates
            radii = self.universe.atoms.radii
            labels = numpy.full(self.atomic.shape[0], "Cage")
            opacity = numpy.full(self.atomic.shape[0], 1.0)

        # Check if cavity is detect. If not, issue a warning.
        if show_cavity:
            if self.cavity is None:
                warnings.warn(
                    "No cavity detected. To visualize cavity, run \
    detect_cavity() first."
                )
            else:
                coordinates = numpy.vstack([coordinates, self.cavity.coordinates])
                radii = numpy.hstack([radii, self.cavity.universe.atoms.radii])
                labels = numpy.hstack(
                    [labels, numpy.full(self.cavity.coordinates.shape[0], "Cavity")]
                )
                opacity = numpy.hstack(
                    [opacity, numpy.full(self.atomic.shape[0], 0.25)]
                )

        # Check if cluster is packed. If not, issue a warning.
        if show_cluster:
            if self.cluster is None:
                warnings.warn(
                    "No cluster packed. To visualize cluster, run pack() \
first."
                )
            else:
                coordinates = numpy.vstack([coordinates, self.cluster.coordinates])
                radii = numpy.hstack([radii, self.cluster.universe.atoms.radii])
                labels = numpy.hstack(
                    [
                        labels,
                        numpy.full(self.cluster.coordinates.shape[0], "Cluster"),
                    ]
                )
                opacity = numpy.hstack([opacity, numpy.full(self.atomic.shape[0], 1.0)])

        # Preview the cage system in a 3D viewer
        x, y, z = coordinates.T
        fig = scatter_3d(
            x=x,
            y=y,
            z=z,
            size=(radii * 2),
            color=labels,
            color_discrete_map={
                "Cage": "rgba(99, 110, 255, 1)",  # blue
                "Cavity": "rgba(102, 102, 102, 0.25)",  # light gray
                "Cluster": "rgba(239, 85, 59, 1)",  # red
            },
            opacity=1.0,
            labels={"color": "Structures"},
            **kwargs,
        )
        # Update layout
        factor = 0.3
        fig.update_layout(
            scene=dict(
                xaxis=dict(
                    showticklabels=False,
                    range=[(x.min() - factor * x.ptp()), (x.max() + factor * x.ptp())],
                ),
                yaxis=dict(
                    showticklabels=False,
                    range=[(y.min() - factor * y.ptp()), (y.max() + factor * y.ptp())],
                ),
                zaxis=dict(
                    showticklabels=False,
                    range=[(z.min() - factor * z.ptp()), (z.max() + factor * z.ptp())],
                ),
            )
        )
        # Update marker
        fig.update_traces(
            marker=dict(
                sizeref=0.002,
            )
        )
        fig.show(renderer)

    def _build_cluster(
        self,
        atom_type: str,
        lattice_type: str,
        lattice_constants: (
            Tuple[float, float, float] | Tuple[float, float] | Tuple[float] | None
        ),
        center: numpy.ndarray,
    ) -> ase.cluster.Cluster:
        """
        Build the cluster of atoms inside cavity.

        Parameters
        ----------
        atom_type : str
            The type of atom in the cluster.
        lattice_type : str
        lattice_constants : Tuple[float, float, float] | Tuple[float, float] | Tuple[float] | None
            The lattice constants `a`, `b`, and `c`. If not specified, the
            experimental value from `ase.data` is used.
        center : numpy.ndarray
            The center of the cluster.

        Returns
        -------
        ase.cluster.Cluster
            The cluster of atoms.
        """
        # Create dummy surfaces
        surfaces = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

        # Based on surfaces=[(1, 0, 0), (0, 1, 0), (0, 0, 1)], get layers
        # for from cage coordinates
        layers = self._get_cluster_layers(atom_type, factor=0.2)

        # Create cluster from `ase.cluster` module
        match lattice_type:
            case "bcc":
                cluster = BodyCenteredCubic(
                    symbols=atom_type,
                    surfaces=surfaces,
                    layers=layers,
                    latticeconstant=lattice_constants,
                )
            case "fcc":
                cluster = FaceCenteredCubic(
                    symbols=atom_type,
                    surfaces=surfaces,
                    layers=layers,
                    latticeconstant=lattice_constants,
                )
            case "hcp":
                cluster = HexagonalClosedPacked(
                    symbols=atom_type,
                    surfaces=surfaces,
                    layers=layers,
                    latticeconstant=lattice_constants,
                )
            case "sc":
                cluster = SimpleCubic(
                    symbols=atom_type,
                    surfaces=surfaces,
                    layers=layers,
                    latticeconstant=lattice_constants,
                )

        # Translate the cluster to the center
        cluster.positions += center

        # Filter atoms inside the cavity
        cluster = self._filter_cluster(cluster)

        return cluster

    def _filter_cluster(
        self,
        cluster: ase.cluster.Cluster,
    ) -> ase.cluster.Cluster:
        """
        Filter atoms in the cluster that are inside the cavity.

        Parameters
        ----------
        cavity_coords : numpy.ndarray
            The coordinates of the cavity.
        cluster_coords : numpy.ndarray
            The coordinates of the cluster.

        Returns
        -------
        numpy.ndarray
            The coordinates of the filtered atoms in the cluster that are
            inside the cavity.
        """
        # Calculate atom indices for the cavity grid
        indices = (
            (cluster.positions - self.cavity._vertices[0]) / self.cavity._step
        ).astype(int)

        # Check if indices are with the cavity grid
        mask = (
            (indices[:, 0] >= 0)
            & (indices[:, 0] < self.cavity.grid.shape[0])
            & (indices[:, 1] >= 0)
            & (indices[:, 1] < self.cavity.grid.shape[1])
            & (indices[:, 2] >= 0)
            & (indices[:, 2] < self.cavity.grid.shape[2])
        )

        # Initialize inside array
        inside = numpy.full(indices.shape[0], False)

        # Recover valid indices
        valid_indices = indices[mask]

        # Select atoms inside the cavity
        inside[mask] = (
            self.cavity.grid[
                valid_indices[:, 0], valid_indices[:, 1], valid_indices[:, 2]
            ]
            > 1
        )

        # Remove atoms outside the cavity. It must removes removing
        # from the end of cluster indices.
        outside = numpy.flip(numpy.invert(inside).nonzero()[0])
        for atom_index in outside:
            cluster.pop(atom_index)

        return cluster

    def _get_cluster_layers(self, atom_type: str, factor: float = 0.2) -> numpy.ndarray:
        """
        Get the number of layers for the cluster.

        Parameters
        ----------
        factor : float, optional
            The factor to multiply the number of layers, by default 0.2.

        Returns
        -------
        numpy.ndarray
            The number of layers for the cluster.
        """
        # Get the length of the xyz coordinates of the cage
        lengths = self.coordinates.ptp(axis=0)

        # Add factor to lengths
        lengths += lengths * factor

        # Get van der Waals radii of atom
        atomic_number = atomic_numbers[atom_type]
        atom_size = covalent_radii[atomic_number] * 2

        # Get number of layers
        layers = numpy.ceil(lengths / atom_size).astype(int)

        return layers

    @property
    def atomic(self) -> numpy.ndarray:
        """
        Get the coordinates of the atoms in the molecular structure.

        Returns
        -------
        numpy.ndarray
            An array containing the atomic information. The array has the
            following columns: resnum, chain, resname, element, x, y, z,
            radius.
        """
        if self.universe is not None:
            return numpy.c_[
                self.universe.atoms.resnums,  # resnums
                self.universe.atoms.chainIDs,  # chains
                self.universe.atoms.resnames,  # resnamess
                self.universe.atoms.elements,  # atom type
                self.universe.atoms.positions,  # x, y, z
                self.universe.atoms.radii,  # radius
            ]

    @property
    def centroid(self) -> numpy.ndarray:
        """
        Get the centroid of the cage structure.

        Returns
        -------
        numpy.ndarray
            An array containing the xyz coordinates of the centroid of the
            cage structure.

        Raises
        ------
        ValueError
            If the cage is not loaded.
        """
        if self.universe is None:
            raise ValueError("No cage loaded. Please run load() first.")
        return self.coordinates.mean(axis=0)

    @property
    def coordinates(self) -> numpy.ndarray:
        """
        Get the coordinates of the atoms in the molecular structure.

        Returns
        -------
        numpy.ndarray
            An array containing the atomic coordinates. The array has the
            following columns: x, y, z.
        """
        if self.universe is not None:
            return self.universe.atoms.positions
