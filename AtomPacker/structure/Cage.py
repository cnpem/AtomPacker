# This source code is part of the AtomPacker package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
A subpackage of *AtomPacker* that contains the `Cage` class.

Macromolecular structure files (PDB, PDBx/mmCIF, XYZ, etc.) can be used to
create a :class:`AtomPacker.structure.Cage` object, which is a container for
the atoms of the structure.
"""

import itertools
import os
import warnings
from copy import deepcopy

import ase.cluster
import numpy
import pandas
from ase.cluster.cubic import BodyCenteredCubic, FaceCenteredCubic, SimpleCubic
from ase.cluster.hexagonal import HexagonalClosedPacked
from ase.data import atomic_numbers, covalent_radii
from joblib import Parallel, delayed
from MDAnalysis import Universe
from plotly.express import scatter_3d
from pyKVFinder import detect, get_vertices
from sklearn.decomposition import PCA
from tqdm import tqdm

from ..core import load_mmcif, load_mol2, load_pdb, load_xyz
from .Cavity import Cavity
from .Cluster import Cluster
from .data import get_lattice_constants


class Cage:
    """
    A class representing a supramolecular cage loaded from a structure file.

    The :class:`AtomPacker.Cage` class is used to store the atoms of a
    macromolecular structure.
    """

    def __init__(self):
        """
        Create a new :class:`AtomPacker.structure.Cage` object.
        """
        self.universe: Universe | None = None
        self.cavity: Cavity | None = None
        self.cluster: Cluster | None = None

    def __repr__(self) -> str:
        return f"<AtomPacker.structure.Cage at {hex(id(self))}>"

    @property
    def atomic(self) -> numpy.ndarray | None:
        """
        Get the coordinates of the atoms in the molecular structure.

        Returns
        -------
        numpy.ndarray | None
            An array containing the atomic information. The array has the
            following columns: resnum, chain, resname, element, x, y, z,
            radius. If the cage is not loaded, returns None.
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
    def centroid(self) -> numpy.ndarray | None:
        """
        Get the centroid of the cage structure.

        Returns
        -------
        numpy.ndarray | None
            An array containing the xyz coordinates of the centroid of the
            cage structure. If the cage is not loaded, returns None.

        Raises
        ------
        ValueError
            If the cage is not loaded.
        """
        if self.universe is not None:
            return self.coordinates.mean(axis=0)

    @property
    def coordinates(self) -> numpy.ndarray | None:
        """
        Get the coordinates of the atoms in the molecular structure.

        Returns
        -------
        numpy.ndarray | None
            An array containing the atomic coordinates. The array has the
            following columns: x, y, z. If the cage is not loaded, returns None.
        """
        if self.universe is not None:
            return self.universe.atoms.positions

    def detect_cavity(
        self,
        step: float = 0.6,
        probe_in: float = 1.4,
        probe_out: float = 10.0,
        removal_distance: float = 2.4,
        volume_cutoff: float = 100.0,
        surface: str = "SES",
        nthreads: int | None = None,
        verbose: bool = False,
        **kwargs: dict[str, object],
    ) -> None:
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
            is `os.cpu_count() - 1`. If -1, uses all available threads.
        verbose : bool, optional
            Print extra information to standard output, by default False.
        kwargs : dict[str, object], optional
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

    def load(self, filename: str, vdw: dict[str, float] | None = None) -> None:
        """
        Load a supramolecular cage structure file into the :class:`MDAnalysis
        .Univese` object.

        Parameters
        ----------
        filename : str
            The filename of the structure file. The file format is determined
            by the suffix.  Supported formats are: .cif, .mol2, .pdb, .xyz.
        vdw : dict[str, float] | None, optional
            A dictionary containing the van der Waals radii for each atom type,
            by default None. If None, the van der Waals radii are looked up
            from the `pyKVFinder` package.

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
        _, ext = os.path.splitext(filename)

        # Match the suffix to the appropriate file format
        match ext:
            case ".cif":
                self.universe = load_mmcif(filename, vdw)
            case ".mol2":
                self.universe = load_mol2(filename, vdw)
            case ".pdb":
                self.universe = load_pdb(filename, vdw)
            case ".xyz":
                self.universe = load_xyz(filename, vdw)
            case _:
                raise ValueError(
                    f"Unsupported file format: {ext}. Supported formats \
are: .cif, .pdb, .xyz, .mol2."
                )

    def pack(
        self,
        atom_type: str,
        lattice_type: str,
        a: float | None = None,
        b: float | None = None,
        c: float | None = None,
        clashing_tolerance: float = 0.0,
        angles: numpy.ndarray | list[float] = [0.0],
        translations: numpy.ndarray | list[float] = [0.0],
        optsave: bool = False,
        optdir: str | None = None,
        nthreads: int | None = None,
        verbose: bool = False,
    ) -> None:
        """
        Pack a cluster of atoms into the cage structure.

        Parameters
        ----------
        atom_type : str
            Atom type for the cluster.
        lattice_type : str
            Lattice type for the cluster. Supported types: 'bcc', 'fcc', 'hcp',
            'sc'.
        a : float | None, optional
            Lattice constant 'a'. If None, uses values from AtomPacker or ASE.
        b : float | None, optional
            Lattice constant 'b'. If None, uses values from AtomPacker or ASE.
        c : float | None, optional
            Lattice constant 'c'. If None, uses values from AtomPacker or ASE.
        clashing_tolerance : float, optional
            Minimum allowed distance (Å) between cluster and cage atoms.
            Default is 0.0.
        angles : numpy.ndarray | list[float], optional
            Rotation angles for cluster optimization. If not specified, no
            optimization is performed. Default is [0.0].
            If specified, angles should be a list or numpy array of angles in
            degrees. Example: [-75, -50, -25, 0, 25, 50, 75].
        translations : numpy.ndarray | list[float], optional
            Translation values for cluster optimization. If not specified, no
            optimization is performed. Default is [0.0].
            If specified, translations should be a list or numpy array of
            translation values in Angstroms. Example: [-0.2, 0.0, 0.2].
        optsave : bool, optional
            If True, saves each optimization step as a PDB file. Default is
            False.
        optdir : str | None, optional
            Directory to save files. If None, uses current working directory.
        nthreads : int | None, optional
            Number of threads to use for parallel processing. If None, uses
            `os.cpu_count() - 1`. Default is None. If -1, uses all available
            threads.
        verbose : bool, optional
            If True, prints detailed information during processing (default is
            False).

        Raises
        ------
        ValueError
            If the cage is not loaded, cavity is not detected, or
            clashing_tolerance < 0.
        """
        if self.universe is None:
            raise ValueError("No cage loaded. Please run load() first.")

        if self.cavity is None:
            raise ValueError("No cavity detected. Please run detect_cavity() first.")

        if nthreads is None:
            nthreads = os.cpu_count() - 1

        # Get lattice constants
        if lattice_type == "hcp":
            if (a is None) and (c is None):
                lattice_constants = get_lattice_constants(atom_type, lattice_type)
                if verbose:
                    print(
                        f"> Using default lattice constants for {atom_type} in {lattice_type} lattice: {lattice_constants}"
                    )
            elif (a is None) or (c is None):
                if a is None:
                    a, _ = get_lattice_constants(atom_type, lattice_type)
                    if verbose:
                        print(
                            f"> Using default lattice constant 'a' for {atom_type} in {lattice_type} lattice: {a}"
                        )
                if c is None:
                    _, c = get_lattice_constants(atom_type, lattice_type)
                    if verbose:
                        print(
                            f"> Using default lattice constant 'c' for {atom_type} in {lattice_type} lattice: {c}"
                        )
                lattice_constants = (a, c)
            else:
                lattice_constants = (a, c)
        elif lattice_type in ["fcc", "bcc", "sc"]:
            if a is None:
                lattice_constants = get_lattice_constants(atom_type, lattice_type)
                if verbose:
                    print(
                        f"> Using default lattice constant for {atom_type} in {lattice_type} lattice: {lattice_constants}"
                    )
            else:
                lattice_constants = a
        else:
            raise ValueError(
                f"Unsupported lattice type: {lattice_type}. Supported \
                lattice types are: 'bcc', 'fcc', 'hcp', and 'sc'."
            )

        # Check if clashing tolerance is greater than or equal to 0
        if clashing_tolerance < 0:
            raise ValueError("Clashing tolerance must be greater than or equal to 0.")

        # Optimize cluster packing
        _cluster, log = self._build_cluster(
            atom_type,
            lattice_type,
            lattice_constants=lattice_constants,
            center=self.centroid,
            clashing_tolerance=clashing_tolerance,
            angles=angles,
            translations=translations,
            optsave=optsave,
            optdir=optdir,
            nthreads=nthreads,
            verbose=verbose,
        )

        # Create `AtomPacker.structure.Cluster` object
        self.cluster = Cluster(cluster=_cluster, cavity=self.cavity, log=log)

    def preview(
        self,
        show_cavity: bool = False,
        show_cluster: bool = False,
        show_openings: bool = False,
        renderer: str = "browser",
        **kwargs: dict[str, object],
    ) -> None:
        """
        Preview the cage system (cage, cavity, and cluster) in a 3D viewer.

        Parameters
        ----------
        show_cavity : bool, optional
            Show the cavity in the 3D viewer, by default False.
        show_cluster : bool, optional
            Show the cluster in the 3D viewer, by default False.
        show_openings : bool, optional
            Show the openings in the cavity in the 3D viewer, by default False.
        renderer : str, optional
            The renderer to use for the 3D viewer. Supported renderers are
            'browser' (default), 'notebook' and 'png'.
        **kwargs : dict[str, object]
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
        UserWarning
            If the openings are not detected, a warning is issued.
        """
        # Check if cage is loaded
        if self.universe is None:
            raise ValueError("No cage loaded. Please run load() first.")

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

        # Check if openings are detected. If not, issue a warning.
        if show_openings:
            if self.cavity.openings is None:
                warnings.warn(
                    "No openings detected. To visualize openings, run \
detect_openings() first."
                )
            else:
                coordinates = numpy.vstack(
                    [coordinates, self.cavity.openings.coordinates]
                )
                radii = numpy.hstack([radii, self.cavity.openings.universe.atoms.radii])
                labels = numpy.hstack(
                    [
                        labels,
                        numpy.full(
                            self.cavity.openings.coordinates.shape[0], "Openings"
                        ),
                    ]
                )
                opacity = numpy.hstack(
                    [
                        opacity,
                        numpy.full(self.cavity.openings.coordinates.shape[0], 1.0),
                    ]
                )

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
                "Openings": "rgba(0, 0, 0, 1)",  # black
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
                sizeref=0.002,
            )
        )
        fig.show(renderer)

    def _build_cluster(
        self,
        atom_type: str,
        lattice_type: str,
        lattice_constants: tuple[float, float] | tuple[float] | None,
        center: numpy.ndarray,
        clashing_tolerance: float = 0.0,
        angles: numpy.ndarray | list[float] = [0.0],
        translations: numpy.ndarray | list[float] = [0.0],
        optsave: bool = False,
        optdir: str | None = None,
        nthreads: int | None = None,
        verbose: bool = False,
    ) -> tuple[ase.cluster.Cluster, pandas.DataFrame]:
        """
        Build the cluster of atoms inside cavity.

        Parameters
        ----------
        atom_type : str
            Atom type for the cluster.
        lattice_type : str
            Lattice type for the cluster. Supported types: 'bcc', 'fcc', 'hcp',
            'sc'.
        lattice_constants : Tuple[float, float] | Tuple[float] | None
            The lattice constants `a`, `b`, and `c`. If not specified,
            the lattice constants will be fetched from `AtomPacker.data
            .lattice_constants` if available. If not, the experimental
            values from `ase.data` will be used.
        center : numpy.ndarray
            The center of the cluster.
        clashing_tolerance : float, optional
            Minimum allowed distance (Å) between cluster and cage atoms.
            Default is 0.0.
        angles : numpy.ndarray | list[float], optional
            Rotation angles for cluster optimization. If not specified, no
            optimization is performed. Default is [0.0].
            If specified, angles should be a list or numpy array of angles in
            degrees. Example: [-75, -50, -25, 0, 25, 50, 75].
        translations : numpy.ndarray | list[float], optional
            Translation values for cluster optimization. If not specified, no
            optimization is performed. Default is [0.0].
            If specified, translations should be a list or numpy array of
            translation values in Angstroms. Example: [-0.2, 0.0, 0.2].
        optsave : bool, optional
            If True, saves each optimization step as a PDB file. Default is
            False.
        optdir : str | None, optional
            Directory to save files. If None, uses current working directory.
        nthreads : int | None, optional
            Number of threads to use for parallel processing. If None, uses
            `os.cpu_count() - 1`. Default is None. If -1, uses all available
            threads.
        verbose : bool, optional
            If True, prints detailed information during processing (default is
            False).

        Returns
        -------
        ase.cluster.Cluster
            The cluster of atoms.
        pandas.DataFrame
            A DataFrame containing the cluster properties, such as the number
            of atoms, the radius, the lattice constants, and the rotation and
            translation angles.
        """
        if optsave:
            if optdir is None:
                optdir = os.getcwd()
            elif not os.path.exists(optdir):
                os.makedirs(optdir)

        if nthreads is None:
            nthreads = os.cpu_count() - 1

        # Create dummy surfaces
        surfaces = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

        # Create Oriented Bounding Box (OBB) for the cage
        obb_axes, obb_extents = self._get_obb()

        # Based on surfaces=[(1, 0, 0), (0, 1, 0), (0, 0, 1)], get layers
        # for from obb extents
        layers = self._get_cluster_layers(atom_type, obb_extents, factor=0.2)

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
            case _:
                raise ValueError(
                    f"Unsupported lattice type: {lattice_type}. Supported \
                    lattice types are: 'bcc', 'fcc', 'hcp', and 'sc'."
                )

        # Rotate the cluster to align with the OBB's orientation
        cluster.set_positions(numpy.dot(cluster.positions, obb_axes))

        # Translate the cluster to the center (cage centroid)
        cluster.positions += center

        # Get radii for nanocluster
        distances = cluster.get_all_distances()
        radii = (distances[distances > 0] / 2).min()
        cluster.info.update({"radii": radii})

        # Get lattice constants for nanocluster
        cluster.info.update({"lattice_constants": lattice_constants})

        def _evaluate_cluster_configuration(self, phi, theta, psi, x, y, z, cluster):
            """
            Helper function to parallelize the optimization process.
            """
            _tmp = deepcopy(cluster)

            # Rotate and translate the cluster
            _tmp.euler_rotate(center="COP", phi=phi, theta=theta, psi=psi)
            _tmp.translate([x, y, z])

            # Filter atoms inside the cavity
            _tmp = self._filter_outside_cavity(_tmp)

            # Filter atoms clashing with the cage
            _tmp = self._filter_clashing_atoms(
                _tmp, clashing_tolerance=clashing_tolerance
            )

            # Save the cluster if requested
            if optsave:
                filename = f"cluster_{x:.2f}_{y:.2f}_{z:.2f}_{phi:.2f}_{theta:.2f}_{psi:.2f}.pdb"
                _tmp.write(os.path.join(optdir, filename))

            # Get the maximum diameter in the cluster
            if len(_tmp) > 1:
                diameter = _tmp.get_all_distances().max()
            elif len(_tmp) == 1:
                diameter = 2 * _tmp.info.get("radii")
            else:
                diameter = 0.0

            # Get maximum number of atoms in the cluster
            maximum_number_of_atoms = int(
                numpy.ceil(
                    self.cavity.volume
                    / ((4 / 3) * numpy.pi * self.universe.atoms.radii[0] ** 3)
                )
            )

            # Get number of atoms in the cluster
            number_of_atoms = len(_tmp)

            return {
                "x": x,
                "y": y,
                "z": z,
                "phi": phi,
                "theta": theta,
                "psi": psi,
                "Atom Type": atom_type,
                "Atom Radius": radii,
                "Cavity Volume (Å³)": self.cavity.volume,
                "Maximum diameter (Å)": diameter,
                "Shape diameter (Å)": (
                    _tmp.get_diameter(method="shape") if len(_tmp) > 0 else 0.0
                ),
                "Volume diameter (Å)": (
                    _tmp.get_diameter(method="volume") if len(_tmp) > 0 else 0.0
                ),
                "Lattice Constants": lattice_constants,
                "Lattice Type": lattice_type,
                "Maximum Number of Atoms": maximum_number_of_atoms,
                "Number of Atoms": number_of_atoms,
            }

        # Iterate over all possible combinations of angles and translations
        combinations = list(
            itertools.product(
                angles, angles, angles, translations, translations, translations
            )
        )

        log = Parallel(n_jobs=nthreads)(
            delayed(_evaluate_cluster_configuration)(
                self,
                phi,
                theta,
                psi,
                x,
                y,
                z,
                cluster,
            )
            for phi, theta, psi, x, y, z in tqdm(
                combinations, desc="> Optimizing cluster"
            )
        )

        # Convert optimization to DataFrame
        log = pandas.DataFrame(log)

        # Get the best cluster based on the volume diameter
        idx = log["Volume diameter (Å)"].idxmax()
        x, y, z, phi, theta, psi = log.loc[idx, ["x", "y", "z", "phi", "theta", "psi"]]

        # Rotate and translate the cluster
        cluster.euler_rotate(center="COP", phi=phi, theta=theta, psi=psi)
        cluster.translate([x, y, z])

        # Filter atoms inside the cavity
        cluster = self._filter_outside_cavity(cluster)

        # Filter atoms clashing with the cage
        cluster = self._filter_clashing_atoms(
            cluster, clashing_tolerance=clashing_tolerance
        )

        return cluster, log

    def _filter_clashing_atoms(
        self,
        cluster: ase.cluster.Cluster,
        clashing_tolerance: float = 0.0,
    ) -> ase.cluster.Cluster:
        """
        Filter atoms in the cluster that are clashing with the cage.

        Parameters
        ----------
        cluster : ase.cluster.Cluster
            The cluster of atoms.
        clashing_tolerance : float, optional
            The clashing tolerance (Å), by default 0.0.
        """
        # Get radii of atoms in the cluster
        # If cluster has more than one atom, get the minimum distance between
        # atoms and use it as the radius of the cluster atom. If cluster has
        # only one atom, get the van der Waals radius of the atom.
        cluster_radii = numpy.full(len(cluster), cluster.info.get("radii"))

        # Get internal limits between cluster atoms and cage atoms
        # (cluster atom radius + cage atom radius)
        limits = (
            cluster_radii[:, numpy.newaxis]
            + self.universe.atoms.radii
            - clashing_tolerance
        )

        # Get distance between cluster positions and cage positions
        distances = numpy.linalg.norm(
            cluster.positions[:, None] - self.coordinates[None], axis=-1
        )

        # Get radii of atoms in the cage
        clashing = (distances < limits).any(axis=1)

        # Remove atoms clashing with the cage
        for atom_index in numpy.flip(clashing.nonzero()[0]):
            cluster.pop(atom_index)

        return cluster

    def _filter_outside_cavity(
        self,
        cluster: ase.cluster.Cluster,
    ) -> ase.cluster.Cluster:
        """
        Filter atoms in the cluster that are outside the cavity.

        Parameters
        ----------
        cluster : ase.cluster.Cluster
            The cluster of atoms.

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

    def _get_cluster_layers(
        self, atom_type: str, obb_extents: numpy.ndarray, factor: float = 0.2
    ) -> numpy.ndarray:
        """
        Get the number of layers for the cluster.

        Parameters
        ----------
        atom_type : str
            The type of atom in the cluster.
        obb_extents : numpy.ndarray
            The extents of the Oriented Bounding Box (OBB).
        factor : float, optional
            The factor to multiply the number of layers, by default 0.2.

        Returns
        -------
        numpy.ndarray
            The number of layers for the cluster.
        """
        # Add factor to lengths
        obb_extents += obb_extents * factor

        # Get van der Waals radii of atom
        atomic_number = atomic_numbers[atom_type]
        atom_size = covalent_radii[atomic_number] * 2

        # Get number of layers
        layers = numpy.ceil(obb_extents / atom_size).astype(int)

        return layers

    def _get_obb(self) -> tuple[numpy.ndarray, numpy.ndarray]:
        """
        Get the Oriented Bounding Box (OBB) of the cage structure.

        Returns
        -------
        numpy.ndarray
            The axes of the OBB.
        numpy.ndarray
            The extents of the OBB.
        """
        # Compute the PCA to get the principal axes
        pca = PCA(n_components=3)
        pca.fit(self.coordinates)

        # Principal axes (directions of the OBB)
        obb_axes = pca.components_

        # Calculate the extents of the OBB
        obb_extents = numpy.ptp(pca.transform(self.coordinates), axis=0)

        return obb_axes, obb_extents
