import math
import os
import pathlib
import subprocess
import warnings
from typing import Dict, List, Optional, Union

import MDAnalysis
import numpy
import pandas
import pyKVFinder
from scipy.spatial.distance import pdist, squareform

__all__ = ["PackmolStructure", "CavityDetector", "AtomPacker"]


class PackmolStructure(object):
    """A Molecule to add to the Packmol system

    Parameters
    ----------
    structure : str
      a single template molecule for Packmol to use
    number : int
      quantity of this molecule to add to new system
    instructions : List[str]
      list of instructions to Packmol for this molecule
      eg 'inside box 0. 0. 0. 40. 40. 40.'
      each item in the list should be a single line instruction
    """

    def __init__(
        self,
        structure: Union[str, pathlib.Path],
        number: int = 1,
        instructions: List[str] = None,
    ):
        self.structure = MDAnalysis.Universe(structure)
        self.number = number
        self.instructions = instructions

    def to_packmol_inp(self, index):
        """Create portion of packmol.inp file from this molecule

        Parameters
        ----------
        index : int
          the index of this template molecule within entire system

        Returns
        -------
        output : str
          text to write to Packmol input for this molecule
        """
        output = f"structure {index}.pdb\n"
        output += f"  number {self.number}\n"
        for instruction in self.instructions:
            output += "  " + instruction + "\n"
        output += "end structure\n\n"

        return output

    def save_structure(self, index):
        """Save this molecule for Packmol to use

        Parameters
        ----------
        index : int
          the index of this template molecule within entire system
        """
        # we mangle Resnames to keep track of which molecule is which
        # so copy the true names, change, write out, change back to real
        old_resnames = self.structure.residues.resnames.copy()
        self.structure.residues.resnames = f"R{index}"
        with MDAnalysis.Writer(f"{index}.pdb") as w:
            w.write(self.structure)
        self.structure.residues.resnames = old_resnames


class CavityDetector(object):
    """Detect cavity within supramolecular cage

    Parameters
    ----------
    vdw : str, pathlib.Path, or dict
        path to vdw radii file or a dictionary with vdw radii (default is None)
    verbose : bool
        print messages (default is False)
    """

    def __init__(
        self,
        step: float = 0.6,
        probe_in: float = 1.4,
        probe_out: float = 4.0,
        removal_distance: float = 2.4,
        volume_cutoff: float = 5.0,
        vdw: Optional[Union[str, pathlib.Path, Dict[str, Dict[str, float]]]] = None,
        verbose: bool = False,
    ):
        # Verbose mode
        self.verbose = verbose

        # Load van der Waals radii dictionary
        self.vdw = self._load_vdw(vdw)

        # Cavity detection parameters
        self.step = step
        self.probe_in = probe_in
        self.probe_out = probe_out
        self.removal_distance = removal_distance
        self.volume_cutoff = volume_cutoff

        # Atomic information
        self.atomic = None

        # Vertices
        self.vertices = None

    def _load_vdw(
        self, vdw: Optional[Union[str, pathlib.Path, Dict[str, Dict[str, float]]]]
    ):
        # Get van der Waals radii dictionary
        if self.verbose:
            print("> Loading atomic dictionary file")

        if vdw is None:
            return pyKVFinder.read_vdw()
        elif type(vdw) in [str, pathlib.Path]:
            return pyKVFinder.read_vdw(vdw)
        elif type(self.vdw) in [dict]:
            return vdw
        else:
            raise TypeError(
                "`vdw` must be a string or a pathlib.Path of a van der Waals radii from .dat file."
            )

    def _get_atomic(self, smc: MDAnalysis.Universe):
        if self.verbose:
            print("> Getting atomic coordinates")

        # Get radius of atoms
        radius = []
        for resname, atom, atom_symbol in zip(
            smc.atoms.resnames,
            smc.atoms.names,
            smc.atoms.types,
        ):
            if resname in self.vdw.keys():
                if atom in self.vdw[resname].keys():
                    radius.append(self.vdw[resname][atom])
            else:
                radius.append(self.vdw["GEN"][atom_symbol.upper()])
        radius = numpy.array(radius)

        # Get atomic coordinates
        atomic = numpy.c_[
            smc.atoms.resnums,
            smc.atoms.chainIDs
            if "chainIDs" in smc.atoms._SETATTR_WHITELIST
            else numpy.full(smc.atoms.resnums.shape, ""),
            smc.atoms.resnames,
            smc.atoms.names,
            smc.atoms.positions,
            radius,
        ]  # shape (n_atoms, 7)

        return atomic

    def _get_vertices(self):
        if self.verbose:
            print("> Getting vertices")

        # Calculate vertices from grid
        vertices = pyKVFinder.get_vertices(self.atomic, self.probe_out, self.step)

        return vertices

    def _detect_cavity(self) -> numpy.ndarray:
        if self.verbose:
            print("> Detecting cavities")

        # Calculate vertices from grid
        self.vertices = self._get_vertices()

        # Detect cavity
        ncav, cavities = pyKVFinder.detect(
            self.atomic,
            self.vertices,
            self.step,
            self.probe_in,
            self.probe_out,
            self.removal_distance,
            self.volume_cutoff,
        )

        # Print number of cavities
        if self.verbose:
            print(f"[==> Number of cavities: {ncav}")

        return cavities

    def _detect_boundary(self, np_atom_radius: float) -> numpy.ndarray:
        if self.verbose:
            print("> Detecting boundaries")

        # Get vertices
        self.vertices = pyKVFinder.get_vertices(self.atomic, self.probe_out, self.step)

        # Calculate new removal distance
        rd = max(
            0.0,
            math.floor((self.removal_distance - np_atom_radius) / self.step)
            * self.step,
        )

        # Detect cavities with a "higher" boundary
        _, cavities = pyKVFinder.detect(
            self.atomic,
            self.vertices,
            self.step,
            self.probe_in,
            self.probe_out,
            rd,
            self.volume_cutoff,
        )

        # Detect boundaries
        nboundaries, boundaries, _ = pyKVFinder.openings(cavities, step=self.step)

        # Print number of boundaries
        if self.verbose:
            print(f"[==> Number of boundaries: {nboundaries}")

        return boundaries


class AtomPacker(object):
    def __init__(
        self,
        smc: PackmolStructure,
        np_atom: PackmolStructure,
        np_atom_radius: float,
        cavity_detector: CavityDetector,
        basedir: Optional[str] = None,
    ):
        # Load Universe of Supramolecular cage
        self.smc = smc

        # Load Universe of nanoparticle atom
        self.np_atom = np_atom

        # Nanoparticle atom radius (Angstroms)
        self.np_atom_radius = np_atom_radius

        # CavityDetector object
        self.cavity_detector = cavity_detector

        # Cavity
        self.cavity = None

        # Boundary
        self.boundary = None

        # AtomPacked structure
        self.packed = None

        # AtomPacking summary
        self.summary = None

        # Set basedir for output files
        if basedir is None:
            self.basedir = os.getcwd()
        else:
            self.basedir = basedir

        # Create base directory
        os.makedirs(self.basedir, exist_ok=True)

    def add_boundary(self) -> None:
        """Add boundary to the system"""
        self.cavity_detector.atomic = self.cavity_detector._get_atomic(
            self.smc.structure
        )

        # Export boundary
        pyKVFinder.export(
            os.path.join(self.basedir, "boundary.pdb"),
            self.cavity_detector._detect_boundary(self.np_atom_radius),
            None,
            self.cavity_detector.vertices,
            self.cavity_detector.step,
        )

        # Load boundary as a PackmolStructure
        self.boundary = PackmolStructure(
            os.path.join(self.basedir, "boundary.pdb"),
            number=1,
            instructions=[
                "center",
                f"radius {self.cavity_detector.step / 2}",
                "fixed 0. 0. 0. 0. 0. 0.",
            ],
        )

    def _make_packmol_input(self):
        if self.boundary is None:
            structures = [self.smc, self.np_atom]
        else:
            structures = [self.smc, self.np_atom, self.boundary]

        # Check if all structures are suitable, fix them if needed
        for structure in structures:
            if not hasattr(structure.structure, "resnames"):
                structure.structure.universe.add_TopologyAttr("resnames")

        # Tolerance is twice the radius of the nanoparticle atom
        tolerance = 2 * self.np_atom_radius

        with open(os.path.join(self.basedir, "packmol.inp"), "w") as out:
            out.write("# autogenerated packmol input\n\n")

            # Output file name
            out.write(f"output {os.path.join(self.basedir, 'output.pdb')}\n\n")

            # Tolerance
            out.write(f"tolerance {tolerance + 0.001}\n\n")

            # Inputs filetype
            out.write("filetype pdb\n\n")

            # Define a random runtime seed
            out.write("seed -1\n\n")

            # Supramolecular cage
            out.write(self.smc.to_packmol_inp(0))
            self.smc.save_structure(0)

            # Nanoparticle atom
            out.write(self.np_atom.to_packmol_inp(1))
            self.np_atom.save_structure(1)

            # Boundary
            if self.boundary is not None:
                out.write(self.boundary.to_packmol_inp(2))
                self.boundary.save_structure(2)

    def _run_packmol(self):
        """Run and check that Packmol worked correctly"""
        try:
            p = subprocess.run(
                "packmol < {}".format(os.path.join(self.basedir, "packmol.inp")),
                check=True,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except subprocess.CalledProcessError as e:
            raise ValueError(
                "Packmol failed with errorcode {}"
                " and stderr: {}".format(e.returncode, e.stderr)
            )
        else:
            with open(os.path.join(self.basedir, "packmol.stdout"), "w") as out:
                out.write(p.stdout.decode())

    def _load_packmol_output(self):
        """Parse the output of Packmol"""
        return MDAnalysis.Universe(
            os.path.join(self.basedir, "output.pdb"), transfer_to_memory=True
        )

    def _reassign_topology(self, replicate: int):
        """Take Packmol created Universe and add old topology features back in

        Attempts to reassign:
        - types
        - names
        - charges
        - masses
        - bonds
        - angles
        - torsions
        - impropers
        - resnames

        Parameters
        ----------
        structures : list
        list of Packmol structures used to create the Packmol Universe
        new : Universe
        the raw output of Packmol

        Returns
        -------
        new : Universe
        the raw output modified to best match the templates given to it
        """
        index = 0

        bonds = []
        angles = []
        dihedrals = []
        impropers = []

        if self.boundary is None:
            structures = [self.smc, self.np_atom]
        else:
            structures = [self.smc, self.np_atom, self.boundary]

        # add required attributes
        for attr in ["types", "names", "charges", "masses"]:
            if any(hasattr(pms.structure, attr) for pms in structures):
                self.packed[replicate].add_TopologyAttr(attr)

                if not all(hasattr(pms.structure, attr) for pms in structures):
                    warnings.warn("added attribute which not all templates had")

        while index < len(self.packed[replicate].atoms):
            # first atom we haven't dealt with yet
            start = self.packed[replicate].atoms[index]
            # the resname was altered to give a hint to what template it was from
            template = structures[int(start.resname[1:])].structure
            # grab atomgroup which matches template
            to_change = self.packed[replicate].atoms[
                index : index + len(template.atoms)
            ]

            # # Update residue names
            # nres = len(template.residues)
            # self.packed[replicate].residues[
            #     start.resindex : start.resindex + nres
            # ].resnames = template.residues.resnames

            # atom attributes
            for attr in ["types", "names", "charges", "masses"]:
                if hasattr(template.atoms, attr):
                    setattr(to_change, attr, getattr(template.atoms, attr))

            # bonds/angles/torsions
            if hasattr(template, "bonds"):
                bonds.extend((template.bonds.to_indices() + index).tolist())
            if hasattr(template, "angles"):
                angles.extend((template.angles.to_indices() + index).tolist())
            if hasattr(template, "dihedrals"):
                dihedrals.extend((template.dihedrals.to_indices() + index).tolist())
            if hasattr(template, "impropers"):
                impropers.extend((template.impropers.to_indices() + index).tolist())

            # update the index pointer to be on next unknown atom
            index += len(template.atoms)

        if bonds:
            # convert to tuples for hashability
            bonds = [tuple(val) for val in bonds]
            self.packed[replicate].add_TopologyAttr("bonds", values=bonds)
        if angles:
            angles = [tuple(val) for val in angles]
            self.packed[replicate].add_TopologyAttr("angles", values=angles)
        if dihedrals:
            dihedrals = [tuple(val) for val in dihedrals]
            self.packed[replicate].add_TopologyAttr("dihedrals", values=dihedrals)
        if impropers:
            impropers = [tuple(val) for val in impropers]
            self.packed[replicate].add_TopologyAttr("impropers", values=impropers)

    def _update_smc(self):
        # Select supramolecular cage in packed structure
        self.smc.structure.atoms = self.packed[0].select_atoms("resname R0")

    def _detect_cavity(self):
        # Detect cavity of packed structure
        self.cavity_detector.atomic = self.cavity_detector._get_atomic(
            self.smc.structure
        )
        self.cavity = self.cavity_detector._detect_cavity()

        # Export cavity
        pyKVFinder.export(
            os.path.join(self.basedir, "cavity.pdb"),
            self.cavity,
            None,
            self.cavity_detector.vertices,
            self.cavity_detector.step,
        )

    def _filter_atoms_inside_cavity(self, replicate: int):
        # Get nanoparticle atoms
        atoms = self.packed[replicate].select_atoms("resname R1").positions

        # Calculate the grid indices for each atom
        indexes = (
            (atoms - self.cavity_detector.vertices[0]) / self.cavity_detector.step
        ).astype(int)

        # Define mask
        mask = self.cavity[indexes[:, 0], indexes[:, 1], indexes[:, 2]] > 1

        # Select the atoms in chain A that are inside the cavities
        inside_cavity = self.packed[replicate].select_atoms("resname R1")[mask]

        # Select the smc in chain X
        smc = self.packed[replicate].select_atoms("resname R0")

        # Merge the two selections
        self.packed[replicate].atoms = inside_cavity.union(smc)
        self.packed[replicate].atoms.write(
            os.path.join(self.basedir, f"packed{replicate}.pdb")
        )

    def _clean_tempfiles(self):
        for f in ["0.pdb", "1.pdb", "2.pdb", os.path.join(self.basedir, "output.pdb")]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass

    def packing(self, replicates: int = 1):
        """Run Packmol"""

        if replicates < 1:
            raise ValueError("replicates must be > 0")

        # Make packmol.inp file
        self._make_packmol_input()

        # Create an empty list of replicates
        self.packed = [None] * replicates

        for replicate in range(replicates):
            # Run packmol: packmol < packmol.inp
            self._run_packmol()

            # Load output.pdb file and reassign topology
            self.packed[replicate] = self._load_packmol_output()
            self._reassign_topology(replicate)

        # Update supramolecular cage positions
        self._update_smc()

        # Detect cavity
        self._detect_cavity()

        for replicate in range(replicates):
            # Filter atoms inside cavity
            self._filter_atoms_inside_cavity(replicate)

        # Clean temporary files
        self._clean_tempfiles()

        # Pandas DataFrame with atoms packed in each replicate
        self.summary = self._summary(replicates)
        self.summary.to_csv(os.path.join(self.basedir, "PackedAtoms.csv"))

    def _summary(self, replicates: int):
        # Packed atoms
        natoms = []
        for replicate in range(replicates):
            natoms.append(len(self.packed[replicate].select_atoms("resname R1")))

        # Number of clashes and clashing atoms
        nclash, nclashing_atoms = [], []
        for replicate in range(replicates):
            # Get positions
            positions = self.packed[replicate].select_atoms("resname R1").atoms.positions

            # Get clashes
            distances = pdist(positions)
            clashes = numpy.sum(distances < self.np_atom_radius * 2)
            nclash.append(clashes)

            # Get clashing atoms
            distances = squareform(distances)
            numpy.fill_diagonal(distances, numpy.inf)
            clashing_atoms = numpy.unique(
                numpy.where(distances < (self.np_atom_radius * 2))[0]
            ).shape[0]
            nclashing_atoms.append(clashing_atoms)

            # Save clashing pairs to file with distance
            clashing_pairs = numpy.argwhere(distances < (self.np_atom_radius * 2))
            clashing_pairs = clashing_pairs[clashing_pairs[:, 0] < clashing_pairs[:, 1]]
            if clashing_pairs.shape[0] > 0:
                with open(os.path.join(self.basedir, f"clashing_pairs{replicate}.err"), 'w') as f:
                    for pair, distance in zip(clashing_pairs, distances[clashing_pairs[:, 0], clashing_pairs[:, 1]]):
                        f.write(f"{pair[0]} {pair[1]} {distance}\n")

        return pandas.DataFrame(
            [natoms, nclash, nclashing_atoms],
            index=["Packed atoms", "Number of clashs", "Clashing atoms"],
            columns=[f"packed{replicate}.pdb" for replicate in range(replicates)],
        ).T
