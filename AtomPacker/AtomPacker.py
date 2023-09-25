import os
import pathlib
import subprocess
import warnings
from typing import Dict, List, Optional, Union, Any

import MDAnalysis
import numpy
import pyKVFinder

__all__ = ["PackmolStructure", "AtomPacker"]


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


class AtomPacker(object):
    def __init__(
        self,
        smc: PackmolStructure,
        np_atom: PackmolStructure,
        np_atom_radius: Optional[float] = None,
        basedir: Optional[str] = None,
    ):
        # Load Universe of Supramolecular cage
        self.smc = smc

        # Load Universe of nanoparticle atom
        self.np_atom = np_atom

        # Nanoparticle atom radius (Angstroms)
        if np_atom_radius is None:
            self.np_atom_radius = 1.0
        else:
            self.np_atom_radius = np_atom_radius

        # Cavity
        self.cavity = None
        self.ncav = None
        self.parameters = None

        # Boundary
        self.boundary = None

        # AtomPacked structure
        self.packed = None

        # Set basedir for output files
        if basedir is None:
            self.basedir = os.getcwd()
        else:
            self.basedir = basedir

        # Create base directory
        os.makedirs(self.basedir, exist_ok=True)

    def add_boundary(self) -> None:
        """Add boundary to the system"""
        # # Calculate vertices from grid
        # vertices = pyKVFinder.get_vertices(atomic, probe_out, step)

        # # Detect boundaries
        # _, self.boundary, _ = pyKVFinder.openings(
        #     self.cavity, step=self.parameters["step"]
        # )
        pass

    def _make_packmol_input(self):
        # Check if all structures are suitable, fix them if needed
        for structure in [self.smc, self.np_atom]:
            if not hasattr(structure.structure, "resnames"):
                structure.structure.universe.add_TopologyAttr("resnames")

        # Tolerance is twice the radius of the nanoparticle atom
        tolerance = 2 * self.np_atom_radius

        with open(os.path.join(self.basedir, "packmol.inp"), "w") as out:
            out.write("# autogenerated packmol input\n\n")

            # Output file name
            out.write(f"output {os.path.join(self.basedir, 'output.pdb')}\n\n")

            # Tolerance
            out.write(f"tolerance {tolerance}\n\n")

            # Inputs filetype
            out.write("filetype pdb\n\n")

            # Supramolecular cage
            out.write(self.smc.to_packmol_inp(0))
            self.smc.save_structure(0)

            # Nanoparticle atom
            out.write(self.np_atom.to_packmol_inp(1))
            self.np_atom.save_structure(1)

            # Boundary
            # TODO: include boundary to inp file

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
            with open("packmol.stdout", "w") as out:
                out.write(p.stdout.decode())

    def _load_packmol_output(self):
        """Parse the output of Packmol"""
        return MDAnalysis.Universe(
            os.path.join(self.basedir, "output.pdb"), transfer_to_memory=True
        )

    def _reassign_topology(self):
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

        # add required attributes
        for attr in ["types", "names", "charges", "masses"]:
            if any(hasattr(pms.structure, attr) for pms in [self.smc, self.np_atom]):
                self.packed.add_TopologyAttr(attr)

                if not all(
                    hasattr(pms.structure, attr) for pms in [self.smc, self.np_atom]
                ):
                    warnings.warn("added attribute which not all templates had")

        while index < len(self.packed.atoms):
            # first atom we haven't dealt with yet
            start = self.packed.atoms[index]
            # the resname was altered to give a hint to what template it was from
            template = [self.smc, self.np_atom][int(start.resname[1:])].structure
            # grab atomgroup which matches template
            to_change = self.packed.atoms[index : index + len(template.atoms)]

            # Update residue names
            nres = len(template.residues)
            self.packed.residues[
                start.resindex : start.resindex + nres
            ].resnames = template.residues.resnames

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
            self.packed.add_TopologyAttr("bonds", values=bonds)
        if angles:
            angles = [tuple(val) for val in angles]
            self.packed.add_TopologyAttr("angles", values=angles)
        if dihedrals:
            dihedrals = [tuple(val) for val in dihedrals]
            self.packed.add_TopologyAttr("dihedrals", values=dihedrals)
        if impropers:
            impropers = [tuple(val) for val in impropers]
            self.packed.add_TopologyAttr("impropers", values=impropers)

    def detect_cavity(
        self,
        parameters: Dict[str, Any],
        vdw: Optional[Union[str, Dict[str, Dict[str, float]]]] = None,
    ):
        """Detect cavity within supramolecular cage

        Parameters
        ----------
        parameters : dict
            cavity detection parameters
        vdw : str, pathlib.Path, or dict
            path to vdw radii file or a dictionary with vdw radii (default is None)
        """
        # Get van der Waals radii dictionary
        if vdw is None:
            vdw = pyKVFinder.read_vdw()
        elif type(vdw) in [str]:
            vdw = pyKVFinder.read_vdw(vdw)
        elif type(vdw) in [dict]:
            vdw = vdw
        else:
            raise TypeError(
                "`vdw` must be a string of a van der Waals radii from .dat file."
            )

        # Get atomic information
        self.smc.strcture.atoms = self.packed.select_atoms("chainID X")
        radius = []
        for resname, atom, atom_symbol in zip(
            self.smc.structure.atoms.resnames,
            self.smc.structure.atoms.names,
            self.smc.structure.atoms.types,
        ):
            if resname in vdw.keys():
                if atom in vdw[resname].keys():
                    radius.append(vdw[resname][atom])
            else:
                radius.append(vdw["GEN"][atom_symbol.upper()])
        radius = numpy.array(radius)

        # Get atomic coordinates
        atomic = numpy.c_[
            self.smc.structure.atoms.resnums,
            self.smc.structure.atoms.chainIDs
            if "chainIDs" in self.smc.structure.atoms._SETATTR_WHITELIST
            else numpy.full(self.smc.structure.atoms.resnums.shape, ""),
            self.smc.structure.atoms.resnames,
            self.smc.structure.atoms.names,
            self.smc.structure.atoms.positions,
            radius,
        ]  # shape (n_atoms, 7)

        # Calculate vertices from grid
        vertices = pyKVFinder.get_vertices(
            atomic,
            parameters["probe_out"],
            parameters["step"],
        )

        # Detect cavity
        self.ncav, self.cavity = pyKVFinder.detect(
            atomic,
            vertices,
            parameters["step"],
            parameters["probe_in"],
            parameters["probe_out"],
            parameters["removal_distance"],
            parameters["volume_cutoff"],
        )

    def filter_inside_cavity(self):
        if self.cavity is None:
            raise ValueError("Cavity is not detected. Run `detect_cavity` first.")
        else:
            # Get nanoparticle atoms
            atoms = self.packed.select_atoms("chainID A").positions

            # Calculate the grid indices for each atom
            indexes = (
                (atoms - self.parameters["vertices"][0]) / self.parameters["step"]
            ).astype(int)

    def packing(
        self,
    ):
        """Run Packmol"""
        # Make packmol.inp file
        self._make_packmol_input()
        # Run packmol: packmol < packmol.inp
        self._run_packmol()
        # Load output.pdb file and reassign topology
        self.packed = self._load_packmol_output()
        self._reassign_topology()
        # Detect cavity

        # Filter nanoparticle atoms inside cavity

        # Write packed structure to file
        self.packed.atoms.write(os.path.join(self.basedir, "packed.pdb"))


if __name__ == "__main__":
    SMC = os.path.join("../", "data", "C1.pdb")
    AU = os.path.join("../", "data", "Au.pdb")

    detection_parameters = {
        "step": 0.25,
        "probe_in": 1.4,
        "probe_out": 10.0,
        "removal_distance": 1.0,
        "volume_cutoff": 5.0,
    }

    # Load Universe of Supramolecular cage (smc)
    smc = PackmolStructure(
        SMC,
        number=1,
        instructions=["center", "fixed 0. 0. 0. 0. 0. 0."],
    )
    # Load Universe of Nanoparticle atom (np_atom)
    np_atom = PackmolStructure(
        AU, number=40, instructions=["inside sphere 0. 0. 0. 7."]
    )

    # Prepare AtomPacker object
    ap = AtomPacker(smc, np_atom, np_atom_radius=1.36, basedir="tests")

    # Run Packmol
    ap.packing()

    # Filter cavity
