import os
import pathlib
import sys
from typing import Dict, List, Optional, Union

if (sys.version_info[0] < 3) or (sys.version_info[1] <= 4):
    import subprocess32 as subprocess
else:
    import subprocess

import warnings

import MDAnalysis
import numpy
import pyKVFinder

__all__ = ["packmol", "PackmolStructure"]

PACKMOL_INP = "packmol.inp"  # name of .inp file given to packmol
PACKMOL_STRUCTURE_FILES = "{}.pdb"
PACKMOL_OUT = "output.pdb"


class PackmolError(Exception):
    # if packmol didn't find solution
    pass


class PackmolStructure(object):
    """A Molecule to add to the Packmol system

    Parameters
    ----------
    structure : MDAnalysis.AtomGroup
      a single template molecule for Packmol to use
    number : int
      quantity of this molecule to add to new system
    instructions : list of strings
      list of instructions to Packmol for this molecule
      eg 'inside box 0. 0. 0. 40. 40. 40.'
      each item in the list should be a single line instruction
    """

    def __init__(
        self, structure: MDAnalysis.AtomGroup, number: int, instructions: List[str]
    ):
        self.structure = structure
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
        output = "structure {}\n".format(PACKMOL_STRUCTURE_FILES.format(index))
        output += "  number {}\n".format(self.number)
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
        self.structure.residues.resnames = "R{}".format(index)
        with MDAnalysis.Writer(PACKMOL_STRUCTURE_FILES.format(index)) as w:
            w.write(self.structure)
        self.structure.residues.resnames = old_resnames


class CavityDetector(object):
    def __init__(
        self,
        smc: PackmolStructure,
        vdw: Optional[Union[str, pathlib.Path, Dict[str, Dict[str, float]]]] = None,
        verbose: bool = False,
    ):
        self.verbose = verbose
        self.vdw = self._load_vdw(vdw)
        self.atomic = self._get_atomic(smc)

    def _load_vdw(self, vdw):
        # Get van der Waals radii dictionary
        if self.verbose:
            print("> Loading atomic dictionary file")

        if vdw is None:
            self.vdw = pyKVFinder.read_vdw()
        elif type(vdw) in [str, pathlib.Path]:
            self.vdw = pyKVFinder.read_vdw(vdw)
        elif type(self.vdw) in [dict]:
            self.vdw = vdw
        else:
            raise TypeError(
                "`vdw` must be a string or a pathlib.Path of a van der Waals radii from .dat file."
            )

    def _get_atomic(self, smc):
        if self.verbose:
            print("> Getting atomic coordinates")

        # Get radius of atoms
        radius = []
        for resname, atom, atom_symbol in zip(
            smc.structure.atoms.resnames,
            smc.structure.atoms.names,
            smc.structure.atoms.types,
        ):
            if resname in self.vdw.keys() and atom in self.vdw[resname].keys():
                radius.append(self.vdw[resname][atom])
            else:
                radius.append(self.vdw["GEN"][atom_symbol.upper()])
        radius = numpy.array(radius)

        # Get atomic coordinates
        atomic = numpy.c_[
            smc.structure.atoms.resnums,
            smc.structure.atoms.chainIDs
            if "chainIDs" in smc.structure.atoms._SETATTR_WHITELIST
            else numpy.full(smc.structure.atoms.resnums.shape, ""),
            smc.structure.atoms.resnames,
            smc.structure.atoms.names,
            smc.structure.atoms.positions,
            radius,
        ]  # shape (n_atoms, 7)

        return atomic

    def detect(
        self,
        step: float = 0.6,
        probe_in: float = 1.4,
        probe_out: float = 4.0,
        removal_distance: float = 2.4,
        volume_cutoff: float = 5.0,
    ):
        if self.verbose:
            print("> Detecting cavities")

        # Calculate vertices from grid
        vertices = pyKVFinder.get_vertices(self.atomic, probe_out, step)

        # Detect cavity
        ncav, cavities = pyKVFinder.detect(
            self.atomic,
            vertices,
            step,
            probe_in,
            probe_out,
            removal_distance,
            volume_cutoff,
        )

        return ncav, cavities, boundaries


def make_packmol_input(
    smc: PackmolStructure, np_atom: PackmolStructure, atom_radius: float = None
):
    """Convert the call into a Packmol usable input file

    Parameters
    ----------
    smc : PackmolStruture
      Supramolecular cage
    np_atom : PackmolStructure
      Nanoparticle atom
    atom_radius : float, optional
      Half of Packmol tolerance, ie., minimum distance between nanoparticle
      atoms, defaults to 1.0
    """
    # Check if all structures are suitable, fix them if needed
    for structure in [smc, np_atom]:
        if not hasattr(structure.structure, "resnames"):
            structure.structure.universe.add_TopologyAttr("resnames")

    # Tolerance is twice the radius of the nanoparticle atom
    if atom_radius is None:
        atom_radius = 1.0
    tolerance = 2 * atom_radius

    with open(PACKMOL_INP, "w") as out:
        out.write("# autogenerated packmol input\n\n")

        # Output file name
        out.write("output {}\n\n".format(PACKMOL_OUT))

        # Tolerance
        out.write("tolerance {}\n\n".format(tolerance))

        # Inputs filetype
        out.write("filetype pdb\n\n")

        # Supramolecular cage
        out.write(smc.to_packmol_inp(0))
        smc.save_structure(0)

        # Nanoparticle atom
        out.write(np_atom.to_packmol_inp(1))
        np_atom.save_structure(1)


def run_packmol():
    """Run and check that Packmol worked correctly"""
    try:
        p = subprocess.run(
            "packmol < {}".format(PACKMOL_INP),
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


def load_packmol_output():
    """Parse the output of Packmol"""
    return MDAnalysis.Universe(PACKMOL_OUT, transfer_to_memory=True)


def clean_tempfiles(structures):
    """Delete files generated by MDAPackmol"""
    structure_files = [
        PACKMOL_STRUCTURE_FILES.format(i) for i in range(len(structures))
    ]

    for f in [PACKMOL_INP, PACKMOL_OUT] + structure_files:
        try:
            os.remove(f)
        except FileNotFoundError:
            pass


def reassign_topology(structures, new):
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
        if any(hasattr(pms.structure, attr) for pms in structures):
            new.add_TopologyAttr(attr)

            if not all(hasattr(pms.structure, attr) for pms in structures):
                warnings.warn("added attribute which not all templates had")

    while index < len(new.atoms):
        # first atom we haven't dealt with yet
        start = new.atoms[index]
        # the resname was altered to give a hint to what template it was from
        template = structures[int(start.resname[1:])].structure
        # grab atomgroup which matches template
        to_change = new.atoms[index : index + len(template.atoms)]

        # Update residue names
        nres = len(template.residues)
        new.residues[
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
        new.add_TopologyAttr("bonds", values=bonds)
    if angles:
        angles = [tuple(val) for val in angles]
        new.add_TopologyAttr("angles", values=angles)
    if dihedrals:
        dihedrals = [tuple(val) for val in dihedrals]
        new.add_TopologyAttr("dihedrals", values=dihedrals)
    if impropers:
        impropers = [tuple(val) for val in impropers]
        new.add_TopologyAttr("impropers", values=impropers)

    return new


def packmol(
    smc: PackmolStructure, np_atom: PackmolStructure, atom_radius: float = None
):
    """Take molecules and settings and create a larger system

    Parameters
    ----------
    smc : PackmolStruture
      Supramolecular cage
    np_atom : PackmolStructure
      Nanoparticle atom
    atom_radius : float, optional
      Half of Packmol tolerance, ie., minimum distance between nanoparticle
      atoms, defaults to 1.0

    Returns
    -------
    new : MDAnalysis.Universe
      Universe object of the system created by Packmol
    """
    try:
        make_packmol_input(smc, np_atom, atom_radius=atom_radius)

        run_packmol()
    except PackmolError:
        # TODO: Deal with error
        new = None
    else:
        new = load_packmol_output()
        reassign_topology(structures, new)
    # finally:
    # clean_tempfiles(structures)

    return None
