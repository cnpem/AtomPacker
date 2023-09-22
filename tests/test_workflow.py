#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `AtomPacker` workflow."""
import os
import pytest
from AtomPacker import *
import MDAnalysis

# HERE = os.path.abspath(os.path.dirname(__file__))
# WATER_PDB = os.path.join(HERE, 'water.pdb')
# UREA_PDB = os.path.join(HERE, 'urea.pdb')

SMC_PDB = os.path.join("data", "C1.pdb")
AU_PDB = os.path.join("data", "Au.pdb")

# Load individual molecule files
smc = PackmolStructure(
    MDAnalysis.Universe(SMC_PDB),
    number=1,
    instructions=["center", "fixed 0. 0. 0. 0. 0. 0."],
)
np_atom = PackmolStructure(
    MDAnalysis.Universe(AU_PDB), number=40, instructions=["inside sphere 0. 0. 0. 7."]
)

# Cavity Detector
ncav, cavities, boundaries = CavityDetector(MDAnalysis.Universe("0.pdb"), verbose=True).detect(
    step=0.25, probe_out=10.0, removal_distance=1.0
)

# Run Packmol
packmol(smc, np_atom, atom_radius=1.36)
