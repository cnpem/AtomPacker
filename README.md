# AtomPacker

A python package for packing nanoparticle atoms into a supramolecular cage.

## System requirements

- [gfortran](https://gcc.gnu.org/wiki/GFortran)

    ```bash
    apt install gfortran
    ```

- [packmol](http://m3g.iqm.unicamp.br/packmol/home.shtml)

    ```bash
    ./build-packmol.sh
    ```

## Python requirements

- [MDAnalysis](https://pypi.org/project/MDAnalysis)
- [pyKVFinder](https://pypi.org/project/pyKVFinder)

## Installation

```bash
git clone https://github.com/jvsguerra/AtomPacker.git
pip install -e AtomPacker
```

## Usage

There are two pipelines available to pack nanoparticle atoms inside a supramolecular cage:

### 1. Packing nanoparticle atoms freely with packmol and filter atoms inside the cavity;

```python
from AtomPacker import *

# Load supramolecular cage into PackmolStructure object
smc = PackmolStructure(
    os.path.join("data", "C1.pdb"),
    number=1,
    instructions=["center", "fixed 0. 0. 0. 0. 0. 0."],
)

# Load nanoparticle atoms into PackmolStructure object
# NOTE: You must set an appropriate number of atoms to fill the cavity and define a sphere that contains the whole supramolecular cage
np_atom = PackmolStructure(
    os.path.join("data", "Au.pdb"), number=40, instructions=["inside sphere 0. 0. 0. 6."]
)

# Create a CavityDetector with detection parameters appropriate for the supramolecular cage
cd = CavityDetector(step=0.25, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0, vdw=None)

# Create the AtomPacked object
# NOTE: For gold, atom radius is 1.36 Å
ap = AtomPacker(smc, np_atom, np_atom_radius=1.36, cavity_detector=cd,
basedir="pipeline-1")

# Run Packing algorithm with 10 replicates
ap.packing(replicates=10)

# Show number of packed atoms
print(ap.summary)
```

### [WIP] 2. Packing nanoparticle atoms with packmol inside a boundary and filter atoms inside the cavity;

```python
from AtomPacker import *

# Load supramolecular cage into PackmolStructure object
smc = PackmolStructure(
    os.path.join("data", "C1.pdb"),
    number=1,
    instructions=["center", "fixed 0. 0. 0. 0. 0. 0."],
)

# Load nanoparticle atoms into PackmolStructure object
# NOTE: You must set an appropriate number of atoms to fill the cavity and define a sphere that contains the whole supramolecular cage
np_atom = PackmolStructure(
    os.path.join("data", "Au.pdb"), number=40, instructions=["inside sphere 0. 0. 0. 6."]
)

# Create a CavityDetector with detection parameters appropriate for the supramolecular cage
cd = CavityDetector(step=0.25, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0, vdw=None)

# Create the AtomPacked object
# NOTE: For gold, atom radius is 1.36 Å
ap = AtomPacker(smc, np_atom, np_atom_radius=1.36, cavity_detector=cd,
basedir="pipeline-2")

# Add a boundary to pack nanoparticle atoms
ap.add_boundary()

# Run Packing algorithm with 10 replicates
ap.packing(replicates=10)

# Show number of packed atoms
print(ap.summary)
```

## Citing

If you find `AtomPacker` useful for you, please cite the following sources:

- L Martinez, R Andrade, E G Birgin, J M Martinez, "Packmol: A package for building initial configurations for molecular dynamics simulations". Journal of Computational Chemistry, 30, 2157-2164, 2009.

- R J Gowers, M Linke, J Barnoud, T J E Reddy, M N Melo, S L Seyler, D L Dotson, J Domanski, S Buchoux, I M Kenney, and O Beckstein. "MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations." In S. Benthall and S. Rostrup, editors, Proceedings of the 15th Python in Science Conference, pages 102-109, Austin, TX, 2016.

- J V S Guerra, H V Ribeiro-Filho, G E Jara, L O Bortot, J G C Pereira and P S Lopes-de-Oliveira. "pyKVFinder: an efficient and integrable Python package for biomolecular cavity detection and characterization in data science." BMC bioinformatics, 22(1), 607, 2021.

## License

The software is licensed under the terms of the GNU General Public License version 3 (GPL3) and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
