# AtomPacker

[![PyPI - Version](https://img.shields.io/pypi/v/AtomPacker)](https://pypi.org/project/AtomPacker/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/AtomPacker)](https://pypi.org/project/AtomPacker/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/AtomPacker)](https://pypi.org/project/AtomPacker/)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/cnpem/AtomPacker/testing.yml?label=testing)
![GitHub](https://img.shields.io/github/license/cnpem/AtomPacker)

A Python package for packing nanoclusters into supramolecular cages.

See also:

- [Documentation](https://cnpem.github.io/AtomPacker/)
- [GitHub repository](https://github.com/cnpem/AtomPacker/)
- [Issues](https://github.com/cnpem/AtomPacker/issues)

## Requirements

- [ASE](https://pypi.org/project/ase)
- [biopython](https://pypi.org/project/biopython)
- [MDAnalysis](https://pypi.org/project/MDAnalysis)
- [numpy](https://pypi.org/project/numpy)
- [plotly](https://pypi.org/project/plotly)
- [pyKVFinder](https://pypi.org/project/pyKVFinder)

## Installation

To install the latest release on [PyPI](https://pypi.org/project/AtomPacker/), run:

```bash
pip install AtomPacker
```

Or, to install the development version, run:

```bash
pip install git+https://github.com/cnpem/AtomPacker.git
```

## Usage

Packing nanoparticle atoms, based on ASE nanocluster, and filter atoms inside a target cavity.

```python
from AtomPacker import Cage
# 1: Load structure from file
cage = Cage()
cage.load("tests/data/ZOCXOH.pdb")
# Uncomment to preview the cage structure.
# cage.preview()
# 2: Detect cavity
cage.detect_cavity(step=0.25, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0)
# Uncomment to preview the cavity structure for detection quality control.
# cage.cavity.preview()
# Show volume
print(f"Cavity volume: {cage.cavity.volume} A^3")
# Uncomment to save the cavity structure.
# cage.cavity.save("tests/cavity.pdb")
# 3: Pack nanocluster into the cavity
cage.pack(atom_type="Au", lattice_type="fcc", a=None, b=None, c=None)
# Uncomment to preview the cluster structure for quality control.
# cage.cavity.preview()
# Uncomment to save the cluster structure.
# cage.cluster.save("tests/cluster.pdb")
# Uncomment to preview the cage, cavity and cluster structure.
# cage.preview(show_cavity=True, show_cluster=True)
# Show summary
print(cage.cluster.summary)
```

## Architecture

The package is organized as follows:

```mermaid
classDiagram
direction LR
Cage "1" o-- "1" Cavity : has
Cage "1" o-- "1..*" Cluster : fits in
Cavity "1" <|-- "1..*" Cluster : needs
namespace AtomPacker {
class Cage {
+ numpy.ndarray atomic
+ Cavity cavity
+ numpy.ndarray centroid
+ Cluster cluster
+ numpy.ndarray coordinates
+ MDAnalysis.Universe universe
+ detect_cavity(float step, float probe_in, float probe_out, float removal_distance, float volume_cutoff, str surface, int nthreads, bool verbose, Dict~str,Any~ **kwargs) Cavity
+ load(filename) MDAnalysis.Universe
+ pack(str lattice_type, str atom_type, float atom_radius, float a, float b, float c) ase.cluster.Cluster
+ preview(bool show_cavity, bool show_cluster, str renderer, Dict~str,Any~ **kwargs) void
# _build_cluster(str atom_type, str lattice_type, Tuple~float~ lattice_constants, numpy.ndarray center) ase.cluster.Cluster
# _filter_cluster(ase.cluster.Cluster cluster) ase.cluster.Cluster
# _get_cluster_layers(str atom_type, float factor) numpy.ndarray            
}
class Cavity {
+ numpy.ndarray coordinates
+ numpy.ndarray grid
+ Universe universe
+ numpy.ndarray volume
# float step
# float probe_in
# float probe_out
# float removal_distance
# numpy.ndarray vertices
# float volume_cutoff
# str surface
+ preview(str renderer, float opacity, Dict~str,Any~ **kwargs) void
+ select_cavity(List~int~ indexes) void
+ save(str filename) void
# _get_universe() Universe
}
class Cluster {
+ str atom_type
+ numpy.ndarray coordinates
+ str lattice_type
+ Tuple~float~ lattice_constants
+ int number_of_atoms
+ int maximum_number_of_atoms
+ pandas.DataFrame summary
+ Universe universe
+ numpy.ndarray volume
# Cavity cavity
# ase.cluster.Cluster cluster
+ diameter(str method) float
+ preview(str renderer, float opacity, Dict~str,Any~ **kwargs) void
+ save(str filename) void
# _get_distances() numpy.ndarray
# _get_lattice_constants() Tuple~float~
# _get_radii() float
# _get_universe() Universe  
}
}
```

## Citing

If you find `AtomPacker` useful for you, please cite the following references:

- Guerra, J. V. S., Ribeiro-Filho, H. V., Jara, G. E., Bortot, L. O., Pereira, J. G. C., & Lopes-de-Oliveira, P. S. (2021). pyKVFinder: an efficient and integrable Python package for biomolecular cavity detection and characterization in data science. BMC bioinformatics, 22(1), 607. https://doi.org/10.1186/s12859-021-04519-4.

- Manuscript in preparation.

## License

The software is licensed under the terms of the GNU General Public License version 3 (GPL3) and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
