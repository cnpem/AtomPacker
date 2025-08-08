===========
Quick start
===========

AtomPacker is a Python package for packing atoms of nanoparticles into supramolecular cages.

It provides tools to load cage structures, detect cavities and their openings, fill them with atoms according to specific lattice parameters, preview the results, and export the generated structures.

.. code-block:: python
    
    >>> import AtomPacker
    >>> # Load the cage structure
    >>> cage = AtomPacker.Cage()
    >>> cage.load("ZOCXOH.pdb")
    >>> # Detect cavity and openings
    >>> cage.detect_cavity(step=0.6, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0)
    >>> cage.cavity.detect_openings()
    >>> # Pack gold atoms into the cavity (FCC lattice)
    >>> cage.pack(atom_type="Au", lattice_type="fcc")
    >>> # Preview all structures together
    >>> cage.preview(show_cavity=True, show_openings=True, show_cluster=True)
    >>> # Show cluster summary
    >>> print(cage.cluster.summary)
    >>> # Save results
    >>> cage.cavity.save("cavity.pdb")
    >>> cage.cavity.openings.save("openings.pdb")
    >>> cage.cluster.save("cluster.pdb")
