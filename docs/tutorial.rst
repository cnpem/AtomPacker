========
Tutorial
========

This section provides a tutorial on how to use the AtomPacker package to pack nanoparticle atoms inside a target supramolecular cage. For detailed reference documentation of the functions and classes contained in the package, see the `API reference <../api.html>`_.

For this tutorial, we will use the `ZOCXOH <https://www.ccdc.cam.ac.uk/structures/Search?Ccdcid=ZOCXOH>`_ cage structure from the `ZOCXOH.pdb <https://github.com/cnpem/AtomPacker/tests/data/ZOCXOH.pdb>`_ file.

First of all, import **AtomPacker** package on Python and create a AtomPacker.Cage object.

.. code-block:: python

    >>> import AtomPacker
    >>> cage = AtomPacker.Cage()

The Cage object will be used to load the cage structure from the PDB file (eg., `ZOCXOH.pdb <https://github.com/cnpem/AtomPacker/tests/data/ZOCXOH.pdb>`_).

.. code-block:: python

    >>> cage.load('ZOCXOH.pdb')
    
If you want to preview the cage structure, you can use the `Cage.preview <../api.html#AtomPacker.Cage.preview>`_ method.

.. code-block:: python

    >>> cage.preview()

.. image:: _static/cage-preview-1.png
    :width: 600
    :align: center

The next step is to detect the cavity inside the cage structure. The `Cage.detect_cavity <../api.html#AtomPacker.Cage.detect_cavity>`_ method will detect the cavity using the `pyKVFinder <https://lbc-lnbio.github.io/pyKVFinder/_api_reference/index.html>`_ detection parameters (``step``, ``probe_in``, ``probe_out``, ``removal_distance``, ``volume_cutoff``, ``surface``). We adjust the parameters to detect the cavity inside the cage structure.

The detected cavity will be stored in the ``cavity`` attribute of the Cage object. 

.. code-block:: python
    
    >>> cage.detect_cavity(step=0.6, probe_in=1.4, probe_out=10.0, removal_distance=1.0, volume_cutoff=5.0)

If you want to preview the cavity structure for quality control in your cavity detection, you can use the `Cavity.preview <../api.html#AtomPacker.Cavity.preview>`_ method.

.. code-block:: python

    >>> cage.cavity.preview()

.. image:: _static/cavity-preview.png
    :width: 600
    :align: center

You can access the cavity coordinates using the ``coordinates`` attribute of the Cavity object.

.. code-block:: python

    >>> print(cage.cavity.coordinates)
    [[-11.998       28.644       13.149     ]
     [-11.998       28.644       13.749     ]
     [-11.998       28.644       14.349     ]
     ...
     [ -1.1980003   30.443998    14.349     ]
     [ -1.1980003   30.443998    14.948999  ]
     [ -0.59800035  29.844       13.749     ]]


You can access the cavity grid using the ``grid`` attribute of the Cavity object.

.. code-block:: python

    >>> print(cage.cavity.grid)
    [[[-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      ...
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]]

     [[-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      ...
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]]

     [[-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      ...
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]]

     ...

     [[-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      ...
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]]

     [[-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      ...
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]]

     [[-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      ...
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]
      [-1 -1 -1 ... -1 -1 -1]]]

You can access the cavity volume using the ``volume`` attribute of the Cavity object.

.. code-block:: python

    >>> print(cage.cavity.volume)
    531.58

You can also save the cavity structure using the `Cavity.save <../api.html#AtomPacker.Cavity.save>`_ method.

.. code-block:: python

    >>> cage.cavity.save("cavity.pdb")

The next step is to pack the nanoparticle atoms inside the cavity. The `Cage.pack <../api.html#AtomPacker.Cage.pack>`_ method will pack the atoms inside the cavity using the specified parameters (``atom_type``, ``lattice_type``, ``a``, ``b``, ``c``). We pack a gold (``Au``) nanoparticle inside the cavity using the face-centered cubic (``fcc``) lattice.

The packed cluster will be stored in the ``cluster`` attribute of the Cage object.

.. code-block:: python

    >>> cage.pack(atom_type="Au", lattice_type="fcc", a=None, b=None, c=None)

If you want to preview the cluster structure for quality control, you can use the `Cluster.preview <../api.html#AtomPacker.Cluster.preview>`_ method.

.. code-block:: python

    >>> cage.cluster.preview()

.. image:: _static/cluster-preview.png
    :width: 600
    :align: center

Also, you can preview the cage, cavity and cluster structure using the `Cage.preview <../api.html#AtomPacker.Cage.preview>`_ method.

.. code-block:: python

    >>> cage.preview(show_cavity=True, show_cluster=True)

.. image:: _static/cage-preview-2.png
    :width: 600
    :align: center

You can access the cluster coordinates using the ``coordinates`` attribute of the Cluster object.

.. code-block:: python

    >>> print(cage.cluster.coordinates)
    [[ -8.56353   25.718132  13.873403]
     [ -8.56353   27.758131  11.833403]
     [-10.60353   27.758131  13.873403]
     [-10.60353   29.798132  15.913403]
     [ -8.56353   27.758131  15.913403]
     [ -8.56353   29.798132  13.873403]
     [ -8.56353   29.798132  17.953403]
     [ -8.56353   31.838133  11.833403]
     [ -8.56353   31.838133  15.913403]
     [ -8.56353   33.87813   13.873403]
     [ -6.52353   25.718132  15.913403]
     [ -4.48353   25.718132  13.873403]
     [ -6.52353   29.798132  11.833403]
     [ -4.48353   27.758131  11.833403]
     [ -4.48353   29.798132   9.793403]
     [ -6.52353   27.758131  13.873403]
     [ -6.52353   29.798132  15.913403]
     [ -4.48353   27.758131  15.913403]
     [ -4.48353   29.798132  13.873403]
     [ -4.48353   29.798132  17.953403]
     [ -6.52353   31.838133   9.793403]
     [ -6.52353   33.87813   11.833403]
     [ -4.48353   31.838133  11.833403]
     [ -6.52353   31.838133  13.873403]
     [ -6.52353   33.87813   15.913403]
     [ -4.48353   31.838133  15.913403]
     [ -4.48353   33.87813   13.873403]
     [ -6.52353   31.838133  17.953403]
     [ -2.44353   29.798132  11.833403]
     [ -2.44353   27.758131  13.873403]
     [ -2.44353   29.798132  15.913403]
     [ -2.44353   31.838133  13.873403]]

You can also save the cavity structure using the `Cluster.save <../api.html#AtomPacker.Cluster.save>`_ method.

.. code-block:: python

    >>> cage.cluster.save("cluster.pdb")

Finally, you can access the summary of the cluster using the ``summary`` attribute of the Cluster object.

.. code-block:: python

    >>> print(cage.cluster.summary)
                                 Au32
    Atom Type                      Au
    Atom Radius              1.442498
    Cavity Volume (Å³)         531.58
    Diameter (maximum)       9.123157
    Diameter (shape)             8.16
    Diameter (volume)        10.12412
    Lattice Constants            4.08
    Lattice Type                  fcc
    Maximum Number of Atoms        43
    Number of Atoms                32
