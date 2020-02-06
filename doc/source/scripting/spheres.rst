.. _spheres:


Spheres setup
-------------

MSTM code requires explicit setup of the positons and sizes of spheres. There are several classes designed in MSTM-studio to help this setup.

Example: particles aggregate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Script to construct spheres with random sizes placed on a regular grid:

.. literalinclude:: spheres_aggr.py
   :lines: 1-8

Sample output::

    Box size estimated as: 77.0 nm
    Desired number of particles: 9
    Number of particles in a box: 8
    Resulted number of particles: 8
    spheres are overlapping, regenerating...
    Box size estimated as: 77.0 nm
    Desired number of particles: 9
    Number of particles in a box: 8
    Resulted number of particles: 8
    [ 9.307773    8.61185299  9.92867988  8.84140858  9.87175352  8.71090184
      9.71505038 12.40459688]

    

Classes
^^^^^^^

.. autoclass:: mstm_studio.mstm_spectrum.Spheres
    :members:

.. autoclass:: mstm_studio.mstm_spectrum.SingleSphere
    :members:

.. autoclass:: mstm_studio.mstm_spectrum.ExplicitSpheres
    :members:

.. autoclass:: mstm_studio.mstm_spectrum.LogNormalSpheres
    :members:
    
    
MSTM run
--------

Example: core-shell particle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Class
^^^^^


.. [MSTM] U. Kreibig, M. Vollmer, "Optical Properties of Metal Clusters" (1995) 553

