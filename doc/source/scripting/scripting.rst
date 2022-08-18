.. _scripting:

Manual: Scripting
=================


Use Python scripts for full control of calculations.


..note:: The easy way to point the MSTM binary in script is the usage of the `os` module: 

    .. code-block:: python

        import os
        os.environ['MSTM_BIN'] = 'your path to mstm binary'
    
    The default path is '~/bin/mstm.x'


.. toctree::
   materials
   contribs
   spheres
   nearfield
   nonspherical
   fitting





