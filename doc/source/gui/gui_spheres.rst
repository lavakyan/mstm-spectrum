.. _gui_spheres:

Spheres
-------

.. image:: interface_mstm.png

Plus button -- add new sphere.

.. image:: gui_spheres_add.png

Important to specify material label for the sphere.


Button with pencil or doulbe-click on a table row -- edit selected sphere. 

Note: double click on a row in material table allow to change the viewed material color.

Circle-arrows -- refresh 3D view.

3D view can be **rotated** with pressed left mouse button and **zoomed** in or out with mouse wheel.


Environment material by default is `m0`. This can be changed using menu:

.. image::  gui_spheres_matrix.png

Cross button deletes the selected sphere.


MSTM run
^^^^^^^^

.. image:: gui_mstm.png

"min" -- minimal wavelength (in nm),

"max" -- maximal wavelength (in nm),

"count" -- number of wavelength points. By default the spacing is 10 nm.

"Calculate" button runs MSTM binary in temporary directory (OS-dependent) and reads the results.

"scale" -- total outer multiplier.

save button -- save extinction to column file.

plot button -- plot without re-calculation (i.e. with new `scale`).

.. image:: gui_mstm_plot.png

The plot controls are rendered by Matplotlib, and can depend on the library version. Generally, it is possible to zoom region of interest and save graphic as a raster or vector image.








