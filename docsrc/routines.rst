Routines and wrappers
#########################

Routines
*************************

prjbuild.py
=========================

Constructs a template and assigns default values from a PNG image.

prjrun.py
=========================

Reads the project file assigns coefficients given that all the required
fields are satisfied then runs the specified 2D forward model. You can
suppress modeling and edit the stiffness and/or permittivity and conductivity
coefficients. Once they are provided in the project file, they won't be
computed or overwritten from the material values. If you would like to change
the material values and recompute the tensor coefficients, you need to delete
the existing tensor coefficients if included in the project file.

im2gif.py
=========================

Create a gif from the model outputs. Currently, this takes some time to run
which you can speed up by increasing the 'write' value in the project file.
This only takes 2D models, and there are bugs with matplotlib that cause
red/blue flashing.

arrayplot.py
=========================

Plot the seismograms or radargrams for the wide angle survey. You can
suppress plotting which will return a .csv file. An auto-controlled gain
function can be called for better visualization. The receiver locations are
given by a text file with the header X,Y,Z. These locations can be given in
meters relative to (0,0,0) or in indices. (0,0,0) is top left when viewing
the image.

rcxdisplay.py
=========================

.. note:: Originally "codisplay" in legacy code

Display the outputs of the common offset survey. This is also called to
display the common midpoint survey. Similar to ``arrayplot.py``, the gain
function can be called.

orientation_tensor.py
=========================

Compute the Euler angles and orientation tensor for a fabric defined by
it's trend and plunge angles. The orientation tensor isn't required by the
program but it provides useful quantitative information describing the
orientation of the fabric.

Wrappers
*************************

common_offset.sh
=========================

This is a wrapper that simulates a common offset survey. The receiver
.xyz file is input to give the points of the survey and the source is
offset from this location given the offsets for the x, y, and z directions.

common_midpoint.sh
=========================

This is similar to the common offset survey but it shifts the source and
reciever away from a common midpoint. The midpoint is specified by the
source location in the project file. By default the source will be to the
viewer's right of the midpoint but to flip the location of the source and
reciever, set the midpoint x-value to negative.

.. note::

    The aspect ratio for the common offset and common midpoint surveys
    determines the axis exaggeration. This will be updated in the future
    to be easier to adjust but to change this value edit the line
    ``ax.set_aspect(aspect=??)`` in ``arrayplot.py`` and ``codisplay.py``
    then run the plotting scripts individually not the wrapper scripts.
