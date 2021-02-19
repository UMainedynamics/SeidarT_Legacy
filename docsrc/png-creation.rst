Creating a png file
#########################

.. |gimp| raw:: html

   <a href="https://www.gimp.org/downloads/install_help.html" target="_blank">GIMP</a>

.. |inkscape| raw:: html

   <a href="http://wiki.inkscape.org/wiki/index.php/Installing_Inkscape" target="_blank">Inkscape</a>


Geometries for the model domain within SeidarT are initiated with a
PNG image. The program identifies unique RGB values, setting material
properties for each. For example, if you wanted to define a geometry
with ice overlying bedrock, you would create a .png image that is one
color for the ice and another for the rock below. Everyone has their
preferences to generate images but |gimp| or |inkscape| provide free
and open software that are more than sufficient.

.. note::

    When creating a PNG, anti-aliasing must be turned off to avoid
    color boundary gradients.

Building images in Inkscape has some advantages other than being free.
Saving a .svg to pdf allows the user to change the number of pixels
and the spatial resolution of the image quite easily. With
ghostscript, the command ::

    gs -q -dBATCH -dNOPAUSE -sDEVICE=png16m -sOutputFile=<file> -r96 <input_file>

will generate a PNG file from a PDF. The resolution :code:`-r` can be
varied to change the pixels. In Inkscape, the image pixels can be set
in Document Properties. When saving the SVG as PDF, you will be
prompted with options, and the value for Resolution for rasterization
(dpi): will determine - in order to get the same pixel setting that
you set in Inkscape - the value for the :code:`-r` (resolution) option
above. If you want to double the resolution, just double this number
(i.e. :code:`-r192`).
