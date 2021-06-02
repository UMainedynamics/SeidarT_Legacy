vectoranim
#########################

This program builds a gif from the set of image output of the FDTD modeling. The images can be in csv or
unformatted Fortran binary, however, the program runs faster to use the latter. To use, ensure the project file is
in the same directory as the output files.

**Usage**

.. code-block:: bash

    vectoranim -p PROJECTFILE -n VALUE


**Inputs**

* ``-p PROJECTFILE``, ``--prjfile PROJECTFILE`` .prj file

    The full file path for the project file

* ``-n VALUE``, ``--num_steps VALUE``

    The time step interval between the images that
    are going to be used. Every time step is written to file which means that
    we can take any equally spaced images to create the gif with an
    appropriate resolution, time to compute, and file size. For example,
    n=20 means that every 20 images will be used thus significantly reducing
    how long it takes to compile the gif.

**Outputs**

    A GIF animation of the wave propagation.



`Back to top ↑ <#top>`_
