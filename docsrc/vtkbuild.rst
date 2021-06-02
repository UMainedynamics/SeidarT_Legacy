vtkbuild
####################

*Creates a series of 3-dimensional vti files that can be imported into*
*Paraview to generate an animation. (Directly building the animation is*
*planned for a future release.)*

**Usage**

.. code-block:: bash

    vtkbuild [-h] -p PROJECTFILE -c [Ex Ey Ez Vx Vy Vz] -n NUM_STEPS

**Inputs**

* ``-p PROJECTFILE``, ``--prjfile PROJECTFILE`` .prj file

    The full file path for the project file.

* ``-c [Ex Ey Ez Vx Vy Vz]``, ``--channel [Ex Ey Ez Vx Vy Vz]``

    Whether to plot the seismic or radar model in the X, Y, or Z direction

* ``-n VALUE``, ``--num_steps VALUE``

    The time step interval between the images that are
    going to be used. n=20 means that 1 out of every 20 images will be used,
    thus significantly reducing how long it takes to compile.

