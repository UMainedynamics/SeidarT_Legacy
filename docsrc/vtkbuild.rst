vtkbuild
####################

*Creates a series of 3-dimensional vti files that can be imported into*
*Paraview to generate an animation. (Directly building the animation is*
*planned for a future release.)*

**Usage**

.. code-block:: bash

    vtkbuild [-h] -c CHANNEL [-f FRAMES_PER_SECOND] -n NUM_STEPS project_file


**Inputs**

* ``-p FILE``, ``--prjfile FILE`` .prj file

    The full file path for the project file, completely filled in for
    the model type used

* ``-c [Ex Ez Vx Vz]``, ``--channel [Ex Ez Vx Vz]``

    Whether to plot the seismic or radar model in the X or Z direction

* ``-n VALUE``, ``--num_steps VALUE``

    The time step interval between the images that are
    going to be used. Every time step is written to file which means that
    we can take any equally spaced images to create the gif with an
    appropriate resolution, time to compute, and file size. For example,
    n=20 means that every 20 images will be used thus significantly reducing
    how long it takes to compile.


`Back to top ↑ <#top>`_
