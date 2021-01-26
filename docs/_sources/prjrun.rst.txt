prjrun
##########################

*Performs single-shot wave propagation with receiver locations specified*
*later in postprocessing.*

**Usage**

.. code-block:: bash

    prjrun -p [project_file].prj -M [n e s]

**Inputs**

* ``-p FILE``, ``--prjfile FILE`` .prj file

    The file path for the project file, completely filled in for the model
    type used except for permittivity and stiffness coefficients, and dt

* ``-i IMAGE_FILE``, ``--image_file IMAGE_FILE``

    The full file path for the .png image

* ``-o PROJECT_FILE``, ``--project_file PROJECT_FILE``

    Name of the output file path with extension, .prj, and excluding
    the full path directory

* ``-M [n e s]`` Which model to run

    * ``n`` calculates only timesteps and material tensors
    * ``e`` electromagnetic propagation
    * ``s`` seismic wave propagation

**Outputs**

.dat files equal in number to the number of time steps specified in the
.prj file for both the x and z direction


`Back to top ↑ <#top>`_
