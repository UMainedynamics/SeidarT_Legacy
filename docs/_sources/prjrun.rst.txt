prjrun
##########################

*Performs single-shot wave propagation with receiver locations specified*
*later in postprocessing.*

**Usage**

.. code-block:: bash

    prjrun -p PRJFILE -m [n e s]

**Inputs**

* ``-p PRJFILE``, ``--prjfile PRJFILE`` .prj file

    The file path for the project file, completely filled in for the model
    type used except for permittivity and stiffness coefficients, and dt

* ``-M [n e s]`` Which model to run

    * ``n`` calculates only timesteps and material tensors, necessary before running sourcefunction
    * ``e`` electromagnetic propagation
    * ``s`` seismic wave propagation

**Outputs**

.dat files equal in number to the number of time steps specified in the .prj file.

`Back to top ↑ <#top>`_
