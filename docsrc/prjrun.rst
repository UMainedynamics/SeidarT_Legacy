prjrun
##########################

*Performs single-shot wave propagation with receiver locations specified*
*later in postprocessing.*

**Usage**

.. code-block:: bash

    prjrun -p PRJFILE -m [n e s] -a [0 1]

**Inputs**

* ``-p PRJFILE``, ``--prjfile PRJFILE`` .prj file

    The file path for the project file, completely filled in for the model
    type used except for permittivity and stiffness coefficients, and dt

* ``-m [n e s]`` Which model to run

    * ``n`` calculates only timesteps and material tensors, necessary before running sourcefunction
    * ``e`` electromagnetic propagation
    * ``s`` seismic wave propagation

* ``-a [0 1] (optional)

    Append/recompute the coefficients to the permittivity and
    stiffness matrices; 1 = yes, 0 = no; default = 1. Do not
    recompute if you have made manual changes to the matrices in the prj file.

**Outputs**

.dat files equal in number to the number of time steps specified in the .prj file.

`Back to top ↑ <#top>`_
