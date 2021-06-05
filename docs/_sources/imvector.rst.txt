imvector
#########################

Plot a snapshot of the vector wavefield in 2D

**Usage**

.. code-block:: bash

    imvector -p PROJECTFILE -v VELOCITYFILE


**Inputs**

* ``-p PROJECTFILE``, ``--prjfile PROJECTFILE`` .prj file

    The full file path for the project file

* ``-v VELOCITYFILE``, ``--velocity VELOCITYFILE``

    The .dat file that corresponds to the velocity in either the
    x-direction or z-direction (e.g. Vx000400.dat). The corresponding
    orthogonal velocity file will be loaded as well.

**Outputs**

    A png image of seismic or radar vector field of the selected timestep overlain on the model geometry



