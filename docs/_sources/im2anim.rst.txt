im2anim
#########################

Create animation of the propagation of waves across a 2-dimensional
study area

**Usage**

.. code-block:: bash

    im2anim -p PRJFILE -c [Ex Ez Vx Vz] -n [steps in time series] -d [delay between frames]


**Inputs**

* ``-p FILE``, ``--prjfile FILE`` .prj file

    The full file path for the project file, completely filled in for
    the model type used

* ``-c [Ex Ez Vx Vz]``, ``--channel [Ex Ez Vx Vz]``

    Whether to plot the seismic or radar model in the X or Z direction

* ``-n VALUE``, ``--num_steps VALUE``

    Number of steps in time series used, 50 works well

* ``-d VALUE``, ``--frames_per_second VALUE``

    Number of frames per second, 10 works well

* ``-a VALUE``, ``--alpha VALUE``

    Change the transparency of the model plotted in the background; default = 0.3.
    Zeros is full transparency, 1 is opaque.

**Outputs**

    A GIF animation of the wave propagation.



`Back to top ↑ <#top>`_
