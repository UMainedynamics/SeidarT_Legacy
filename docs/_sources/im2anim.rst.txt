im2anim
#########################

Create animation of the propagation of waves across a 2-dimensional
study area

**Usage**

.. code-block:: bash

    im2anim -p PROJECTFILE -c [Ex Ez Vx Vz] -n [steps in time series] -d [delay between frames]


**Inputs**

* ``-p PROJECTFILE``, ``--prjfile PROJECTFILE`` .prj file

    The full file path for the project file

* ``-c [Ex Ez Vx Vz]``, ``--channel [Ex Ez Vx Vz]``

    Whether to plot the seismic or radar model in the X or Z direction

* ``-n VALUE``, ``--num_steps VALUE``

    Interval between time series steps used to create the animation

* ``-d VALUE``, ``--delay VALUE`` OPTIONAL

    Delay between frames, default=1

* ``-a VALUE``, ``--alpha VALUE`` OPTIONAL

    Change the transparency of the model plotted in the background; default = 0.3.
    Zero is full transparency, 1 is opaque.

**Outputs**

    A GIF animation of the wave propagation.



`Back to top ↑ <#top>`_
