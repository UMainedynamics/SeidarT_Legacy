im2anim
#########################

Create animation of the propagation of waves across a 2-dimensional
study area

**Usage**

.. code-block:: bash

    im2anim -p [project_file].prj -c [Ex Ez Vx Vz] -n [steps in time series] -f [frames per second]


**Inputs**

* ``-p FILE``, ``--prjfile FILE`` .prj file

    The full file path for the project file, completely filled in for
    the model type used

* ``-c [Ex Ez Vx Vz]``, ``--channel [Ex Ez Vx Vz]``

    Whether to plot the seismic or radar model in the X or Z direction

* ``-n VALUE``, ``--num_steps VALUE``

    Number of steps in time series used, 50 works well

* ``-f VALUE``, ``--frames_per_second VALUE``

    Number of frames per second, 10 works well

* ``-t VALUE``, ``--threshold VALUE``

    (OPTIONAL) Set array values to zero when they go below a specific
    threshold

**Outputs**

``-o OUTPUT``, ``--output OUTPUT``

    Specify the output format. 0 - GIF (default), 1 - MP4

    Outputs an animation containing an image of the wave propagation
    at each frame



`Back to top ↑ <#top>`_
