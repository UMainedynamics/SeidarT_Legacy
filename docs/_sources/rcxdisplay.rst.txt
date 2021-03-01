rcxdisplay
#########################

*Plot the amplitude timeseries based on a csv matrix*

**Usage**

.. code-block:: bash

    rcxdisplay -p PROJECTFILE -f DATAFILE -g VALUE -e VALUE -s [OPTIONAL]


**Inputs**

* ``-p PROJECTFILE``, ``--prjfile PROJECTFILE`` .prj file

    The file path for the project file.

* ``-f DATAFILE``, ``--file DATAFILE``

    Path to the csv file with receiver timeseries data,
    commonly receiver_array.csv or [Vx Vy Vz Ex Ey Ez].co.csv. See arraybuild for how to generate this file for single shots.

* ``-g VALUE``, ``--gain VALUE``

    Gain (smoothing length). Possible values are 1 through the number of timesteps.

* ``-e VALUE``, ``--exaggeration VALUE``

    Vertical exaggeration. Set the aspect ratio between the x and z axes for
    plotting. Default is 0.5

* ``-s VALUE``, ``--seismic VALUE``

    (OPTIONAL) Whether the model is seismic rather than electromagnetic

    * A nonzero value sets this flag to True (i.e. mark as a seismic model)
    * A value of 0 sets this (explicitly) to False (i.e. mark as an electromagnetic
      model, this is the default)


**Outputs**

Plot of amplitude across area surveyed.


`Back to top ↑ <#top>`_
