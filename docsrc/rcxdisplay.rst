rcxdisplay
#########################

*Plot the amplitude timeseries based on a csv matrix*

**Usage**

.. code-block:: bash

    rcxdisplay -p [projectfile] -f [csv file with plottable output] -r [receiver location text file] -g [gain] -e [vertical exaggeration] -s [use if seismic]


**Inputs**

* ``-p FILE``, ``--prjfile FILE`` .prj file

    The file path for the project file, completely filled in for the model
    type used except for permittivity and stiffness coefficients, and dt

* ``-f FILE``, ``--file FILE``

    Path to the csv file with receiver timeseries data,
    commonly receiver_array.csv or Ex.co.csv

* ``-g VALUE``, ``--gain VALUE``

    Gain (smoothing length)

* ``-e value``, ``--exaggeration VALUE``

    Vertical exaggeration. Set the aspect ratio between the x and y axes for
    plotting. Default is 0.5

* ``-s VALUE``, ``--seismic VALUE``

    (OPTIONAL) Whether the model is seismic rather than electromagnetic

    * A nonzero value sets this flag to True (i.e. mark as a seismic model)
    * A value of 0 sets this (explicitly) to False (i.e. mark as an electromagnetic
      model, this is the default)


**Outputs**

Plot of amplitude value across area surveyed


`Back to top ↑ <#top>`_
