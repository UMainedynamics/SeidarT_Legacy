wiggleplot
########################

*Plot single, and select handful of amplitudes at specified X*
*locations and all Z locations*

**Usage**

.. code-block:: bash

    wiggleplot -p PROJECTFILE -r DATAFILE -g VALUE -x VALUE -d VALUE -n VALUE \
      -c [Ex Ey Ez Vx Vy Vz]


**Input**

* ``-p PROJECTFILE``, ``--prjfile PROJECTFILE`` .prj file

    The file path for the project file.

* ``-r receiver_array.csv``

    The file path for the receiver array data, typically receiver_array.csv

* ``-g VALUE``, ``--gain VALUE``

    The linear horizontal exaggeration of the
    amplitude values from the receiver array file.

* ``-x VALUE``, ``--receiver_spacing``

    The horizontal distance between receivers, in meters, to lable the x-axis properly.

* ``-d``, ``--columns``

    The frequency at which columns are pulled for
    plotting from the csv file

* ``-n VALUE``, ``--single_plot_dist``

    The receiver number (column) indicating amplitude values will
    be plotted in a single-wiggle plot

* ``-c [Ex Ey Ez Vx Vy Vz]``, ``--channel [Ex Ey Ez Vx Vy Vz]``

    The channel to query
    (``Vx`` for the seismic x direction, for example)


**Output**

* Plot of the amplitude value of each pixel in a column
* Plot of the amplitude value of every pixel in every nth column
  aligned across a graph


