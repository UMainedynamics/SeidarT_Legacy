wiggleplot
########################

*Plot single, and select handful of amplitudes at specified X*
*locations and all Z locations*

**Usage**

.. code-block:: bash

    wiggleplot -r receiver_array.csv -p [project_file].prj -g [gain] -d [amplitude column plotting frequency] -n [singular amplitude column to plot] -c [Ex Ez Vx Vz]


**Input**

* ``-r receiver_array.csv``

    The file path for the receiver array data

* ``-p FILE``, ``--prjfile FILE`` .prj file

    The file path for the project file, completely filled in for the model
    type used except for permittivity and stiffness coefficients, and dt

* ``-g VALUE``, ``--gain VALUE``

    The linear horizontal exaggeration of the
    amplitude values from the receiver array file

* ``-x VALUE``, ``--receiver_spacing``

    The horizontal distance between receivers, in meters

* ``-d``, ``--columns``

    The frequency at which columns are pulled for
    plotting from the csv file

* ``-n VALUE``, ``--single_plot_dist``

    The distance along the image surface where the amplitude values will
    be plotted

* ``-c [Ex Ez Vx Vz]``, ``--channel [Ex Ez Vx Vz]``

    The channel to query
    (``Vx`` for the seismic x direction, for example)


**Output**

* Plot of the amplitude value of each pixel in a column
* Plot of the amplitude value of every pixel in every nth column
  aligned across a graph


.. |top| raw:: html

   <a href="#top">Back to top ↑</a>

|top|
