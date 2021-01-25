wiggleplot
########################

*Plot single, and select handful of amplitudes at specified X*
*locations and all Z locations*

.. code-block:: bash

    wiggleplot -r receiver_array.csv -p [project_file].prj -g [gain] -d [amplitude column plotting frequency] -n [singular amplitude column to plot] -c [Ex Ez Vx Vz]


**Input**

* ``-r receiver_array.csv``

    Amplitude values for all pixels

* ``-p FILE``: .prj file

    Completely filled in for the model being run

* ``-g VALUE``: Gain, horizontal exaggeration

    1000-5000 usually works, but smaller values can too if the amplitude
    values are large

* ``-d VALUE``: Every nth column to pull to plot the amplitude values of

    5 or 10 usually works

* ``-n VALUE``

    Which singular column to plot the amplitude values of

* ``-c [one of Ex Ez Vx Vz]``

    Whether to plot the seismic or radar model in the X or Z direction
    (``Vx`` for the seismic x direction, for example)


**Output**

* Plot of the amplitude value of each pixel in a column
* Plot of the amplitude value of every pixel in every nth column
  aligned across a graph
