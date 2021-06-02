common_offset
#########################

*A wrapper script to create a common offset survey. The survey is*
*along the x-direction, but can be extended to other directions.*

**Usage**

.. code-block:: bash

    common_offset -p PROJECTFILE -o X Y Z -r RCXFILE -s [OPTIONAL] -c VALUE [optional]


**Inputs**

* ``-p PROJECTFILE``, ``--project PROJECTFILE``

    Project file path

* ``-o X Y Z``, ``--offset X Y Z``

    Source-receiver offset distance for all three directions (meters). Even in 2D calculations, you
    must enter three space-separated X Y Z values.

* ``-r RCXFILE``, ``--receivers RCXFILE``

    The coordinate locations of every survey point (meters). The file is comma-delimited
    in each row, with the first row required to by 'X,Y,Z'

* ``-s``, ``--seismic``

    (OPTIONAL) Specifier to run seismic common offset. Default is electromagnetic.

* ``-c``, ``--cores``

    For parallel computation, specify the number of cores.


**Outputs**

Three CSV files, containing Ex, Ey, and Ez (or Vx, Vy, Vz) values, with time series in columns of the source locations.




.. |top| raw:: html

   <a href="#top">Back to top ↑</a>

|top|
