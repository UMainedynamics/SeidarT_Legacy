common_offset
#########################

*A wrapper script to create a common offset survey. The survey is*
*along the x-direction, but can be extended to other directions.*


.. code-block:: bash

    common_offset -p [project_file].prj -o [offset in X Y Z] -r [receivers.xyz, typically uniform spacing]


**Inputs**

* ``-p FILE``, ``--project FILE``

    Project file path

* ``-o X Y Z``, ``--offset X Y Z``

    Source-reciever offset distance for all three directions (meters)

* ``-r FILE``, ``--receivers FILE``

    The coordinate locations of every survey point (meters)

* ``-s``, ``--seismic``

    (OPTIONAL) Specifier to run seismic common offset. Default is EM.

* ``-c``, ``--cores``

    For parallel computation, specify the number of cores.


**Outputs**

Three CSV files, containing Ex, Ey, and Ez values
