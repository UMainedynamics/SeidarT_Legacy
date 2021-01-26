common_midpoint
########################

*creates receiver setup by defining final and middle locations,*
*and receiver spacing across surface*

**Usage**

.. code-block:: bash

    common_midpoint -p [project_file].prj -t [X] -o [X] -d [distance between each receiver] [-s] [-p]


**Inputs**

* ``-p FILE``, ``--prjfile FILE`` .prj file

    The file path for the project file, completely filled in for the model
    type used except for permittivity and stiffness coefficients, and dt

* ``-t VALUE``, ``--total VALUE``

    The terminal distance between the source and reciever

* ``-o VALUE``, ``--offset VALUE``

    The initial source and reciever offset from the midpoint
    given in (+/- meters). A negative value means that the
    source is on the lookers left of the midpoint. The total
    source and reciever distance is 2*offset.

* ``-d VALUE``, ``--delta VALUE``

    Source and reciever step length (meters); total distance
    between the source and reciever is 2*delta*i + 2*offset.

* ``-s``, ``--seismic``

    (OPTIONAL) Specifier to run seismic common offset

* ``-p``, ``--plot``

    (OPTIONAL) Show plot. Default is none


**Outputs**



`Back to top ↑ <#top>`_
