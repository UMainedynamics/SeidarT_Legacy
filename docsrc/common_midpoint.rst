common_midpoint
########################

*creates receiver setup by defining final and middle locations,*
*and receiver spacing across surface*

*UNDER CONSTRUCTION*

**Usage**

.. code-block:: bash

    common_midpoint -p PRJFILE -t [X] -o [X] -d [distance between each receiver] [-s] [-p]


**Inputs**

* ``-p PRJFILE``, ``--prjfile PRJFILE``

    The project file path

* ``-t VALUE``, ``--total VALUE``

    The terminal distance between the source and receiver

* ``-o VALUE``, ``--offset VALUE``

    The initial source and receiver offset from the midpoint
    given in (+/- meters). A negative value means that the
    source is left of the midpoint. The total
    source and reciever distance is 2*offset.

* ``-d VALUE``, ``--delta VALUE``

    Source and receiver step length (meters); total distance
    between the source and receiver is 2*delta*i + 2*offset.

* ``-s``, ``--seismic``

    (OPTIONAL) Specifier to run seismic common offset

* ``-p``, ``--plot``

    (OPTIONAL) Show plot. Default is none


**Outputs**

`Back to top ↑ <#top>`_
