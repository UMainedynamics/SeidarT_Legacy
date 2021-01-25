arraybuild
###################

*generates a csv file (receiver_array.csv) with the timeseries in columns*
*r each receiver for the given channel*

**Usage**

.. code-block:: bash

    arraybuild [-h] -p PRJFILE -r RCXFILE [-i INDEX] -c CHANNEL


**Inputs**

* ``-h``, ``--help``

    show this help message and exit

* ``-p PRJFILE``, ``--prjfile PRJFILE``

    The project file path.

* ``-r RCXFILE``, ``--rcxfile RCXFILE``

    the file path for the text file of receiver locations, in comma separated rows

* ``-i INDEX``, ``--index INDEX``

    Indicate whether the receiver file contains coordinate
    indices or if these are the locations in meters.
    Default (0 - meters)

* ``-c CHANNEL``, ``--channel CHANNEL``

    The channel to query.

    | For radar: Ex, Ey, or Ez
    | For seismic: Vx, Vy, Vz, S1, S2, S3, S4, S5, S6
    | S# represents the stress tensor:

    * S1 - sigma_11
    * S2 - sigma_22
    * S3 - sigma_33
    * S4 - sigma_23
    * S5 - sigma_13
    * S6 - sigma_12


**Outputs**

CSV of amplitude values for all pixels surveyed (*receiver_array.csv*)



`Back to top ↑ <#top>`_
