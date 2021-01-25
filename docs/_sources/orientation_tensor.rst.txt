orientation_tensor
###########################

*Wrapper to generate the euler angles for the plunge and trend*
*then plot the results*

**Usage**

.. code-block:: bash

    common_offset -o [project_file].prj -n [VALUE] -P [VALUE] -t [VALUE] -a [VALUE] -A [VALUE]

**Inputs**

* ``-o``, ``--outputfile``

    Specify the file to save the outputs

* ``-n VALUE``, ``--npts VALUE``

    Total number of grains in synthetic sample

* ``-P VALUE``, ``--plunge VALUE``

    Plunge angle in degrees

* ``-t VALUE``, ``--trend VALUE``

    Trend angle in degrees

* ``-a VALUE``, ``--anglemin VALUE``

    Minimum angle deviation

* ``-A VALUE``, ``--anglemax VALUE``

    Maximum angle deviation

* ``-S``, ``--suppress_plotting``

    (OPTIONAL) Suppress all plotting

**Outputs**

* ``-o FILE``, ``--outputfile FILE``


`Back to top ↑ <#top>`_
