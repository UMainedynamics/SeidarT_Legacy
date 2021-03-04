orientation_tensor
###########################

*Wrapper to generate the Euler angles for the plunge and trend*
*then plot the results*

**Usage**

.. code-block:: bash

    orientation_tensor -o ANGFILE -n VALUE -P VALUE -t VALUE -a VALUE -A VALUE -S

**Inputs**

* ``-o``, ``--outputfile``

    Specify the file to save the outputs

* ``-n VALUE``, ``--npts VALUE``

    Total number of grains in synthetic sample

* ``-P VALUE``, ``--plunge VALUE``

    Plunge angle in degrees for center of single pole

* ``-t VALUE``, ``--trend VALUE``

    Trend angle in degrees for center of single pole

* ``-a VALUE``, ``--anglemin VALUE``

    Minimum angle deviation from center of single pole

* ``-A VALUE``, ``--anglemax VALUE``

    Maximum angle deviation from center of single pole

* ``-S``, ``--suppress_plotting``

    (OPTIONAL) Suppress all plotting

**Outputs**

* ``-o FILE``, ``--outputfile FILE``

    This will produce a delimited file, the name of which you can enter on the materials line(s) in the prj file.

`Back to top ↑ <#top>`_
