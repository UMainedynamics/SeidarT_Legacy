common_midpoint
########################

*creates receiver setup by defining final and middle locations,*
*and receiver spacing across surface*

**Usage**

.. code-block:: bash

    common_midpoint -f [project_file].prj -t [X] -o [X] -d [distance between each receiver]


**Inputs**

* .prj file

    Filled in, except for permittivity and stiffness coefficients, and dt

* ``-t``: Location of the last receiver

    X pixel location ONLY, on the png file

* ``-o``: Location of half of the distance between the first and last receiver 

    X pixel location ONLY, on the png file

* ``-d``: Spacing between each receiver

OPTIONAL: If wanting to run seismic instead of radar, add component ``-s`` to end of code line


**Outputs**




`Back to top ↑ <#top>`_
