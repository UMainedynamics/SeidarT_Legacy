rcxdisplay
#########################

*Plot the amplitude timeseries based on a csv matrix*

**Usage**

.. code-block:: bash

    rcxdisplay -p [projectfile] -f [csv file with plottable output] -r [receiver location text file] -g [gain] -e [vertical exaggeration] -s [use if seismic]


**Inputs**

* ``-p`` project file (with the .prj extension)
* ``-f`` commonly receiver_array.csv or Ex.co.csv
* ``-r`` commonly receivers.xyz
* ``-g`` gain
* ``-e`` vertical exaggeration, commonly 0.1-0.05
* ``-s`` use if the model is seismic rather than electromagnetic


**Outputs**

Plot of amplitude value across area surveyed


`Back to top ↑ <#top>`_
