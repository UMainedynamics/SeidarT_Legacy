sourcefunction
#########################

*Creates .dat files needed to define the source impulse.*
*Does not need to be run if files are already present.*

.. code-block:: bash

    sourcefunction [-h] -p PROJECTFILE [-S SOURCETYPE] [-m MODELTYPE][-a AMPLITUDE]

**Inputs**


* ``-p``: Project file

    Name of .prj file

* ``-S``: Source type

    Specify the source type. Available wavelets are: gaus0, gaus1,
    gaus2 (gaussian n-th derivative), chirp, chirplet. (Default = gaus0)

* ``-m``: Model type

    Specify whether to construct the source for an electromagnetic or
    seismic model. e: electromagnetic, s: seismic, b: both

* ``-a``: Amplitude

    Input the scalar actor for source amplification. (Default = 1.0)

* ``-h``, ``--help``

    show this help message and exit

* ``-p PROJECTFILE``, ``--projectfile PROJECTFILE``

    The path to the project file

* ``-S SOURCETYPE``, ``--sourcetype SOURCETYPE``

    Specify the source type. Available wavelets are: gaus0, gaus1,
    gaus2 (gaussian n-th derivative), chirp, chirplet. (Default = gaus0)

* ``-m MODELTYPE``, ``--modeltype MODELTYPE``

    Specify whether to construct the source for an em or seismic model.
    s-seismic, e-electromagnetic, b-both

* ``-a AMPLITUDE``, ``--amplitude AMPLITUDE``

    Input the scalar factor for source amplification (Default = 1.0)


**Outputs**

.dat files in x, y, z for each model type
