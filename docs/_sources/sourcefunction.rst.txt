sourcefunction
#########################

*Creates .dat files needed to define the source impulse.*
*Does not need to be run if files are already present.*

**Usage**

.. code-block:: bash

    sourcefunction [-h] -p PRJFILE -S SOURCETYPE -m [s e b] -a AMPLITUDE


**Inputs**

* ``-p``: PRJFILE

    Path to the .prj file, which must have dt and permittivity and stiffness tensors completed
    (usually via the ``-m n`` option of ``prjrun``).

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

**Outputs**

.dat files in x, y, z for each model type


`Back to top ↑ <#top>`_
