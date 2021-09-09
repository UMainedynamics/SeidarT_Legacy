SeidarT |version|
#####################################

| `https://github.com/UMainedynamics/SeidarT <https://github.com/umainedynamics/SeidarT>`_
| **Geophysical wave modeling through complex media**
| **Primary creator:** *Steven Bernsen*
| **with:** *Christopher Gerbi, Ian Nesbitt, Ann Hill, Senthil Vel, Knut Christianson, Seth Campbell, Ben Hills*

| *Support from this work comes from* `NSF awards <https://nsf.gov/awardsearch/showAward?AWD_ID=1643301&HistoricalAwards=false>`_ *1643301 and 1643353.*

The Seismic and Radar Toolbox (SeidarT) is a collaboration between
the Universities of Maine and Washington to provide
an open source platform for forward modeling elastic (seismic) and
electromagnetic (radar) wave propagation. The major objective of the
project is to easily and quickly implement isotropic and anisotropic
complex geometries and/or velocity structures to investigate,
and plan field campaigns to image local and regional subsurface structure,
particularly for the cryosphere. Larger problems require the
curvature of the Earth to be taken into consideration.

Much of this code has been adopted from the SEISMIC_CPML software
provided by `Computational Infrastucture for Geophysics (CIG) <https://geodynamics.org/cig/software/>`_.
Further details to the backend numerical code can be found in
the :doc:`references` section.


.. _tutorial:

.. toctree::
    :numbered:
    :maxdepth: 2
    :caption: Tutorial guide

    about
    installation
    getting_started
    examples
    routines

.. _modules:

.. toctree::
    :maxdepth: 2
    :caption: Command reference

    prjbuild
    sourcefunction
    prjrun
    common_offset
    common_midpoint
    orientation_tensor
    arraybuild
    rcxdisplay
    wiggleplot
    im2anim
    vtkbuild
    implot
    imvector
    vectoranim

.. _extras:

.. toctree::
    :maxdepth: 2
    :caption: Extra material

    png-creation
    build-docs
    materials
    units
    references

Search documentation
======================

Need to look something up?

* :ref:`search`
