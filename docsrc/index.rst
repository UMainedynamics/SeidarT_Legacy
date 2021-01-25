SeidarT |version|
#####################################

| |github|
| **Multi-layer Fortran-based geophysical wave modeling**
| *by Steven Bernsen and Christopher Gerbi*

.. only:: comment

    the following are definitions for raw html to insert in | substitution-tags |


.. |github| raw:: html

   <a href="https://github.com/sbernsen/SeidarT" target="_blank">https://github.com/sbernsen/SeidarT</a>

.. |CIG| raw:: html

   <a href="https://geodynamics.org/cig/software/" target="_blank">Computational Infrastucture for Geophysics (CIG)</a>


The Seismic and Radar Toolbox (SeidarT) is a collaboration between
researchers at the Universities of Maine and Washington to provide
an open source platform for forward modeling mechanical and
electromagnetic wave propagation. The major objective of the
project is to easily and quickly implement isotropic and anisotropic
complex geometries and/or velocity structures to develop prior
constraints for - not limited to - investigating, estimating,
and imaging englacial ice structure, sub-glacial boundary conditions
on the sub-regional scale. Larger problems would require the
curvature of the Earth to be taken into consideration, but many
glacier seismic and radar experiments do not expand into regional
parameter estimation and velocity modeling.

Much of this code has been adopted from the SEISMIC_CPML software
provided by |CIG|.
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
    routines
    examples

.. _extras:

.. toctree::
    :maxdepth: 2
    :caption: Extra material

    appendix
    references

.. _modules:

Code documentation
========================

The modules available in SeidarT are organized by type below.

------------

.. toctree::
    :maxdepth: 2
    :caption: Initialize

    prjbuild

.. toctree::
    :maxdepth: 2
    :caption: Run

    sourcefunction
    prjrun
    common_offset
    common_midpoint

.. toctree::
    :maxdepth: 2
    :caption: Plot and visualize

    arraybuild
    rcxdisplay
    wiggleplot
    im2anim
    vtkbuild


Function index
==================

Need to look something up?

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
