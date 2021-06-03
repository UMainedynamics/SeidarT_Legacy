Building documentation
#############################

This documentation is built with Sphinx 3.5.4 and rst2pdf 0.98.

To build the documentation, ensure you are in the SeidarT conda
environment, then install the requirements::

    pip install sphinx==3.5.4 sphinx_rtd_theme rst2pdf

Then, change directory into ``docsrc/``::

    cd SeidarT/docsrc

Then, use the ``make`` command to direct Sphinx to build
the documentation::

    make github

This will build documentation (both HTML and PDF) and move these
items to the folder where GitHub will look to render them (``docs/``).
Additionally, you can specify which documentation is built by using
the ``make html`` or ``make pdf`` in order to preview these.
The outputs of these intermediary commands will be in
``docsrc/_build/html`` or ``docsrc/_build/pdf``.

Additional information can be found
`here <https://github.com/sbernsen/SeidarT/blob/master/docsrc/README.md>`_.


.. note::

    Sphinx 4.0 and above depreciates a function that rst2pdf relies
    on to build the PDF file. At this time it is necessary to downgrade
    to Sphinx 3.5.4 in order to build documentation in PDF format.
