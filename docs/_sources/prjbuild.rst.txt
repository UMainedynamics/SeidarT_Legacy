prjbuild
#######################

*generates project text file template*

.. code-block:: bash

    prjbuild -i [input_image].png -o [output_project_file].prj

**Inputs**

* ``-i``: .png file of cross sectional view of study site

    * No antialiasing
    * Use pixels dimensions proportional to unit dimensions
    * Start with small dimensions (150x150, 200x500), as larger
      dimensions quickly become more time intensive

**Outputs**

* ``-o``: .prj filename to be filled in based on site characteristics