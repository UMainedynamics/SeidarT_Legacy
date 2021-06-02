prjbuild
#######################

*generates project text file template*

**Usage**

.. code-block:: bash

    prjbuild -i IMAGEFILE -p PROJECTFILE

**Inputs**

* ``-i``: .png file of cross sectional view of study site

    * No antialiasing
    * Use pixels dimensions proportional to unit dimensions
    * Start with small dimensions (150x150, 200x500), as larger
      dimensions quickly become more time intensive

* ``-p``: filename, including .prj extension


**Outputs**

* Project file to fill in based on site characteristics

