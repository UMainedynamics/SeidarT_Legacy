Installation
###################

The dynamical programming language of Python3 is used as a front
end to run the more computationally extensive modeling schemes in
Fortran. There are a number of ways to install this software
depending on your desired level of control over the process. Most
users should be fine with the ":ref:`auto-install`" in the
section below.

If you are an advanced user and looking to take more in-depth
control of the installation process, you can use the
:ref:`manual-install`.

.. note::

    **Mac OS X users**: you will have to install
    `Homebrew <https://brew.sh/>`_ prior to installing this software
    if you do not already have it.


.. _auto-install:

Auto-installation
*********************************

First, some background for the uninitiated:

This method will install SeidarT on most Unix-based systems using a
combination of Anaconda and Fortran compilation. It's not necessary
to know or do much more than execute a few command line entries.
However, it will help to know that the installer will create what's
called a "virtual environment" or more specifically an "Anaconda
environment". Without getting into details, this conda environment
is basically a way to create a virtual Python configuration that
will let this software run but not affect your system's Python
configuration, which is a fancy way of saying that conda will avoid
causing cascading bad effects on your computer. In short: conda is
nice because it plays nicely with various systems. Note well that
you will have to activate this conda environment whenever you wish
to run SeidarT. This is not hard but it is critical. Read on for
instructions on how to both install and activate your SeidarT conda
environment.


#. First, make sure you have the proper prerequisites
    On Mac OS X systems with `Homebrew <https://brew.sh/>`_,
    copy/paste/enter the following two lines::

        brew update
        brew install gcc git

    **or** on Debian/Ubuntu systems that use the :code:`aptitude`
    package manager, do the following to get those prerequisites::

        sudo apt update
        sudo apt install gcc-9 gfortran git

#. Get the software
    Change to the directory where you'd like SeidarT to go,
    clone it from GitHub, and change into the SeidarT directory
    created by Git::

        cd /path/to/parent/directory
        git clone https://github.com/sbernsen/SeidarT
        cd SeidarT

#. Install
    Run the auto-installer script and follow the prompts::

        bash install.sh

    This script will first check if Anaconda or Miniconda already
    exists on your system. It will download and install Miniconda
    if necessary. It will also add conda to your :code:`$PATH`,
    which is to say that it will tell your computer where to look
    to use the conda software. Then, it will use conda to create a
    SeidarT environment and install the rest of the requirements into
    that environment. Finally, it will compile the Fortran code in
    this repository.

You should now be able to move on to the :doc:`getting_started`
section.

`Back to top ↑ <#top>`_

.. _manual-install:

Manual installation instructions
***********************************

This is the more involved and less sure-fire installation method,
probably best left for intermediate level Unix users and up. The
following system dependencies are required:

**Python3, Fortran95, GCC, pip, git, dos2unix, ghostscript,**
**imageMagick**

and additionally, these Python dependencies are also required:
*numpy, scipy, matplotlib, mplstereonet*
(optional for viewing fabric distributions).

#. Get system prerequisites
    First, install what you will need to compile the Fortran code. This
    can be with Homebrew (on OS X) using the commands::

        brew update
        brew install gcc git dos2unix ghostscript imagemagick numpy vtk python pip

    and via **apt** (on Linux) with::

        sudo apt update
        sudo apt install gcc-10 git dos2unix ghostscript imagemagick python3.8 python3-numpy python3-vtk python3-pip

#. Install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
    `Anaconda <https://www.anaconda.com/products/individual#Downloads>`_
    will work as well, but miniconda is a smaller initial installation,
    and will only install what you need.

#. Get Python prerequisites
    From a Terminal window in which the :code:`conda` command is accessible,
    run the following commands::

        conda create -n seidart python=3 pip git ghostscript imagemagick numpy matplotlib scipy pyevtk vtk
        conda activate seidart
        pip install mplstereonet
    
#. Get the software
    ::

        cd /path/to/parent/directory
        git clone git@github.com:sbernsen/SeidarT.git
        cd SeidarT

#. Run the installer
    ::

        bash manual_install.sh

#. Update PATH
    When the compilation is finished, we can add the folder to the path
    directory and the python path directory. Currently, this software is
    supported with bash so append the following lines to the
    :code:`~/.bashrc` file if using Ubuntu::

        export PATH=$PATH:/path/to/SeidarT/bin

        export PYTHONPATH=$PYTHONPATH:/path/to/SeidarT/bin

    and if Python 2 is the default version, create an alias by adding this
    line to your aliases (either in :code:`~/.bashrc` or
    :code:`~/.bash_aliases`) ::

        alias python=python3

    .. note::
        Notes for inexperienced users:

        Depending on the OS release (El Capitan, High Sierra, Mojave, etc.)
        and whether you have Anaconda installed appending a path might be
        different. Anaconda may set aliases so troubleshooting on a Mac can
        be cumbersome. Before editing the :code:`/etc/path`,
        :code:`.bash_profile`, :code:`.profile`, or :code:`.bashrc` files,
        it is a good idea to create a backup especially if you are not
        familiar with either or any of those files. To do this copy the
        original to a new name. For example, ::

            cp <location/of/path/definitions> <location/of/path/definitions>_original

        that way you can always revert back to the working script.

        There are a variety of ways to edit the documents but for simplicity
        change directories to the home folder::

            cd ~

        and input into the command line::

            sudo nano .bashrc
        
        and append the :code:`export PATH=...` lines at the bottom.
        Save and close the file (*CTRL+X*, then *Y* and enter) then check
        to make sure it is included in the path::

            . ~/.bashrc
            echo $PATH
            echo $PYTHONPATH


................

* :ref:`genindex`
* :ref:`search`


`Back to top ↑ <#top>`_

