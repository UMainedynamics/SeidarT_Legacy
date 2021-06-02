# SeidarT

### Table of Contents
[Introduction](#introduction)  
[Installation](#install)  


Full documentation appears in the docs folder. 

[comment]: ======================================================================

## Introduction <a name="introduction"></a>

The Seismic and Radar Toolbox (SeidarT) is a collaboration between researchers at the Universities of Maine and Washington to provide an open source platform for forward modeling mechanical and electromagnetic wave propagation. The major objective of the project is to easily and quickly implement isotropic and anisotropic complex geometries and/or velocity structures to develop prior constraints for - not limited too - investigating, estimating, and imaging englacial ice structure, sub-glacial boundary conditions on the sub-regional scale. Larger problems would require the curvature of the Earth to be taken into consideration, but many glacier seismic and radar experiments do not expand into regional parameter estimation and velocity modeling.

Much of this code has been adopted from the *SEISMIC_CPML* software provided by [Computational Infrastucture for Geophysics (CIG)](https://geodynamics.org/cig/software/). Further details to the backend numerical code can be found in the [References](#references) section below.

## Installation <a name="install"></a>

The dynamical programming language of **Python3** is used as a front end to run the more computationally extensive modeling schemes in **Fortran**. There are a number of ways to install this software depending on your desired level of control over the process. Most users should be fine with the "automatic installation" in the [section below](#auto-install).

### "Auto" Installation <a name="auto-install"></a>

First, some background:

This method will install SeidarT on most Unix-based systems using a combination of Anaconda and Fortran compilation. It's not necessary to know or do much more than execute a few command line entries. However, it will help to know that the installer will create what's called a "virtual environment" or more specifically an "Anaconda environment". Without getting into details, this conda environment is basically a way to create a virtual Python configuration that will let this software run but not affect your system's Python configuration, which is a fancy way of saying that conda will avoid causing cascading bad effects on your computer. In short: conda is nice because it plays nicely with various systems. Note well that **you will have to activate this conda environment whenever you wish to run SeidarT.** This is not hard but it is critical. Read on for instructions on how to both install and activate your SeidarT conda environment.

1. First, make sure you have the proper prerequisites. On Debian/Ubuntu systems, open a Terminal window and do

    ~~~
    sudo apt update
    sudo apt install gcc-9 gfortran git
    ~~~

    **or** in a Mac OS X Terminal window with [Homebrew](https://brew.sh), copy/paste/enter the following

    ~~~
    brew update
    brew install gcc git
    ~~~

    *Note to OS X users: you may have to install [Homebrew](https://brew.sh) if you do not have it.*

2. Now, change to the directory where you'd like SeidarT to go, clone it from GitHub, and change into the SeidarT directory created by Git

    ~~~
    cd /path/to/parent/directory
    git clone https://github.com/sbernsen/SeidarT
    cd SeidarT
    ~~~

3. Run the auto-installer script and follow the prompts

    ~~~
    bash install.sh
    ~~~

    This script will first check if Anaconda or Miniconda already exists on your system. It will download and install Miniconda if necessary. It will also add conda to your `$PATH`, which is to say that it will tell your computer where to look to use the conda software. Then, it will use conda to create a SeidarT environment and install the rest of the requirements into that environment. Finally, it will compile the Fortran code in this repository.

4. Activate the SeidarT environment (you will have to do this whenever you want to use SeidarT)

    ~~~
    conda activate SeidarT
    ~~~

    *Note: you may have to open a new Terminal window for your system to recognize `conda` as a valid command.*

5. You should now be ready to run SeidarT scripts!

    Remember that whenever you would like to use SeidarT in the future, you will have to activate the SeidarT conda environment first by using the command in step 4.

[comment]: ======================================================================
### Manual Installation Instructions <a name="manual-install"></a>

This is the more involved and less sure-fire installation method, probably best left for intermediate level Unix users and up. The following system dependencies are required:

**Python3, Fortran95, GCC, pip, git, dos2unix, ghostscript, imageMagick**

and additionally, the **Python** dependencies: *numpy*, *scipy*, *matplotlib*, *mplstereonet* (optional for viewing fabric distributions).

First, install what you will need to compile the Fortran code. This can be done via **apt** with
~~~
sudo apt update
sudo apt install gcc-10 git dos2unix ghostscript imagemagick python3.8 python3-numpy python3-vtk python3-pip
~~~

and with homebrew using the commands

~~~
brew update
brew install gcc git dos2unix ghostscript imagemagick numpy vtk python pip
~~~

Then, install the rest of the dependencies using pip

~~~
pip3 install -U matplotlib scipy pyevtk vtk mplstereonet
~~~

Anaconda installation:

~~~
conda create -n seidart python=3 pip git ghostscript imagemagick numpy matplotlib scipy pyevtk vtk
conda activate seidart
pip install mplstereonet
~~~

To download the software, open a terminal (*Ctrl-Alt-t* for Ubuntu and *Ctrl-Opt-Shift-t* for Mac) and change directories

~~~
cd /path/to/parent/directory
~~~

to where you would like SeidarT to be located then clone the files from *github*

~~~
git clone git@github.com:sbernsen/SeidarT.git
~~~

This will create a folder named **SeidarT** that contains all executables and modules. Change directories to the **SeidarT** folder

~~~
cd SeidarT
~~~

and run the install file

~~~
bash noconda_install.sh
~~~

When the compilation is finished, we can add the folder to the path directory and the python path directory. Currently, this software is supported with **bash** so append the following lines to the *.bashrc* if using Ubuntu

~~~
export PATH=$PATH:/path/to/SeidarT/bin

export PYTHONPATH=$PYTHONPATH:/path/to/SeidarT/bin
~~~

and if Python 2 is the default version create an alias in the *.bashrc* file by adding the line in the alias section

~~~
alias python=python3
~~~

Depending on the OS release (El Capitan, High Sierra, Mojave, etc.) and whether you have *Anaconda* installed appending a path might be different. Anaconda may set aliases so troubleshooting on a Mac can be cumbersome. Before editing the */etc/path*, *.bash_profile*, *.profile* or *.bashrc* file it is a good idea to create a backup especially if you are not familiar with either or any of those files. To do this copy the original to a new name. For example,

~~~
cp <location/of/path/definitions> <location/of/path/definitions>_original
~~~

that way you can always revert back to the working script.


There are a variety of ways to edit the documents but for simplicity change directories to the home folder

~~~
cd ~
~~~

and input into the command line

~~~
sudo nano .bashrc
~~~

or either

~~~
sudo nano /etc/paths
sudo nano ~/.profile
~~~

then scroll down to the bottom of the file and append the path. Save and close (*Ctrl-x* then *Y* and enter) the file then check to make sure it is included in the path

~~~
echo $PATH
echo $PYTHONPATH
~~~
