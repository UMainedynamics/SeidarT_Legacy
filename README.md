# SeidarT

Seismic and radar modeling using staggered grid finite differences (FDTD) 


### Table of Contents
[Introduction](#introduction)  
[Installation](#install)  
[Getting Started](#getting_started)  
    [Examples](#examples)  
[Appendix](#appendix)  
[References](#references)  

[comment]: ======================================================================

## Introduction <a name="introduction"></a>

The Seismic and Radar Toolbox (SeidarT) is a collaboration between researchers at the Universities of Maine and Washington to provide an open source platform for forward modeling mechanical and electromagnetic wave propagation. The major objective of the project is to easily and quickly implement isotropic and anisotropic complex geometries and/or velocity structures to develop prior constraints for - not limited too - investigating, estimating, and imaging englacial ice structure, sub-glacial boundary conditions on the sub-regional scale. Larger problems would require the curvature of the Earth to be taken into consideration, but many glacier seismic and radar experiments do not expand into regional parameter estimation and velocity modeling. 

Much of this code has been adopted from the *SEISMIC_CPML* software provided by [Computational Infrastucture for Geophysics (CIG)](https://geodynamics.org/cig/software/). Further details to the backend numerical code can be found in the [References](#references) section below.


[comment]: ======================================================================
## Installation <a name="install"></a>

The dynamical programming language of **Python3** is used as a front end to run the more computationally extensive modeling schemes in **Fortran**. The following system dependencies are required:

**Python3, Fortran95, GCC, pip, git, ghostscript, imageMagick**

and additionally, the **Python** dependencies:

*numpy*, *scipy*, *glob3*, *matplotlib*.

Installation of these packages is different between MacOSx and Linux distributions. [*Anaconda*](https://docs.continuum.io/anaconda/) provides a convenient **Python** environment for installing and developing programs, however, Anaconda doesn't provide *GCC-7* or later so compilation of the **Fortran** code will return an error. This can be remedied by uninstalling **GCC** from *Anaconda* and upgrading via *apt* to *GCC-7* or greater. 

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
bash install.sh
~~~

When the compilation is finished, we can add the folder to the path directory and the python path directory. Currently, this software is supported with **bash** so append the following lines to the *.bashrc* if using Ubuntu or in the */etc/path* for Mac.

~~~
placeholder/for/now
~~~

There are a variety of ways to edit the documents but for simplicity change directories to the home folder 

~~~
cd ~
~~~

and input into the command line

~~~
sudo nano .bashrc
~~~

or 

~~~ 
sudo nano /etc/paths
~~~

then scroll down to the bottom of the file and append the path. Save and close (*Ctrl-x* then *Y* and enter) the file then check to make sure it is included in the path 

~~~
echo $PATH
echo $PYTHONPATH
~~~

[comment]: ======================================================================

## Getting Started <a name="getting_started"></a>

Geometries are initiated with a PNG image and the program identifies unique RGB values. Everyone has their preferences to generate images but [*GIMP*](https://www.gimp.org/downloads/install_help.html) and [*Inkscape*](http://wiki.inkscape.org/wiki/index.php/Installing_Inkscape) provide free and open software that are more than sufficient. When creating a PNG anti-aliasing must be turned off to avoid color boundary gradients. 

To get started on a new project create a new folder and save the image to the folder. From the command line, change directories to the project folder then enter into the command line

~~~
prjbuild -i /path/to/geometry/image.png -o project_filename.prj
~~~

The project filename is optional if the -o option is omitted and the default file *jordan_downs.prj* will be generated. For any of the executables the -h or --help option will provide 

[comment]: ======================================================================




[comment]: ======================================================================




