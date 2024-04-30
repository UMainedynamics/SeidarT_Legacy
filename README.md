# SeidarT

<!-- ### Table of Contents -->
<!-- [Introduction](#introduction)  
[Installation](#install)  
[Auto-Install](#auto-installation) 
[Manual Install](#manual-installation)  
[Hardware Requirements](#hardware-requirements)  
[Operating System Requirements](#operating-system-requirements)   -->

Full documentation appears in the docs folder. 

[comment]: ======================================================================

## Introduction <a name="introduction"></a>

The Seismic and Radar Toolbox (SeidarT) is a collaboration between researchers at the Universities of Maine and Washington to provide an open source platform for forward modeling mechanical and electromagnetic wave propagation. The major objective of the project is to easily and quickly implement isotropic and anisotropic complex geometries and/or velocity structures to develop prior constraints for - not limited too - investigating, estimating, and imaging englacial ice structure, sub-glacial boundary conditions on the sub-regional scale. Larger problems would require the curvature of the Earth to be taken into consideration, but many glacier seismic and radar experiments do not expand into regional parameter estimation and velocity modeling.

Much of the Staggered Grid FDTD code has been adopted from the *SEISMIC_CPML* software provided by [Computational Infrastucture for Geophysics (CIG)](https://geodynamics.org/cig/software/). Further details to the backend numerical code can be found in the [References](#references) section below.

## Installation <a name="install"></a>

The dynamical programming language of **Python3** is used as a command line interface to run the more computationally extensive modeling schemes in **Fortran**. There are a number of ways to install this software depending on your desired level of control over the process. Most users should be fine with the "automatic installation" in the [section below](#auto-install).

SeidarT package binaries are publicly availble on the [PyPi repository](https://pypi.org/project/seidart/) and source code can be found at github (not yet public). 

### "Auto" installation <a name="auto-install"></a>

Extract the *install* directory from the *install.tar.gz* or *install.bz* for Unix/Linux or Windows users, respectively

There are 2 install scripts, *full_install.sh* and *full_install.bat*, which will cover Linux, MacOS 13 and 14, and Windows , respectively. It's not necessary to know or do much more than execute a few command line entries via a bash terminal or powershell terminal. A virtual environment is created using the Miniconda/Anaconda package manager. This will avoid causing system incompatibilities and complicated software dependencies. Documentation for managing conda environments with Miniconda or Anaconda can be found [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). The full Anaconda release has a GUI called Navigator for managing environments. This can be found [here](https://docs.anaconda.com/free/navigator/tutorials/manage-environments/).  

The *full_install.?* checks for an Anaconda/Miniconda and will install it if it isn't found. If Anaconda/Miniconda is not defined in your 'PATH' variable than it will be installed using the default install location. After installing the Conda package, the *seidart* environment is built using pre-defined dependencies in a YAML (Yet Another Markup Language) file. Both *Bash* executables and Python modules are built during install. In order to use either, the environment must be active. This can be easily done from the *Bash* command line interface (CLI) using the command
```
conda activate seidart
```

[comment]: ======================================================================
### Manual installation<a name="manual-install-pip"></a>

The full repo can be found on GitHub and hosted on PyPi. *SeidarT* has been tested on Python 3.11 and is not yet supported with Python 3.12. For users that prefer building virtual environments with Anaconda, the install folder contains the *seidart-environment.yml* or it can be found in the root directory of the GitHub repo. The command 
```
conda env create -f seidart-environment.yml
```
will install all dependencies and the *seidart* package. For users who prefer more control in their installation, their is a small list of dependencies that must be met. These are:  *gcc*>10, *gfortran*, *ghostscript*, *imagemagick*, *numpy*, *pandas*, *matplotlib*, *scipy*, *glob2*, *pyevtk*, *mplstereonet*. Following install of all dependencies, 
```
pip install seidart
```
will pull and install the package from PyPi. 

### Hardware Requirements <a name="hardware-requirement"></a>

*SeidarT* was tested and developed on a quad core 5th gen i7 processor with 16 Gb of RAM without any burden on the system so a typical modern laptop is sufficient for many application. When running models with large domains or a high number of time steps, the computational load is obviously increased, however the storage requirements become more significant. It can be easy to fill up 10's of Gb of storage, but an external drive can resolve that problem. The Apple M-chips may have compatability issues with particular types of software and Python packages, but we have maintained a relatively simple design along with leveraging some of the most commonly used Python packages which should help to mitigate any issues with computing on an M-chip. 

### Operating System Requirements <a name="operating-system"></a>

All of the development was carried out on a Linux operating system and limited to Debian, Ubuntu, Solus 2, and Fedora. No compatibility issues between Linux flavors arose. The binaries are built on Github Actions for Windows 10 and 11 (latest), MacOS 13 and 14 (latest), and most flavors of Linux. Cross-platform usability is one of the core tenets in the development of the software and needs to be maintained in future development. 