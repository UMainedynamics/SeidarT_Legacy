# SeidarT

Seismic and radar modeling using staggered grid finite differences (FDTD) and convolutional perfectly matched boundaries (CPML)

### Table of Contents
[Introduction](#introduction)  
[Installation](#install)  
[Getting Started](#getting_started)  
[Routines](#routines)  
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

Installation of these packages is different between MacOSx and Linux distributions. [*Anaconda*](https://docs.continuum.io/anaconda/) provides a convenient **Python** environment for installing and developing programs, however, Anaconda doesn't provide *GCC-7* or later so compilation of the **Fortran** code will return an error. This can be remedied by uninstalling **GCC** from *Anaconda* and upgrading via *apt* to *GCC-7* or greater. For *Anaconda* users, it's necessary to 

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

When the compilation is finished, we can add the folder to the path directory and the python path directory. Currently, this software is supported with **bash** so append the following lines to the *.bashrc* if using Ubuntu 
~~~
PATH=$PATH:/path/to/SeidarT/
export $PATH

PYTHONPATH=$PYTHONPATH:/path/to/SeidarT
export PYTHONPATH
~~~

or in the */etc/path* for Mac.
~~~

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

# Getting Started <a name="getting_started"></a>

Geometries are initiated with a PNG image and the program identifies unique RGB values. Everyone has their preferences to generate images but [*GIMP*](https://www.gimp.org/downloads/install_help.html) and [*Inkscape*](http://wiki.inkscape.org/wiki/index.php/Installing_Inkscape) provide free and open software that are more than sufficient. When creating a PNG anti-aliasing must be turned off to avoid color boundary gradients. 

To get started on a new project create a new folder and save the image to the folder. From the command line, change directories to the project folder then enter into the command line

~~~
prjbuild -i /path/to/geometry/image.png -o project_filename.prj
~~~

The project filename is optional if the -o option is omitted and the default text file *jordan_downs.prj* will be generated. For any of the executables the -h or --help option will provide a quick reference to any positional, optional or required arguments.

Using the text editor of choice, you can edit the .prj file. 


[comment]: ======================================================================




[comment]: ======================================================================

## Routines <a name="routines"></a>

*prjbuild.py*  
constructs a template and assigns default values from a PNG image.  

*prjrun.py*  
reads the project file assigns coefficients given that all the required fields are satisfied then runs the specified 2D forward model. You can suppress modeling and edit the stiffness and/or permittivity and conductivity coefficients. Once they are provided in the project file, they won't be computed or overwritten from the material values. If you would like to change the material values and recompute the tensor coefficients, you need to delete the existing tensor coefficients if included in the project file. 

*im2gif.py*  
Create a gif from the model outputs. Currently, this takes some time to run which you can speed up by increasing the 'write' value in the project file. 

*arrayplot.py*  
Plot the seismograms or radargrams for the wide angle survey. You can suppress plotting which will return a .csv file. A simple exponential gain function can be called and applied to each of the time series. 


*codisplay.py*  
Display the outputs of the common offset survey. This is also called to display the common midpoint survey. Similar to arrayplot.py, the gain function can be called.

<u>Wrappers</u>

*wide_angle.sh*  
Script to model and display an equally spaced reciever array. If the .dat files are not deleted, any node within the model domain can be called. The time series for the recievers are saved as 'reciever_array.csv'. A default gain is set in the script file but can be varied by editing the appropriate line. 

*common_offset.sh*  
This is a wrapper that shifts the source and queries the respective reciever location over a given distance and assembles the subsequent time series into a common offset image. The reciever trails the source so leave space for the reciever location. The initial source location is given in the corresponding .prj file. During the routine this value is changed each iteration then restored to the original value when it is complete. A default gain is set in the script file but can be varied by editing the appropriate line. 

*common_midpoint.sh*  
This is similar to the common offset survey but it shifts the source and reciever away from a common midpoint. The midpoint is specified by the source location in the project file. By default the source will be to the viewer's right of the midpoint but to flip the location of the source and reciever, set the midpoint x-value to negative. A default gain is set in the script file but can be varied by editing the appropriate line. 

#### *Note: The aspect ratio for the wide angle, common offset and common midpoint surveys determines the axis exaggeration. This will be updated in the future to be easier to adjust but to change this value edit the line 'ax.set_aspect(aspect=??)' in arrayplot.py and codisplay.py then run the plotting scripts individually not the wrapper scripts. 



[comment]: ======================================================================

## Examples <a name="examples"></a>

#### Dipping Bed Model  

This is a small simple model to introduce the routines. After running the install script without errors then adding the program file parent directy to the path, open a terminal and change directories into the EXAMPLES folder. Input into the command line

~~~
python3 -m prjbuild -i dipping_bed.png -o dipping_bed.prj
~~~

Fill in the domain values and material values as








