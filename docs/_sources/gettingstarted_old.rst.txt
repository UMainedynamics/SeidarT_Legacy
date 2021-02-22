[comment]: ======================================================================

# Getting Started <a name="getting_started"></a>

Geometries for the model domain within SeiDarT are initiated with a PNG image. The program identifies unique RGB values, setting material properties for each. For example, if you wanted to define a geometry with ice overlying bedrock, you would create a .png image that is one color for the ice and another for the rock below. Everyone has their preferences to generate images but [*GIMP*](https://www.gimp.org/downloads/install_help.html) and [*Inkscape*](http://wiki.inkscape.org/wiki/index.php/Installing_Inkscape) provide free and open software that are more than sufficient. When creating a PNG anti-aliasing must be turned off to avoid color boundary gradients.

Building images in *Inkscape* has some advantages other than being free. Saving a .svg to pdf allows the user to change the number of pixels and the spatial resolution of the image quite easy. With ghostscript the command

~~~
gs -q -dBATCH -dNOPAUSE -sDEVICE=png16m -sOutputFile=<file> -r96 <input_file>
~~~

will generate a PNG file from a PDF. The resolution -r can be varied to change the pixels. In *Inkscape* the image pixels can be set in *Document Properties*. When saving the SVG as PDF you will be prompted with options, and the value for *Resolution for rasterization (dpi):* will determine - in order to get the same pixel setting that you set in *Inkscape* - the value for the *-r* (resolution) option above. If you want to double the resolution just double this number (i.e. -r192 ).


To get started on a new project create a new folder and save the image to the folder. From the command line, change directories to the project folder then enter into the command line

~~~
prjbuild -i /path/to/geometry/image.png -o project_filename.prj
~~~

The project filename is optional if the -o option is omitted then the default text file *jordan_downs.prj* will be generated. For any of the executables the -h or --help option will provide a quick reference to any positional, optional or required arguments.

Using the text editor of choice, you can edit the .prj file. TODO: Describe what all is in this file?


[comment]: ======================================================================




[comment]: ======================================================================

## Routines <a name="routines"></a>

*prjbuild.py*
constructs a template and assigns default values from a PNG image.

*prjrun.py*
reads the project file assigns coefficients given that all the required fields are satisfied then runs the specified 2D forward model. You can suppress modeling and edit the stiffness and/or permittivity and conductivity coefficients. Once they are provided in the project file, they won't be computed or overwritten from the material values. If you would like to change the material values and recompute the tensor coefficients, you need to delete the existing tensor coefficients if included in the project file.

*im2gif.py*
Create a gif from the model outputs.

*arraybuild.py*
Create a timeseries of signals for each specified receiver

*rcxdisplay.py*
(originally codisplay in legacy code) Display survey outputs of the common offset survey.

*orientation_tensor.py*
Compute the Euler angles and orientation tensor for a fabric defined by its trend and plunge angles.


<u>Wrappers</u>

*common_offset.sh*
This is a wrapper that simulates a common offset survey. The receiver .xyz file is input to give the points of the survey and the source is offset from this location given the offsets for the x, y, and z directions.

*common_midpoint.sh*
This is similar to the common offset survey but it shifts the source and receiver away from a common midpoint. The midpoint is specified by the source location in the project file. By default the source will be to the viewer's right of the midpoint but to flip the location of the source and receiver, set the midpoint x-value to negative. 
