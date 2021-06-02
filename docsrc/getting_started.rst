Getting Started
########################

General workflow
******************************

Using SeidarT follows a relatively simple workflow.

#. You need two or three files to start:

  * A 2D image saved in *png* format.
  * A *csv* file listing the X,Y,Z coordinates of receivers for your survey
  * If your material is anisotropic, you need a file in the format delimited file specifying the

  Euler angles for a number of crystals, with one triplet per line. See an example orientation
  file and/or generate one using the ``orientation_tensor`` function.

#. Generate a project file (using ``prjbuild``) and edit that text file to set up your survey.
#. Create files describing the radar or seismic source (``sourcefunction``).
#. Choose the style of survey you want to do [single shot, common offset, common midpoint, or (in development) polarimetric] and run the calculations.
#. For single shot, you can create an animation of the wave propagation (``im2anim`` for 2D or ``vtkbuild`` for 2.5D).
#. Display your results as radar- or seismograms, or wiggle plots. You can also save the timeseries-receiver data in a *csv* file for further processing in different software.

Output from the seismic model is m/s and from the radar model is 

Files to generate or edit
******************************

* *PNG image (.png)*

    This defines the geometry of your model system. A good starting size is
    50 to 500 pixels for each direction. Each RGB color represents a different
    material, so **the file must be saved with no antialiasing**. Typically each pixel represents the same distance in x and z (in meters).
    To get started on a new project create a new folder and save the image
    to the folder. From the command line, change directories to the
    project folder then enter the following::

        prjbuild -i /path/to/geometry/image.png -p project_filename.prj

    Below, we describe the *prj* file structure and how to edit it.

* *receiver locations (text file, commonly receivers.xyz)*

    A comma separated list of X,Y,Z coordinates (one set per line,
    with X,Y,Z as the first line) for receiver locations. Can use pixels, but
    more typically meters as the units.


* *project file (.prj)*

    This file is the heart of the software. It defines domain values, material properties, and survey conditions for
    electromagnetic and seismic runs. Here, we identify what each line means and which to edit.
    All lines with # are comments. Bold text indicates a line the user should edit.


    ``I,fill2.png`` The image file associated with this *.prj* file.

    ``D,dim,2`` **Choose either 2D or 2.5D.** 2.5D is the 2D image extruded in the y-direction.

    ``D,nx,240`` Read from the image file. Do not change.

    ``D,ny,1`` **Number of pixels in the extruded direction if using 2.5D.**

    ``D,nz,50`` Read from the image file. Do not change.

    ``D,dx,1`` **Number of meters each pixel represents in the x direction.**

    ``D,dy,1`` **Number of meters each pixel represents in the y direction.**

    ``D,dz,1`` **Number of meters each pixel represents in the z direction.**

    ``D,cpml,20`` **Thickness of absorbing boundary layer. A typical value is 20.**

    ``D,nmats,3`` Read from the image file. Do not change.

    ``D,tfile,`` An attenuation processing value that is not yet implemented.

    ``M,1,ice1h,98/197/178,-10,2,910,0,0,TRUE,test.ang`` One comma-separated-values line per material (per color in the model image). **User should change/add
    the material name (see list), temperature (in degrees Celsius), density, porosity, water content,
    whether the material is anisotropic (TRUE or FALSE), and if anistropic, the name of the anisotropy file. Use a dummy value of 2 for attenuation,
    recognizing that attenuation is not yet incorporated in the calculations.** User should not change the material ID or R/G/B values. Note: Since large density
    gradients cause numerical instabilities, the
    density for air must be increased. A value of 400.0 works until a
    better formulation of the air-rock interface is implemented.

    ``S,dt,`` Timestep will be calculated automatically.

    ``S,time_steps,500`` **Decide how many timesteps you want the model to run.**

    ``S,x,100`` **x-coordinate of the seismic source**

    ``S,y,0`` **y-coordinate of the seismic source**

    ``S,z,0`` **z-coordinate of the seismic source**

    ``S,f0,60`` **Frequencyt of hte seismic source**

    ``S,theta,0`` **Inclination of the seismic source (+ is down)**

    ``S,phi,0`` **Angle of seismic source from x-axis in the x-y plane (+ is counterclockwise when viewed from above)**

    ``C,0.0,`` Stiffness tensor for each material. User *can* enter or change this manually if desired. If blank, calculated from materials information
    in the earlier section.

    ``E,dt,`` Timestep will be calculated automatically.

    ``E,time_steps,500`` **Decide how many timesteps you want the model to run.**

    ``E,x,100`` **x-coordinate of the radar source**

    ``E,y,0`` **y-coordinate of the radar source**

    ``E,z,0`` **z-coordinate of the radar source**

    ``E,f0,1e8`` **Frequencyt of hte radar source. 10-100MHz is a good range to start.**

    ``E,theta,0`` **Inclination of the radar source (+ is down)**

    ``E,phi,0`` **Angle of radar source from x-axis in the x-y plane (+ is counterclockwise when viewed from above)**

    ``P,0.0,`` Permittivity tensor for each material. User *can* enter or change this manually if desired. If blank, calculated from materials information
    in the earlier section.

* *orientation file*

    A delimited file of one entry of Bunge notation Euler angles per line.
    A typical number of entries is 500 to ensure a smooth data field.


