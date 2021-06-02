#### Dipping Bed Model

*Seismic and Electromagnetic Wavefields*

These are small, simple models to introduce the routines.

![Alt text](EXAMPLES/dipping_bed.png)
After running the install script without errors then adding the program file parent directy to the path, open a terminal and change directories into the EXAMPLES folder. Input into the command line

~~~
prjbuild -i dipping_bed.png -o dipping_bed.prj
~~~

Fill in the missing domain values for 1 meter spacing in both directions and an absorbing boundary of 30 pixels as

>D,dx,1
>D,dy,n/a
>D,dz,1
>D,cpml,30
>D,write,10

then add the following material values

>M,0,granite,0/0/0,-5,1,2540,1,10,False,
>M,1,water,0/100/200,-1,2,1000,0,0,False,
>M,2,ice1h,50/50/200,-2,3,910,2,30,False,
>M,3,wet_sand,80/120/120,-2,2,1500,4,40,False,
>M,4,ice1h,120/152/200,-10,2,910,2,5,False,
>M,5,basalt,200/0/0,-10,2,1200,0,0,False,
>M,6,air,200/240/255,0,0,1.3,0,0,False,

the source parameters for a seismic source just below the surface on the right side of the domain

>S,dt,
>S,time_steps,2500
>S,x,280
>S,y,0
>S,z,21
>S,f0,80
>S,theta,0
>S,phi,0

and the source parameters for a transverse radar pulse source in the same location as the seismic source except at the surface

>E,dt,
>E,time_steps,1600
>E,x,150
>E,y,0
>E,z,20
>E,f0,1e7
>E,theta,90

All lines that start with a '#' are commented lines but without the correct identifiers 'D', 'M', 'S', 'C', 'E', and 'P' the line will not be read by the computer so it is possible to make notes in the project file. If in doubt, and for good habit, use a hashtag.

For reference, you can use the RGB values to identify the materials so it is important to keep track of these when creating a geometry. In this case, the top very light blue/green is air (M,6), the blue dot (M,1) is water, the red line (M,5) is a dirty horizon of ice, the dark blue (M,2) is saturated porous ice (2% porosity, 30% water content) while the overlying lavender (M,4) is porous colder ice, the green (M,3) is a moderately saturated till layer, and all underlain by granite bedrock (M,0).

After filling in the domain and material values save then run the command

~~~
prjrun dipping_bed.prj -M n
~~~

Even though we have all the required fields entered the model didn't run because we used the '-M n' option. This is because basalt was inserted as an arbitrary material type in the project file as a placeholder and to output values. We can just as well use 'ice1h' or any other material in the python dictionary. The project file is now updated so we can edit the 'S,5' and 'P,5' coefficients to be between the 'ice1h' values and basalt values. Values closer to ice are reasonable with conductivity values closer to granite or basalt. I chose

>C,5.0,1.58e10, 8.44e9, 8.44e9, 1.58e10, 8.44e9, 1.58e10, 2.67e10, 2.67e10, 2.67e10, 1050.0

>P,5.0,3.8,3.8,3.8,1e-05,1e-05,1e-05

Since large density gradients cause numerical instabilities, the density for air must be increased. A value of 400.0 works until a better formulation of the air-rock interface is implemented. Now you can run the model (both seismic and electromagnetic)

~~~
prjrun dipping_bed.prj -M b
~~~

When that is finished let's build the GIF to view the wavefield. Starting with the seismic wavefield, enter

~~~
im2anim dipping_bed.prj -c -n 10 Vx -f 30
~~~

then repeat when the command line returns but change the '-c' option to 'Vz', 'Ex', or 'Ez'. Alternatively, you could change the frame rate (-f) to a higher or smaller number to adjust the speed of the GIF. We specified the 'D,write' parameter to 10 which is a little undersampled to view seismograms or radargrams but creating the GIF takes some time and we don`t need that much resolution. If you would like to create the seismograms decrease the write parameter between 2-4.


*Common Offset Survey*
