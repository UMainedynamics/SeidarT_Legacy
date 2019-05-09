#!/bin/bash

# Wide angle seismic and radar profiles using the SeidarT program


# Define some constants
imfile=./negis1_crop2m.png
prjfile=./negis1.prj

# Reciever values
rcxi='70 0 70'
rcxf='370 0 70'
dr=4 #Large array for plotting purposes

# -----------------------------------------------------------------------------
# We'll do this one model at a time just to demonstrate the work flow

# Read in the project file and run the seismic model
python3 -m prjrun $prjfile --model s

# Get the seismograms. A CSV is created called 'reciever_array.csv' that 
# includes the timeseries for each reciever. Time is in the row direction.

# To view the array uncomment the following line
# python3 -m arrayplot $prjfile -c Vx -i $rcxi -f $rcxf -d $dr -L 1

# To see the t-x plot
python3 -m arrayplot $prjfile -c Vx -i $rcxi -f $rcxf -d $dr -g 1
python3 -m arrayplot $prjfile -c Vz -i $rcxi -f $rcxf -d $dr -g 1

# Create the gif from each snapshot. This may take a long time
python3 -m im2gif $prjfile -c Vx -f 20
python3 -m im2gif $prjfile -c Vz -f 20


# Convert the Fortran unformatted binary to CSV. This takes quite a while
# python3 -m seidart_io -o csv -d 1 -f $prjfile


# -----------------------------------------------------------------------------
# Now run the radar model, create the gifs and 

# python3 -m prjrun $prjfile --model s

# python3 -m im2gif $prjfile -c Ex -f 1
# python3 -m im2gif $prjfile -c Ez -f 1

# Get the radargrams

# Do the common offset



# Clear up some stuff 
# rm *.csv
