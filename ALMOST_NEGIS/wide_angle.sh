#!/bin/bash

# Wide angle seismic and radar profiles using the SeidarT program


# Define some constants
imfile=./negis1_crop2m.png
prjfile=./negis1.prj

# -----------------------------------------------------------------------------
# We'll do this one model at a time just to demonstrate the work flow

# Read in the project file and run the seismic model
python3 -m prjrun $prjfile --model s

# Create the gif from each snapshot
python3 -m im2gif $prjfile -c Vx -f 32
python3 -m im2gif $prjfile -c Vz -f 32

# Get the seismograms


# Convert the Fortran unformatted binary to CSV
python3 -m seidart_io -o csv -d 1 -f $prjfile


# Now run the radar model, create the gifs but convert to csv first just for 
# proof of concept
python3 -m prjrun $prjfile --model s

python3 -m seidart_io -o csv -d 1 -f $prjfile

python3 -m im2gif $prjfile -c Ex -f 1
python3 -m im2gif $prjfile -c Ez -f 1

# Get the radargrams

# Do the common offset



# Clear up some stuff 
# rm *.csv
