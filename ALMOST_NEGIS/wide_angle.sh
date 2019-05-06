#!/bin/bash

# Wide angle seismic and radar profiles using the SeidarT program


# Define some constants
imfile=./negis1_crop2m.png
prjfile=./negis1.prj

# -----------------------------------------------------------------------------
# We'll do this one model at a time just to demonstrate the work flow

# Read in the project file and run the seismic model
python3 -m prjrun $prjfile --model s

# Convert the Fortran unformatted binary to CSV
python3 -m seidart_io -o csv -d 1 -f $prjfile

# Create the gif from each snapshot


# Create the reciever array

