#!/bin/bash

:'
POLARIMETRIC_POINT.SH wrapper script for conducting a polarimetric point estimate
for a reciever in plane with the source and one cross polar reciever. We can
rotate the source, reciever, or both.

INPUT   
    dphi - the change in rotation angle in the x-y plane
    
'

#! Define inputs - This will change to command line arguments


# ============================== Print Help Menu ============================== 
#! Add options as command line arguments are incorporated. 
# For now we can add a description of what we're doing

# =========================== Create the output File ==========================
# We want to save the time series for each rotation. There will be two files: 
# one for in plane and the other for cross polar. 

# Create the files. Time is in the i-direction, Angle is in the j-direction


# ================================ Get to work ================================

# Is it easier to rotate the tensors relative to the source/reciever or rotate 
# the source/reciever?

# Save the initial values

# --------------------------------------------------------
# If we keep the source stationary and rotate the reciever

# --------------------------------------------------------
# If we rotate the source we must also rotate the reciever

# Loop through each rotation angle 

    # Run the model 

    # Get the Ex and Ey values then rotate the timeseries 

    # Append the values to the csv file

    # Change the source angle


# Reset to the initial values


# Make a quick and dirty plot 


# Say goodbye
echo "Look! A pair of boobs!"
echo "(.Y.)"