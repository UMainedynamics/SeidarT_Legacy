#!/usr/bin/env python3

# From the set of image outputs, we can build a gif. The images will be in csv 
# or fortran unformatted binary


import numpy as np
import argparse

import display_functions as df


# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""This program builds a gif from 
	the set of image output of the FDTD modeling. The images can be in csv or 
	unformatted Fortran binary, however, the program runs faster to use the 
	latter. """ )

parser.add_argument( 'project_file', nargs=1, type=str, 
						help='the full file path for the project file', default=None)

parser.add_argument( '-c', '--channel', nargs = 1, type = str, required = False,
	help = """Specify whether a particular channel is going to be used. The 
	available channels are Ex, Ez, Vx, and Vz for the electric field and 
	seismic velocities, respectively.""")



# Get the arguments
args = parser.parse_args()
project_file = ''.join(args.project_file)
channel = ''.join(args.channel)

# ----------------------- Definitions ----------------------- 