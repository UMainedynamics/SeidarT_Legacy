#!/usr/bin/env python3

# Create a receiver array and plot the timeseries

import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from imdefinitions import *
from class_definitions import *

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""This program creates a csv file of time series for each receiver location
	listed in the the specified receiver file.""" )

parser.add_argument(
    '-p', '--prjfile',
    nargs = 1, type = str, required = True,
    help = 'The project file path'
)

parser.add_argument(
    '-r', '--rcxfile',
    nargs=1, type=str, required = True,
    help='the file path for the text file of receiver locations'
)

parser.add_argument(
    '-i', '--index',
    nargs = 1, type = int, default = [0], required = False,
    help = """Indicate whether the receiver file contains coordinate indices or
    if these are the locations in meters. Default (0 - meters)"""
)

parser.add_argument(
    '-c', '--channel',
    nargs = 1, type = str, required = True,
    help = """The channel to query. """
)

# Get the arguments
args = parser.parse_args()
project_file = ''.join(args.prjfile)
receiver_file = ''.join(args.rcxfile)
channel = ''.join(args.channel)
rind = args.index[0] == 0

# ==================== Create the object and assign inputs ====================
# We don't need the material values
domain, material, seismic, electromag = loadproject(
    project_file,
    Domain(),
    Material(),
    Model(),
    Model()
)

xyz = np.genfromtxt(
    receiver_file,
    delimiter = ',',
    names = True,
    dtype = float
)

# We need to make sure the recievers are ordered correctly and the absorbing
# boundary is corrected for
# First check to see if the inputs are indices or
domain.dim = domain.dim[0]
cpml = int(domain.cpml[0])
if rind:
    receiver_locations = np.vstack(
        [
            xyz['X']/float(domain.dx[0]) + cpml,
            xyz['Y']/float(domain.dy[0]) + cpml,
            xyz['Z']/float(domain.dz[0]) + cpml
        ]
    ).T
else:
    receiver_locations = np.vstack(
        [
            xyz['X'] + cpml,
            xyz['Y'] + cpml,
            xyz['Z'] + cpml
        ]
    ).T

# Adjust the object fields relative to the cpml
domain.nx = int(domain.nx[0]) + 2*cpml
domain.ny = int(domain.ny[0]) + 2*cpml
domain.nz = int(domain.nz[0]) + 2*cpml

if channel == 'Ex' or channel == 'Ey' or channel == 'Ez':
	source_location = np.array(
		[
			float(electromag.x[0])/float(domain.dx[0]) + cpml,
			float(electromag.y[0])/float(domain.dy[0]) + cpml,
			float(electromag.z[0])/float(domain.dz[0]) + cpml
		]
	)
else:
	source_location = np.array(
		[
			float(seismic.x[0])/float(domain.dx[0]) + cpml,
			float(seismic.y[0])/float(domain.dy[0]) + cpml,
			float(seismic.z[0])/float(domain.dz[0]) + cpml
		]
	)

# Get the timeseries for each receiver. It will be saved in a file called receiver_array.csv
getrcx(channel, receiver_locations, domain)
