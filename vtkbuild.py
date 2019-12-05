#!/usr/bin/env python3

#

import numpy as np
import glob
import argparse
from scipy.io import FortranFile

from evtk.hl import imageToVTK

# -------------------------- Command Line Arguments ---------------------------
# parser = argparse.ArgumentParser(description="""This program builds a gif from
# 	the set of image output of the FDTD modeling. The images can be in csv or
# 	unformatted Fortran binary, however, the program runs faster to use the
# 	latter. """ )

# parser.add_argument( 'project_file', nargs=1, type=str,
# 						help='the full file path for the project file')

# parser.add_argument( '-c', '--channel', nargs = 1, type = str, required = True,
# 	help = """Specify whether a particular channel is going to be used. The
# 	available channels are Ex, Ez, Vx, and Vz for the electric field and
# 	seismic velocities, respectively.""")

# parser.add_argument( '-f', '--frames_per_second', nargs = 1, type = int,
# 	required = False, default = 1, help = """The number of frames per second
# 	to build the GIF.""")

# parser.add_argument( '-n', '--num_steps', nargs = 1, type = int,
# 	required = True, help = """The time step interval between the images that
# 	are going to be used. Every time step is written to file which means that
# 	we can take any equally spaced images to create the gif with an
# 	appropriate resolution, time to compute, and file size. For example,
# 	n=20 means that every 20 images will be used thus significantly reducing
# 	how long it takes to compile the gif.""")

# parser.add_argument( '-t', '--threshold', nargs = 1, type = float,
# 	required = False, default=[0.0001], help = """Set values to zero when they
# 	below a specific threshold. Default = 0.0001""")

# parser.add_argument( '-o', '--output', nargs = 1, type = int, required = False,
# 	default = [0], help = """Specify the output format. 0 - GIF (default), 1 - MP4 """)

#-- Get the arguments
# args = parser.parse_args()
# project_file = ''.join(args.project_file)
# channel = ''.join(args.channel)
# frame_rate = args.frames_per_second[0]
# num_steps = args.num_steps[0]
# threshold = args.threshold[0]
# output_format = args.output[0]


project_file = "dipping_bed.prj"
channel = 'Ex'
frame_rate = 20
num_steps = 10

# ------------------------------ Run the program ------------------------------

# Get the values we need
f = open(project_file)

for line in f:

	# Get the image file
	if line[0] == 'I':
		# There's a trailing new line value
		imfile = line[2:-1]

	# All domain inputs must be input except for nz and dy
	if line[0] == 'D':
		temp = line.split(',')

		if temp[1] == 'nx':
			nx = int( temp[2].rsplit()[0] )
		if temp[1] == 'ny':
			ny = int( temp[2].rsplit()[0] )
		if temp[1] == 'nz':
			nz = int( temp[2].rsplit()[0] )
		if temp[1] == 'dx':
			dx = float(temp[2].rsplit()[0] )
		if temp[1] == 'dy':
			dy = float(temp[2].rsplit()[0] )
		if temp[1] == 'dz':
			dz = float(temp[2].rsplit()[0] )
		if temp[1] == 'cpml':
			cpml = int( temp[2].rsplit()[0])

	if channel == 'Ex' or channel == 'Ey' or channel == 'Ez':
		if line[0] == 'E':
			temp = line.split(',')
			if temp[1] == 'dt':
				edt = float(temp[2].rsplit()[0])
			if temp[1] == 'x':
				ex = float(temp[2].rsplit()[0])
			if temp[1] == 'y':
				ey = float(temp[2].rsplit()[0])
			if temp[1] == 'z':
				ez = float(temp[2].rsplit()[0])
	else:
		if line[0] == 'S':
			temp = line.split(',')
			if temp[1] == 'dt':
				sdt = float(temp[2].rsplit()[0])
			if temp[1] == 'x':
				sx = float(temp[2].rsplit()[0])
			if temp[1] == 'y':
				sy = float(temp[2].rsplit()[0])
			if temp[1] == 'z':
				sz = float(temp[2].rsplit()[0])

f.close()

# Define some plotting inputs
nx = nx + 2*cpml
ny = ny + 2*cpml
nz = nz + 2*cpml
ncells = nx*ny*nz



# Create the coordinate system
X = np.linspace(1, nx, num = nx)*dx
Y = np.linspace(1, nx, num = ny)*dy
Z = np.linspace(nz, 1, num = nz)*dz
x = np.zeros([nx, ny, nz])
y = np.zeros([nx, ny, nz])
z = np.zeros([nx, ny, nz])

for i in range(0, nx):
    for j in range(0, ny):
        for k in range(0, nz):
            x[i,j,k] = X[i]
            y[i,j,k] = Y[j]
            z[i,j,k] = Z[k]


# Add the source location to plot
if channel == 'Ex' or channel == 'Ey' or channel == 'Ez':
    ex = (ex/dx + cpml+1)
    ey = (ey/dy + cpml+1)
    ez = (ez/dz + cpml+1)
    source_location = np.array([ex, ey, ez])
    dt = edt
else:
    sx = (sx/dx + cpml+1)
    sy = (sy/dy + cpml+1)
    sz = (sz/dz + cpml+1)
    source_location = np.array([sx, sy, sz])
    dt = sdt

# Proceed accordingly to the channel flag

# Check if the .dat files are still around
files = glob.glob(channel + '*.dat')
ind = 0
files.sort()

# We'll start counting with the first frame
n=num_steps
for fn in files:
    if n == num_steps:
        f = FortranFile(fn, 'r')
#        dat = f.read_reals('float32').reshape( (nx, ny, nz), order = "F" )
        dat = f.read_reals(dtype = 'float32')
#        print(dat.min())
        # dat = dat.reshape(nz, ny, nx) # This is what works with Ey and Ez 
        dat = dat.reshape(nz, ny, nx)
        
        # Zero out any values below our given threshold
        duration = dt*ind
        
#        if channel == 'Vx' or channel == 'Vz':
#            time_label = 'Time (s): ' + str(np.round(duration, 5) )
#        else:
#            time_label = 'Time (s): ' + str(np.round(duration, 8) )
#			
		# Reset the countern = 1 
        ind = ind + 1
        n = 1
        
        vtkfilename = "./image" + channel + "." + str(ind)
        imageToVTK(vtkfilename, cellData = {"Displacement" : dat} )
        
    else:
        ind = ind + 1
        n = n + 1

# Dimensions 

# Variables 

