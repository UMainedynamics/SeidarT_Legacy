#!/usr/bin/env python3

#

import numpy as np
import glob
import argparse
from scipy.io import FortranFile

from evtk.hl import imageToVTK

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

	if channel == 'Ex' or channel == 'Ez':
		if line[0] == 'E':
			temp = line.split(',')
			if temp[1] == 'dt':
				edt = float(temp[2].rsplit()[0])
			if temp[1] == 'x':
				ex = float(temp[2].rsplit()[0])
			if temp[1] == 'z':
				ez = float(temp[2].rsplit()[0])
	else:
		if line[0] == 'S':
			temp = line.split(',')
			if temp[1] == 'dt':
				sdt = float(temp[2].rsplit()[0])
			if temp[1] == 'x':
				sx = float(temp[2].rsplit()[0])
			if temp[1] == 'z':
				sz = float(temp[2].rsplit()[0])


f.close()

# Define some plotting inputs
nx = nx + 2*cpml
ny = ny + 2*cpml
nz = nz + 2*cpml
npoints = nx*ny*nz
x = np.linspace(1, nx, num = nx)*dx
y = np.linspace(1, nx, num = ny)*dy
z = np.linspace(nz, 1, num = nz)*dz

# Add the source location to plot
if channel == 'Ex' or channel == 'Ez':
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

print('Creating GIF.')

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
		dat = f.read_reals('float32').reshape( (nx, ny, nz), order = "F" )

		# Normalize the values
		max_amplitude = np.abs(dat).max()
		
		# Zero out any values below our given threshold
		dat_normalize[np.abs(dat_normalize) < (max_amplitude*threshold) ] = 0.0
		duration = dt*ind

		if channel == 'Vx' or channel == 'Vz':
			time_label = 'Time (s): ' + str(np.round(duration, 5) )
		else:
			time_label = 'Time (s): ' + str(np.round(duration, 8) )

		# Reset the counter
		n = 1
		ind = ind + 1

	else:
		ind = ind + 1
		n = n + 1


# Dimensions 

ncells = nx * ny * nz 

npoints = (nx + 1) * (ny + 1) * (nz + 1) 

# Variables 

pressure = np.random.rand(ncells).reshape( (nx, ny, nz), order = 'C') 
temp = np.random.rand(npoints).reshape( (nx + 1, ny + 1, nz + 1)) 

imageToVTK("./image", cellData = {"Displacement" : dat} )