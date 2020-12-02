#!/usr/bin/env python3

# 

import numpy as np
import glob
import argparse
from scipy.io import FortranFile

from pyevtk.hl import imageToVTK

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""This program builds .VTI 
    (Visualization Toolkit Image) files from the 3d array outputs of the FDTD 
    modeling. These files can be displayed using Paraview.""" )

parser.add_argument(
    '-p', '--prjfile', 
    nargs=1, type=str, required = True,
    help='the full file path for the project file'
)

parser.add_argument(
    '-c', '--channel', 
    nargs = 1, type = str, required = True,
	help = """Specify whether a particular channel is going to be used. The
	available channels are Ex, Ez, Vx, and Vz for the electric field and
	seismic velocities, respectively."""
)

parser.add_argument(
    '-n', '--num_steps', 
    nargs = 1, type = int, required = True, 
    help = """The time step interval between the images that
	are going to be used. Every time step is written to file which means that
	we can take any equally spaced images to create the gif with an
	appropriate resolution, time to compute, and file size. For example,
	n=20 means that every 20 images will be used thus significantly reducing
	how long it takes to compile."""
)

#-- Get the arguments
args = parser.parse_args()
project_file = ''.join(args.prjfile)
channel = ''.join(args.channel)
num_steps = args.num_steps[0]

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

# There are slight differences between the seismic and radar domains in
# terms of staggered grid geometry 
if channel == 'Ex':
	NX = nz
	NY = ny
	NZ = nx-1
elif channel == 'Ey':
	NX = nz
	NY = ny-1
	NZ = nx
elif channel == 'Ez':
	NX = nz-1
	NY = ny
	NZ = nx
else:
    NX = nz
    NY = ny 
    NZ = nx 

# We'll start counting with the first frame
n=num_steps
for fn in files:
    if n == num_steps:
        f = FortranFile(fn, 'r')
        dat = f.read_reals(dtype = 'float32')
        dat = dat.reshape(NX, NY, NZ)
        
        # Zero out any values below our given threshold
        duration = dt*ind
        		
		# Reset the countern = 1 
        ind = ind + 1
        n = 1
        
        vtkfilename = "./image" + channel + "." + str(ind)
        imageToVTK(vtkfilename, cellData = {"Displacement" : dat} )
        
    else:
        ind = ind + 1
        n = n + 1


