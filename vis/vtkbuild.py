#!/usr/bin/env python3

# 

import numpy as np
import glob
import argparse
from scipy.io import FortranFile
from definitions import *
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

domain, material, seismic, electromag = loadproject(
    project_file,
    Domain(), 
    Material(),
    Model(),
    Model()
)

# Define some plotting inputs
domain.cpml = int(domain.cpml[0])
nx = int(domain.nx[0]) + 2*domain.cpml
ny = int(domain.ny[0]) + 2*domain.cpml
nz = int(domain.nz[0]) + 2*cpml
ncells = nx*ny*nz

# Create the coordinate system
domain.dx = float(domain.dx[0])
domain.dy = float(domain.dy[0])
domain.dz = float(domain.dz[0])

X = np.linspace(1, nx, num = nx)*domain.dx
Y = np.linspace(1, nx, num = ny)*domain.dy
Z = np.linspace(nz, 1, num = nz)*domain.dz
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
    electromag.x = float(electromag.x[0])
	electromag.y = float(electromag.y[0])
	electromag.z = float(electromag.z[0])
    ex = electromag.x/domain.dx + domain.cpml+1
    ey = electromag.y/domain.dy + domain.cpml+1
    ez = electromag.z/domain.dz + domain.cpml+1
    source_location = np.array([ex, ey, ez])
    dt = float(electromag.dt[0])
else:
    seismic.x = float(seismic.x[0])
	seismic.y = float(seismic.y[0])
	seismic.z = float(seismic.z[0])
    sx = seismic.x/domain.dx + domain.cpml+1
    sy = seismic.y/domain.dy + domain.cpml+1
    sz = seismic.z/domain.dz + domain.cpml+1
    source_location = np.array([sx, sy, sz])
    dt = float(seismic.dt[0])

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


