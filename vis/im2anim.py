#!/usr/bin/env python3

# From the set of image outputs, we can build a gif. The images will be in csv
# or fortran unformatted binary


import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.io import FortranFile
from definitions import *

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""This program builds a gif from
	the set of image output of the FDTD modeling. The images can be in csv or
	unformatted Fortran binary, however, the program runs faster to use the
	latter. """ )

parser.add_argument(
    '-p','--prjfile', 
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
    '-f', '--frames_per_second', 
    nargs = 1, type = int, required = False, default = 1, 
    help = """The number of frames per second
	to build the GIF."""
)

parser.add_argument(
    '-n', '--num_steps', 
    nargs = 1, type = int, required = True, 
    help = """The time step interval between the images that
	are going to be used. Every time step is written to file which means that
	we can take any equally spaced images to create the gif with an
	appropriate resolution, time to compute, and file size. For example,
	n=20 means that every 20 images will be used thus significantly reducing
	how long it takes to compile the gif."""
)

# parser.add_argument(
#     '-t', '--threshold', 
#     nargs = 1, type = float, required = False, default=[0.0001], 
#     help = """Set values to zero when they
# 	below a specific threshold. Default = 0.0001"""
# )

parser.add_argument(
    '-o', '--output', 
    nargs = 1, type = int, required = False, default = [0], 
    help = """Specify the output format. 0 - GIF (default), 1 - MP4 """
)

# Get the arguments
args = parser.parse_args()
project_file = ''.join(args.prjfile)
channel = ''.join(args.channel)
frame_rate = args.frames_per_second[0]
num_steps = args.num_steps[0]
threshold = args.threshold[0]
output_format = args.output[0]
# ===

# 03/16/2020 matplotlib got COVID-19 and now there are issues saving images to
# GIF so we're  trying some things 
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

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
nz = int(domain.nz[0]) + 2*domain.cpml
domain.dx = float(domain.dx[0])
domain.dz = float(domain.dz[0])

if channel == 'Ex':
    nx = nx-1

if channel == 'Ez':
	nz = nz-1

x = np.linspace(1, nx, num = nx)*domain.dx
z = np.linspace(nz, 1, num = nz)*domain.dz
extent = (
    domain.cpml, 
    (nx-domain.cpml), 
    (nz-domain.cpml), 
    domain.cpml
)

# Create the gif object
animated_gif = AnimatedGif( size=(nx, nz) )
animated_gif.output_format = output_format
# Load the model image
# animated_gif.background = np.zeros( [nz, nx, 3] )
animated_gif.background = mpimg.imread(domain.imfile)
# animated_gif.background[cpml:(nz-cpml), cpml:(nx-cpml)] = mpimg.imread(imfile)

# Add the source location to plot
if channel == 'Ex' or channel == 'Ez':
    electromag.x = float(electromag.x[0])
    electromag.z = float(electromag.z[0])
    ex = electromag.x/domain.dx + domain.cpml+1
    ez = electromag.z/domain.dz + domain.cpml+1
    animated_gif.source_location = np.array([ex, ez])
    dt = float(electromag.dt[0])
else:
    seismic.x = float(seismic.x[0])
    seismic.z = float(seismic.z[0])
    sx = seismic.x/domain.dx + domain.cpml+1
    sz = seismic.z/domain.dz + domain.cpml+1
    animated_gif.source_location = np.array([sx, sz])
    dt = float(seismic.dt[0])

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
		dat = f.read_reals(dtype = 'float32')
		dat = dat.reshape(nz, nx)
		# Normalize the values
		max_amplitude = np.abs(dat).max()
		dat_normalize = dat#/np.max([max_amplitude,1])
		# dat_normalize[ dat_normalize < -1.0 ] = -1.0
		# dat_normalize[ dat_normalize > 1.0 ] = 1.0
		# Zero out any values below our given threshold
		dat_normalize[np.abs(dat_normalize) < (max_amplitude*threshold) ] = 0.0
		duration = dt*ind
		if channel == 'Vx' or channel == 'Vz':
			time_label = 'Time (s): ' + str(np.round(duration, 5) )
		else:
			time_label = 'Time (s): ' + str(np.round(duration, 8) )
		animated_gif.add( dat_normalize, time_label, extent)
		# Reset the counter
		n = 1
		ind = ind + 1
	else:
		ind = ind + 1
		n = n + 1


if output_format == 1:
	imfilename = channel + '.mp4'
else:
	imfilename = channel + '.gif'


animated_gif.save(imfilename, frame_rate = frame_rate)


