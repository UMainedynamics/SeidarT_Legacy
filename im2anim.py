#!/usr/bin/env python3

# From the set of image outputs, we can build a gif. The images will be in csv 
# or fortran unformatted binary


import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as anim
from scipy.io import FortranFile

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""This program builds a gif from 
	the set of image output of the FDTD modeling. The images can be in csv or 
	unformatted Fortran binary, however, the program runs faster to use the 
	latter. """ )

parser.add_argument( 'project_file', nargs=1, type=str,
						help='the full file path for the project file')

parser.add_argument( '-c', '--channel', nargs = 1, type = str, required = True,
	help = """Specify whether a particular channel is going to be used. The 
	available channels are Ex, Ez, Vx, and Vz for the electric field and 
	seismic velocities, respectively.""")

parser.add_argument( '-f', '--frames_per_second', nargs = 1, type = int, 
	required = False, default = 1, help = """The number of frames per second
	to build the GIF.""")

parser.add_argument( '-n', '--num_steps', nargs = 1, type = int, 
	required = True, help = """The time step interval between the images that
	are going to be used. Every time step is written to file which means that
	we can take any equally spaced images to create the gif with an 
	appropriate resolution, time to compute, and file size. For example, 
	n=20 means that every 20 images will be used thus significantly reducing 
	how long it takes to compile the gif.""")

parser.add_argument( '-t', '--threshold', nargs = 1, type = float,
	required = False, default=[0.0001], help = """Set values to zero when they 
	below a specific threshold. Default = 0.0001""")

parser.add_argument( '-o', '--output', nargs = 1, type = int, required = False,
	default = [0], help = """Specify the output format. 0 - GIF (default), 1 - MP4 """)

# Get the arguments
args = parser.parse_args()
project_file = ''.join(args.project_file)
channel = ''.join(args.channel)
frame_rate = args.frames_per_second[0]
num_steps = args.num_steps[0]
threshold = args.threshold[0]
output_format = args.output[0]
# ===
# ----------------------- Definitions ----------------------- 

class AnimatedGif:
	def __init__(self, size=(640,480) ):
		self.fig = plt.figure()
		self.fig.set_size_inches(size[0]/100, size[1]/100)
		ax = self.fig.add_axes([0, 0, 1, 1], frameon=False, aspect=1)
		ax.set_xticks([])
		ax.set_yticks([])
		self.images = []
		self.background = []
		self.source_location = []

	def add(self, image, label='', extent=None ):

		plt_im = plt.imshow(image,cmap='seismic', animated=True, extent=(0, (nx), (nz), 0))
		plt_bg = plt.imshow(self.background,alpha = 0.3, extent=extent, animated = True)
		plt.scatter(self.source_location[0], self.source_location[1], 
			marker = '*', s = 30, linewidths = 1, 
			edgecolor = (0.2, 0.2, 0.2, 1 ) )
		
		plt_txt = plt.text(extent[0] + 20, extent[2] + 20, label, color='red') # Lower left corner
		self.images.append([plt_im, plt_bg, plt_txt])

	def save(self, filename, frame_rate = 50):
		animation = anim.ArtistAnimation(self.fig, self.images, interval = frame_rate, blit = True)
		
		if output_format == 1:
			animation.save(filename, dpi = 200) 
		else:
			animation.save(filename, dpi = 200, writer = 'imagemagick')



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
		if temp[1] == 'nz':
			nz = int( temp[2].rsplit()[0] )
		if temp[1] == 'dx':
			dx = float(temp[2].rsplit()[0] )
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
nz = nz + 2*cpml
x = np.linspace(1, nx, num = nx)*dx
z = np.linspace(nz, 1, num = nz)*dz
extent = (cpml, (nx-cpml), (nz-cpml), cpml)




# Create the gif object
animated_gif = AnimatedGif( size=(nx, nz) ) 

# Load the model image
# animated_gif.background = np.zeros( [nz, nx, 3] )
animated_gif.background = mpimg.imread(imfile)
# animated_gif.background[cpml:(nz-cpml), cpml:(nx-cpml)] = mpimg.imread(imfile)

# Add the source location to plot
if channel == 'Ex' or channel == 'Ez':
	ex = (ex/dx + cpml+1)
	ez = (ez/dz + cpml+1)
	animated_gif.source_location = np.array([ex, ez]) 
	dt = edt
else:
	sx = (sx/dx + cpml+1)
	sz = (sz/dz + cpml+1)
	animated_gif.source_location = np.array([sx, sz])
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
		dat = f.read_reals(dtype = 'float32')
		dat = dat.reshape(nz, nx)


		# Normalize the values
		max_amplitude = np.abs(dat).max()
		dat_normalize = dat#/max_amplitude
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


