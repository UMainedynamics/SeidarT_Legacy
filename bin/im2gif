#!/usr/bin/python3

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
						help='the full file path for the project file', 
						default=None)

parser.add_argument( '-c', '--channel', nargs = 1, type = str, required = False,
	help = """Specify whether a particular channel is going to be used. The 
	available channels are Ex, Ez, Vx, and Vz for the electric field and 
	seismic velocities, respectively.""", default = None)

parser.add_argument( '-f', '--frames_per_second', nargs = 1, type = int, 
	required = False, default = 1, help = """The number of frames per second
	to build the GIF.""")

# Get the arguments
args = parser.parse_args()
project_file = ''.join(args.project_file)
channel = ''.join(args.channel)
frame_rate = args.frames_per_second[0]


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

		plt_im = plt.imshow(image,cmap='seismic', animated=True)
		plt_bg = plt.imshow(self.background,alpha = 0.3, extent=extent)
		plt.scatter(self.source_location[0], self.source_location[1], 
			marker = '*', s = 30, linewidths = 1, 
			edgecolor = (0.2, 0.2, 0.2, 1 ) )
		
		plt_txt = plt.text(extent[0] + 20, extent[2] + 20, label, color='red') # Lower left corner
		self.images.append([plt_im, plt_bg, plt_txt])

	def save(self, filename, frame_rate = 50):
		animation = anim.ArtistAnimation(self.fig, self.images, interval = frame_rate)
		animation.save(filename, writer='imagemagick')


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
		if temp[1] == 'write':
			write = int(temp[2].rsplit()[0])

	if line[0] == 'S':
		temp = line.split(',')
		if temp[1] == 'dt':
			sdt = float(temp[2].rsplit()[0])
		if temp[1] == 'x':
			sx = float(temp[2].rsplit()[0])
		if temp[1] == 'z':
			sz = float(temp[2].rsplit()[0])

	if line[0] == 'E':
		temp = line.split(',')
		if temp[1] == 'dt':
			edt = float(temp[2].rsplit()[0])
		if temp[1] == 'x':
			ex = float(temp[2].rsplit()[0])
		if temp[1] == 'z':
			ez = float(temp[2].rsplit()[0])


f.close()

# Define some plotting inputs
nx = nx + 2*cpml
nz = nz + 2*cpml
x = np.linspace(1, nx, num = nx)*dx
z = np.linspace(nz, 1, num = nz)*dz
extent = (cpml, nx-cpml, nz-cpml, cpml)

ex = (ex + cpml+1)/dx
ez = (ez + cpml+1)/dz
sx = (sx + cpml+1)/dx
sz = (sz + cpml+1)/dz

# Create the gif object
animated_gif = AnimatedGif( size=(nx, nz) ) 

# Load the model image
animated_gif.background = mpimg.imread(imfile)

# Add the source location to plot
if channel == 'Ex' or channel == 'Ez':
	animated_gif.source_location = np.array([ex, ez]) 
	dt = edt
else:
	animated_gif.source_location = np.array([sx, sz])
	dt = sdt

print('Creating GIF.')
# Proceed accordingly to the channel flag
if channel:
	
	# Check if the .dat files are still around
	files = glob.glob(channel + '*.dat')

	# if they aren't then get the csv files
	if not files:
		pass
		# files = glob.glob(channel+'*.csv')
		# files.sort()
		# ind = 0
		# for fn in files:
		# 	dat = np.genfromtxt(fn, delimiter = ',')
			
		# 	duration = dt*ind
		# 	time_label = 'Time (s): ' + str(np.round(duration, 5) )
		# 	animated_gif.add(dat, time_label, extent)
			
		# 	ind = ind + write

		# animated_gif.save(channel + '.gif', fps = frame_rate)

		
	else:
		ind = 0
		files.sort()
		for fn in files:
			f = FortranFile(fn, 'r')
			dat = f.read_reals(dtype = 'float32')
			dat = dat.reshape(nz, nx)

			duration = dt*ind
			time_label = 'Time (s): ' + str(np.round(duration, 5) )
			animated_gif.add(dat, time_label, extent)

			ind = ind + write

		animated_gif.save(channel + '.gif', frame_rate = frame_rate)

else:

	# First do the E-field
	animated_gif.source_location = np.array([ex, ez])
	Ex_files = glob.glob('Ex*.dat')
	animated_gif.save('Ez.gif')

	Ez_files = glob.glob('Ez*.dat')
	animated_gif.save('Ez.gif')

	# Now do the seismic
	animated_gif.source_location = np.array([sx, sz])
	Vx_files = glob.glob('Vx*.dat')
	animated_gif.save('Vx.gif')

	Vz_files = glob.glob('Vz*.dat')
	animated_gif.save('Vz.gif')	
