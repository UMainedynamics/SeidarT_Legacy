


import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt

from scipy.io import FortranFile



# -------------------------- Command Line Arguments ---------------------------

# parser = argparse.ArgumentParser(description="""This program builds a gif or 
# 	mp4 video from the set of image output of the FDTD modeling. The images are
# 	in unformatted Fortran binary, and each time step is saved as a PNG then
# 	the video/gif is compiled using a more efficient memory usage program (e.g.
# 	ImageMagick). """ )

# parser.add_argument( 'project_file', nargs=1, type=str,
# 						help='the full file path for the project file')

# parser.add_argument( '-c', '--channel', nargs = 1, type = str, required = True,
# 	help = """Specify whether a particular channel is going to be used. The 
# 	available channels are Ex, Ey, Ez, Vx, Vy, and Vz for the electric field 
# 	and seismic velocities, respectively.""")

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

# parser.add_argument( '-p', '--position', nargs = 2 type = float, 
# 	required = False, default = [10, 10], help = """Give the camera position in
# 	terms of elevation and azimuth, respectively.""" )


# # Get the arguments
# args = parser.parse_args()
# project_file = ''.join(args.project_file)
# channel = ''.join(args.channel)
# frame_rate = args.frames_per_second[0]
# num_steps = args.num_steps[0]
# threshold = args.threshold[0]
# output_format = args.output[0]
# ele_camera = arg.position[0]
# axi_camera = arg.position[1]


# Get the arguments
args = parser.parse_args()
project_file = 'dipping_bed.prj'
channel = 'Ex'
frame_rate = 40
num_steps = 40
threshold = 0.01
output_format = 'mp4'
ele_camera = 10
axi_camera = 10




# ===
# ----------------------- Definitions ----------------------- 


class Frame:
	def __init__(self, size=(640,480) ):
		self.fig = plt.figure()
		self.fig.set_size_inches(size[0]/100, size[1]/100)
		ax = self.fig.add_axes([0, 0, 1, 1], frameon=False, aspect=1)
		ax.set_xticks([])
		ax.set_yticks([])
		self.images = []
		self.background = []
		self.source_location = []

	def MakeBackground25(self, imfile):

		# If this is a 2.5D image we must extrude the 2D model into the third dimension
		if dim == '2.5':
			# Read the image
			im = mpimg.imread(imfile)
			self.background = np.zeros( (nz, nx, ny ) )
			
			for i in range(0, ny):
				self.background[:,:,i] = im 

		if dim == '3': # We'll construct 3D backgrounds later
			pass



	def plot(self, image, label='', extent=None ):

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

		it temp[1] == 'dim':
			dim = str(temp[2].rsplit()[0] )
		if temp[1] == 'nx':
			nx = int(temp[2].rsplit()[0] )
		if temp[1] == 'ny':
			ny = int(temp[2].rsplit()[0] )
		if temp[1] == 'nz':
			nz = int(temp[2].rsplit()[0] )
		if temp[1] == 'dx':
			dx = float(temp[2].rsplit()[0] )
		if temp[1] == 'dy':
			dy = float(temp[2].rsplit()[0] )
		if temp[1] == 'dz':
			dz = float(temp[2].rsplit()[0] )
		if temp[1] == 'cpml':
			cpml = int(temp[2].rsplit()[0])

	if channel == 'Ex' or channel == 'Ez':
		if line[0] == 'E':
			temp = line.split(',')
			if temp[1] == 'dt':
				edt = float(temp[2].rsplit()[0])
			if temp[1] == 'x':
				ex = float(temp[2].rsplit()[0])
			it temp[1] == 'y':
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

# --------------------------------- Read model --------------------------------

Frame.MakeBackground25(imfile)


# Plot the background


# ------------------------------ Plotting Inputs ------------------------------

# Define some plotting inputs
nx = nx + 2*cpml
nz = nz + 2*cpml
x = np.linspace(1, nx, num = nx)*dx
z = np.linspace(nz, 1, num = nz)*dz
extent = (cpml, (nx-cpml), (nz-cpml), cpml)

# Add the source location to plot
if channel == 'Ex' or channel == 'Ez':
	ex = (ex/dx + cpml+1)
	ey = (ey/dy + cpml+1)
	ez = (ez/dz + cpml+1)
	Frame.source_location = np.array([ex, ey, ez]) 
	dt = edt
else:
	sx = (sx/dx + cpml+1)
	sy = (sy/dy + cpml+1)
	sz = (sz/dz + cpml+1)
	Frame.source_location = np.array([sx, sy, sz])
	dt = sdt


