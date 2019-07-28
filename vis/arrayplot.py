#!/usr/bin/env python3

# Create a reciever array and plot the timeseries

import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.io import FortranFile

# import io_functions as iof

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""This program creates an equally
	spaced array of recievers given the x, y, and z coordinates. If the model 
	isn't specified as 2.5d in the project file then """ )

parser.add_argument( 'meta_file', nargs=1, type=str, 
						help='the file path for the survey metadata', 
						default=None)

parser.add_argument( '-c', '--channel', nargs = 1, type = str, required = True,
	help = """The channel to query. """)

parser.add_argument('-g', '--gain', nargs = 1, type = float, required = False,
	help = "The exponential value for 2^m of the gain function (default=None)", default = None)

parser.add_argument('-L', '--layout', nargs = 1, type = int, required = False,
	help = """Plot the reciever layout (1/0). This suppresses generating the 
	t-x timeseries plot.""", default = [0])

parser.add_argument('-S', '--suppress_plotting', nargs = 1, type = int, 
	help = """Suppress all plotting (1/0). The reciever outputs will be located in 
	the reciever_array.csv file.""", required = False, default = [0])

parser.add_argument('-e', '--exaggeration', nargs=1, type = float, required = False,
	help = """Set the aspect ratio between the x and y axes for 
	plotting. Default is 0.5""", default = [0.5])


# ---
# ---
# ======================= Create the class variables ==========================

# The timeseries plot
class Array:
	def __init__(self):
		super().__init__()
		self.build()

	def build(self):
		# Initialize variables
		self.dt = None
		self.dx = None
		self.dy = None
		self.dz = None
		self.source_location = None # Must be m-by-2 array in indices [x_ind, y_ind]
		self.initial_position = None
		self.final_position = None
		self.reciever_locations = None
		self.modelfile = None
		self.dr = None
		self.timeseries = None
		self.t = None
		self.channel = None
		self.color = ['red', 'blue'] # The default positive and negative fill
		self.nx = None
		self.ny = None
		self.nz = None

	# -------------------------- Function Definitions -------------------------
	def getrcx(self):
	# input rcx as an n-by-2 array integer values for their indices. 

		all_files = glob.glob(self.channel + '*.dat')
		all_files.sort()

		m = len(all_files)
		n = len(self.reciever_locations[:,1])

		self.timeseries = np.zeros([m,n])

		i = 0
		for fn in all_files:
			npdat = read_dat(fn, self.nx, self.nz)
			for j in range(0, n):
				pass
				# Don't forget x is columns and z is rows
				self.timeseries[i,j] = npdat[ self.reciever_locations[j,1], self.reciever_locations[j,0] ]	
			
			i = i + 1

		# Save the array as csv for other types of processing
		np.savetxt("reciever_array.csv", self.timeseries, delimiter = ",")

	# -------------------------------------------------------------------------
	def tsplot(self):
		
		fig, ax = plt.subplots()
		# fig.set_size_inches(8, 10)
		m,n = self.timeseries.shape 
		self.t = np.linspace(1, m, m)*self.dt 

		# extent = (1, n, self.t[-1], self.t[0])
		if gain:
			gain_function = np.zeros([m,n])
			for j in range(0, m):
				gain_function[j,:] = np.exp(-j/(2**gain))
			im = ax.imshow(self.timeseries/gain_function, cmap = 'Greys' )
		else:
			im = ax.imshow(self.timeseries, cmap = 'Greys')

		# Label the x-axis
		ax.xaxis.tick_top()
		ax.xaxis.set_label_position('top')
		ax.set_xlabel(r"Reciever #")
		
		# Label the y-axis
		if self.channel == 'Vx' or self.channel == 'Vz':
			mult = 1e2
		else:
			mult = 1e6

		time_locations = np.linspace(1, m, 10) 
		time_labels = np.round( time_locations*self.dt*mult, 4)
		ax.set_ylabel(r"Two way travel time (s)")
		plt.yticks(ticks=time_locations, labels = time_labels.astype(str) )

		# Other figure operations
		if self.channel == 'Vx' or self.channel == 'Vz':
			plt.figtext(0.30, 0.07, 'x $10^{-2}$')	
		else:
			plt.figtext(0.30, 0.07, 'x $10^{-6}$')

		ax.set_aspect(aspect = exaggeration	)
		plt.show()

	# -------------------------------------------------------------------------
	def arraybuild(self):
		
		# Calculate the length
		rcx_len = np.sqrt( sum( (self.final_position - self.initial_position)**2) )
		
		if rcx_len == 0:
			rxc_len = 1 # This is for common offset
			nrcx = 1
			xz = np.zeros([nrcx, 2])
		else:
			nrcx = np.floor(rcx_len/self.dr)
			xz = np.zeros([nrcx.astype(int),2])
		
		 

		xz[:,0] = np.linspace(self.initial_position[0], 
			self.final_position[0], nrcx)/self.dx
		xz[:,1] = np.linspace(self.initial_position[2], 
			self.final_position[2], nrcx)/self.dz

		self.reciever_locations = xz.astype(int) + self.cpml

	# -------------------------------------------------------------------------
	def plot_layout(self):

		figmod, axmod = plt.subplots()

		img = mpimg.imread(self.modelfile)
		axes_extent = [0, img.shape[1],img.shape[0], 0 ]

		implot = axmod.imshow(img, extent = extent)

		# add white upside down triangle recievers
		axmod.scatter(self.reciever_locations[:,0]- self.cpml, 
			self.reciever_locations[:,1] - self.cpml, 
			marker = 'v', s = 30, c = (0.8, 0.8, 0.8, 1), 
		    linewidths = 0.5, edgecolor = (0.2, 0.2, 0.2, 1 ) )

		# add white star source
		axmod.scatter(self.source_location[0], self.source_location[1], 
			marker = '*', s = 30, c = (0.8, 0.8, 0.8, 1 ), 
		    linewidths = 1, edgecolor = (0.2, 0.2, 0.2, 1 ) )

		axmod.set_ylabel("z-indice")
		axmod.set_xlabel("x-indice")

		plt.show()

# ------------- One more quick global definition

def read_dat(fn, nx, ny):
	f = FortranFile(fn, 'r')
	dat = f.read_reals(dtype = 'float32')
	dat = dat.reshape(ny, nx)
	f.close()

	return(dat)


# ==================== Create the object and assign inputs ====================
array = Array()

# Fill in object fields 

# Get the arguments
args = parser.parse_args()
meta_file = ''.join(args.meta_file)
array.channel = ''.join(args.channel)


# Optional inputs
gain = args.gain[0]
layout = args.layout[0] == 1
showplots = args.suppress_plotting[0] == 1
exaggeration = args.exaggeration[0]

# --------------------- Get the values from the meta file ---------------------

project_file = meta_file.split('.')
project_file = '.'.join( project_file[:-3] ) + '.prj'

f = open(meta_file)

# If running the wrapper functions, the meta file will save the same each time,
# but for whatever reason we'll assume that this isn't the case
for line in f:

	 temp = line.rstrip().rsplit()

	 if temp[0] == 'delta:':
	 	array.dr = float(temp[1])
	 
	 if temp[0] == 'initial_position:':
	 	array.initial_position = np.asarray(temp[1:]).astype(float)

	 if temp[0] == 'final_position:':
	 	array.final_position = np.asarray( temp[1:]).astype(float)

f.close()


# --------------------------- Query the project file --------------------------
#!!!!!! Add the 2.5D components later

# Get the values we need
f = open(project_file)

for line in f:

	# Get the image file
	if line[0] == 'I':
		# There's a trailing new line value
		array.modelfile = line[2:-1]

	# All domain inputs must be input except for nz and dy
	if line[0] == 'D':
		temp = line.split(',')

		if temp[1] == 'nx':
			array.nx = int( temp[2].rsplit()[0] )
		if temp[1] == 'nz':
			array.nz = int( temp[2].rsplit()[0] )
		if temp[1] == 'dx':
			array.dx = float(temp[2].rsplit()[0] )
		if temp[1] == 'dz':
			array.dz = float(temp[2].rsplit()[0] )
		if temp[1] == 'cpml':
			array.cpml = int( temp[2].rsplit()[0])
		if temp[1] == 'write':
			array.write = int(temp[2].rsplit()[0])

	# Source
	
	if array.channel == 'Ex' or array.channel == 'Ez':
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

# ----------------------------- Create the object -----------------------------


# Adjust the object fields relative to the cpml
array.nx = array.nx + 2*array.cpml
array.nz = array.nz + 2*array.cpml
array.dx = array.dx
array.dz = array.dz


extent = (array.cpml, array.nx-array.cpml, array.nz-array.cpml, array.cpml)
if array.channel == 'Ex' or array.channel == 'Ez':
	array.source_location = np.array([ (ex + array.cpml)/array.dx, 
		(ez + array.cpml)/array.dz])
	array.dt = edt
else:
	array.source_location = np.array([(sx + array.cpml)/array.dx,
		(sz + array.cpml)/array.dz]) 
	array.dt = sdt

# Construct the array
array.arraybuild()

# Get the timeseries for each reciever
array.getrcx()

if showplots:
	quit()

if layout:
	array.plot_layout()
else:
	array.tsplot()

