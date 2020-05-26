#!/usr/bin/env python3

# Create a receiver array and plot the timeseries

import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.io import FortranFile

# import io_functions as iof
from imdefinitions import * 

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""This program creates an equally
	spaced array of receivers given the x, y, and z coordinates. If the model 
	isn't specified as 2.5d in the project file then """ )

parser.add_argument(
    '-p', '--prjfile',
    nargs = 1, type = str, required = True,
    help = 'The project file path.'
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

parser.add_argument(
    '-g', '--gain',
    nargs = 1, type = float, required = False,
    help = "The exponential value for 2^m of the gain function (default=None)", default = None
)

parser.add_argument(
    '-L', '--layout', 
    nargs = 1, type = int, required = False,
    help = """Plot the receiver layout (1/0). This suppresses generating the
    t-x timeseries plot.""", default = [0]
)

parser.add_argument(
    '-S', '--suppress_plotting',
    nargs = 1, type = int,
    help = """Suppress all plotting (1/0). The receiver outputs will be located in
    the receiver_array.csv file.""", required = False, default = [0]
)

parser.add_argument(
    '-e', '--exaggeration',
    nargs=1, type = float, required = False,
    help = """Set the aspect ratio between the x and y axes for
    plotting. Default is 0.5""", default = [0.5]
)


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
        self.receiver_locations = None
        self.modelfile = None
        # self.dr = None 
        self.timeseries = None
        self.t = None
        self.channel = None
        self.color = ['red', 'blue'] # The default positive and negative fill
        self.nx = None
        self.ny = None
        self.nz = None
        self.dim = None
        
    # -------------------------- Function Definitions -------------------------
    def getrcx(self):
        # input rcx as an n-by-2 array integer values for their indices. 
        all_files = glob.glob(self.channel + '*.dat')
        all_files.sort()
        
        m = len(all_files)
        n = len(self.receiver_locations[:,1])
        self.timeseries = np.zeros([m,n])
        i = 0
        if self.dim == '2':
        	self.ny = None
        
        if self.dim == '2.5':
            for fn in all_files:
                npdat = self.read_dat(fn)
                for j in range(0, n):
                    # Don't forget x is columns and z is rows
                    self.timeseries[i,j] = npdat[
                        int(self.receiver_locations[j,2]),
                        int(self.receiver_locations[j,1]),
                        int(self.receiver_locations[j,0])
					]
                    
                i = i + 1
        else:
            for fn in all_files:
                npdat = self.read_dat(fn)
                for j in range(0, n):
                    # Don't forget x is columns and z is rows
                    self.timeseries[i,j] = npdat[
                        int(self.receiver_locations[j,2]),
                        int(self.receiver_locations[j,0])
                    ]
                    
                i = i + 1
        
        # Save the array as csv for other types of processing
        np.savetxt("receiver_array.csv", self.timeseries, delimiter = ",")
    
    # ---------------------------------------
    def tsplot(self):
        m, n = self.timeseries.shape
        # if the gain is 0, then the window is 1
        if self.gain == 0:
            self.gain = int(1)
        # The gain can't exceed the length of the time series
        if self.gain > m:
            self.gain = m 

        for ind in range(0, n):
            self.timeseries[:,ind] = agc(
                self.timeseries[:,ind], self.gain, "mean"
            )
        
        # Create the values for the y-axis
        timelocs = np.arange(0,m, int(m/10) ) #10 tick marks
        timevals = timelocs*self.dt
        
        # Create the reciever location
        xlocs = np.arange(0, n, int(n/5) ) #5 tick marks 
        
        # Create the figure
        fig = plt.figure(figsize =(n/2,m/2) )
        ax1 = plt.gca()
        # ax2 = ax1.twinx() 
        ax1.imshow(self.timeseries, cmap = 'Greys', aspect = 'auto')
        ax1.set_xlabel(r'Receiver #')
        ax1.xaxis.tick_top()
        ax1.xaxis.set_label_position('top')
        ax1.set_xticks(xlocs)
        ax1.set_xticklabels(xlocs)
        ax1.set_ylabel(r'Two-way Travel Time (s)')
        ax1.set_yticks(timelocs)
        ax1.set_yticklabels(timevals)
        
        # ax2.set_ylabel('Depth (m)')
        # ax2.set_yticks(timelocs) 
        # ax2.set_yticklabels(twt)
        
        # # Label the x-axis
        # ax.xaxis.tick_top()
        # ax.xaxis.set_label_position('top')
        # ax.set_xlabel(r"Receiver #")
        
        # Label the y-axis
        if self.channel == 'Vx' or self.channel == 'Vz':
            mult = 1e2
        else:
            mult = 1e6
        
        time_locations = np.linspace(1, m, 10)
        time_labels = np.round( time_locations*self.dt*mult, 4)
        # ax.set_ylabel(r"Two way travel time (s)")
        plt.yticks(ticks=time_locations, labels = time_labels.astype(str) )
        
        # # Other figure operations
        # if self.channel == 'Vx' or self.channel == 'Vz':
        #     plt.figtext(0.30, 0.07, 'x $10^{-2}$')	
        # else:
        #     plt.figtext(0.30, 0.07, 'x $10^{-6}$')
        
        ax1.set_aspect(aspect = exaggeration	)
        plt.show()
    
    def plot_layout(self):
        figmod, axmod = plt.subplots()
        img = mpimg.imread(self.modelfile)
        axes_extent = [0, img.shape[2],img.shape[0], 0 ]
        
        implot = axmod.imshow(img, extent = extent)
        
        # add white upside down triangle receivers
        axmod.scatter(
            self.receiver_locations[:,0]- self.cpml,
            self.receiver_locations[:,1] - self.cpml,
            marker = 'v', s = 30, c = (0.8, 0.8, 0.8, 1),
            linewidths = 0.5, edgecolor = (0.2, 0.2, 0.2, 1 )
        )
        
        # add white star source
        axmod.scatter(
            self.source_location[0], self.source_location[1],
            marker = '*', s = 30, c = (0.8, 0.8, 0.8, 1 ),
            linewidths = 1, edgecolor = (0.2, 0.2, 0.2, 1 ) 
        )
        
        axmod.set_ylabel("z-indice")
        axmod.set_xlabel("x-indice")
        plt.show()
    
    # ------------- 
    def read_dat(self, fn):
        if self.ny:
            if self.channel == 'Ex':
                NX = self.nz
                NY = ny
                NZ = nx-1
            elif self.channel == 'Ey':
                NX = self.nz
                NY = self.ny-1
                NZ = self.nx
            elif self.channel == 'Ez':
                NX = self.nz-1
                NY = self.ny
                NZ = self.nx
            else:
                NX = self.nz
                NY = self.ny 
                NZ = self.nx
            
        f = FortranFile(fn, 'r')
        dat = f.read_reals(dtype = 'float32')
        
        if self.ny:
            dat = dat.reshape(NX, NY, NZ)
        else:
            dat = dat.reshape(self.nx, self.nz)
        
        f.close()
        return(dat)

# ==================== Create the object and assign inputs ====================
array = Array()

# Fill in object fields 

# Get the arguments
args = parser.parse_args()
project_file = ''.join(args.prjfile)
receiver_file = ''.join(args.rcxfile)
array.channel = ''.join(args.channel)
rind = args.index[0] == 0

# Optional inputs
array.gain = args.gain[0]
layout = args.layout[0] == 1
showplots = args.suppress_plotting[0] == 1
exaggeration = args.exaggeration[0]

# --------------------------- Query the project file --------------------------
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

		if temp[1] == 'dim':
			array.dim = str( temp[2].rsplit()[0] )
		if temp[1] == 'nx':
			array.nx = int( temp[2].rsplit()[0] )
		if temp[1] == 'ny':
			array.ny = int( temp[2].rsplit()[0] )
		if temp[1] == 'nz':
			array.nz = int( temp[2].rsplit()[0] )
		if temp[1] == 'dx':
			array.dx = float(temp[2].rsplit()[0] )
		if temp[1] == 'dy':
			array.dy = float(temp[2].rsplit()[0] )
		if temp[1] == 'dz':
			array.dz = float(temp[2].rsplit()[0] )
		if temp[1] == 'cpml':
			array.cpml = int( temp[2].rsplit()[0])

	# Source
	
	if array.channel == 'Ex' or array.channel == 'Ey' or array.channel == 'Ez':
		if line[0] == 'E':
			temp = line.split(',')
			if temp[1] == 'dt':
				array.dt = float(temp[2].rsplit()[0])
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
				array.dt = float(temp[2].rsplit()[0])
			if temp[1] == 'x':
				sx = float(temp[2].rsplit()[0])
			if temp[1] == 'y':
				sy = float(temp[2].rsplit()[0])
			if temp[1] == 'z':
				sz = float(temp[2].rsplit()[0])
				

f.close()

# ------------------------------ Build the object -----------------------------
xyz = np.genfromtxt(
    receiver_file, 
    delimiter = ',', 
    names = True, 
    dtype = float
)

# We need to make sure the recievers are ordered correctly and the absorbing 
# boundary is corrected for
# First check to see if the inputs are indices or 
if rind:
    array.receiver_locations = np.vstack(
        [
            xyz['X']/array.dx + array.cpml,
            xyz['Y']/array.dy + array.cpml,
            xyz['Z']/array.dz + array.cpml
        ]
    ).T
else:
    array.receiver_locations = np.vstack(
        [
            xyz['X'] + array.cpml,
            xyz['Y'] + array.cpml,
            xyz['Z'] + array.cpml
        ]
    ).T

# Adjust the object fields relative to the cpml
array.nx = array.nx + 2*array.cpml
array.ny = array.ny + 2*array.cpml
array.nz = array.nz + 2*array.cpml

# extent = (array.cpml, array.nx-array.cpml, array.nz-array.cpml, array.cpml)
if array.channel == 'Ex' or array.channel == 'Ey' or array.channel == 'Ez':
	array.source_location = np.array(
		[ 
			ex/array.dx + array.cpml,
			ey/array.dy + array.cpml,
			ez/array.dz + array.cpml
		]
	)
else:
	array.source_location = np.array(
		[
			sx/array.dx + array.cpml,
			sy/array.dy + array.cpml,
			sz/array.dz + array.cpml
		]
	) 

# Get the timeseries for each receiver
array.getrcx()

if showplots:
	quit()

if layout:
	array.plot_layout()
else:
	array.tsplot()