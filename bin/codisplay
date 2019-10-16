#!/usr/bin/env python3

# Plot the common offset survey

import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""CODISPLAY reads the common 
	offset or midpoint survey file created from common_offset.sh. The """ )

parser.add_argument( 'meta_file', nargs=1, type=str,
	help='File path that contains the survey metadata.', )

parser.add_argument( '-c', '--channel', nargs=1, type=str, required=True,
	help="""The channel that you would like to plot. Valid options 
	are Ex, Ez, Vx, Vz.""")

parser.add_argument('-g', '--gain', nargs = 1, type = float, required = False,
	help = "The exponential value for 2^m of the gain function (default=None)", default = None)

parser.add_argument('-e', '--exaggeration', nargs=1, type = float, required = False,
	help = """Set the aspect ratio between the x and y axes for 
	plotting. Default is 0.5""", default = [0.5])

# ----------------- Delete after debug -------------------


args = parser.parse_args()
meta_file = ''.join(args.meta_file)
gain = args.gain
exaggeration = args.exaggeration[0]
channel = ''.join(args.channel)

if gain:
	gain = gain[0]

# --------------------- Get the values from the meta file ---------------------

f = open(meta_file)

# If running the wrapper functions, the meta file will save the same each time,
# but for whatever reason we'll assume that this isn't the case
for line in f:

	 temp = line.rstrip().rsplit()

	 if temp[0] == 'offset:':
	 	offset = float(temp[1])
	 if temp[0] == 'delta:':
	 	ds = float(temp[1])
	 if temp[0] == 'survey_type:':
	 	survey_type = temp[1]
	 	
f.close()

cofile = '.'.join( meta_file.split('.')[:-2] )+ '.' + channel + '.' + survey_type + '.csv'
project_file = '.'.join( meta_file.split('.')[:-3] )+ '.prj'

# -----------------------------------------------------------------------------
# Get the values we need
f = open(project_file)

for line in f:
	if channel == 'Vx' or channel == 'Vz':
		# Source
		if line[0] == 'S':
			temp = line.split(',')
			if temp[1] == 'dt':
				dt = float(temp[2].rsplit()[0])
	else:
		if line[0] == 'E':
			temp = line.split(',')
			if temp[1] == 'dt':
				dt = float(temp[2].rsplit()[0])

f.close()


dat = np.genfromtxt(cofile, delimiter = ',')
m,n = dat.shape

if channel == 'Vx' or channel == 'Vz':
	mult = 1e2
else:
	mult = 1e6

time_locations = np.linspace(1, m, 10) 
time_labels = np.round( time_locations*dt*mult, 4)

if survey_type == 'cmp':
	dist_locations = np.round( np.linspace(1, n, 7) )
	dist_labels = 2*dist_locations*ds + offset*2
	dist_labels = dist_labels.astype(int)
else:
	dist_locations = np.round(np.linspace(0, n-1, 7) )
	dist_labels = dist_locations*ds
	dist_labels = dist_labels.astype(int)

# Create the figure object using subplots
fig, ax = plt.subplots()

if gain:
	gain_function = np.zeros([m,n])
	for j in range(0, m):
		gain_function[j,:] = np.exp(-j/(2**gain))
	im = ax.imshow(dat/gain_function, cmap = 'Greys' )
else:
	im = ax.imshow(dat, cmap = 'Greys')

# Label the x axis
plt.xticks(dist_locations, dist_labels.astype(str) )
ax.set_xlabel(r"Source-Reciever Distance (m)")
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

# Label the y axis
ax.set_ylabel(r"Two way travel time (s)")
plt.yticks(time_locations, time_labels.astype(str) )

# Other figure handle operations
ax.set_aspect(aspect = exaggeration)

if channel == 'Vx' or channel == 'Vz':
	plt.figtext(0.30, 0.07, 'x $10^{-2}$')	
else:
	plt.figtext(0.30, 0.07, 'x $10^{-6}$')


plt.show()

