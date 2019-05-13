#!/usr/bin/env python3

# Plot the common offset survey

import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""CODISPLAY reads the common 
	survey file created from common_offset.sh. The """ )

parser.add_argument( 'project_file', nargs=1, type=str,
	help=""" The project file used to model the common offset survey """)

parser.add_argument( '-s', '--survey_file', nargs=1, type=str, required = True,
						help='the survey .csv file', 
						default=None)
parser.add_argument( '-d', '--delta', nargs=1, type=str, required = False,
	help="""The change in source distance in meters along the profile. The 
	change in reciever distance is the same.""", default=[1])

parser.add_argument('-g', '--gain', nargs = 1, type = float, required = False,
	help = "The exponential value for 2^m of the gain function (default=None)", default = None)

parser.add_argument('-m', '--model_type', nargs=1, type= str, required = False,
	help = "Specify the type of model; s - seismic; e - electromag (default))")

args = parser.parse_args()
project_file = ''.join(args.project_file)
cofile = ''.join(args.survey_file)
ds = args.delta[0]
gain = args.gain
model=args.model_type

if gain:
	gain = gain[0]

# -----------------------------------------------------------------------------
# Get the values we need
f = open(project_file)

for line in f:
	if model== 's':
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


dat = np.genfromtxt(cofile, delimiter = ' ')
m,n = dat.shape

fig, ax = plt.subplots()

if gain:
	gain_function = np.zeros([m,n])
	for j in range(0, m):
		gain_function[j,:] = np.exp(-j/(2**gain))
	im = ax.imshow(dat/gain_function, cmap = 'Greys' )
else:
	im = ax.imshow(dat, cmap = 'Greys')

ax.set_xlabel(r"Normalized Distance (m)")
ax.set_ylabel(r"Normalized Time")
ax.set_aspect(aspect = 0.01)
plt.show()

