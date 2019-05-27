#!/usr/bin/python3

# Plot the common offset survey

import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""CODISPLAY reads the common 
	survey file created from common_offset.sh. The """ )

# parser.add_argument( 'project_file', nargs=1, type=str,
# 	help=""" The project file used to model the common offset survey """)

parser.add_argument( 'survey_file', nargs=1, type=str, required = True,
						help='the survey .csv file', 
						default=None)

parser.add_argument('-g', '--gain', nargs = 1, type = float, required = False,
	help = "The exponential value for 2^m of the gain function (default=None)", default = None)

parser.add_argument('-e', '--exaggeration', nargs=1, type = float, required = False,
	help = """Set the aspect ratio between the x and y axes for 
	plotting. Default is 0.5""", default = [0.5])

# ----------------- Delete after debug -------------------
# parser.add_argument('-m', '--model_type', nargs=1, type= str, required = False,
# 	help = "Specify the type of model; s - seismic; e - electromag (default))")


# parser.add_argument( '-d', '--delta', nargs=1, type=float, required = False,
# 	help="""The change in source distance in meters along the profile. The 
# 	change in reciever distance is the same.""", default=[1])

# parser.add_argument( '-o', '--offset', nargs =1, type = float, required = False,
# 	help=""" The initial offset of the source and the reciever from the 
# 	midpoint""", default = [5])
# --------------------------------------------------------


args = parser.parse_args()
cofile = ''.join(args.survey_file)
gain = args.gain
exaggeration = args.exaggeration[0]

if gain:
	gain = gain[0]


# The filename contains information about the survey
survey_info = cofile.split('.')

# The file was saved as <basename>.<delta>.<offset>.<channel>.<co/cmp>.csv but 
# basename could be <yada_yada>.<blah_blah> 

project_file = survey_info[0:-5] + '.prj'
ds = survey_info[-5]
offset = survey_info[-4]

if survey_info[-3] == 'Vx' or survey_info[-3] == 'Vz':
	model = 's'

survey_type = survey_info[-2]

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


dat = np.genfromtxt(cofile, delimiter = ',')
m,n = dat.shape

time_locations = np.linspace(1, m, 10) 
time_labels = np.round( time_locations*dt*1e6, 4)


if survey_type == 'cmp':
	dist_locations = np.round( np.linspace(1, n, 7) )
	dist_labels = 2*dist_locations*ds + offset*2
	dist_labels = dist_labels.astype(int)
else:
	dist_location = np.round(0, n-1, 7) * ds 
	dist_labels = dist_locations.astype(int)



fig, ax = plt.subplots()

if gain:
	gain_function = np.zeros([m,n])
	for j in range(0, m):
		gain_function[j,:] = np.exp(-j/(2**gain))
	im = ax.imshow(dat/gain_function, cmap = 'Greys' )
else:
	im = ax.imshow(dat, cmap = 'Greys')

ax.set_xlabel(r"Source-Reciever Distance (m)")
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

ax.set_ylabel(r"Two way travel time (s)")

ax.set_aspect(aspect = exaggeration)

# Label the y axis
plt.yticks(ticks=time_locations, labels = time_labels.astype(str) )
plt.figtext(0.30, 0.07, 'x $10^{-6}$')

# Label the x axis
plt.xticks(ticks = dist_locations, labels = dist_labels.astype(str) )

plt.show()

