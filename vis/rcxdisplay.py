#!/usr/bin/env python3

# Plot the common offset survey
#
# Available backends are:

# ['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg',
# 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg',
# 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']

# Set matplotlib backend
import matplotlib
# Set matplotlib backend 
import matplotlib 
try:
    matplotlib.use('Qt4Agg')
except Exception as e:
    matplotlib.use('TkAgg')
    print('Qt4 import error. Now using Tk instead of Qt. Details:')
    print(e)

import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from class_definitions import *

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(
    description="""Plots a csv file of timeseries data (columns) for targeted receivers """ )

parser.add_argument(
    '-p', '--prjfile',
    nargs = 1, type = str, required = True,
    help = """Path to the project file."""
)

parser.add_argument(
    '-f', '--file',
    nargs=1, type=str, required=True,
    help="""Path to the csv file with receiver timeseries data"""
)

parser.add_argument(
    '-g', '--gain',
    nargs = 1, type = float, required = False, default = 100,
    help = "The smoothing length"
)

parser.add_argument(
    '-e', '--exaggeration',
    nargs=1, type = float, required = False, default = [0.5],
	  help = """Set the aspect ratio between the x and y axes for
    plotting. Default is 0.5"""
)

#parser.add_argument(
#    '-r', '--receiver',
#    nargs = 1, type = str, required = True,
#    help = """The csv file path that contains the receiver locations"""
#)
parser.add_argument(
    '-s', '--seismic',
    nargs = 1, type = int, required = False, default = [1],
    help = """Flag whether this is a seismic or electromagnetic model"""
)

#parser.add_argument(
#    '-t', '--surveytype',
#    nargs = 1, type = str, required = True,
#    help = """The type of survey that you are plotting. The available options
#    are: war, co, cmp. The correspond to wide angle, common offset, and common
#    midpoint, respectively."""
#)


args = parser.parse_args()
project_file = ''.join(args.prjfile)
cofile = ''.join(args.file)
gain = args.gain[0]
exaggeration = args.exaggeration[0]
seismic = args.seismic[0] == 1
survey_type = ''.join(args.surveytype)


# project_file = 'easy_greenland.prj'
# cofile = 'receiver_array.csv'
# gain = 101
# exaggeration = 0.25
# seismic = True
# -----------------------------------------------------------------------------
# Load the values from the project file
domain, material, seismic, electromag = loadproject(
    project_file,
    Domain(),
    Material(),
    Model(),
    Model()
)

# Open the csv containing time series values
dat = np.genfromtxt(cofile, delimiter = ',')
m,n = dat.shape

if seismic:
    mult = 1e2
else:
    mult = 1e6

timelocs = np.arange(0, m, int(m/10) ) # 10 tick marks along y-axis
rcxlocs = np.arange(0, n, int(n/5) ) # 5 tick marks along x-axis 


if seismic:
    timevals = np.round(timelocs*float(seismic.dt[0]) * mult, 2)
else:
    timevals = np.round(timelocs*float(electromag.dt[0]) * mult, 2)

if gain == 0:
    gain = 1
elif gain < m:
    for j in range(0, n):
        # Subtract the mean value
        dat[:,j] = dat[:,j] - np.mean(dat[:,j])
        dat[:,j] = agc(dat[:,j], gain, "mean")
else:
    gain = m


fig = plt.figure(figsize =(n/2,m/2) )
ax1 = plt.gca()
# ax2 = ax1.twinx() 

ax1.imshow(dat, cmap = 'Greys', aspect = 'auto')
ax1.set_xlabel(r'Receiver #')
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top')
ax1.set_xticks(rcxlocs)
ax1.set_xticklabels(rcxlocs)
ax1.set_ylabel(r'Two-way Travel Time (s)')
ax1.set_yticks(timelocs)
ax1.set_yticklabels(timevals)

# Other figure handle operations
ax1.set_aspect(aspect = exaggeration)

if seismic:
	ax1.text(0, m + 0.03*m, 'x $10^{-2}$')	
else:
	ax1.text(0, m + 0.03*m, 'x $10^{-6}$')

ax1.update_datalim( ((0,0),(m, n)))
plt.show()

