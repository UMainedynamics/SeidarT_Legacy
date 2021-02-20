#!/usr/bin/env python3

# From the set of image outputs, we can build a gif. The images will be in csv
# or fortran unformatted binary


import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.io import FortranFile
# from definitions import *
from imgen import *
import os
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
    '-n', '--num_steps', 
    nargs = 1, type = int, required = True, 
    help = """The time step interval between the images that
	are going to be used. Every time step is written to file which means that
	we can take any equally spaced images to create the gif with an
	appropriate resolution, time to compute, and file size. For example,
	n=20 means that every 20 images will be used thus significantly reducing
	how long it takes to compile the gif."""
)

parser.add_argument(
    '-c', '--channel', 
    nargs = 1, type = str, required = True,
    help = """Specify whether a particular channel is going to be used. The
	available channels are Ex, Ez, Vx, and Vz for the electric field and
	seismic velocities, respectively."""
)

parser.add_argument(
    '-d', '--delay', 
    nargs = 1, type = int, required = False, default = [1], 
    help = """The amount of delay between two frames"""
)

parser.add_argument(
    '-a', '--alpha',
    nargs = 1, type = float, required = False, default = [0.3],
    help = """(OPTIONAL FLOAT [0,1]) Change the transparency of the model 
    plotted in the background; default = 0.3. Zeros is full transparency, 1 is 
    CIA transparency."""
)

# Get the arguments
args = parser.parse_args()
prjfile = ''.join(args.prjfile)
channel = ''.join(args.channel)
delay = str(args.delay[0])
num_steps = args.num_steps[0]
alpha = min([1, args.alpha[0]] ) 


# =============================================================================

# Check if the .dat files are still around
files = glob.glob(channel + '*.dat')
files.sort()

# We'll start counting with the first frame
n=num_steps

for fn in files:
    if n == num_steps:
        mag = FDTDImage(prjfile, fn)
        mag.getprjvals()
        mag.magnitudeplot(alpha = alpha)
        mag.addlabels()
        mag.plotfile = 'magnitude.' + fn[:-3] + '.png'
        plt.savefig(mag.plotfile)
        plt.close()
        n = 1
    else:
        n = n + 1


print('Creating the GIF')
# Use imagemagick via shell command to create the gif
shellcommand = 'convert -delay ' + \
    delay + ' -loop 0 magnitude.' + channel + '*.png ' + \
        channel + '.gif'
call(shellcommand, shell = True)

# Remove the png files 
for filepath in glob.glob('magnitude.' + channel + '*.png'):
    os.remove(filepath)
        
