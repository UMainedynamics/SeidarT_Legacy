#!/usr/bin/env python3

import glob
import argparse
import matplotlib.pyplot as plt
from imgen import *
from subprocess import call
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


# Get the arguments
args = parser.parse_args()
prjfile = ''.join(args.prjfile)
num_steps = args.num_steps[0]

# This is setup for seismic
files = glob.glob('Vx*.dat')
files.sort()

print('Creating PNG snapshots')
# We want to do the first file; initial condition
n=num_steps
for fn in files:
    if n == num_steps:
        vi = FDTDImage(prjfile, fn)
        vi.getprjvals()
        vi.quiverplot(papercolumnwidth=10)
        vi.addlabels()
        plt.savefig(vi.plotfile)
        plt.close()
        n = 1
    else:
        n = n + 1


# Create the gif using Imagemagick
print('Creating the GIF')
call('convert -delay 20 -loop 0 vector.*.png vector.gif', shell = True)
    
