#!/usr/bin/env python3

"""
Plot a snapshot of the wavefield with axes
"""

import numpy as np 
import matplotlib.pyplot as plt
import argparse
from definitions import *
import matplotlib.image as mpimg
from imgen import *


# =============================================================================
parser = argparse.ArgumentParser(
    description="""Plot a snapshop of the vector wavefield in 2D"""
)

parser.add_argument(
    '-p', '--prjfile',
    nargs = 1, type = str, required = True,
    help = 'The full file path to the project file'
)

parser.add_argument(
    '-v', '--velocity',
    nargs = 1, type = str, required = True,
    help = """The .dat file that corresponds to the velocity in either the 
    x-direction or z-direction (e.g. Vx000400.dat). The corresponding 
    orthogonal velocity file will be loaded as well."""
) 

args = parser.parse_args()
prjfile = ''.join(args.prjfile)
velocityfile = ''.join(args.velocity)


# ============================ Create the objects =============================

mag = FDTDImage(prjfile, velocityfile)
mag.getprjvals()
mag.magnitudeplot()
mag.addlabels()
mag.plotfile = velocityfile[:-3] + 'plot.png'
plt.savefig(mag.plotfile)
plt.close()
