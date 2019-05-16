#!/usr/bin/env python3

# Wrapper to generate the euler angles for the plunge and trend 
# then plot the results

import numpy as np
import argparse
import matplotlib.pyplot as plt

import orientsynth as ot

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description=""" """ )

parser.add_argument( 'output_file', nargs=1, type=str, required = True
						help='the full file path for the project file', 
						default=None)

parser.add_argument( '-n', '--npts', nargs = 1, type = int, required = False,
	help = """ Total number of grains in synthetic sample.""", default = [100])

parser.add_argument( '-p', '--plunge', nargs = 3, type = float, 
	required = True, help = """ Plunge angle in degrees.""")

parser.add_argument( '-t', '--trend', nargs =1, type = float, required = True,
	help = """ Trend angle in degrees.""")

parser.add_argument( '-a', '--anglemin', nargs = 1, type = float, 
	required = True, help = """ Minimum angle deviation.""")

parser.add_argument( '-A', '--anglemax', nargs = 1, type = float, 
	required = True, help = """Maximum angle deviation. """)

parser.add_argument('-S', '--suppress_plotting', nargs = 1, type = int, 
	help = """Suppress all plotting (1/0)""", required = False, default = [0])

# Get the arguments
args = parser.parse_args()
npts = args.npts[0]
plunge=args.plunge[0]
trend=args.trend[0]
amin=args.anglemin[0]
amax=args.anglemax[0]

euler_list= np.array([])
end_orientation = np.array([])


euler_list, end_orientation = ot.fabric_ortosynth(trend, plunge, amin, amax, npts, &
	euler_list, end_orientation)