#!/usr/bin/env python3
#
# Convert data to SEG-Y format
# Either a shot gather with an array of receivers
# or a receiver gather with multiple transmits (can be common offset of cmp)
#
# Benjamin Hills
# University of Washington
# bhills@uw.edu
# 9/23/19
#
# -----------------------------------------------------------------------------

import segyio
import numpy as np
import argparse

# -------------------------- Command-Line Parser ---------------------------

parser = argparse.ArgumentParser(description="""Write data to SEG-Y format. Data are in the form of traces from
                    multiple receivers (shot-gather) or from one receiver and multiple sources (receiver-gather). Multiple
                    transmits can be for common offset, cmp, etc.""" )

parser.add_argument( '-d', '--data_file', nargs=1, type=str, required = True,
		    help="""Path to the file with data to be converted""", default=None)

parser.add_argument( '-p', '--project_file', nargs=1, type=str, required = True,
                    help = """Path to project file from which the time step will be taken""", default=None)

parser.add_argument( '-o', '--segy_file', nargs=1, type=str, required = False,
                    help = """Path to location where the segy file should be saved. If left as None, the data
                    file name will be used with extension .sgy.""", default=None)

parser.add_argument( '-s', '--seismic', nargs=1, type=bool, required = False,
                    help = """Boolean for radar data (default) or seismic. This will control which time step to take from
                    the .prj file.""", default=False)

parser.add_argument( '-m', '--meta_file', nargs=1, type = str, required = False,
	            help = """Path to metadata file.""", default = None)

# Get the arguments
args = parser.parse_args()
data_file = ''.join(args.data_file)
project_file = ''.join(args.project_file)
if args.segy_file is None:
    segy_file = data_file[:-4]+'.sgy'
else:
    segy_file = ''.join(args.segy_file)
seismic = args.seismic
if args.meta_file is not None:
    meta_file = ''.join(args.meta_file)

# -------------------------- Do the Conversion With segyio ---------------------------

# Get the data from a .csv file
data = np.transpose(np.genfromtxt(data_file,delimiter=','))
# Get the time step from the project file
# can be seismic or electromagnetic and would have a different time step
with open(project_file,'r') as fid:
    prj_contents = fid.read()
if seismic:
    dt_start = prj_contents.find("S,dt,") + 5
    dt_factor = 1e6
else:
    dt_start = prj_contents.find("E,dt,") + 5
    dt_factor = 1e12
dt_end = prj_contents[dt_start:].find("\n") + dt_start
dt = prj_contents[dt_start:dt_end]
# dt needs to be an integer for the segy file, so multiply by a large number
dt = int(float(dt)*dt_factor)

# save as segy
segyio.tools.from_array2D(segy_file,data,dt=dt)

