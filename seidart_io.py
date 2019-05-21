#!/usr/bin/python3

# Some basic I/O operations for the SeidarT programs to convert the unformatted
# Fortran 

import glob 
import numpy as np
import scipy.io as scio
from scipy.io import FortranFile
import argparse
from subprocess import call


# ========================== Command Line Arguments ==========================
parser = argparse.ArgumentParser(description="""This contains a set of 
	functions to read the all of the .dat outputs of the SeisarT modeling 
	routines.""" )

parser.add_argument( '-o', '--output', nargs=1, type=str, required=True,
					help="""The data type to convert to; csv (default), 
					mat. """, default = 'csv')

parser.add_argument( '-f', '--project_file', nargs=1, type=str, required=True, 
					help="""The project file associated with the .dat files""", 
					default=None)

parser.add_argument( '-d', '--delete', nargs=1, type=int, required=False,
					help = """Specify whether to delete (1) the .dat files 
					following the routines""", default = 0)

# Get the arguments
args = parser.parse_args()
project_file = ''.join(args.project_file)
fmt = ''.join(args.output)
rm = args.delete[0]

# =========================== Function Definitions ============================

def dat2csv(nx, ny):

	for fn in glob.glob('*.dat'):
		npdat = read_dat(fn, nx, ny)
		# replace the .dat extension with .csv
		sfn = fn[0:-3] + 'csv'
		np.savetxt(sfn, npdat, delimiter=",")

# -----------------------------------------------------------------------------

def dat2mat(nx, ny):
	# similar to above but saving to matlab specific file
	for fn in glob.glob(chan + '*.dat'):
		npdat = read_dat(fn, nx, ny)
		sfn = fn[0:-3] + 'mat'
		scio.savemat(sfn, npdat)


# -----------------------------------------------------------------------------

def read_dat(fn, nx, ny):
	f = FortranFile(fn, 'r')
	dat = f.read_reals(dtype = 'float64')
	dat = dat.reshape(ny, nx)
	f.close()

	return(dat)

# ============================== Run the program ==============================

# Read the project file and get the nx and ny values
f = open(project_file)

for line in f:

	# All domain inputs must be input except for ny and dy
	if line[0] == 'D':
		temp = line.split(',')

		if temp[1] == 'nx':
			nx = int( temp[2].rsplit()[0] )
		if temp[1] == 'nz':
			nz = int( temp[2].rsplit()[0] )
		if temp[1] == 'cpml':
			cpml = int( temp[2].rsplit()[0])

f.close()

nx = nx + 2*cpml
nz = nz + 2*cpml

if fmt == 'mat':
	print('Converting .dat files to .mat. This may take a moment.\n')
	dat2mat(nx, nz)
else:
	print('Converting .dat files to .csv. This may take a moment.\n')
	dat2csv(nx, nz)


if rm == 1:
	print('Deleting all .dat files.\n')
	call('rm *.dat', shell=True)
