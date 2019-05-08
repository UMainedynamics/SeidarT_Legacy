#!/usr/bin/env

# Some io functions that will be useful

import glob 
import numpy as np
import scipy.io as scio
from scipy.io import FortranFile


# =============================================================================
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
