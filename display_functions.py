#!/usr/bin/env python3

# This is the suite of display functions for the seismicfdtd2d and emfdtd2d 
# outputs. The double precision data is saved as a Fortran unformatted binary 
# but fortunately numpy and scipy make this an easy process

import glob 
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg

# ======================= Create the class variables ==========================

# The timeseries plot
class ts:
	def __init__(self):
		super().__init__()
		self.build()

	def build(self):
		# Initialize variables
		self.dt = None
		self.location = None # Must be m-by-2 array in indices [x_ind, y_ind]
		self.timeseries = None
		self.t = None
		self.channel = None
		self.color = ['red', 'blue'] # The default positive and negative fill
		self.nx = None
		self.nz = None

	# -------------------------- Function Definitions -------------------------



	def getrcx(self):
	# input rcx as an n-by-2 array integer values for their indices. 

		all_files = glob.glob(self.channel + '*.dat')
		m = len(all_files)
		n = len(self.location[:,1])

		self.timeseries = np.zeros([m,n])
		for fn in all_files:
			npdat = read_dat(fn, self.nx, self.ny)

			i = 0
			for j in range(0, n):
				self.timeseries[i,j] = npdat[ self.location[j,0], self.location[j,1] ]

		

	def tsplot(self):
		
		fig, ax = plt.subplots()

		m,n = self.timeseries.shape 
		self.t = np.linspace(1, m, 1)*self.dt 

		if vertical:
			for offset in range(0, n):
				y = self.timeseries[:,offset] + offset
				ax.plot(y, self.t, 'k-')
				ax.fill_between(self.t, offset, y, where = (y > offset), color = 'k')
		else:

			for offset in range(0, n):
				y = self.timeseries[:,offset] + offset

				ax.plot(self.t, y, 'k-')
				ax.fill_between(self.t, offset, y, where = (y > offset), color = 'k')

		plt.show()



# =========================== General Use Functions ===========================
def read_rcx(rcxfile):
	# Import the reciever locations from a text list
	pass

def modelplt(imfile, rcx=False):
	# Formatting for plotting the image
	
	# Load the image

	if rcx:
		pass
	else:
		pass


