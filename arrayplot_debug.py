#!/usr/bin/env python3

# Create a reciever array and plot the timeseries

import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.io import FortranFile

import io_functions as iof

# ------------
project_file = 'negis1.prj'
channel = 'Vx'

initial_coords = np.array([70, 0, 70])
final_coords = np.array([370, 0, 70])
delta = 2


# ---
# ---
# ======================= Create the class variables ==========================

# The timeseries plot
class Array:
	def __init__(self):
		super().__init__()
		self.build()

	def build(self):
		# Initialize variables
		self.dt = None
		self.dx = None
		self.dy = None
		self.dz = None
		self.source_location = None # Must be m-by-2 array in indices [x_ind, y_ind]
		self.initial_position = None
		self.final_position = None
		self.reciever_locations = None
		self.dr = None
		self.timeseries = None
		self.t = None
		self.channel = None
		self.color = ['red', 'blue'] # The default positive and negative fill
		self.nx = None
		self.ny = None
		self.nz = None
		self.vertical = False

	# -------------------------- Function Definitions -------------------------
	def getrcx(self):
	# input rcx as an n-by-2 array integer values for their indices. 

		all_files = glob.glob(self.channel + '*.dat')
		all_files.sort()

		m = len(all_files)
		n = len(self.reciever_locations[:,1])

		self.timeseries = np.zeros([m,n])

		i = 0
		for fn in all_files:
			npdat = iof.read_dat(fn, self.nx, self.nz)

			for j in range(0, n):
				# Don't forget x is columns and z is rows
				self.timeseries[i,j] = npdat[ self.reciever_locations[j,1], self.reciever_locations[j,0] ]	
			
			i = i + 1

		# Save the array as csv for other types of processing
		np.savetxt("reciever_array.csv", self.timeseries, delimiter = ",")

	# -------------------------------------------------------------------------
	def tsplot(self):
		
		fig, ax = plt.subplots()

		m,n = self.timeseries.shape 
		self.t = np.linspace(1, m, m)*self.dt 

		im = ax.imshow(self.timeseries, cmap = 'Greys')

		# if self.vertical:
		# 	for offset in range(0, n):
		# 		y = self.timeseries[:,offset] + offset
		# 		ax.plot(y, self.t, 'k-')
		# 		ax.fill_between(self.t, offset, y, where = (y > offset), color = 'k')
		# else:

		# 	for offset in range(0, n):
		# 		y = self.timeseries[:,offset] + offset

		# 		ax.plot(self.t, y, 'k-')
		# 		ax.fill_between(self.t, offset, y, where = (y > offset), color = 'k')

		plt.show()


	def arraybuild(self):
		
		# Calculate the length
		rcx_len = np.sqrt( sum( (self.final_position - self.initial_position)**2) )
		nrcx = np.floor(rcx_len/self.dr)


		xz = np.zeros([nrcx.astype(int),2]) 
		xz[:,0] = np.linspace(self.initial_position[0], 
			self.final_position[0], nrcx)/self.dx
		xz[:,1] = np.linspace(self.initial_position[2], 
			self.final_position[2], nrcx)/self.dz

		self.reciever_locations = xz.astype(int) + self.cpml




# --------------------------- Query the project file --------------------------
array = Array()

# Fill in object fields 
array.channel = channel

#!!!!!! Add the 2.5D components later

# Get the values we need
f = open(project_file)

for line in f:

	# Get the image file
	if line[0] == 'I':
		# There's a trailing new line value
		imfile = line[2:-1]

	# All domain inputs must be input except for nz and dy
	if line[0] == 'D':
		temp = line.split(',')

		if temp[1] == 'nx':
			array.nx = int( temp[2].rsplit()[0] )
		if temp[1] == 'nz':
			array.nz = int( temp[2].rsplit()[0] )
		if temp[1] == 'dx':
			array.dx = float(temp[2].rsplit()[0] )
		if temp[1] == 'dz':
			array.dz = float(temp[2].rsplit()[0] )
		if temp[1] == 'cpml':
			array.cpml = int( temp[2].rsplit()[0])
		if temp[1] == 'write':
			array.write = int(temp[2].rsplit()[0])

	# Source
	if line[0] == 'S':
		temp = line.split(',')
		if temp[1] == 'dt':
			sdt = float(temp[2].rsplit()[0])
		if temp[1] == 'x':
			sx = float(temp[2].rsplit()[0])
		if temp[1] == 'z':
			sz = float(temp[2].rsplit()[0])

	if line[0] == 'E':
		temp = line.split(',')
		if temp[1] == 'dt':
			edt = float(temp[2].rsplit()[0])
		if temp[1] == 'x':
			ex = float(temp[2].rsplit()[0])
		if temp[1] == 'z':
			ez = float(temp[2].rsplit()[0])


f.close()

# ----------------------------- Create the object -----------------------------


# Adjust the object fields relative to the cpml
array.nx = array.nx + 2*array.cpml
array.nz = array.nz + 2*array.cpml
array.dx = array.dx
array.dz = array.dz


extent = (array.cpml, array.nx-array.cpml, array.nz-array.cpml, array.cpml)
if array.channel == 'Ex' or array.channel == 'Ez':
	array.source_location = np.array([ (ex + array.cpml)/array.dx, 
		(ez + array.cpml)/array.dz])
	array.dt = edt
else:
	array.source_location = np.array([(sx + array.cpml)/array.dx,
		(sz + array.cpml)/array.dz]) 
	array.dt = sdt

array.initial_position = initial_coords
array.final_position = final_coords
array.dr = delta

# Construct the array
array.arraybuild()

# Get the timeseries for each reciever
array.getrcx()

array.tsplot()



