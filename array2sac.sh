#!/bin/bash

# Write individual reciever timeseries as SAC timeseries 

# Define some useful variables 


# Get information from the project file 


# Tally up the number of the recievers
num_cols=`head -1 $rcxfile | sed 's/[^,]//g' | wc -c`


for col in {1..$num_cols}
do
	echo $col 



# Get the reciever 


# Get the location of the reciever


# Write to Seismic Analysis Code format



done








# Some basic I/O operations to allow for converting the reciever timeseries to 
# # different outputs

# import glob 
# import numpy as np
# import scipy.io as scio
# from scipy.io import FortranFile
# import argparse
# from subprocess import call


# # ========================== Command Line Arguments ==========================
# parser = argparse.ArgumentParser(description="""This contains a set of 
# 	functions to read the all of the .dat outputs of the SeisarT modeling 
# 	routines.""" )

# parser.add_argument( '-o', '--output', nargs=1, type=str, required=True,
# 					help="""The data type to convert to. Available options 
# 					are sac (default), segy, and mat. """, default = 'csv')

# parser.add_argument( '-f', '--project_file', nargs=1, type=str, required=True, 
# 					help="""The project file associated with the .dat files""", 
# 					default=None)

# parser.add_argument( '-m', '--meta_file', nargs = 1, type = str, required=True,
# 					help=""" """, default=None)


# # Get the arguments
# args = parser.parse_args()
# project_file = ''.join(args.project_file)
# meta_file = ''.join(args.meta_file)
# fmt = ''.join(args.output)

# # ============================= Create the Object =============================

# class Recievers:
# 	def __init__(self):
# 		super().__init__()
# 		self.build()

# 	def build(self):
# 		# Initialize variables
# 		self.project_file
# 		self.meta_file
# 		self.dt
# 		self.source_xyz
# 		self.reciever_locations 
# 		self.channel 
# 		self.time_vector
# 		self.distance 
# 		self.


# # =========================== Function Definitions ============================


# def file_read():


# def array2csv(nx, ny):

# 	for fn in glob.glob('*.dat'):
# 		npdat = read_dat(fn, nx, ny)
# 		# replace the .dat extension with .csv
# 		sfn = fn[0:-3] + 'csv'
# 		np.savetxt(sfn, npdat, delimiter=",")

# # -----------------------------------------------------------------------------

# def dat2mat(nx, ny):
# 	# similar to above but saving to matlab specific file
# 	for fn in glob.glob(chan + '*.dat'):
# 		npdat = read_dat(fn, nx, ny)
# 		sfn = fn[0:-3] + 'mat'
# 		scio.savemat(sfn, npdat)


# # -----------------------------------------------------------------------------

# def read_dat(fn, nx, ny):
# 	f = FortranFile(fn, 'r')
# 	dat = f.read_reals(dtype = 'float64')
# 	dat = dat.reshape(ny, nx)
# 	f.close()

# 	return(dat)

# # ============================== Run the program ==============================

# # Read the project file and get the nx and ny values
# f = open(project_file)

# for line in f:

# 	# All domain inputs must be input except for ny and dy
# 	if line[0] == 'D':
# 		temp = line.split(',')

# 		if temp[1] == 'nx':
# 			nx = int( temp[2].rsplit()[0] )
# 		if temp[1] == 'nz':
# 			nz = int( temp[2].rsplit()[0] )
# 		if temp[1] == 'cpml':
# 			cpml = int( temp[2].rsplit()[0])

# f.close()

# nx = nx + 2*cpml
# nz = nz + 2*cpml

# if fmt == 'mat':
# 	print('Converting .dat files to .mat. This may take a moment.\n')
# 	dat2mat(nx, nz)
# else:
# 	print('Converting .dat files to .csv. This may take a moment.\n')
# 	dat2csv(nx, nz)


# if rm == 1:
# 	print('Deleting all .dat files.\n')
# 	call('rm *.dat', shell=True)
