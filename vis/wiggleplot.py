#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 13:23:36 2020

@author: Ann Hill
University of Maine
"""
#CG: updated for receiver distance and gain2

import numpy as np
import matplotlib.pyplot as plt
import argparse


#__________________________SET UP COMMAND LINE ARGUMENTS_______________________


parser = argparse.ArgumentParser(description='run wiggle plots from the \
	 command line')

parser.add_argument(
    '-r', '--receiver_data',
    nargs=1, type=str,
    help="""The file path for the receiver array data""",
    required=True
)

parser.add_argument(
    '-p', '--prjfile',
    nargs=1, type=str,
    help='Project file path',
    required=True
)

parser.add_argument(
    '-g', '--gain', nargs = 1, type = int,
    required = True, help = """The linear horizontal exaggeration of the
    amplitude values from the receiver array file""",
    default = None
)

parser.add_argument(
    '-x', '--receiver_spacing', nargs = 1, type = int,
    required = True, help = """The horizontal distance between receivers, in meters""",
    default = None
)

parser.add_argument(
    '-d', '--columns', nargs = 1, type = int,
    required = True,
    help = """The frequency at which columns are pulled for
    plotting from the csv file""",
    default = None
)

parser.add_argument(
    '-n', '--single_plot_dist',
    nargs = 1, type = int,
    help="""the distance along the image surface where the amplitude values will
    be plotted""",
    required=True
)

parser.add_argument(
    '-c', '--channel',
    nargs = 1, type = str, required = True,
    help = """The channel to query. """
)



#____________________ASSIGN ARGUMENTS__________________________________________
args = parser.parse_args()
gain=args.gain[0]
prjfile = args.prjfile[0]
rcxfile=args.receiver_data[0]
spacing=args.receiver_spacing[0]
d = args.columns[0]
spd = args.single_plot_dist[0]
channel=args.channel[0]

#____________________PULL NECESSARY INFO FROM THE PRJ FILE_____________________
class Array:
    def __init__(self):
        super().__init__()
        self.build()

    def build(self):
      
# Initialize variables
        self.dt = None
        self.dr = None
        self.timeseries = None
        self.t = None
        self.channel= None

array = Array()
array.channel = channel

#Get the values
f=open(prjfile)

for line in f:
	#this section wrong, I think (CG, Sept 24 2020): reading the pixels, rather than the receivers. Assumes receivers are spaced every meter.
	if line[0] == 'D':
		temp = line.split(',')
		if temp[1] == 'dx':
			array.dx = float(temp[2].rsplit()[0] )

	# Source
	if array.channel == 'Ex' or array.channel == 'Ez':
		if line[0] == 'E':
			temp = line.split(',')
			if temp[1] == 'dt':
				array.dt = float(temp[2].rsplit()[0] )
	else:
		if line[0] == 'S':
			temp = line.split(',')
			if temp[1] == 'dt':
				array.dt = float(temp[2].rsplit()[0] )
f.close()



#___________________THE SCRIPT_________________________________________________

#Get the data from the file
A1 = np.genfromtxt(rcxfile, delimiter = ",")
#Convert horizontal meters to column number
n=np.around(spd/spacing)
n=n.astype(int)
#Extract that chosen column as an array from the matrix of all the data.
B1=A1[:,n]
#Define the number of rows in the column. This will be plotted on the vertical
#axis
C1a=len(B1)

#Add one to the number of rows as Python does not include the last number when
#counting
C1a=C1a+1
#Make the number of rows into a vector
C1=np.arange(1,C1a)
#Transpose C1 to make it more useable
C1=np.transpose(C1)
#Multiply C1 by the time step factor to get the two way travel time
C1=C1*array.dt

#Define a figure
plt.figure(1)
#Plot the number of rows vs. the amplitude values
plt.plot(B1,C1,c='red',ls='-',lw=0.5)
#Format the axes
plt.gca().invert_yaxis()
plt.gca().xaxis.set_ticks_position('top')
plt.gca().xaxis.set_label_position('top')
plt.xlabel('Amplitude')
# plt.gca().yaxis.set_ticks([])
plt.ylabel('Two-Way Travel Time (s)')
#Display the figure
plt.show()

#Define an empty matrix the dimensions of csv file. This will be filled in
#through the upcoming for loop
A3=np.nan*np.ones((A1.shape), np.float)

#Multipy the values in the csv file by the gain. This will make the variations in
#the amplitude visible in the upcoming plot. Base gain is so amplitude reaches across columns (gain2; spacing)
max_A1 = np.max(A1)
gain2=spacing/max_A1
A1=A1*gain2
A1=A1*gain

#Create a for loop to loop over the amplitude values by column. Add the column number
#to the amplitude values so when the values are plotted they are spaced out
#corresponding to the distance along the surface. Put these new values in the
#previously created empty matrix
for aa in range(0, A1.shape[1] ):
    A3[:,aa]=A1[:,aa]+aa

#Create another empty matrix. This will contain each of the selected columns
#that will be plotted
B3=np.nan*np.ones((A3.shape), np.float)
#Create a for loop to identify every nth column
for bb in range(0, np.size(A1, axis=1), d):

    #Add the identified column to the previously created empty matrix, and
	#multiply by the pixels to meters conversion factor, spacing
    B3[:,bb]=(A3[:,bb])*spacing
    # print(B3)
    #Create an accompanying vector of the number of rows corresponding to the
	#column by first creating an empty matrix
    C3=[]
	#Loop through the number of rows the column values will be plotted against
    for i in range(0, np.size(A1, axis=0)):
		#Multiply the values by the time step factor and append to the empty
		#matrix
        C3.append(i*array.dt)
    # print(C3)
	# C3=range(1, np.size(A1, axis=0)+1)
	#Multiply C3 by the timestep to get the two way travel time
    #Create a figure to plot the data
plt.figure(2)
    #Add one column to the figure of amplitude vs. "depth" for every run through
    #the loop
plt.plot(B3,C3,c='red',linewidth=0.5)
    #Format the axes
plt.gca().invert_yaxis()
plt.ylabel('Two-Way Travel Time (s)')
plt.gca().xaxis.set_ticks_position('top')
plt.gca().xaxis.set_label_position('top')
plt.xlabel('Distance Along Surface (m)')
#plt.gca().invert_yaxis()

#Display the figure
plt.show()

f.close()
