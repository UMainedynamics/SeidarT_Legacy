#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 13:23:36 2020

@author: Ann
"""

import numpy as np
#!!!!! For ion import pyplot from matplotlib 
#!!!!! If you want to import matplotlib then call
#!!!!!      import matplotlib as mpl
#!!!!! It will be easier to understand what everyone else is saying
import matplotlib.pyplot as plt
import argparse


#__________________________SET UP COMMAND LINE ARGUMENTS_______________________

    
parser = argparse.ArgumentParser(description='run wiggle plots from the \
	 command line')

parser.add_argument( '-r', '--reciever_array_file', nargs=1, type=str, 
	dest='input', help='the file path for the receiver array data', 
	required=True)

parser.add_argument( '-p', '--project_file', nargs=1, type=str, dest='prj',
 	help='the file path for the project file', required=True)

parser.add_argument( '-g', '--gain', nargs = 1, type = int, 
	required = True, dest='gain', help = "The horizontal exaggeration of the \
		 amplitude values from the receiver array file", 
		    default = None)

parser.add_argument( '-d', '--columns', nargs = 1, type = int, 
 	required = True, dest='cols', help = "The frequency at which columns are \
		pulled for plotting from the csv file", 
		    default = None)

parser.add_argument( '-n', '--column_number', nargs = 1, type = int, 
					dest='n', help='the distance along the imagesurface \
	 where the amplitude values will be plotted', required=True)

parser.add_argument( '-c', '--channel', nargs = 1, type = str, required = True,
 	dest='channel', help = """The channel to query. """)



#____________________ASSIGN ARGUMENTS__________________________________________
args = parser.parse_args()    
gain=args.gain[0]
prjfile = args.prj[0] 
rcxfile=args.input[0]
d = args.cols[0]
n = args.n[0]
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
channel=array.channel

#Get the values
f=open(prjfile)

for line in f:
	
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
plt.plot(B1,C1,ls='-',lw=0.5)
#Format the axies
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
#the amplitude visible in the upcoming plot
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
	#multiply by the pixels to meters conversion factor, array.dx
    B3[:,bb]=(A3[:,bb])*array.dx
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
    plt.plot(B3,C3,c='blue',linewidth=0.5)
    #Format the axies
    plt.gca().invert_yaxis()
    plt.ylabel('Two-Way Travel Time (s)')
    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    plt.xlabel('Distance Along Surface (m)')

#Display the figure
plt.show()

f.close()