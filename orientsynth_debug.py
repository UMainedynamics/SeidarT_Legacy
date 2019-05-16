#!/usr/bin/env python3

# Wrapper to generate the euler angles for the plunge and trend 
# then plot the results

import numpy as np
import argparse
import matplotlib.pyplot as plt

import orientsynth as ot


# Get the arguments
npts = 1000
plunge= 0
trend= 0
amin= 0
amax = 15

euler_list = np.array([])
end_orientation = np.array([])

# The order that the tuples are returned is the same order that they are 
# passed to the function
euler_list, end_orientation = ot.orientsynth(trend, plunge, amin, amax, 
	npts)


plot_orientation = np.zeros([npts,2])
plot_orientation[:,0] = -end_orientation[:,0]/(1-end_orientation[:,2])
plot_orientation[:,1] = end_orientation[:,1]/(1-end_orientation[:,2])

fig, ax = plt.subplots()

ax.scatter(plot_orientation[:,0], plot_orientation[:,1])

plt.show()