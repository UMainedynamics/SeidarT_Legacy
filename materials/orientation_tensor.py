#!/usr/bin/env python3

# Wrapper to generate the euler angles for the plunge and trend 
# then plot the results

import numpy as np
import argparse
import matplotlib.pyplot as plt
import mplstereonet

import orientsynth as ot

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description=""" """ )

parser.add_argument( 'output_file', nargs=1, type=str,
						help='Specify the file to save the outputs', 
						default=None)

parser.add_argument( '-n', '--npts', nargs = 1, type = int, required = False,
	help = """ Total number of grains in synthetic sample.""", default = [100])

parser.add_argument( '-p', '--plunge', nargs = 1, type = float, 
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
plot_bool = args.suppress_plotting[0] == 1
output_file = ''.join(args.output_file)

print(plot_bool)

# -----------------------------------------------------------------------------
euler_list= np.array([])
orten = np.array([])


euler_list, orten = ot.orientsynth(trend, plunge, amin, amax, npts)

# Save the euler angles
print(type(euler_list) )
np.savetxt(output_file, euler_list, delimiter = " ")

if plot_bool:
	fig = plt.figure( figsize = (7, 5) )
	ax = fig.add_subplot(111,projection='stereonet')
	ax.pole(euler_list[:,0]*180/np.pi, euler_list[:,1]*180/np.pi, 'g^', markersize=6)
	ax.grid()

	# plt.rc('text', usetex=True)
	plt.figtext(0.01, 0.95, 'trend = ' + str(trend)   )
	plt.figtext(0.01, 0.9, 'plunge = ' + str(plunge) )
	plt.figtext(0.01, 0.85, 'angle min/max = ' + str(amin) + '/' + str(amax) )
	plt.figtext(0.01, 0.8, 'N = ' + str(npts) )

	# plt.figtext(0.01, 0.7, r'\underline{Tensor Coefficients}')
	plt.figtext(0.01, 0.7, 'Tensor Coefficients')
	plt.figtext(0.01, 0.65, '$a_{11} = $' + str( round( orten[0,0], 5) ) )
	plt.figtext(0.01, 0.6, '$a_{22} = $' + str( round( orten[1,1], 5) ) )
	plt.figtext(0.01, 0.55, '$a_{33} = $' + str( round( orten[2,2], 5) ) )

	plt.figtext(0.01, 0.5, '$a_{12} = $' +str( round( orten[0,1],5 ) ) )
	plt.figtext(0.01, 0.45, '$a_{13} = $' +str( round( orten[0,2], 5) ) )
	plt.figtext(0.01, 0.4, '$a_{23} = $' +str( round( orten[1,2], 5) ) )

	plt.show()
