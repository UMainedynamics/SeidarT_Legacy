#!/usr/bin/env python3
#
# This script will read an image and build the template project file template
# to be used in the seisarT program
#
# -----------------------------------------------------------------------------

import argparse
import numpy as np
import matplotlib.image as mpimg

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""The SeidarT software requires a
	.PNG image that is used to construct the model domain for seismic and
	electromagnetic wave propagation. Given the image file, a project file
	will be constructed which contains all the necessary parameters to be
	read in to the finite differences time domain modeling schemes.""" )


parser.add_argument(
    '-i','--imagefile', 
    nargs=1, type=str, required = True,
    help='the full file path for the image', default=None
)

parser.add_argument(
    '-p', '--prjfile',
    nargs=1, type=str, required = False, default = 'jordan_downs.prj',
    help = """name of output file path with extension .prj and excluding
    the full path directory"""
)

# Get the arguments
args = parser.parse_args()
image_file = ''.join(args.imagefile)
project_file = ''.join(args.prjfile)

new_line = '\n'
# ------------------------ Some Necessary Definitions -------------------------

def image2int(imfilename):
    # read the image
    img = mpimg.imread(imfilename)

    # Convert RGB to a single value
    rgb_int = np.array(65536*img[:,:,0] +  255*img[:,:,1] + img[:,:,2])

    # Get the unique values of the image
    rgb_uni = np.unique(rgb_int)

    # We want the unique rgb values too
    rgb = np.zeros( [len(rgb_uni), 3] )

    # reshape the image. We know it's three channels
    img_vect = np.zeros( [np.prod(rgb_int.shape), 3] )
    img_vect[:,0] = np.reshape(img[:, :, 0], np.prod(np.shape(img[:, :, 0]) ) )
    img_vect[:,1] =	np.reshape(img[:, :, 1], np.prod(np.shape(img[:, :, 1]) ) )
    img_vect[:,2] =	np.reshape(img[:, :, 2], np.prod(np.shape(img[:, :, 2]) ) )

    # mat_id = np.array( range(0, len(rgb_uni) ) )
    for ind in range(0, len(rgb_uni) ):
    	rgb_ind = np.reshape(rgb_int == rgb_uni[ind], [np.prod(rgb_int.shape)])
    	rgb[ind,:] = (img_vect[rgb_ind,:])[0,:]
    	rgb_int[ rgb_int == rgb_uni[ind] ] = ind


    if np.max(rgb) <= 1.0:
    	rgb = rgb * 255
    	rgb = rgb.astype(int)

    return rgb_int.astype(int), rgb


# -------------------------------- Add a header -------------------------------
header_comment = """
# This is a project file template for the SeidarT software. In order to run the
# model for seismic, electromagnetic or both, the required inputs must be
#
# Domain Input Values:
#	dim 		- STR; either '2' or '2.5'; default is '2'
#	nx,ny,nz 	- INT; the dimensions of the image. If dim = 2.5, and ny is
#			  empty then default ny=1
#	dx,dy,dz	- REAL; the spatial step size for each dimension in meters. If
#			  dim = 2.5 and dy is empty then default dy=min(dx,dz)
#
# Material Input Values:
#	id 		- INT; the identifier given to each unique rgb value as it
#			  is read into the computer. It's recommended to use this
#			  script to make sure it's sorted correctly.
#	R/G/B 		- STR; the 0-255 values for each color code.
#	Temperature 	- REAL; temperature in Celsius.
#	Attenuation 	- REAL; (placeholder) will be attenuation length soon.
#	Density 	- REAL; density in kg/m^3
#	Porosity 	- REAL; percent porosity
#	Water_Content 	- REAL; percent of pores that contain water
#	Anisotropic 	- BOOL; whether the material is anisotropic (True) or
#			  isotropic (False).
#	ANG_File 	- STR; if Anisotrpic is True then the full path to the
#			  .ang file is supplied. The .ang file is a delimited text
#			  file that contains the 3-by-n array of euler rotation
#			  angles in radians.
#
#		or alternatively...
#	C11-C66 	- REAL; the stiffness coefficients with the appropriate id
#	E11-E33,S11-S33	- REAL; the permittivity and conductivity coefficients and
#			  'id' value corresponding to the coefficients along the diagonal
#			  of their respective tensors.
#
#
# Source Input Values:
#	dt 		- REAL; dx/(2*maxvelcity)
#	steps 		- INT; the total number of time steps
#	x,y,z 		- REAL; locations in meters, +x is to the right, +z is down, +y is into the screen
#	f0 		- REAL; center frequency for the guassian pulse function if
#			  'source_file' isn't supplied
#	theta 		- REAL; source orientation in the x-z plane,
#	phi 		- REAL; source orientation in the x-y plane for 2.5/3D only,
#	source_file	- STR; the pointer to the text file that contains the source
#			  timeseries as a steps-by-1 vector.
#
# 	**phi and theta are the rotation angles for spherical coordinates so
#		x = r sin(theta)cos(phi)
#		y = r sin(theta)sin(phi)
#		z = r cos(theta)
#
#	Theta is the angle from the z-axis (+ down relative to image), phi is the
#	angle from x-axis in the x-y plane (+ counterclockwise when viewed from above)
#
# Written by Steven Bernsen
# University of Maine
# -----------------------------------------------------------------------------

"""

# ---------------------------- Read the Image File ----------------------------
im, rgb = image2int(image_file)
im = im.transpose()
mat_id = np.unique(im)
# Start writing the project file. To allow for headers we will start all
# pertinant information after

with open(project_file, 'w') as prj:
	prj.write(header_comment)
	prj.write(new_line)
	prj.write('I,'+image_file+new_line)
	prj.write(new_line)



# ------------------------- Write Domain Parameters ---------------------------
dim = 'D,dim,2'
nx = 'D,nx,' + str(np.shape(im)[0])
ny = 'D,ny,n/a'
nz = 'D,nz,' + str(np.shape(im)[1])
dx = 'D,dx,'
dy = 'D,dy,n/a'
dz = 'D,dz,'
cpml = 'D,cpml,20'
nmat = 'D,nmats,' + str(len( np.unique(im) ))
tfile = 'D,tfile,'
with open(project_file, 'a') as prj:
	prj.write(dim+new_line)
	prj.write(nx+new_line)
	prj.write(ny+new_line)
	prj.write(nz+new_line)
	prj.write(dx+new_line)
	prj.write(dy+new_line)
	prj.write(dz+new_line)
	prj.write(cpml+new_line)
	prj.write(nmat+new_line)
	prj.write(tfile + new_line)
	prj.write(new_line)




# ------------------------- Write Material Parameters -------------------------

header = ("# number, id, R/G/B, Temperature, Attenuation, Density, Porosity, "
				"Water_Content, Anisotropic, ANG_File")


with open(project_file, 'a') as prj:
	i = 0

	prj.write(header + new_line )
	for x in mat_id:
		ln = ('M,' + str(x) + ',,' + str(rgb[x,0])  + '/' +
			str(rgb[x,1]) + '/' + str(rgb[x,2]) +
			',,,,,,,')
		prj.write( ln + new_line)

	prj.write(new_line)


# ------------------------- Write Seismic Parameters --------------------------
dt = 'dt,'
steps = 'time_steps,'
x = 'x,'
y = 'y,'
z = 'z,'
f0 = 'f0,'
theta = 'theta,0'
phi = 'phi,0'
source_file='source_file,'

comm = '# The source parameters for the seismic model'
header = '# id, C11, C12, C13, C22, C23, C33, C44, C55, C66, rho'

with open(project_file, 'a') as prj:
	i = 0
	prj.write(comm + new_line)
	prj.write('S,' + dt + new_line)
	prj.write('S,' + steps + new_line)
	prj.write('S,' + x + new_line)
	prj.write('S,' + y + new_line)
	prj.write('S,' + z + new_line)
	prj.write('S,' + f0 + new_line)
	prj.write('S,' + theta + new_line)
	prj.write('S,' + phi + new_line)
	prj.write(new_line)

	prj.write(header + new_line )
	for ind in mat_id:
		prj.write( 'C,' + str(ind) + ',,,,,,,,,,' + new_line)

	prj.write(new_line)



# -------------------------- Write Radar Parameters ---------------------------

comm = '# The source parameters for the electromagnetic model'
header = '# id, e11, e22, e33, s11, s22, s33'

with open(project_file, 'a') as prj:
	i = 0

	prj.write(comm + new_line)
	prj.write('E,' + dt + new_line)
	prj.write('E,' + steps + new_line)
	prj.write('E,' + x + new_line)
	prj.write('E,' + y + new_line)
	prj.write('E,' + z + new_line)
	prj.write('E,' + f0 + new_line)
	prj.write('E,' + theta + new_line)
	prj.write('E,' + phi + new_line)

	prj.write(new_line)

	prj.write(header + new_line )
	for ind in mat_id:
		prj.write( 'P,' + str(ind) + ',,,,,,,,,,' + new_line)

	prj.write(new_line)



# ------------------- Write Additional Survey File Tempates -------------------
# meta_header = """
# # Options for the following fields
# # project_file - (STRING) the file path to the project file
# # survey_type - (STRING) the type of survey you would like to model. Available
# #				options are 'co' = common offset, 'cmp' = common midpoint,
# #				'wa' = wide angle.
# #
# # The following inputs change given the survey type. There are additional
# # values that need to be passed in the wrapper
# #
# # delta (FLOAT)
# #					'wa' the spacing between each reciever
# #					'cmp' the change in the source and the reciever distance from
# #						the common midpoint (given below)
# #					'co' the shift in the same direction of the source and
# #						reciever. The spacing between the source and reciever
# #						remains constant so they are moved in the same direction
# #
# # initial_position (FLOAT)
# #					'wa' the initial reciever location along the array in meters
# #					'cmp' the reciever location
# #					'co'  the reciever location
# #
# # final_position (FLOAT)
# #					'wa' the final reciever location along the array in meters
# #					'cmp' this is the same value as initial position; moot
# #					'co' 			'' ''			''	''
# #
# """

# if args.meta_file is not None:
# 	with open(meta_file, 'w') as meta:
# 		meta.write(meta_header + new_line)
# 		meta.write('project_file: ' + project_file + new_line)
# 		meta.write('survey_type: wa' + new_line)
# 		meta.write('delta: 1' + new_line)
# 		meta.write('initial_position: 0  0  0' + new_line)
# 		meta.write('final_position: ' + str(np.shape(im)[0]) + '0 ' + str(np.shape(im)[0]) + new_line )
# 		meta.write('reciever_file: None' + new_line)
# 		meta.write('source_file: None' )
