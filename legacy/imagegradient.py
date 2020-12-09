#!/usr/bin/env python3 

# Development script to smooth gradients (i.e. high density contrasts at air interface)


import numpy as np 
from class_definitions import * 
from scipy.io import FortranFile
# Globals 
prjfile = "dipping_bed.prj"

# Read the png and rho.dat images 
domain, material = loadproject(prjfile, Domain(), Material(), Model(), Model() )[0:2]
# NX = int(domain.nx[0]) + 2*int(domain.cpml[0])
# NZ = int(domain.nz[0]) + 2*int(domain.cpml[0])


# zeroshear
# imagegradient
airnum = int(
    material.material_list[
        material.material_list[:,1] == 'air', 0
    ][0]
)

gradmatrix = (domain.geometry != airnum).astype(int)
# Take the gradient in both directions
gradz = np.diff(gradmatrix, axis = 0)
gradx = np.diff(gradmatrix, axis = 1)

# For grady we will append a column of zeros at the beginning so that the value
# 1 is located at the interface but on the non-air side 
gradzpos = np.row_stack([np.zeros([gradz.shape[1] ]),gradz])
# For gradx we will append a row of zeros at the beginning 
gradxpos = np.column_stack([np.zeros([gradx.shape[0] ]),gradx]) 

# -1 also means that there is a air interface. We will need to append to the
# end of the array, then we can just flip the sign
gradzneg = -(np.row_stack( [gradz, np.zeros([gradz.shape[1]]) ] ) )
gradxneg = -(np.column_stack( [gradx, np.zeros([gradx.shape[0]]) ] ) )

# At the surface we want to have 15% of density 
grad15 = gradzpos + gradxpos + gradzneg + gradxneg
grad15[grad15>0] = 0.15 

# The next node down will be 50%. Do the same as above but append opposite
gradzpos = np.row_stack( [np.zeros([gradz.shape[1] ]),gradz] )[:-1,:]
gradxpos = np.column_stack( [np.zeros([gradx.shape[0] ] ),gradx])[:,:-1]
gradzneg = np.row_stack( [ gradz, np.zeros( [gradz.shape[1] ]) ] )[1:,:]
gradxneg = np.column_stack( [ gradx, np.zeros( [gradx.shape[0] ]) ] )[:,1:]

grad50 = gradzpos + gradxpos + gradzneg + gradxneg 
grad50[grad50 > 0] = 0.5

# Put grad50 and grad15 together but take the lowest; where gradient is both
# 0.15 and 0.5, we want 0.15
grad = (grad50 > 0).astype(int) - (grad15 > 0).astype(int)  
grad = 0.5 * (grad > 0).astype(int) 
grad = grad + grad15 

