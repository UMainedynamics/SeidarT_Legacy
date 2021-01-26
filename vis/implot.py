#!/usr/bin/env python3

"""
Plot a snapshot of the wavefield with axes
"""

import numpy as np 
import matplotlib.pyplot as plt
import argparse
from definitions import *
import matplotlib.image as mpimg



project_file = 'twolayer.prj'
fn = 'Vx000500.dat'
channel = fn[0:2]
rcx = 'receivers.xyz'
plotfile = fn[:-3] + 'plot.png'


# ============================ Create the objects =============================
# Let's initiate the domain
domain, material, seismic, electromag = loadproject(
    project_file,
    Domain(), 
    Material(),
    Model(),
    Model()
)

# The prj file is loaded as string lists so we need to put values in the 
# correct formats
domain.cpml = int(domain.cpml[0])
nx = int(domain.nx[0]) + 2*domain.cpml
nz = int(domain.nz[0]) + 2*domain.cpml
domain.dx = float(domain.dx[0])
domain.dz = float(domain.dz[0])
extent = (
    domain.cpml, 
    (nx-domain.cpml), 
    (nz-domain.cpml), 
    domain.cpml
)

if channel == 'Ex':
    nx = nx-1

if channel == 'Ez':
	nz = nz-1

# Add the source location to plot
if channel == 'Ex' or channel == 'Ez':
    electromag.x = float(electromag.x[0])
    electromag.z = float(electromag.z[0])
    srcx = electromag.x/domain.dx + domain.cpml+1
    srcz = electromag.z/domain.dz + domain.cpml+1
    dt = float(electromag.dt[0])
else:
    seismic.x = float(seismic.x[0])
    seismic.z = float(seismic.z[0])
    srcx = seismic.x/domain.dx + domain.cpml+1
    srcz = seismic.z/domain.dz + domain.cpml+1
    dt = float(seismic.dt[0])

f = FortranFile(fn, 'r')
dat = f.read_reals(dtype = 'float32')
dat = dat.reshape(nz, nx)

background = mpimg.imread(domain.imfile)
bound = np.max([abs(np.min(dat)),abs(np.max(dat))])

# Set the figure size to be for a full two column width
papercolumnwidth = 7.2
fig = plt.figure(
    figsize = [papercolumnwidth, nz*papercolumnwidth/nx]
)
ax = plt.gca()
ax.imshow(
    dat,cmap='seismic', 
    extent=(0, (nx), (nz), 0),
    vmin=-bound,vmax=bound
)
ax.imshow(
    background,
    alpha = 0.3, 
    extent=extent, 
)

# Source location
ax.scatter(
    srcx, 
    srcz,
    marker = '*', 
    s = 30, 
    linewidths = 1,
    edgecolor = (0.2, 0.2, 0.2, 1 ) 
)


# Reciever locations if applicable
# if
xticklocs = np.array(
    [domain.cpml, nx/4, int(nx/2), 3*nx/4, nx - domain.cpml]
)
yticklocs = np.array(
    [domain.cpml, nz/4, int(nz/2), 3*nz/4, nz - domain.cpml]
)
# Set axes labels
xticklabels = (xticklocs-domain.cpml)*domain.dx
yticklabels = (yticklocs-domain.cpml)*domain.dz


ax.set_xlabel(r'X (m)')
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.set_xticks(xticklocs)
ax.set_xticklabels(xticklabels)
ax.set_ylabel(r'Z (m)')
ax.set_yticks(yticklocs)
ax.set_yticklabels(yticklabels)

plt.show()