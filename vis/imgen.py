#!/usr/bin/env python3

"""
Functions for plotting a 2D vector/quiver plot over the model image
"""

import numpy as np 
import matplotlib.pyplot as plt
import argparse
from definitions import *
import matplotlib.image as mpimg


# ============================ Create the objects =============================
class FDTDImage:
    def __init__(self, prjfile, inputdat):
        self.prjfile = prjfile
        self.x = None
        self.z = None
        self.velocityfile = inputdat
        self.xvelocityfile = 'Vx' + inputdat[2:]
        self.zvelocityfile = 'Vz' + inputdat[2:]
        self.channel = inputdat[0:2]
        self.srcx = None
        self.srcz = None
        self.extent = None
        self.background = None
        self.ax = None 
        self.fig = None
        self.nx = None 
        self.nz = None 
        self.dx = None 
        self.dz = None 
        self.cpml = None
        self.plotfile = 'vector.' + inputdat[2:-3] + '.png'
    
    def getprjvals(self):
        # Let's initiate the domain
        domain, material, seismic, electromag = loadproject(
            self.prjfile,
            Domain(), 
            Material(),
            Model(),
            Model()
        )
        cpml = int(domain.cpml[0])
        nx = int(domain.nx[0]) + 2*cpml
        nz = int(domain.nz[0]) + 2*cpml
        dx = float(domain.dx[0])
        dz = float(domain.dz[0])
        # The EM code uses a slightly different mesh
        if self.channel == 'Ex':
            nx = nx - 1
        if self.channel == 'Ez':
            nz = nz - 1
        
        self.extent = (
            -cpml, 
            (nx-cpml), 
            (nz-cpml), 
            -cpml
        )

        if self.channel == 'Ex' or self.channel == 'Ez':
            self.srcx = float(electromag.x[0])/dx + cpml + 1
            self.srcz = float(electromag.z[0])/dz + cpml + 1
        else:    
            self.srcx = float(seismic.x[0])/dx + cpml + 1
            self.srcz = float(seismic.z[0])/dz + cpml + 1
        
        # Define tick locations for plotting
        self.xticklocs = np.array(
            [
                0, 
                (nx-2*cpml)/4, 
                int((nx-2*cpml)/2), 
                3*(nx-2*cpml)/4, 
                nx - 2*cpml
            ]
        )
        self.yticklocs = np.array(
            [
                0, 
                (nz-2*cpml)/4, 
                int((nz-2*cpml)/2), 
                3*(nz-2*cpml)/4, 
                nz - 2*cpml
            ]
        )
        # Load the model image and assign variables
        self.background = mpimg.imread(domain.imfile)
        self.cpml = cpml
        self.nx = nx 
        self.nz = nz
        self.dx = dx 
        self.dz = dz
        
    # -------------------------------------------------------------------------
    def quiverplot(self, papercolumnwidth = 7.2):
        # buildmesh
        # Get axes values
        x = (np.arange(-self.cpml, self.nx-self.cpml))
        z = (np.arange(-self.cpml, self.nz-self.cpml))
        # Append the cpml values 
        x, z = np.meshgrid(x,z)

        u = datinput(self.xvelocityfile, self.nx, self.nz)
        v = datinput(self.zvelocityfile, self.nx, self.nz)
        # Set the figure size to be for a full two column width
        
        # Create the figure and axes objects
        self.fig = plt.figure(
            figsize = [papercolumnwidth, self.nz*papercolumnwidth/self.nx]
        )
        self.ax = plt.gca()
        
        # add the model 
        self.ax.imshow(
            self.background,
            alpha = 0.7, 
            extent=[0, self.nx-2*self.cpml, self.nz-2*self.cpml, 0]
        )
        
        # add quiver
        q = self.ax.quiver(
            x, z, u, v, 
            headwidth = 0.5, 
            headlength = 1,
            headaxislength = 1, 
            scale = 16, 
            minlength = 0.1
        )
    
    def magnitudeplot(self, papercolumnwidth = 7.2):
        dat = datinput(self.velocityfile, self.nx, self.nz)
        
        self.fig = plt.figure(
            figsize = [papercolumnwidth, self.nz*papercolumnwidth/self.nx]
        )
        self.ax = plt.gca()
        self.ax.imshow(
            dat, 
            cmap = 'seismic',
            extent = [0, self.nx, self.nz, 0]
        )
        self.ax.imshow(
            self.background,
            alpha = 0.3,
            extent = [0, self.nx-2*self.cpml,self.nz-2*self.cpml, 0]
        )
        
    
    # -------------------------------------------------------------------------
    def addlabels(self):
        # Set axes labels
        xticklabels = (self.xticklocs)*self.dx
        yticklabels = (self.yticklocs)*self.dz
        self.ax.set_xlabel(r'X (m)')
        self.ax.xaxis.tick_top()
        self.ax.xaxis.set_label_position('top')
        self.ax.set_xticks(self.xticklocs)
        self.ax.set_xticklabels(xticklabels)
        self.ax.set_ylabel(r'Z (m)')
        self.ax.set_yticks(self.yticklocs)
        self.ax.set_yticklabels(yticklabels)
        
    def addsource(self):
        # Source location
        self.ax.scatter(
            self.srcx, 
            self.srcz,
            marker = '*', 
            s = 30, 
            linewidths = 1,
            edgecolor = (0.2, 0.2, 0.2, 1 ) 
        )

    def addrcx(self):
        pass 
    
    
# =============================================================================
def datinput(fn, nx, nz):
    # Get the U,V components; the velocity values
    f = FortranFile(fn, 'r')
    dat = f.read_reals(dtype = 'float32')
    dat = dat.reshape(nz, nx)
    return(dat)
   