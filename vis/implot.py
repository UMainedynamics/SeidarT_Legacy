#!/usr/bin/env python3

"""
Plot a snapshot of the wavefield with axes
"""

import numpy as np 
import matplotlib.pyplot as plt
import argparse
from definitions import *



project_file = 'twolayer.prj'

# ============================ Create the objects =============================
# Let's initiate the domain
domain, material, seismic, electromag = loadproject(
    project_file,
    Domain(), 
    Material(),
    Model(),
    Model()
)

# Get axes values
x,y,z,t = indvar( electromag, domain)

domain.cpml = int(domain.cpml[0])
nx = domain.nx + 2*domain.cpml
nz = domain.nz + 2*domain.cpml
domain.dx = float(domain.dx[0])
domain.dz = float(domain.dz[0])
extent = 