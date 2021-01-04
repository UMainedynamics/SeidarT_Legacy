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


x,y,z,t = indvar( electromag, domain)

extent = 