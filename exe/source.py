#!/usr/bin/env python3

"""Build the source function. The default is an impulse function with center 
frequency f0, but we can also generate other wavelets and chirplets."""

import argparse
import numpy as np 

from class_definitions import *


# =========================== Command Line Arguments ==========================


project_file = "dipping_bed.prj"
emsourcefile = 'em_source.dat'
seissourcefile = 'seismic_source.dat'

type = "impulse"

# ================================ Definitions ================================

def gaus0(t, dim, f):
    pass 

def gaus1(t, dim, f):
    pass 

def gaus2(t, dim, f):
    pass 

def 



# =============================================================================


# Load the project file 
# Let's initiate the domain
domain = Domain()
material = Material()
seismic = Model()
electromag = Model()


domain, material, seismic, electromag = loadproject(
    project_file, 
    domain, material, 
    seismic, 
    electromag
)
# Create the time vector 

# Create the source function

# Rotate 
