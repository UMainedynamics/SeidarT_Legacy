#!/bin/bash

# Install script for SEISART toolbox

# Compile the fortran code
f2py3 -c --fcompiler=gnu95 -m emfdtd2d emFDTD2d.f95
f2py3 -c --fcompiler=gnu95 -m seismicfdtd2d seismicFDTD2d.f95

# Compile the python modules
python3 -m py_compile material_functions.py dispaly_functions.py

# Create the executables
chmod +x prjrun.py
chmod +x prjbuild.py
