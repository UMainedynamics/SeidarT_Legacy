#!/bin/bash

# Install script for SEISART toolbox

# Compile the fortran code
f2py3 -c --fcompiler=gnu95 -m emfdtd2d emFDTD2D.f95
f2py3 -c --fcompiler=gnu95 -m seismicfdtd2d seismicFDTD2D.f95


