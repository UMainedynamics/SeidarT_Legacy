#!/bin/bash

# Install script for SEISART toolbox
if [ ! -d "bin" ]; then
	mkdir bin	
fi

# Compile the fortran code
f2py3 -c --fcompiler=gnu95 -m emfdtd2d emFDTD2D.f95
f2py3 -c --fcompiler=gnu95 -m seismicfdtd2d seismicFDTD2D.f95
f2py3 -c --fcompiler=gnu95 -m orientsynth orientsynth.f95

mv *.so bin

# --------------------------- Create the executables --------------------------
# Start with the python scripts
cp prjbuild.py bin/prjbuild
cp prjrun.py bin/prjrun
cp arrayplot.py bin/arrayplot
cp codisplay.py bin/codisplay
cp im2gif.py bin/im2gif
cp orientation_tensor.py bin/orientation_tensor

# Change them to executables
chmod +x bin/prjbuild \
        bin/prjrun \
        bin/arrayplot \
        bin/codisplay \
        bin/im2gif \
        bin/orientation_tensor


# Now do the bash scripts
cp wide_angle bin/wide_angle
cp common_offset bin/common_offset
cp common_midpoint bin/common_midpoint

chmod +x bin/wide_angle bin/common_offset bin/common_midpoint


# ---------------- Move all other required files to bin folder ----------------
mv material_functions.py bin/material_functions.py