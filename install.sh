#!/bin/bash

# Install script for SEISART toolbox
if [ ! -d "bin" ]; then
	mkdir bin	
fi

# Compile the fortran code
#2D
cd fdtd
f2py3 -c --fcompiler=gnu95 -m emfdtd2d emFDTD2D.f95
f2py3 -c --fcompiler=gnu95 -m seismicfdtd2d seismicFDTD2D.f95
mv *.so ../bin
cd ..

# 2.5D
f2py3 -c --fcompiler=gnu95 -m seismicfdtd25d seismicFDTD25D.f95
f2py3 -c --fcompiler=gnu95 -m emfdtd25d emFDTD25D.f95
mv *.so bin

# Synthetic microstructure
cd materials
f2py3 -c --fcompiler=gnu95 -m orientsynth orientsynth.f95
mv *.so ../bin
cd ..

# --------------------------- Create the executables --------------------------
# Start with the python scripts
cp exe/prjbuild.py bin/prjbuild
cp exe/prjrun.py bin/prjrun
cp materials/orientation_tensor.py bin/orientation_tensor

# Move the visualization tools 
cp vis/arrayplot.py bin/arrayplot
cp vis/codisplay.py bin/codisplay
cp vis/im2anim.py bin/im2anim

# Change them to executables
chmod +x bin/prjbuild \
        bin/prjrun \
        bin/arrayplot \
        bin/codisplay \
        bin/im2anim \
        bin/orientation_tensor


# Now do the bash scripts
cp survey_wrappers/wide_angle bin/wide_angle
cp survey_wrappers/common_offset bin/common_offset
cp survey_wrappers/common_midpoint bin/common_midpoint

cp io/array2sac bin/array2sac

chmod +x bin/wide_angle bin/common_offset bin/common_midpoint bin/array2sac


# ---------------- Move all other required files to bin folder ----------------
cp materials/material_functions.py bin/material_functions.py