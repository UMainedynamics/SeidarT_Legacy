#!/bin/bash

# Install script for SEIDART toolbox

# Make sure we have a folder to put everything in
if [ ! -d "bin" ]; then
	mkdir bin
fi

# Compile the fortran code
#2D
cd fdtd
f2py3 -c --fcompiler=gnu95 -m emfdtd2d emFDTD2D.f95
f2py3 -c --fcompiler=gnu95 -m seismicfdtd2d seismicFDTD2D.f95
f2py3 -c --fcompiler=gnu95 -m emfdtd25d emFDTD25D.f95
f2py3 -c --fcompiler=gnu95 -m seismicfdtd25d seismicFDTD25D.f95
mv *.so ../bin
cd ..

# Synthetic microstructure
cd materials
f2py3 -c --fcompiler=gnu95 -m orientsynth orientsynth.f95
mv *.so ../bin
cd ..

# --------------------------- Create the executables --------------------------
# Start with the python scripts
cp exe/prjbuild.py bin/prjbuild
cp exe/prjrun.py bin/prjrun
cp exe/sourcefunction.py bin/sourcefunction
cp materials/orientation_tensor.py bin/orientation_tensor
# cp materials/class_definitions.py bin/class_definitions

# Move the visualization tools
cp vis/arrayplot.py bin/arrayplot
cp vis/codisplay.py bin/codisplay
cp vis/im2anim.py bin/im2anim
cp vis/vtkbuild.py bin/vtkbuild 
cp vis/wiggleplot.py bin/wiggleplot 

# move the conversion scripts
cp io/array2segy.py bin/array2segy

# Change them to executables
chmod +x bin/prjbuild \
        bin/prjrun \
        bin/sourcefunction \
        bin/arrayplot \
        bin/codisplay \
        bin/wiggleplot \
        bin/im2anim \
        bin/orientation_tensor \
        bin/array2segy \
	bin/vtkbuild 


# Now do the bash scripts
cp survey_wrappers/wide_angle bin/wide_angle
cp survey_wrappers/common_offset bin/common_offset
cp survey_wrappers/common_midpoint bin/common_midpoint

cp io/array2sac bin/array2sac

chmod +x bin/wide_angle bin/common_offset bin/common_midpoint bin/array2sac


# ---------------- Move all other required files to bin folder ----------------
cp materials/material_functions.py bin/material_functions.py
cp materials/class_definitions.py bin/class_definitions.py
cp vis/imdefinitions.py bin/imdefinitions.py 