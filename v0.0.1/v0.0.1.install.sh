#!/bin/bash

# Install script for SEIDART toolbox

# Make sure we have a folder to put everything in
if [ ! -d "../bin" ]; then
	mkdir ../bin
fi

# --------------------------------- Clean Up ----------------------------------
# Typically, install will overwrite everything but when there is a compile 
# error this can make debugging difficult if it goes unnoticed. Not everyone 
# can read the f2py or gfortran outputs

# Just clean up if specified 
clean=${1:-"none"}
if [[ $clean == "clean" ]] ; then 
        # Clear the bin folder
        rm -rf  ../bin/* 
        # Remove .mod and .o (if any) files generated during the fortran compile
        rm source/*.mod
        rm source/*.o
        exit 1
fi

# -----------------------------------------------------------------------------
# Make sure everything is in unix format 
cd source
dos2unix * 

# Compile the fortran code
f2py3 -c --fcompiler=gnu95 -m emfdtd2d emFDTD2D.f95
f2py3 -c --fcompiler=gnu95 -m seismicfdtd2d seismicFDTD2D.f95
f2py3 -c --fcompiler=gnu95 -m emfdtd25d emFDTD25D.f95
f2py3 -c --fcompiler=gnu95 -m seismicfdtd25d seismicFDTD25D.f95
mv *.so ../../bin


# Synthetic microstructure
f2py3 -c --fcompiler=gnu95 -m orientsynth orientsynth.f95
mv *.so ../../bin

# --------------------------- Create the executables --------------------------
# Start with the python scripts
cp ./prjbuild.py ../../bin/prjbuild
cp ./prjrun.py ../../bin/prjrun
cp ./sourcefunction.py ../../bin/sourcefunction
cp ./orientation_tensor.py ../../bin/orientation_tensor

# Move the visualization tools
cp ./arraybuild.py ../../bin/arraybuild
cp ./rcxdisplay.py ../../bin/rcxdisplay
cp ./im2anim.py ../../bin/im2anim
cp ./vtkbuild.py ../../bin/vtkbuild 
cp ./wiggleplot.py ../../bin/wiggleplot 

# move the conversion scripts
cp ./array2segy.py ../../bin/array2segy

# Change them to executables
chmod +x ../../bin/prjbuild \
        ../../bin/prjrun \
        ../../bin/sourcefunction \
        ../../bin/arraybuild \
        ../../bin/rcxdisplay \
        ../../bin/wiggleplot \
        ../../bin/im2anim \
        ../../bin/orientation_tensor \
        ../../bin/array2segy \
	../../bin/vtkbuild 


# Now do the bash scripts
cp ./common_offset ../../bin/common_offset
cp ./common_midpoint ../../bin/common_midpoint

cp ./array2sac ../../bin/array2sac

chmod +x ../../bin/common_offset \
        ../../bin/common_midpoint \
        ../../bin/array2sac


# ---------------- Move all other required files to bin folder ----------------
cp ./material_functions.py ../../bin/material_functions.py
cp ./class_definitions.py ../../bin/class_definitions.py
cp ./imdefinitions.py ../../bin/imdefinitions.py 

cd ..