#!/bin/bash

# Install script for SEIDART toolbox

# ----------------------------- Anaconda install ------------------------------
ver="v0.2"
echo "--------------------------------------------"
echo "SeidarT Anaconda-based installer $ver"
echo "Univ. of Maine / Univ. of Washington, 2020"
echo "--------------------------------------------"
echo "This installer will check for Anaconda/Miniconda and install a SeidarT environment prior to compiling."
echo "You will have the option to install Miniconda if no existing conda is found."
echo ""
bash conda_deps.sh ||
echo "Conda installation failed. Try installing dependencies the run the noconda_install script." ||
exit 1

`grep etc/profile.d/conda.sh ~/.bashrc`
conda activate SeidarT &&
echo "Successfully activated SeidarT environment for compiling" ||
echo "Could not find SeidarT conda environment. Exiting." ||
exit 1

echo ""
echo "Starting compiling/installation process..."

# -----------------------------------------------------------------------------
# Make sure we have a folder to put everything in
if [ ! -d "bin" ]; then
	mkdir bin
fi

# --------------------------------- Clean Up ----------------------------------
# Typically, install will overwrite everything but when there is a compile 
# error this can make debugging difficult if it goes unnoticed. Not everyone 
# can read the f2py or gfortran outputs

# Just clean up if specified 
clean=${1:-"none"}
if [[ $clean == "clean" ]] ; then 
        # Clear the bin folder
        rm -rf bin/* 
        # Remove .mod and .o (if any) files generated during the fortran compile
        rm fdtd/*.mod
        rm fdtd/*.o
        exit 1
fi

# -----------------------------------------------------------------------------
# Make sure everything is in unix format 
dos2unix vis/* 
dos2unix survey_wrappers/*
dos2unix exe/* 
dos2unix materials/*

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
cp vis/arraybuild.py bin/arraybuild
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
        bin/arraybuild \
        bin/codisplay \
        bin/wiggleplot \
        bin/im2anim \
        bin/orientation_tensor \
        bin/array2segy \
	bin/vtkbuild 


# Now do the bash scripts
cp survey_wrappers/common_offset bin/common_offset
cp survey_wrappers/common_midpoint bin/common_midpoint

cp io/array2sac bin/array2sac

chmod +x bin/common_offset bin/common_midpoint bin/array2sac


# ---------------- Move all other required files to bin folder ----------------
cp materials/material_functions.py bin/material_functions.py
cp materials/class_definitions.py bin/class_definitions.py
cp vis/imdefinitions.py bin/imdefinitions.py 

echo ""
echo "Done."
echo 'Type "conda activate SeidarT" to access the SeidarT environment.'