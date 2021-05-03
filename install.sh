#!/bin/bash

# Install script for SEIDART toolbox

# ----------------------------- Anaconda install ------------------------------
VER="v0.2"
echo "--------------------------------------------
SeidarT Anaconda-based installer $VER
Univ. of Maine / Univ. of Washington, 2020
--------------------------------------------
This installer will check for Anaconda/Miniconda
and install a SeidarT environment prior to compiling from
source.
You will have the option to install Miniconda
if no existing conda is found.

Please follow instructions in prompts.
"
read -rp $'Press Enter to continue...\n'


echo '
Do you wish to delete source files after installation?
This will not affect how the program runs, but you will
not be able to edit source code.
'
# ask if we should enable pure end-user mode and delete the source scripts
read -rp $'Type "yes" and press Enter to delete source files. (default: no)\n' EUMODE

if [[ "$EUMODE" == "yes" || "$EUMODE" == "Yes" || "$EUMODE" == "YES" ]] ; then
        # move instead of copy (end-user mode)
        FX=mv
else
        # straight copy (for developer mode)
        echo "Developer mode enabled. Source scripts will not be deleted."
        FX=cp
fi

bash conda_deps.sh ||
echo "Conda installation failed. Try installing dependencies the run the noconda_install script." ||
exit 1

`grep etc/profile.d/conda.sh ~/.bashrc`
conda activate SeidarT &&
echo "conda activate SeidarT : Successfully activated SeidarT environment." ||
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
CLEAN=${1:-"none"}
if [[ $CLEAN == "clean" ]] ; then 
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
dos2unix fdtd/*

# Compile the fortran code
#2D
cd fdtd
f2py3 -c --fcompiler=gnu95 -DNPY_NO_DEPRECATED_API -m emfdtd2d emFDTD2D.f95
f2py3 -c --fcompiler=gnu95 -DNPY_NO_DEPRECATED_API -m seismicfdtd2d seismicFDTD2D.f95
f2py3 -c --fcompiler=gnu95 -DNPY_NO_DEPRECATED_API -m emfdtd25d emFDTD25D.f95
f2py3 -c --fcompiler=gnu95 -DNPY_NO_DEPRECATED_API -m seismicfdtd25d seismicFDTD25D.f95
mv *.so ../bin
cd ..

# Synthetic microstructure
cd materials
f2py3 -c --fcompiler=gnu95 -m orientsynth orientsynth.f95
mv *.so ../bin
cd ..

# --------------------------- Create the executables --------------------------
# Start with the python scripts
$FX exe/prjbuild.py bin/prjbuild
$FX exe/prjrun.py bin/prjrun
$FX exe/sourcefunction.py bin/sourcefunction
$FX materials/orientation_tensor.py bin/orientation_tensor

# Move the visualization tools
$FX vis/arraybuild.py bin/arraybuild
$FX vis/rcxdisplay.py bin/rcxdisplay
$FX vis/im2anim.py bin/im2anim
$FX vis/vtkbuild.py bin/vtkbuild 
$FX vis/wiggleplot.py bin/wiggleplot 
$FX vis/imgen.py bin/imgen.py # The generalized image functions module
$FX vis/imvector.py bin/imvector
$FX vis/vectoranim.py bin/vectoranim
$FX vis/implot.py bin/implot

# move the conversion scripts
$FX io/array2segy.py bin/array2segy

# Change them to executables
chmod +x bin/prjbuild \
        bin/prjrun \
        bin/sourcefunction \
        bin/arraybuild \
        bin/rcxdisplay \
        bin/wiggleplot \
        bin/im2anim \
        bin/orientation_tensor \
        bin/array2segy \
	bin/vtkbuild \
        bin/imvector \
        bin/vectoranim \
        bin/implot



# Now do the bash scripts
$FX survey_wrappers/common_offset bin/common_offset
$FX survey_wrappers/common_midpoint bin/common_midpoint

$FX io/array2sac bin/array2sac

chmod +x bin/common_offset bin/common_midpoint bin/array2sac


# ---------------- Move all other required files to bin folder ----------------
$FX materials/material_functions.py bin/material_functions.py
$FX materials/definitions.py bin/definitions.py

echo ""
echo "Done."
echo 'Type "conda activate SeidarT" to access the SeidarT environment.'
