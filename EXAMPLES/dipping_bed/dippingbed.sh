#!/bin/bash 

prjfile='dipping_bed.prj'
mod='s' 


# Run a 2.5D model and output VTI files to be viewed in Paraview 

# Create the source function
sourcefunction -p $prjfile -S gaus0 -m $mod -a 1000

# Run the model 
prjrun $prjfile -M $mod 

# Create the .vti files 
vtkbuild $prjfile -c Vx -f 10 -n 10 
