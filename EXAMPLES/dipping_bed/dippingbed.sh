#!/bin/bash 

prjfile='dipping_bed.prj'
mod='s' 
rcxfile='receivers.xyz'

# Run a 2.5D model and output VTI files to be viewed in Paraview 

# Create the source function
sourcefunction -p $prjfile -S gaus0 -m $mod -a 1000

# Run the model 
prjrun $prjfile -M $mod 

# Create the .vti files 
vtkbuild $prjfile -c Vx -f 10 -n 10 

# plot the array of recievers
# The receiver file is in indices
# arrayplot -p $prjfile -c Vx -r $rcxfile -g 101 -e 0.1 -S 0 -i 1
# mv reciever_array.csv reciever_array_Vx.csv

# arrayplot -p $prjfile -c Vy -r $rcxfile -g 101 -e 0.1 -S 0 -i 1
# mv reciever_array.csv reciever_array_Vy.csv

# arrayplot -p $prjfile -c Vz -r $rcxfile -g 101 -e 0.025 -S 0 -i 1
# mv reciever_array.csv reciever_array_Vz.csv
