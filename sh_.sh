#!/bin/bash

#######################################
#Shell script to help with bookkeeping SeidarT runs
#Each analysis is copied to its own timestamped folder
#######################################

#ENTER RUN INFO
projectfolder='/EXAMPLES/fill'
prjfile='fill.prj'
img='fill2.png'
mod='e'
rx='receivers.xyz'
chan='Ex'
proj='co' #co or shot
desc='2m' #comment out or give short description to be part of folder name
dte="$(date +%Y%m%d_%H%M%S)"
stamp=$desc'_'$dte
clean=0 #1 means remove .dat files (i.e., from a single shot)
anim=0 #1 means create animation

offset='1 0 0'
gain=200
exag=0.04

#NO NEED TO EDIT PAST THIS POINT
###########################################################

#COPY FILES AND ENTER THE NEW DIRECTORY
cd $projectfolder
mkdir $stamp
cp $rx $stamp'/'$rx
mv $stamp'/'$rx $stamp'/'receivers.xyz
cp $prjfile $stamp'/'$prjfile
cp $img $stamp'/'$img
cd $stamp

#RUN THE MODEL
prjrun -p $prjfile -m n
sourcefunction -p $prjfile -S gaus0 -m $mod -a 1

if [ $proj == 'co' ]
then
  if [ $mod == 's' ]
  then
  common_offset -p $prjfile -o $offset -r receivers.xyz -s
  else
  common_offset -p $prjfile -o $offset -r receivers.xyz
  fi

else #single shot
prjrun -p $prjfile -m $mod
#  if grep -q 'D,dim,2.5' $prjfile
#  then
#    if [ $mod == 's' ]
#      arraybuild
#    then
#    else
arraybuild -p $prjfile -r receivers.xyz -c $chan

#    fi
fi

# ANIMATIONS
if [ $anim == '1' ]
then
#FOR 2.5D
  if grep -q 'D,dim,2.5' $prjfile
  then
    vtkbuild -p $prjfile -c $chan -n 1
  else
    im2anim  -p $prjfile -c $chan -n 1 -d 1
  fi
fi

#PLOT PROFILE
if [ $proj == 'co' ]
then
  rcxdisplay -p $prjfile -f $chan'.co.csv' -g $gain -e $exag # -r receivers.xyz -t $proj
else
  rcxdisplay -p $prjfile -f 'receiver_array.csv' -g $gain -e $exag # -r receivers.xyz -t $proj
fi

#REMOVE FILES
if [ $clean == '1' ]
then
  rm E*.dat
  rm S*.dat
  rm V*.dat
fi
