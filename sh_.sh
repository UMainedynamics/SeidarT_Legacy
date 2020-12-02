#!/bin/bash

#######################################
#Shell script to help with bookkeeping SeidarT runs
#Each analysis is copied to its own timestamped folder
#######################################

#ENTER RUN INFO
projectfolder='/Users/chris/SD/fill'
prjfile='fill.prj'
img='fill2.png'
mod='s'
rx='receivers10.xyz'
chan='Vx'
proj='shot' #co or shot
stamp="$(date +%Y%m%d_%H%M%S)"
clean=0 #1 means remove .dat files (i.e., from a single shot)
anim=1 #1 means create animation

offset='1 0 0'
gain=100
exag=0.01

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
sourcefunction -p $prjfile -S gaus0 -m $mod -a 1

if [ $proj == 'co' ]
then
  if [ $mod == 's' ]
  then
  common_offset -p $prjfile -o $offset -r receivers.xyz -s
  else
  common_offset -p $prjfile -o $offset -r receivers.xyz
  fi
else
prjrun $prjfile -M $mod
arraybuild -p $prjfile -r receivers.xyz -c $chan
fi

# ANIMATIONS
if [ $anim == '1' ]
then
#FOR 2.5D
  if grep -q 'D,dim,2.5' $prjfile;
  then
    vtkbuild $prjfile -c $chan -n 1
  else
    im2anim  $prjfile -c $chan -n 10 -f 10
  fi
fi

#PLOT PROFILE
if [ $proj == 'co' ]
then
  codisplay -p $prjfile -f $chan'.co.csv' -g $gain -e $exag # -r receivers.xyz -t $proj
else
  codisplay -p $prjfile -f 'receiver_array.csv' -g $gain -e $exag # -r receivers.xyz -t $proj
fi

#REMOVE FILES
if [ $clean == '1' ]
then
  rm E*.dat
  rm S*.dat
  rm V*.dat
fi
