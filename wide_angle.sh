#!/bin/bash

: '
WIDE_ANGLE.SH is a wrapper script to create a common offset survey using the
SeidarT programs. The survey is along the x-direction, but can be extended to 
other directions.

INPUT
	-f, --project 	Project file path
	-I, --initial 	The inital coordinates of the recievers. Requires 3 
					arguments x y z
	-F, --final		The final coordinates of the source. Requires 3
					arguments  x y z
	-d, --delta		Reciever spacing (meters)
	-m, --model 	(OPTIONAL) Specifier to run model; s-seismic, e-electromag,
					b-both seismic and electromag, n-none (default) 

'

# ------------------------------ Input Arguments ------------------------------
while [[ $# -gt 0 ]]; do
	key="$1"

	case $key in
	    -f|--project_file)
	    prjfile="$2"
	    shift # past argument
	    shift # past value
	    ;;
	    -F|--final) 
	    xf="$2"
	    yf="$3"
	    zf="$4"
	    shift;shift;shift;shift # past value
	    ;;
	    -I|--initial) 
	    xi="$2"
	    yi="$3"
	    zi="$4"
	    shift;shift;shift;shift # past value
	    ;;
		-d|--delta)
		ds="$2"
		shift;shift
	    ;;
	    -g|--gif)
		gif=true
		shift
		;;
	    -s|--seismic)
		seismic=true
		shift
	esac
done


# ------------------------------- Check Inputs --------------------------------

if [[ -z $prjfile ]]; then
	echo ERROR: Missing project file
	exit 1
elif [[ -z $xi || -z $yi || -z $zi ]]; then
	echo ERROR: Missing initial reciever location
	exit 1
elif [[ -z $xf || -z $yf || -z $zf ]]; then
	echo ERROR: Missing final reciever location
	exit 1
elif [[ -z $ds ]]; then
	echo ERROR: Missing reciever spacing
	exit 1
else
	echo
fi

echo $gif

# ------------------------------- Run the Model -------------------------------
if [[ $seismic ]]; then
	c1=Vx
	c2=Vz
	mod=s
	gain=7
	echo Modeling seismic wave propagation
else
	c1=Ex
	c2=Ez
	mod=e
	gain=7
	echo Modeling electromagnetic wave propagation
fi


python3 -m prjrun $prjfile --model $mod

echo Creating reciever array
python3 -m arrayplot $prjfile -c $c1 -I $xi $yi $zi -F $xf $yf $zf -d $ds -g $gain
python3 -m arrayplot $prjfile -c $c2 -I $xi $yi $zi -F $xf $yf $zf -d $ds -g $gain


if [[ $gif ]]; then
	python3 -m im2gif $prjfile -c $c1 -f 30
	python3 -m im2gif $prjfile -c $c2 -f 30
fi