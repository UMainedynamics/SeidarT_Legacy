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
	-o, --offset 	Source and reciever offset distance (meters)
	-d, --delta		Source and reciever step length (meters)
	-m, --model 	(OPTIONAL) Specifier to run model; s-seismic, e-electromag,
					b-both seismic and electromag, n-none (default) 


'

cofile='common_offset.csv' # The output file
wild='*.*'


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
	    -o|--offset)
		offset="$2"
		shift;shift
		;;
		-d|--delta)
		ds="$2"
		shift;shift
	    ;;
	    -s|--seismic)
		seismic=true
		shift
	esac
done


python3 -m prjrun $prjfile --model s
python3 -m arrayplot $prjfile -c Vx -i $rcxi -f $rcxf -d $dr -g 1