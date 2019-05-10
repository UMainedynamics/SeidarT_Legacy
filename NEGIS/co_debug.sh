#!/bin/bash

# The final position (-F --final) require 3 positional arguments x y z


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

# ================================ Error Checks ===============================
# A quick check to make sure required arguments are given
if [ -z $prjfile ]; then
	echo ERROR: Missing project file
	exit 1
elif [ -z $offset ]; then
	echo ERROR: Missing source to reciever offset
	exit 1
elif [ -z $ds ]; then
	echo ERROR: Missing step value for source
	exit 1
elif [[ -z $xf || $xy || $xz ]]; then
	echo ERROR: Missing final coordinates
	exit 1
else
	echo All input arguments satisfied.
fi



# ================================ Get to work ================================

if [ $seismic ]; then
	xstring='S,x,'
	zstring='S,z,'
else
	xstring='E,x,'
	zstring='E,z,'
fi

# Get the initial reciever location
x=`grep -F "S,x" $prjfile`
x=$(echo $x | cut -f3 -d,)

# sed -i 's/S,x,*.*/S,x,replacedline/' negis1.prj

# Read the line that is 

echo $x

# sed -i 's/$xstring$wild/$xstring$newval/' negis1.prj