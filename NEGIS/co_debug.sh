#!/bin/bash

: '
COMMON_OFFSET.SH is a wrapper script to create a common offset survey using the
SeidarT programs. The survey is along the x-direction, but can be extended to 
other directions.

INPUT
	-f, --project 	Project file path
	-F, --final		The final coordinates of the source. Requires 3
					arguments  x y z
	-o, --offset 	Source and reciever offset distance (meters)
	-d, --delta		Source and reciever step length (meters)
	-s, --seismic 	(OPTIONAL) Specifier to run seismic common offset 


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
	echo 
fi

# ================================ Get to work ================================

if [ $seismic ]; then
	xstring='S,x,'
	ystring='S,y,'
	zstring='S,z,'
	c1='Vx'
	c2='Vz'
else
	xstring='E,x,'
	ystring='E,y,'
	zstring='E,z,'
	c1='Ex'
	c2='Ez'
fi

# Get the initial reciever location
x=`grep -F $xstring $prjfile`
x=$(echo $x | cut -f3 -d,)

y=`grep -F $ystring $prjfile`
y=$(echo $x | cut -f3 -d,)

z=`grep -F $zstring $prjfile`
z=$(echo $x | cut -f3 -d,)


# Create the initial reciever location
rx=`echo "$x + $offset" | bc`
ry=$y
rz=$z

# Create the common offset file
touch $cofile

while [ $rx -lt $xf ]; do
	# Run the first simulation
	python3 -m prjrun $prjfile --model s

	# Get the reciever timeseries
	python3 -m arrayplot $prjfile -c $c1 -i $rx $ry $rz -f $rx $ry $rz -d 1 -g 0

	# append the timeseries to the others
	paste -d' ' $cofile reciever_array.csv > temp.csv
	mv temp.csv $cofile

	# Shift the source
	newsrc=`echo "$x + $ds" | bc`
	sed -i 's/$xstring$wild/$xstring$newsrc/' negis1.prj

	# Shift the reciever
	rx=`echo "$rx + $ds" | bc`
	
done

