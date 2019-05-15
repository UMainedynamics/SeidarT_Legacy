#!/bin/bash

: '
COMMON_MIDPOINT.SH is a wrapper script to create a common midpoint survey using 
the SeidarT programs. The survey is along the x-direction, but can be extended 
to other directions. The midpoint is initialized from the x, y, and z locations
provided from the source location in the project file.

INPUT
	-f, --project 	Project file path
	-t, --total		The terminal distance between the source and reciever. 
	-o, --offset 	The initial source and reciever offset from the midpoint
					given in (+/- meters). A negative value means that the 
					source is on the lookers left of the midpoint. The total 
					source and reciever distance is 2*offset. 
	-d, --delta		Source and reciever step length (meters); total distance
					between the source and reciever is 2*delta*i + 2*offset.
	-s, --seismic 	(OPTIONAL) Specifier to run seismic common offset 

'

cmpfilex='common_midpoint_x.csv' # The output file
cmpfilez='common_midpoint_z.csv' # The output file
wild='.*'


while [[ $# -gt 0 ]]; do
	key="$1"

	case $key in
	    -f|--project_file)
	    prjfile="$2"
	    shift # past argument
	    shift # past value
	    ;;
	    -t|--total)
		total=$2
		shift;shift
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
elif [[ -z $total ]]; then
	echo ERROR: Missing terminal offset distance
	exit 1
elif [ -z $offset ]; then
	echo ERROR: Missing source to reciever offset
	exit 1
elif [ -z $ds ]; then
	echo ERROR: Missing step value for source
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
	mod=s
else
	xstring='E,x,'
	ystring='E,y,'
	zstring='E,z,'
	c1='Ex'
	c2='Ez'
	mod=e
fi

# Create the common midpoint file
touch $cmpfilex
touch $cmpfilez

# Get the midpoint from the project file
xs=`grep -F $xstring $prjfile`
xs=$(echo $xs | cut -f3 -d,)

ys=`grep -F $ystring $prjfile`
ys=$(echo $ys | cut -f3 -d,)

zs=`grep -F $zstring $prjfile`
zs=$(echo $zs | cut -f3 -d,)


# For now we aren't varying the other dimensions
yr=$ys
zr=$zs

# We need to save the original value so that we can rewrite it at the end
xorig=$xs

# Check if the midpoint is given as negative or positive
if [[ $xs -lt 0 ]]; then
	echo Source is on left of midpoint
	xs=`echo "$xs * -1" | bc `
	ds=`echo "$ds * -1" | bc `
	offset=`echo "$offset * -1" | bc `
	sign=-1
else
	echo Source is on right of midpoint
	sign=1
fi

# Initialize the source and reciever locations
xr=`echo "$xs - $offset" | bc` 
xs=`echo "$xs + $offset" | bc`
total_offset=`echo "$sign*($xs - $xr)" | bc`

# Rewrite the midpoint with the source location because prjrun requires the 
# source location
# sed -i "s/$xstring$wild/$xstring$xs/" $prjfile


while [[ total_offset -lt total ]]; do

	python3 -m prjrun $prjfile --model $mod

	# Get the reciever timeseries for the x-direction
	python3 -m arrayplot $prjfile -c $c1 -I $xr $yr $zr -F $xr $yr $zr -d 1 -g 0 -S 1

	# append the timeseries to the others
	paste -d' ' $cmpfilex reciever_array.csv > temp.csv
	mv temp.csv $cmpfilex
	rm reciever_array.csv

	# Get the reciever timeseries for the z-direction
	python3 -m arrayplot $prjfile -c $c2 -I $xr $yr $zr -F $xr $yr $zr -d 1 -g 0 -S 1

	# append the timeseries to the others
	paste -d' ' $cmpfilez reciever_array.csv > temp.csv
	mv temp.csv $cmpfilez
	rm reciever_array.csv

	# Shift the source
	xs=`echo "$xs + $ds" | bc`
	sed -i "s/$xstring$wild/$xstring$xs/" $prjfile
	# Shift the reciever
	xr=`echo "$xr - $ds" | bc`

	total_offset=`echo "$sign*($xs - $xr)" | bc`
	echo $xr $xs $total_offset
done

sed -i "s/$xstring$wild/$xstring$xorig/" $prjfile


# Display the results for the Ex field. We can do the same for Ez if we like
python3 -m codisplay $prjfile -s $cmpfilex -d $ds -m e

