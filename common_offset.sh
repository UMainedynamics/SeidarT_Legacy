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

cofilex='common_offset_x.csv' # The output file
cofilez='common_offset_z.csv' # The output file
wild='.*'


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
echo $seismic

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

# Create the common offset file
touch $cofilex
touch $cofilez

# Get the initial reciever location
x=`grep -F $xstring $prjfile`
x=$(echo $x | cut -f3 -d,)
xorig=$x

y=`grep -F $ystring $prjfile`
y=$(echo $y | cut -f3 -d,)

z=`grep -F $zstring $prjfile`
z=$(echo $z | cut -f3 -d,)


ry=$y
rz=$z	

# The reciever will trail the source relative to the direction the source is 
# moving. Negative direction is viewers left
if [ $x -lt $xf ]; then # Moving right
	rx=`echo "$x - $offset" | bc` 

	while [ $x -lt $xf ]; do

		python3 -m prjrun $prjfile --model $mod

		# Get the reciever timeseries for the x-direction
		python3 -m arrayplot $prjfile -c $c1 -I $rx $ry $rz -F $rx $ry $rz -d 1 -g 0 -S 1

		# append the timeseries to the others
		paste -d' ' $cofilex reciever_array.csv > temp.csv
		mv temp.csv $cofilex
		rm reciever_array.csv

		# Get the reciever timeseries for the z-direction
		python3 -m arrayplot $prjfile -c $c2 -I $rx $ry $rz -F $rx $ry $rz -d 1 -g 0 -S 1

		# append the timeseries to the others
		paste -d' ' $cofilez reciever_array.csv > temp.csv
		mv temp.csv $cofilez
		rm reciever_array.csv

		# Shift the source
		x=`echo "$x + $ds" | bc`
		sed -i "s/$xstring$wild/$xstring$x/" $prjfile
		# Shift the reciever
		rx=`echo "$rx + $ds" | bc`

	done

else # Moving left
	rx=`echo "$x + $offset" | bc` 
	
	while [ $x -gt $xf ]; do

		python3 -m prjrun $prjfile --model $mod

		# Get the reciever timeseries for the x-direction
		python3 -m arrayplot $prjfile -c $c1 -I $rx $ry $rz -F $rx $ry $rz -d 1 -g 0 -S 1

		# append the timeseries to the others
		paste -d' ' $cofilex reciever_array.csv > temp.csv
		mv temp.csv $cofilex
		rm reciever_array.csv

		# Get the reciever timeseries for the z-direction
		python3 -m arrayplot $prjfile -c $c2 -I $rx $ry $rz -F $rx $ry $rz -d 1 -g 0 -S 1

		# append the timeseries to the others
		paste -d' ' $cofilez reciever_array.csv > temp.csv
		mv temp.csv $cofilez
		rm reciever_array.csv

		# Shift the source
		x=`echo "$x - $ds" | bc`
		sed -i "s/$xstring.*/$xstring$x/" $prjfile
		
		# Shift the reciever
		rx=`echo "$rx - $ds" | bc`

		# remove the .dat files
		rm *.dat
	done
fi


# Display the results for the Ex field. We can do the same for Ez if we like
python3 -m codisplay $prjfile -s $cofilex -d $ds -m e

