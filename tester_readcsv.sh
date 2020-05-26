
#!/bin/bash

# Tester to read the csv file 

INPUT=receivers.xyz

sed 1d receivers.xyz | while IFS=, read X Y Z; 
do 
    echo "$Z and $X $Y"; 
done
