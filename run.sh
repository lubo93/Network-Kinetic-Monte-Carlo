#!/bin/bash

TIME=50
NMEAS=100
P=0.1
R=10.0
M=10
QInt=1.0
QExt=1.0

INPUT=facebook_combined.txt

for i in {1..1}
do 
	echo "Running with TIME=$TIME NMEAS=$NMEAS P=$P R=$R M=$M QInt=$QInt QExt=$QExt..."
	./dynamicNet $TIME $NMEAS $P $R $M $QInt $QExt $INPUT
done
