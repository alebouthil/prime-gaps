#!/bin/bash

echo 'Compiling'
rm output.txt > /dev/null # delete old values
touch output.txt

# Compile program from source
mpicc -o prime-gaps prime-gaps.c -lgmp

# Test with all needed processor numbers
for i in 2 3 4 5 6 7 8
do
	echo "running with $i processors"
	# Output to textfile for graphing
	mpirun -np $i prime-gaps >> output.txt
done

echo 'Complete, testing output available in output.txt'

python3 plot.py
