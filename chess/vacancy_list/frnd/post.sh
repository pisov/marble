#!/bin/bash


#check if the output directories are made if not make them
#mkdir -p img
mkdir -p postdata

rm nuclu.dat
#run the program multiple times for the different steps 
#for "$(ls marble.ch*)" in $currentdir
for   file in marble.ch*
do
#	echo $file
	./dump.x $file >> nuclu.dat
#	./dump.x $file 
done;


