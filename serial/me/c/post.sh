#!/bin/bash

currentdir="$(pwd)"
subdir="/data/*"

#check if the output directories are made if not make them
mkdir -p img
mkdir -p postdata

#run the program multiple times for the different steps 
for filename in $currentdir$subdir
do
	./post.x $filename
done;


