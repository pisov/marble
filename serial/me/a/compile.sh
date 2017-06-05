#!/bin/bash

mkdir -p img
echo -n Compiling ... 

gfortran -c utils.f90
gfortran -g npv.f utils.f90 -o npv.x

echo Done !!!
