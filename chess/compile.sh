#!/bin/bash

echo -n Compiling ... 

gfortran -c utils.f90
mpif90 npv.f90 utils.f90 -o npv.x

echo Done !!!
