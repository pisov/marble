#!/bin/bash

echo -n Compiling ... 

gfortran -O3 -c utils.f90
gfortran -O3 -g npv.f -o npv.x
gfortran -O3 post.f utils.f90 -o post.x


echo Done !!!
