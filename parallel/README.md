#Readme

How to build? You need parallel Fortran 90 compiler: mpif90

1. Compile the code
```
make clean
make
```
2. Run
```
mpirun -np numproc marble.x
```
