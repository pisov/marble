# Chess parallelization without vacancy list

Parallel version of the program, where vacancies are NOT in seperate array.
 This implementation uses the fortran function random_number for generating random numbers.

## Compilation and execution

The defaule compiler is ```mpif90``` and it is set in Makefile with ```FC=mpif90```

To compile

```
make clean
make
```
To execute the program 

```
mpirun -np N marble.x
```

## Configuration the initial state of the system

The configuration the initial state of the system is hardcoded.
