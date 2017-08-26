# chess paralelization of marble

In this folder there are diferent subfolders for the diferent implementations of the program. In this implementation the post processing is made after the simulation.

The program must be compiled with a paralel version of the ```intel``` compiler. The defaule compiler is ```mpif90``` and it is set in Makefile with ```FC=mpif90```.

Due to the restrictions of the algorithm the program must be runned on ```(n+2)(m+2)``` number of cores. The size of the grid must also be able to divide properly into the subgrids for each processor.

-------------------------------------------------------------------------------------
## Folder frnd

Parallel version of the program, where vacancies are in seperate array.
 This implementation uses the fortran function random_number for generating random numbers.

-------------------------------------------------------------------------------------
## Folder SPRNG

Parallel version of the program, where vacancies are in seperate array.
 This implementation uses the The Scalable Parallel Random Number Generators Library (SPRNG) for generating random numbers.