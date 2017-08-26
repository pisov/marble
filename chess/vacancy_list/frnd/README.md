# Chess parallelization with vacancy list

Parallel version of the program, where vacancies are in seperate array.
 This implementation uses the fortran function random_number for generating random numbers.

## Compilation and execution

The program must be compiled with a paralel version of the ```intel``` compiler. The defaule compiler is ```mpif90``` and it is set in Makefile with ```FC=mpif90```.

To compile

```
make clean
make
```
To execute the program 

```
mpirun -np N marble.x

```

Where N is the number of processors.

Due to the restrictions of the algorithm the program must be runned on ```(n+2)(m+2)``` number of cores. The size of the grid must also be able to divide properly into the subgrids for each processor.

## Configuration the initial state of the system

This can be made by edditing the ```marble.in``` file. 

```nit``` - number of iterations <br />
```nout``` - the number of iterations on which a checkpoint file is made <br />
```nsize``` - size of the grid <br />
```pbvac``` - the percentage content of the vacancies in the system <br />
```pbratio``` - the ratio of the two substances in the system <br />
```pbfrmv``` -  the probability of exclusion; ```(1-pbfrmv)``` - probability of diffusion <br />

## Postprocessing

The post processing is made via the program ```dump.x```. The program takes as an argument the checkpoint file of the simulation. The program reads the ```marble.in``` file in order to get the initial state of the system, so it must be in the folder.

```
./dump.x marble.chkp-XXXXXX
```


## Other 

```post.sh``` - example script for making the post processing

```sim.sh``` - example job file for making a simulation where the percentage content of the vacancies in the system is parametrized.