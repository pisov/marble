OBJ=marble.o mc.o utils.o globvars.o
CFLAGS=-O3 
INCLUDES=-I/opt/lx/u16.04/sprng/5/openmpi-1.10.3/intel/2017/include
LIBPATH=-L/opt/lx/u16.04/sprng/5/openmpi-1.10.3/intel/2017/lib
LIBS=-lsprng -lstdc++
CFLAGS=-g -traceback

#CFLAGS = -lc -DINTEL -I/opt/lx/u16.04/sprng/5/openmpi-1.10.3/intel/2017/include -L/opt/lx/u16.04/sprng/5/openmpi-1.10.3/intel/2017/lib -lsprng -lstdc++

FC=mpif90
EXE=marble.x

all: marbel.x dump.x

dump.x: dump.o utils.o globvars.o
	$(FC) -o dump.x $(LD_CFLAGS) dump.o utils.o globvars.o
marbel.x: $(OBJ)
	$(FC) -o $(EXE) $(LD_CFLAGS) $(OBJ) -lc -DINTEL $(INCLUDES) $(LIBPATH) $(LIBS)

marble.o: mc.o utils.o globvars.o marble.F90
	$(FC) $(CFLAGS) -c marble.F90 -lc -DINTEL $(INCLUDES) $(LIBPATH) $(LIBS)

dump.o: utils.o globvars.o dump.f90
	$(FC) $(CFLAGS) -c dump.f90

mc.o: globvars.o mc.F90
	$(FC) $(CFLAGS) -c mc.F90 globvars.o -lc -DINTEL $(INCLUDES) $(LIBPATH) $(LIBS)

utils.o: globvars.o utils.f90
	$(FC) $(CFLAGS) -c utils.f90 globvars.o

globvars.o: globvars.f90
	$(FC) $(CFLAGS) -c globvars.f90
	
clean:
	rm -f dump.x *.o *.mod $(EXE)
