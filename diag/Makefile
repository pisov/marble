OBJ=diag.o mc.o utils.o globvars.o
CFLAGS=-O3
#CFLAGS=-g -traceback

FC=mpif90
EXE=diag.x

all: diag.x dump.x

dump.x: dump.o utils.o globvars.o
	$(FC) -o dump.x $(LD_CFLAGS) dump.o utils.o globvars.o
diag.x: $(OBJ)
	$(FC) -o $(EXE) $(LD_CFLAGS) $(OBJ)

diag.o: mc.o utils.o globvars.o diag.f90
	$(FC) $(CFLAGS) -c diag.f90

dump.o: utils.o globvars.o dump.f90
	$(FC) $(CFLAGS) -c dump.f90

mc.o: globvars.o mc.f90
	$(FC) $(CFLAGS) -c mc.f90 globvars.o

utils.o: globvars.o utils.f90
	$(FC) $(CFLAGS) -c utils.f90 globvars.o

globvars.o: globvars.f90
	$(FC) $(CFLAGS) -c globvars.f90
	
clean:
	rm -f dump.x *.o *.mod $(EXE)
