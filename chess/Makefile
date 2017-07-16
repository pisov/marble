OBJ=npv.o  utils.o 
CFLAGS=-O3
#CFLAGS=-g -traceback

FC=mpif90
EXE=npv.x

all: marbel.x 

dump.x:  utils.o 
	$(FC) $(LD_CFLAGS) utils.o
marbel.x: $(OBJ)
	$(FC) -o $(EXE) $(LD_CFLAGS) $(OBJ)

utils.o: utils.f90
	$(FC) $(CFLAGS) -c utils.f90

npv.o: utils.o  npv.f90
	$(FC) $(CFLAGS) -c npv.f90


#mc.o: globvars.o mc.f90
#	$(FC) $(CFLAGS) -c mc.f90 globvars.o

#utils.o: globvars.o utils.f90
#	$(FC) $(CFLAGS) -c utils.f90 globvars.o

#globvars.o: globvars.f90
#	$(FC) $(CFLAGS) -c globvars.f90
	
clean:
	rm -f *.o *.mod $(EXE)
