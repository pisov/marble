OBJ=npv.o  utils.o 
CFLAGS=-O3
#CFLAGS=-g -traceback

FC=mpif90
EXE=npv.x

all: npv.x 

dump.x:  utils.o 
	$(FC) $(LD_CFLAGS) utils.o
npv.x: $(OBJ)
	$(FC) -o $(EXE) $(LD_CFLAGS) $(OBJ)

utils.o: utils.f90
	$(FC) $(CFLAGS) -c utils.f90

npv.o: utils.o  npv.f90
	$(FC) $(CFLAGS) -c npv.f90

	
clean:
	rm -f dump.x *.o *.mod $(EXE)
