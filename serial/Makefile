OBJ=npv1.o
CFLAGS=-O3

FC=gfortran
EXE=npv.x

all: npv.x

npv.x: npv1.o
	$(FC) -o npv.x $(LD_CFLAGS) npv1.o

npv1.o: npv1.f
	$(FC) $(CFLAGS) -c npv1.f
	
clean:
	rm -f *.x *.o *.mod $(EXE)
