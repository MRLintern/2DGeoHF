FC=gfortran
FFLAGS=-std=f95 -fall-intrinsics -Wall -Wpedantic -O2 -fopenmp
EXECUTABLE=main
LINKER_OPTIONS=-lm

$(EXECUTABLE): main.o 
	$(FC) main.o -o main

main.o: main.f90
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) main.f90

.PHONY clean:

clean:
	rm -rf *.o *.mod $(EXECUTABLE)
