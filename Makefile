GC=gfortran
#GC=ifort
CFLAGS=-Wall -Wextra -pedantic -std=f2008 -O3 -march=native -fopenmp
#CFLAGS=-warn all -std08 -02 -xHost
LFLAGS=-fopenmp
#LFLAGS=-fopenmp
OUT=simul
OUTG=gen

all: constants.o physics.o io_tools.o main.o
	$(GC) $(LFLAGS) $^ -o $(OUT)

gen: constants.o main_gen.o
	$(GC) $(LFLAGS) $^ -o $(OUTG)

%.o: %.f90
	$(GC) $(CFLAGS) -c $<

clean:
	rm -f *.o *.mod $(OUT) $(OUTG)

watch:
	bash autocompile.sh
