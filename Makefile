GC=mpif90
CFLAGS=-warn all -std08 -02 -xHost -openmp
LFLAGS=-openmp
MPIRUN=mpirun
MPIPPN=-npernode

OUT=simul
OUTG=gen

all: gen simul

simul: constants.o io_tools.o mpi_tools.o physics.o main.o
	$(GC) $(LFLAGS) $^ -o $(OUT)

gen: constants.o main_gen.o
	$(GC) $(LFLAGS) $^ -o $(OUTG)

%.o: %.f90
	$(GC) $(CFLAGS) -c $<

clean:
	rm -f *.o *.mod $(OUT) $(OUTG)

sync:
	rsync -ahz --progress *.f90 Makefile mesopsl1.obspm.fr:"~/n-body"

