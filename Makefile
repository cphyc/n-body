# Configuration for gfortran
#GC=gfortran
#CFLAGS=-Wall -Wextra -pedantic -std=f2008 -O2 -march=native -fopenmp
#LFLAGS=-fopenmp

# Configuration for mpifort
GC=mpifort
CFLAGS=-Wall -Wextra -pedantic -std=f2008 -O3 -march=native -fopenmp
LFLAGS=-fopenmp
MPIRUN=mpirun

# Configuration for mpif90
#GC=mpif90
#CFLAGS=-warn all -std08 -02 -xHost -openmp
#LFLAGS=-openmp
#MPIRUN=mpirun

# Configuration for ifort
#GC=ifort
#CFLAGS=-warn all -std08 -02 -xHost -openmp
#LFLAGS=-openmp

# Configuration for mpiifort
#GC=mpiifort
#CFLAGS=-warn all -std08 -02 -xHost -openmp
#LFLAGS=-openmp
#MPIRUN=mpiexec.hydra

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

watch:
	bash autocompile.sh

sync:
	rsync -ahz --progress *{.f90,.sh,.slurm} mesopsl1.obspm.fr:"~/n-body"

run: all
	./$(OUTG)
	/usr/bin/time -f "Total: %es User: %Us System: %Ss CPU: %P" ./$(OUT)

runmpi: all
	./$(OUTG)
	/usr/bin/time -f "Total: %es User: %Us System: %Ss CPU: %P" \
		$(MPIRUN) -n $(MPI_PROC) --npernode $(MPI_PROC_PER_NODE) ./$(OUT)
