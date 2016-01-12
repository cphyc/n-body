GC=gfortran
CFLAGS=-Wall -Wextra -pedantic -std=f2008 -march=native -O3
LFLAGS=
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
