GC=gfortran
CFLAGS=-Wall -Wextra -pedantic -std=f2008 -g -O2
LFLAGS=
OUT=simul
OUTG=gen

all: constants.o physics.o initial_conditions.o main.o
	$(GC) $(LFLAGS) $^ -o $(OUT)

gen: constants.o main_gen.o
	$(GC) $(LFLAGS) $^ -o $(OUTG)

%.o: %.f90
	$(GC) $(CFLAGS) -c $<

clean:
	rm *.o *.mod $(OUT) $(OUTG)
