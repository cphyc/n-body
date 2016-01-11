GC=gfortran
CFLAGS=-Wall -Wextra -pedantic -std=f2008 -g -O2
LFLAGS=
OUT=simul

all: constants.o initial_conditions.o main.o
	$(GC) $^ $(LFLAGS) -o $(OUT)

main.o: main.f90
	$(GC) -c $^ $(CFLAGS)

initial_conditions.o: initial_conditions.f90
	$(GC) -c $^ $(CFLAGS)

constants.o: constants.f90
	$(GC) -c $^ $(CFLAGS)

