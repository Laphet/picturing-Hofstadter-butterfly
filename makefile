default: test

CFLAGS=-g -Wall -fopenmp
LDFLAGS=
LDLIBS=-lm

TriMatEigen.o: TriMatEigen.c

TriMatEigen_omp.o: TriMatEigen_omp.c

test: TriMatEigen.o TriMatEigen_omp.o test.c

clean:
	rm -f *.o
