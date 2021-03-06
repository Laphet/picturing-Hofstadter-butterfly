default: PicHofButterfly

MKLPATH=$(MKLROOT)/lib/intel64_lin
MKLINCLUDE=$(MKLROOT)/include
CC=icc
CFLAGS=-Wall -O3 -I$(MKLINCLUDE) -qopenmp
LDFLAGS=-L$(MKLPATH)
LDLIBS=-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

TriMatEigen.o: TriMatEigen.c

TriMatEigen_groups.o: TriMatEigen_groups.c

test: TriMatEigen.o TriMatEigen_groups.o test.c

lapacke_test: lapacke_test.c

PicHofButterfly: PicHofButterfly.c

clean:
	rm -f *.o
