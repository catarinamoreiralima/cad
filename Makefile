CC := gcc
MPICC ?= mpicc
CFLAGS := -O3 -Wall -Wextra -std=c11 -fopenmp
LDLIBS := -fopenmp

.PHONY: all clean

all: fire_mpi_omp $(patsubst %.c,%,$(wildcard fire_seq.c fire_omp.c))

fire_seq: fire_seq.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

fire_omp: fire_omp.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

fire_mpi_omp: fire_mpi_omp.c
	$(MPICC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f fire_seq fire_omp fire_mpi_omp
