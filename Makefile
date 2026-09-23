CC := gcc
CFLAGS := -O3 -Wall -Wextra -std=c11 -fopenmp
LDLIBS := -fopenmp

.PHONY: all clean

all: fire_seq fire_omp

fire_seq: fire_seq.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

fire_omp: fire_omp.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f fire_seq fire_omp
