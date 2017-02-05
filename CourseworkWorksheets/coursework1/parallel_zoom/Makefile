NUM_NODES?=4
HOSTFILE?=../../../machines.txt
PROGNAME=parallel_zoom

default: build

build:
	mpicc mandelbrot_zoom_cyclic.c -I. -o $(PROGNAME)

run: build
	mpiexec --machinefile $(HOSTFILE) -n $(NUM_NODES) $(PROGNAME)

clean:
	@rm -f $(PROGNAME) *.o

.PHONY: clean run