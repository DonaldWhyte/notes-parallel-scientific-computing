NUM_NODES?=4
HOSTFILE?=../../../machines.txt
PROGNAME=hello

default: build

build:
	mpicc series_sum.c -I. -o $(PROGNAME)

run: build
	mpiexec --machinefile $(HOSTFILE) -n $(NUM_NODES) $(PROGNAME)

clean:
	@rm -f $(PROGNAME) *.o

.PHONY: clean run