NUM_NODES?=4
HOSTFILE?=../../../machines.txt
PROGNAME?=specify_program_using_PROGNAME_variable

default: build

build:
	mpicc parallel_scalar_product.c -I. -o parallel_scalar_product
	mpicc serial_scalar_product.c -I. -o serial_scalar_product
	mpicc vector_addition.c -I. -o vector_addition
	mpicc vector_addition_remainder.c -I. -o vector_addition_remainder

run: build
	mpiexec --machinefile $(HOSTFILE) -n $(NUM_NODES) $(PROGNAME)

clean:
	@rm -f $(PROGNAME) *.o

.PHONY: clean run