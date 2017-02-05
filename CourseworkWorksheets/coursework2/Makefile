NUM_NODES?=4
HOSTFILE?=../../../machines.txt
PROGNAME?=specify_program_using_PROGNAME_variable
DOC_NAME=parallel-cw2

default: build

build:
	mpicc heat_blocks.c -I. -o heat_blocks
	mpicc heat_blocks_fixed_iter.c -I. -o heat_blocks_fixed_iter
	mpicc heat_rows.c -I. -o heat_rows
	mpicc heat_rows_fixed_iter.c -I. -o heat_rows_fixed_iter
	mpicc jacobi.c -I. -o jacobi
	mpicc lexicographic_gauss.c -I. -o lexicographic_gauss
	mpicc redblack_gauss.c -I. -o redblack_gauss
	mpicc sor.c -I. -o sor

run: build
	mpiexec --machinefile $(HOSTFILE) -n $(NUM_NODES) $(PROGNAME)

clean:
	@rm -f *.o heat_blocks heat_blocks_fixed_iter heat_rows \
		heat_rows_fixed_iter jacobi lexicographic_gauss redblack_gauss sor

# Build document twice (first time to build TOC, second time to use it).
# Skip the first build if the toc index files have already been generated.
docs:
	@if [ ! -f $(DOC_NAME).toc ] ; \
	then \
	    pdflatex $(DOC_NAME).tex ; \
	fi;
	@pdflatex $(DOC_NAME).tex

clean-docs:
	@rm -f $(DOC_NAME).pdf *.aux *.lof *.log *.lot *.fls *.out *.toc *.fmt \
		*.fot *.cb *.cb2

.PHONY: clean run clean-docs