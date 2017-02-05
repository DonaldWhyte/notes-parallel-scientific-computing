DOC_NAME = parallel_and_scientific_computing

default: all

all: programs docs

clean: clean-programs clean-docs

programs:
	$(MAKE) -C CourseworkWorksheets

clean-programs:
	$(MAKE) -C CourseworkWorksheets clean

# Build document twice (first time to build TOC, second time to use it).
# Skip the first build if the toc index files have already been generated.
docs:
	@if [ ! -f $(DOC_NAME).toc ] ; \
	then \
	    pdflatex $(DOC_NAME).tex ; \
	fi;
	@pdflatex $(DOC_NAME).tex
	$(MAKE) -C CourseworkWorksheets docs

clean-docs:
	@rm -f *.pdf *.aux *.lof *.log *.lot *.fls *.out *.toc *.fmt *.fot *.cb *.cb2
	$(MAKE) -C CourseworkWorksheets clean-docs

.PHONY: clean clean-programs clean-docs