## Parallel Scientific Computing Study Notes

Scientific Computing is the study of methods and procedures to approximate
solutions to mathematical problems. Often, we have to approximate solutions
because computing the optimal solution takes too long or is not possible to
due machine floating point inaccuracy.

The notes in this repo explain how to model physical phenomena as mathematical
models. It then models explains how thwy are used to derive numerical models,
which are solved to discover some information about the physical phenomena.

Solving these models takes an enourmous amount of numerical computation. To
speed the computation, the models are often solved in parallel, across
multiple machines. Methods of distributing work across multiple machines is
also covered in these nodes.

### Get the Source

```
git clone https://github.com/DonaldWhyte/notes-parallel-scientific-computing.git
```

### Building the Docs

The notes are written in LaTeX. To build the notes, first install LaTeX and
pdflatex. [**This document**](https://en.wikibooks.org/wiki/LaTeX/Installation)
explains how to install LaTeX on Windows, Mac OS X and various Linux
distributions.

Afterwards, navigate into the cloned and execute the build script:

```
cd notes-parallel-scientific-computing
make docs
```

This will build two PDF documents in the current working directory:

* `parallel-scientific-computing.pdf` -- main document that contains almost all the notes
* `parallel-cw2.pdf` -- document comparing the performance of various numeric solvers (e.g. Jacobi, Gauss-Seidel, SOR)

### Buidling the Test Programs

The test programs are written in C. All the test programs are located in
directories within `CourseworkWorksheets`. To build them, run `make programs`.

### Building Everything

Simply run `make` or `make all` to build all test programs and documents.

### Executing the Test Programs

TODO

