## Parallel Scientific Computing Study Notes

Scientific Computing is the study of methods and procedures to approximate
solutions to mathematical problems. Often, we have to approximate solutions
because computing the optimal solution takes too long or is not possible to
due machine floating point inaccuracy.

The notes in this repo explain how to model physical phenomena as mathematical
models. It then models explains how they are used to derive numerical models,
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

The programs use OpenMPI for parallelism. This means they use `mpicc` for
building and `mpiexec` to execute the programs across multiple machines.

You'll need to install [OpenMPI](https://www.open-mpi.org/) to build and run
the programs. Instructions on how to install OpenMPI for some platforms is
specified below. If your platform is not lsited, then check out the
[OpenMPI](https://www.open-mpi.org/) website for instructions on how to
install OpenMPI on your platform.

#### Mac OS X

```
brew install open-mpi
```

#### CentOS

```
yum install openmpi openmpi-devel environment-modules
module add openmpi-x86_64
```

### Building Everything

Simply run `make` or `make all` to build all test programs and documents.

To just build the test programs, run `make programs`.

### Executing the Test Programs

The distributed test programs use `mpiexec`, since several processes are
executed across multiple machines. `mpiexec` takes a text file that contains
the hostnames of all hosts that can be used to run the distributed programs.
The text file has the format:

```
hostname1
hostname2
...
hostname4
```

This repo's makefiles pass `mpiexec` the path to a host list file. This is
`machines.txt` at the root fo this repository. You must add the hostnames of
all your worker machinres to `machines.txt` to run any of the test programs.

After filling in `machines.txt`, navigate to the directory that contains the program you want to run (e.g. `CourseworkWorksheets/worksheet0`).

To run that program with the default number of nodes (4), then simply type
`make run`. To specify the number of nodes/processes to use, type:

```
make run NUM_NODES=<numberOfNodes>
```

##### Overriding Machine List File

You can also override the default machines list file path to use a file at a
different location. This can be done using the `HOSTFILE` argument, like so:

```
make run HOSTFILE=<pathToMachinesListTextFile>
```

### Acknowledgements

Credit to University of Leeds and Matthew Hubbard for doing a great job
teaching me the fundamentals of parallel scientific computation.
