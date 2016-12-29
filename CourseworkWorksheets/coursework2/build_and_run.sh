#!/bin/sh

mpicc -lm -o $1 "$1.c"
mpiexec -n $2 "./$1"
