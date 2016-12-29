// MODIFICATION OF vector_addittion.c
// Created to handle sitations where the problem size
// (length of vectors to sum) is not divisble by the
// number of available processes

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main( int argc, char **argv )
{
  int noprocs, nid, i, n;
  int size, remainder, startOfLeftover;
  float *a, *b, *c;
  FILE *fp;

  MPI_Status status;

  // Get ID of current process and the number of processes overall
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &nid );
  MPI_Comm_size( MPI_COMM_WORLD, &noprocs );

  // NOTE: remainder method has a problem is # processes larger
  // than the problem size. It just computes EVERYTHING on
  // the master process the remainder is the N (problem size)

  if ( nid==0 ) { // master process -- all IO here
    fp = fopen( "vecs_large.dat", "rt" );
    fscanf( fp, "%d", &n );

    remainder = n % noprocs; // guaranteed to be <= even sized chunk
    startOfLeftover = (n - remainder);

    printf("%d %d %d %d", n, noprocs, remainder, startOfLeftover);

    // Send problem size (WITHOUT REMAINDER) to all processes
    MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD );

    size = startOfLeftover / noprocs; // guaranteed to be integer

    // Reads vectors in from the file
    a = (float *) calloc( n, sizeof(float) );
    b = (float *) calloc( n, sizeof(float) );
    for ( i=0; i<n; i++ )
      fscanf( fp, "%f %f", &a[i], &b[i] );
    fclose( fp );

    // Send a chunk of the array to each worker process
    for ( i=1; i<noprocs; i++ ) {
      MPI_Send( &a[size*i], size, MPI_FLOAT, i, 10, MPI_COMM_WORLD );
      MPI_Send( &b[size*i], size, MPI_FLOAT, i, 20, MPI_COMM_WORLD );
    }

    // Compute our chunk of the problem (master process also does work)
    c = (float *) calloc( n, sizeof(float) );
    for ( i=0; i<size; i++ )
      c[i] = a[i]+b[i];
    // Also compute the leftover chunk that could not be distributed evenly
    for ( i=startOfLeftover; i<n; i++ )
      c[i] = a[i]+b[i];
    // Receive solutions to remaining chunks from worker preocesses
    for ( i=1; i<noprocs; i++ ) {
      MPI_Recv( &c[size*i], size, MPI_FLOAT, i, 30, MPI_COMM_WORLD, &status );
    }
  }
  else { // worker process (as far as they know, the remainder is NOT part of the problem size!
    MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD ); // receive process size from master (DOESN'T INCLUDE REMAINDER)

    size = n/noprocs;
    a = (float *) calloc( size, sizeof(float) );
    b = (float *) calloc( size, sizeof(float) );

    MPI_Recv( &a[0], size, MPI_FLOAT, 0, 10, MPI_COMM_WORLD, &status );
    MPI_Recv( &b[0], size, MPI_FLOAT, 0, 20, MPI_COMM_WORLD, &status );

    // Compute solution to chunk of problem
    c = (float *) calloc( size, sizeof(float) );
    for ( i=0; i<size; i++ )
      c[i] = a[i]+b[i];
    // Send sub-solution back to master process
    MPI_Send( &c[0], size, MPI_FLOAT, 0, 30, MPI_COMM_WORLD );
  }

  // Master process displays results to user
  if ( nid==0 ) {
    printf( "The first 12 elements of the sum of the two vectors is \n" );
    for ( i=0; i<12; i++ ) {
      printf( "%6.2f", c[i] );
    }
    printf( "\n" );
  }

  free( a ); free( b ); free( c );

  MPI_Finalize();

  return 0;
}
