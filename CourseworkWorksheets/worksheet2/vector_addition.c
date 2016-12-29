#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main( int argc, char **argv )
{
  int noprocs, nid, i, n, size;
  float *a, *b, *c;
  FILE *fp;

  MPI_Status status;

  // Get ID of current process and the number of processes overall
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &nid );
  MPI_Comm_size( MPI_COMM_WORLD, &noprocs );

  if ( nid==0 ) { // master process -- all IO here
    fp = fopen( "vecs_small.dat", "rt" );
    fscanf( fp, "%d", &n );

    MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD ); // send problem size to all processes

    if ( n%noprocs ) { // double check problem size is divisible by number of processors
      printf( "Number of processes is not a factor of n.\n" );
      MPI_Abort( MPI_COMM_WORLD, -1 );
    }

    // Reads vectors in from the file
    a = (float *) calloc( n, sizeof(float) );
    b = (float *) calloc( n, sizeof(float) );
    for ( i=0; i<n; i++ )
      fscanf( fp, "%f %f", &a[i], &b[i] );
    fclose( fp );

    // Send a chunk of the array to each worker process
    for ( i=1, size=n/noprocs; i<noprocs; i++ ) {
      MPI_Send( &a[size*i], size, MPI_FLOAT, i, 10, MPI_COMM_WORLD );
      MPI_Send( &b[size*i], size, MPI_FLOAT, i, 20, MPI_COMM_WORLD );
    }
 
    // Compute our chunk of the problem (master process also does work)
    c = (float *) calloc( n, sizeof(float) );
    for ( i=0; i<size; i++ )
      c[i] = a[i]+b[i];
    // Receive solutions to remaining chunks from worker preocesses
    for ( i=1; i<noprocs; i++ ) {
      MPI_Recv( &c[size*i], size, MPI_FLOAT, i, 30, MPI_COMM_WORLD, &status );
    }
  }
  else { // worker process
    MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD ); // receive process size from master

    if ( n%noprocs ) {
      printf( "Number of processes is not a factor of n.\n" );
      MPI_Abort( MPI_COMM_WORLD, -1 );
    }

    // Receive chunk of problem to compute from master
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
