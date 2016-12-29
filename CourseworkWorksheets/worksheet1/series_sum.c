#include <stdio.h>
#include <mpi.h>

// ADDITIONAL COMMENTS ADDED BY DONALD WHYTE ON 10/10/13

int main( int argc, char **argv )
{
  int i, N, noprocs, nid;
  float sum=0, Gsum;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &nid );
  MPI_Comm_size( MPI_COMM_WORLD, &noprocs );

  // In this case, the process with ID 0 is the MASTER process
  // which asks the user for the number of terms to use.
  // Processing will not start until this is complete.
  if ( nid==0 ) {
    printf( "Please enter the number of terms N -> " );
    scanf( "%d", &N ); // master process blocks
  }

  /* Broadcast the problem size (N) to all processes */
  // Tbhe 0 in this function call states WHICH process does
  // the broadcast and which processes LISTEN for the broadcast.
  // If nid==0, then the process sends N to the other processes
  // If nid!=0, then we BLOCK HERE and wait until the broadcast
  // has been sent!
  // Therefore, broadcast  is a synchronised communication
  // mechanism.
  MPI_Bcast( &N, 1, MPI_INT, 0, MPI_COMM_WORLD );

  // Start from the ID of the process and go up in steps
  // matching the number of processes for a sum interleaved
  // throughout all the processes
  // (NOTE: REVERSED LOOP TO GO BACKWARDS SO SMALLEST TERMS
  // ADDED FIRST)
  for ( i=(N-nid-1); i>=0; i-=noprocs )
    if ( i%2==0 )
      sum += (float) 1 / (i+1);
    else
      sum -= (float) 1 / (i+1);

  /* Add up the sums from each process and store the result in Gsum */
  MPI_Reduce( &sum, &Gsum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD );

  // Master process is also in charge of displaying the result to the user
  // NOTE HOW ALL IO-BOUND PROCESSES ARE USED BY ONE PROCESS TO PREVENT
  // RACE CONDITIONS!!!!
  if ( nid==0 )
    printf( "An estimate of ln(2) is %f \n", Gsum );

  MPI_Finalize();

  return 0;
}
