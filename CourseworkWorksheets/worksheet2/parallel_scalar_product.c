#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

void masterProcess(int processID, int numProcesses)
{
  FILE *fp;
  int i, n;
  int partSize;
  float* a;
  float* b;
  float chunkSum = 0.0f;
  float globalSum = 0.0f;

  // Read two vectors from file
  fp = fopen( "vecs_small.dat", "rt" );
  fscanf( fp, "%d", &n );
  a = (float *) calloc( n, sizeof(float) );
  b = (float *) calloc( n, sizeof(float) );
  for ( i=0; i<n; i++ )
    fscanf( fp, "%f %f", &a[i], &b[i] );
  fclose( fp );
  // Broadcast problem size to workers
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  partSize = n / numProcesses;
  // Scatter problem to worker processes
  MPI_Scatter(&a, partSize, MPI_FLOAT, &chunkSum, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  // Solve the master's own chunk
  for ( i = 0; (i < partSize); i++)
      chunkSum += a[i] * b[i];
  // Reduce the worker's solutions into a single sum
  MPI_Reduce(&chunkSum, &globalSum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  // Output solution
  printf("The scalar product is %f\n", globalSum);
  // Free resources
  free(a);
  free(b);
}

void workerProcess(int processID, int numProcesses)
{
    int i, n;
    int partSize;
    float* a;
    float* b;
    float chunkSum = 0.0f;
    float globalSum = 0.0f;
    MPI_Status status;

    // Get problem size from master process
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    partSize = n / numProcesses;
    // Wait for the chunk of the problem this process deals with
    a = (float*)calloc(partSize, sizeof(float));
    b = (float*)calloc(partSize, sizeof(float));
    MPI_Scatter(&a[0], partSize, MPI_FLOAT, &chunkSum, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(&b[0], partSize, MPI_FLOAT, &chunkSum, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // Compute solution to chunk
    for (i = 0; (i < partSize); i++)
	chunkSum = a[i] * b[i];
    // Send solution back to master process
    MPI_Reduce(&chunkSum, &globalSum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    // Free resources
    free(a);
    free(b);
}

int main( int argc, char **argv )
{
  int processID;
  int numProcesses;

  // Initialise MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

  // Decide which type of process is currently running 
  if (processID == 0) {
    masterProcess(processID, numProcesses);
  } else {
    workerProcess(processID, numProcesses);
  }

  MPI_Finalize();

  return 0;
}
