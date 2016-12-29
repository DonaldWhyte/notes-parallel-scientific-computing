/*
 * Parallel implementation of 2D Laplace equation solver
 *   -- Jacobi iteration
 *   -- partitioning done by blocks
 *   -- non-constant boundary conditions specified
 *   -- exact solution compared with
 *
 * Code supplied by Peter Jimack and modified by Matthew Hubbard */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

/* Define parameter values */
#define N        961
#define MAX_ITER 8000

/* Define the geometric bounds of the domain */
#define xmin -2.0
#define xmax  2.0
#define ymin -2.0
#define ymax  2.0


/* Allocate memory for a two-dimensional array (pointers to pointers) */
double **matrix( int m, int n )
{
  int i;
  double **ptr;

  ptr = (double **) calloc( m, sizeof(double *) );
  for ( i=0; i<m; i++ )
    ptr[i] = (double *) calloc( n, sizeof(double) );

  return (ptr);
}



/* Carry out a single Jacobi iteration on the specified range of mesh nodes */
double iteration( double **old, double **new, int first_row, int last_row, int first_col, int last_col )
{
  double diff, maxdiff=0.0;
  int i, j;

  for ( i=first_row; i<=last_row; i++ )
    for ( j=first_col; j<=last_col; j++ ) {
      /* The Jacobi update for node (i,j) */
      new[i][j] = 0.25 * ( old[i+1][j] + old[i-1][j] + old[i][j+1] + old[i][j-1] );

      /* Update the max norm of  x^(k+1) - x^(k) */
      diff = new[i][j] - old[i][j];
      if ( diff<0 )
        diff = -diff;
      if ( maxdiff<diff )
        maxdiff = diff;
    }

  /* Return the max norm of  x^(k+1) - x^(k)  for the specified range of nodes */
  return (maxdiff);
}



/* Calculate the exact solution at point (x,y) */
double exact( double x, double y )
{
  double solution;

  solution = 1.0;

  return (solution);
}



/* Calculate the max norm error in the final approximation on the specified range of
 * mesh nodes                                                                        */
double final_error( double **new, int start_row, int nrows, int start_col, int ncols, double dx, double dy )
{
  double x, y, diff, maxerr=0.0;
  int i, j;

  for ( i=1; i<=nrows; i++ ) {
    for ( j=1; j<=ncols; j++ ) {
      /* Calculate the geometric coordinates of point (i,j) */
      x = xmin + dx * (double) (start_row+i-1);
      y = ymin + dy * (double) (start_col+j-1);

      /* Update the max norm of  approximate - exact */
      diff = new[i][j] - exact( x, y );
      if ( diff<0 )
        diff = -diff;
      if ( maxerr<diff )
        maxerr = diff;
    }
  }

  /* Return the max norm of  approximate - exact  for the specified range of nodes */
  return (maxerr);
}



/* Write the approximate solution out to a single data file from process 0
 *   -- the format is row-by-row                                           */
void write_file( int nid, double **new, int rowmin, int rowmax, int colmin, int colmax, int first_row, int first_col, int noprocs )
{
  double **full;
  double *send_buffer, *recv_buffer;
  int row, col, counter, i, j, ip;
  int imin, imax, jmin, jmax;
  int send_indices[4], recv_indices[4];
  int send_size, recv_size;
  FILE *fp;
  MPI_Status status;


  /* Fix local row and column indcies */
  imin = rowmin - first_row + 1;
  imax = rowmax - first_row + 1;
  jmin = colmin - first_col + 1;
  jmax = colmax - first_col + 1;


  if ( nid==0 ) {
    /* Allocate memory on process 0 to store the full solution */
    full = matrix( N+1, N+1 );

    /* Insert the solution values calculated on process 0 in to the full array */
    for ( row=rowmin, i=imin; row<rowmax+1; row++, i++ )
      for ( col=colmin, j=jmin; col<colmax+1; col++, j++ )
        full[row][col] = new[i][j];

    /* Receive solution values from the other processes and enter them in to the full array */
    for ( ip=1; ip<noprocs; ip++ ) {
      /* Receive the indices indicating the solution values calculated on process ip */
      MPI_Recv( &recv_indices[0], 4, MPI_INT, ip, 25, MPI_COMM_WORLD, &status );

      rowmin = recv_indices[0]; rowmax = recv_indices[1];
      colmin = recv_indices[2]; colmax = recv_indices[3];

      /* Calculate the number of solution values to be sent by process ip */
      recv_size   = (rowmax-rowmin+1)*(colmax-colmin+1);
      recv_buffer = (double *) calloc( recv_size, sizeof(double) );

      /* Receive the solution values calculated on process ip */
      MPI_Recv( &recv_buffer[0], recv_size, MPI_DOUBLE, ip, 15, MPI_COMM_WORLD, &status );

      /* Place the solution values in the message in to the full array, in order */
      counter = 0;
      for ( row=rowmin; row<rowmax+1; row++ )
        for ( col=colmin; col<colmax+1; col++ ) {
          full[row][col] = recv_buffer[counter];
          counter++;
        }

      free( recv_buffer );
    }

    /* Write out the full array to  heat_solution.dat  row by row */
    fp = fopen( "heat_solution.dat", "wt" );

    for ( row=0; row<N+1; row++ )
      for ( col=0; col<N+1; col++ )
        fprintf( fp, "%6d %6d %8.6f\n", row, col, full[row][col] );

    fclose( fp );

    free( full );
  }
  else {
    send_indices[0] = rowmin; send_indices[1] = rowmax;
    send_indices[2] = colmin; send_indices[3] = colmax;

    /* Send the indices indicating the solution values to be sent to process 0 */
    MPI_Send( &send_indices[0], 4, MPI_INT, 0, 25, MPI_COMM_WORLD );

    /* Calculate the number of solution values to be sent to process 0 */
    send_size   = (imax-imin+1)*(jmax-jmin+1);
    send_buffer = (double *) calloc( send_size, sizeof(double) );

    /* Place the solution values to be sent in to the message in a specific order */
    counter = 0;
    for ( i=imin; i<imax+1; i++ )
      for ( j=jmin; j<jmax+1; j++ ) {
        send_buffer[counter] = new[i][j];
        counter++;
      }

    /* Send the solution values calculated on this process to process 0 */
    MPI_Send( &send_buffer[0], send_size, MPI_DOUBLE, 0, 15, MPI_COMM_WORLD );

    free( send_buffer );
  }
}



/* Main program */
int main( int argc, char **argv )
{
  double **new, **old, **tmp;
  double *send_columnS, *recv_columnS;
  double *send_columnN, *recv_columnN;
  double diff, MaxDiff, MaxDiffG, MaxErr, MaxErrG;
  double dx, dy, x, y;
  double start_time=0.0, end_time=0.0;
  int first_row, nrows, first_col, ncols;
  int noprocs, nid, i, j, iter;
  int northid, southid, westid, eastid;
  int sqrtprocs, errorcode=0;
  int rowmin, rowmax, colmin, colmax;
  MPI_Status status;
  MPI_Request req_sendN, req_sendS, req_sendW, req_sendE;
  MPI_Request req_recvN, req_recvS, req_recvW, req_recvE;


  /* Initialise for MPI */
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &nid );
  MPI_Comm_size( MPI_COMM_WORLD, &noprocs );


  /* Calculate the mesh size */
  dx = ( xmax - xmin ) / (double) N;
  dy = ( ymax - ymin ) / (double) N;


  sqrtprocs = (int) sqrt(noprocs);
  if ( nid==0 ) {
    /* Check that the number of processors is a perfect square */
    if ( noprocs!=sqrtprocs*sqrtprocs ) {
      fprintf( stdout, "Quitting...number of MPI tasks not a perfect square.\n" );
      MPI_Abort( MPI_COMM_WORLD, errorcode );
      exit( 0 );
    }
    /* Check that N-1 is divisible by sqrtprocs */
    if ( (N-1)%sqrtprocs!=0 ) {
      fprintf( stdout, "Quitting...N-1 not divisible by sqrtprocs.\n" );
      MPI_Abort( MPI_COMM_WORLD, errorcode );
      exit( 0 );
    }
  }


  /* Calculate the first row and the number of rows allocated to this process */
  nrows     = (N-1)/sqrtprocs;
  first_row = 1 + (nid%sqrtprocs)*nrows;

  /* Calculate the first column and the number of columns allocated to this process */
  ncols     = (N-1)/sqrtprocs;
  first_col = 1 + (nid/sqrtprocs)*ncols;


  /* Allocate memory for the local work arrays
   *   -- note that an extra layer of nodes is required all round the updated region */
  new = matrix( nrows+2, ncols+2 );
  old = matrix( nrows+2, ncols+2 );

  /* Allocate memory for send and recv buffers */
  send_columnS = calloc( nrows, sizeof(double) );
  recv_columnS = calloc( nrows, sizeof(double) );
  send_columnN = calloc( nrows, sizeof(double) );
  recv_columnN = calloc( nrows, sizeof(double) );


  /* Set up ID numbers for the neighbouring processors (-1 at a boundary) */
  westid  = nid - 1;
  if ( nid%sqrtprocs==0 )
    westid  = -1;
  eastid  = nid + 1;
  if ( (nid+1)%sqrtprocs==0 )
    eastid  = -1;
  southid = nid - sqrtprocs;
  if ( nid<sqrtprocs )
    southid = -1;
  northid = nid + sqrtprocs;
  if ( nid>=noprocs-sqrtprocs )
    northid = -1;

  if ( nid==0 )
    start_time = MPI_Wtime();


  /* Apply boundary conditions at the xmin boundary */
  if ( westid==-1 ) {
    for ( j=0; j<=ncols+1; j++ ) {
      x = xmin;
      y = ymin + dy * (double) (first_col+j-1);
      new[0][j]       = old[0][j]       = exact( x, y );
    }
  }

  /* Apply boundary conditions at the xmax boundary */
  if ( eastid==-1 )
    for ( j=0; j<=ncols+1; j++ ) {
      x = xmax;
      y = ymin + dy * (double) (first_col+j-1);
      new[nrows+1][j] = old[nrows+1][j] = exact( x, y );
    }

  /* Apply boundary conditions at the ymin boundary */
  if ( southid==-1 )
    for ( i=0; i<=nrows+1; i++ ) {
      x = xmin + dx * (double) (first_row+i-1);
      y = ymin;
      new[i][0]       = old[i][0]       = exact( x, y );
    }

  /* Apply boundary conditions at the ymax boundary */
  if ( northid==-1 )
    for ( i=0; i<=nrows+1; i++ ) {
      x = xmin + dx * (double) (first_row+i-1);
      y = ymax;
      new[i][ncols+1] = old[i][ncols+1] = exact( x, y );
    }


  /* Carry out a single iteration to initialise MaxDiffG */
  MaxDiff = iteration( old, new, 1, nrows, 1, ncols );
  MPI_Allreduce( &MaxDiff, &MaxDiffG, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );


  iter = 1;
  if ( nid==0 )
    fprintf( stdout, "Iteration = %i, MaxDiffG = %12.5e\n", iter, MaxDiffG );


  while ( iter<MAX_ITER ) { /* The error estimate is still too large */
    /* Replace the old values with the new values */
    tmp = new;
    new = old;
    old = tmp;


    /* Deal with the north boundary of the block (needs buffer) */
    req_sendN = req_recvS = MPI_REQUEST_NULL;
    if ( northid!=-1 ) {
      for ( i=1; i<=nrows; i++ )
        send_columnN[i-1] = old[i][ncols];

      MPI_Isend( &send_columnN[0], nrows, MPI_DOUBLE, northid, 10, MPI_COMM_WORLD, &req_sendN );
      MPI_Irecv( &recv_columnN[0], nrows, MPI_DOUBLE, northid, 20, MPI_COMM_WORLD, &req_recvS );
    }

    /* Deal with the south boundary of the block (needs buffer) */
    req_sendS = req_recvN = MPI_REQUEST_NULL;
    if ( southid!=-1 ) {
      for ( i=1; i<=nrows; i++ )
        send_columnS[i-1] = old[i][1];

      MPI_Isend( &send_columnS[0], nrows, MPI_DOUBLE, southid, 20, MPI_COMM_WORLD, &req_sendS );
      MPI_Irecv( &recv_columnS[0], nrows, MPI_DOUBLE, southid, 10, MPI_COMM_WORLD, &req_recvN );
    }

    /* Deal with the west boundary of the block (doesn't need buffer) */
    req_sendW = req_recvE = MPI_REQUEST_NULL;
    if ( westid!=-1 ) {
      MPI_Isend( &old[1][1], ncols, MPI_DOUBLE, westid, 30, MPI_COMM_WORLD, &req_sendW );
      MPI_Irecv( &old[0][1], ncols, MPI_DOUBLE, westid, 40, MPI_COMM_WORLD, &req_recvE );
    }

    /* Deal with the east boundary of the block (doesn't need buffer) */
    req_sendE = req_recvW = MPI_REQUEST_NULL;
    if ( eastid!=-1 ) {
      MPI_Isend( &old[nrows][1],   ncols, MPI_DOUBLE, eastid, 40, MPI_COMM_WORLD, &req_sendE );
      MPI_Irecv( &old[nrows+1][1], ncols, MPI_DOUBLE, eastid, 30, MPI_COMM_WORLD, &req_recvW );
    }

    /* Update the solution values which don't require communication with other processes */
    MaxDiff = iteration( old, new, 2, nrows-1, 2, ncols-1 );


    /* Wait until solution values adjacent to the north boundary have been received... */
    if ( northid!=-1 ) {
      MPI_Wait( &req_recvS, &status );
      for ( i=1; i<=nrows; i++ )
        old[i][ncols+1] = recv_columnN[i-1];
    }
    /* ...before updating the solution values that require this information */
    diff = iteration( old, new, 2, nrows-1, ncols, ncols );
    if ( diff>MaxDiff )
      MaxDiff = diff;


    /* Wait until solution values adjacent to the south boundary have been received... */
    if ( southid!=-1 ) {
      MPI_Wait( &req_recvN, &status );
      for ( i=1; i<=nrows; i++ )
        old[i][0] = recv_columnS[i-1];
    }
    /* ...before updating the solution values that require this information */
    diff = iteration( old, new, 2, nrows-1, 1, 1 );
    if ( diff>MaxDiff )
      MaxDiff = diff;


    /* Wait until solution values adjacent to the west boundary have been received... */
    if ( westid!=-1 )
      MPI_Wait( &req_recvE, &status );
    /* ...before updating the solution values that require this information */
    diff = iteration( old, new, 1, 1, 1, ncols );
    if ( diff>MaxDiff )
      MaxDiff = diff;


    /* Wait until solution values adjacent to the east boundary have been received... */
    if ( eastid!=-1 )
      MPI_Wait( &req_recvW, &status );
    /* ...before updating the solution values that require this information */
    diff = iteration( old, new, nrows, nrows, 1, ncols );
    if ( diff>MaxDiff )
      MaxDiff = diff;
   

    /* Estimate the error for the complete update */
    MPI_Allreduce( &MaxDiff, &MaxDiffG, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );


    /* Write out the error estimate every 1000 iterations */
    iter++;
    if ( nid==0 && iter%1000==0 )
      fprintf( stdout, "Iteration = %i, MaxDiffG = %12.5e\n", iter, MaxDiffG );


    /* Make sure that the Isends are tidied up by MPI */
    if ( northid!=-1 )
      MPI_Wait( &req_sendN, &status );
    if ( southid!=-1 )
      MPI_Wait( &req_sendS, &status );
    if ( westid!=-1 )
      MPI_Wait( &req_sendW, &status );
    if ( eastid!=-1 )
      MPI_Wait( &req_sendE, &status );
  }


  /* Compute the error in the final approximation */
  MaxErr = final_error( new, first_row, nrows, first_col, ncols, dx, dy );
  MPI_Allreduce( &MaxErr, &MaxErrG, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );


  if ( nid==0 ) {
    fprintf( stdout, "Iteration = %i, MaxDiffG = %12.5e\n", iter, MaxDiffG );
    fprintf( stdout, "Error = %12.5e\n", MaxErrG );

    end_time = MPI_Wtime();
    fprintf( stdout, "Execution time in seconds = %12.5e\n", end_time-start_time );
  }


  /* Assign bounds for the local and global indices of the rows and columns to indicate
   * the nodes of the mesh that are updated by this process, purely for the purpose of
   * ending all of the solution to process 0 for writing to a single file               */
  rowmin = first_row;
  rowmax = first_row + nrows - 1;
  colmin = first_col;
  colmax = first_col + ncols - 1;

  /* The data values to be written to file should also include the boundary values, so
   * adjust the row and column bounds on the appropriate processors to include these   */
  if ( westid==-1 )
    rowmin--;
  if ( eastid==-1 )
    rowmax++;
  if ( southid==-1 )
    colmin--;
  if ( northid==-1 )
    colmax++;

  /* Write the data to file */
  //write_file( nid, new, rowmin, rowmax, colmin, colmax, first_row, first_col, noprocs );


  /* Free allocated memory */
  free( new );
  free( old );


  /* Finish up */
  MPI_Finalize();

  return (0);
}
