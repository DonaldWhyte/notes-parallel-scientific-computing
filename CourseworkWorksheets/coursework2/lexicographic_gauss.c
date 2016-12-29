/*
 * Serial implementation of 2D Laplace equation solver
 *   -- Lexicographic Gauss-Seidel
 *   -- non-constant boundary conditions specified
 *   -- exact solution compared with
 *
 * Code supplied by Peter Jimack and modified by Matthew Hubbard
 * and Donald Whyte */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

/* Define parameter values */
#define N   241
#define Tol 0.000000000001
// 10^-5 = 0.00001
// 10^-6 = 0.000001
// 10^-8 = 0.00000001
// 10^-10 = 0.0000000001
// 10^-12 = 0.000000000001

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



/* Carry out a single Gauss Seidel iteration on the specified range of mesh nodes */
double iteration( double **old, double **new, int first_row, int last_row, int first_col, int last_col )
{
  double diff, maxdiff=0.0;
  int i, j;

  for ( i=first_row; i<=last_row; i++ )
    for ( j=first_col; j<=last_col; j++ ) {
      /* The Gauss-Seidel update for node (i,j) */
      new[i][j] = 0.25 * ( old[i+1][j] + new[i-1][j] + old[i][j+1] + new[i][j-1] );

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
  double diff, MaxDiff, MaxDiffG, MaxErr, MaxErrG;
  double dx, dy, x, y;
  double start_time=0.0, end_time=0.0;
  int first_row, nrows, first_col, ncols;
  int noprocs, nid, remainder, i, j, iter;
  int rowmin, rowmax, colmin, colmax;
  MPI_Status status;
  MPI_Request req_send10, req_send20, req_recv10, req_recv20;


  /* Initialise for MPI */
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &nid );
  MPI_Comm_size( MPI_COMM_WORLD, &noprocs );


  /* Calculate the mesh size */
  dx = ( xmax - xmin ) / (double) N;
  dy = ( ymax - ymin ) / (double) N;


  /* Calculate the first row and the number of rows allocated to this process */
  first_row = 1;
  remainder = (N-1)%noprocs;
  nrows     = (N-1-remainder)/noprocs;
  /* Care is taken with any left over rows after dividing the number of rows by the number of
   * processes: they are allocated one-by-one to the first few processes to avoid having one
   * process taking much more work than all of the others                                     */
  for ( i=0; i<nid; i++ ) {
    if ( i<remainder )
      first_row = first_row + nrows + 1;
    else
      first_row = first_row + nrows;
  }
  if ( nid<remainder )
    nrows++;

  /* Calculate the first column and the number of columns allocated to this process */
  first_col = 1;
  ncols     = N-1;


  /* Allocate memory for the local work arrays
   *   -- note that an extra layer of nodes is required all round the updated region */
  new = matrix( nrows+2, ncols+2 );
  old = matrix( nrows+2, ncols+2 );


  if ( nid==0 )
    start_time = MPI_Wtime();


  /* Apply boundary conditions at the ymin and ymax boundaries */
  for ( i=0; i<=nrows+1; i++ ) {
    x = xmin + dx * (double) (first_row+i-1);
    y = ymin;
    new[i][0]       = old[i][0]       = exact( x, y );
    y = ymax;
    new[i][ncols+1] = old[i][ncols+1] = exact( x, y );
  }

  /* Apply boundary conditions at the xmin and xmax boundaries */
  if ( nid==0 )
    for ( j=1; j<=ncols; j++ ) {
      x = xmin;
      y = ymin + dy * (double) (first_col+j-1);
      new[0][j]       = old[0][j]       = exact( x, y );
    }
  if ( nid==noprocs-1 )
    for ( j=1; j<=ncols; j++ ) {
      x = xmax;
      y = ymin + dy * (double) (first_col+j-1);
      new[nrows+1][j] = old[nrows+1][j] = exact( x, y );
    }


  /* Carry out a single iteration to initialise MaxDiffG */
  MaxDiff = iteration( old, new, 1, nrows, 1, N-1 );
  MPI_Allreduce( &MaxDiff, &MaxDiffG, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );


  iter = 1;
  if ( nid==0 )
    fprintf( stdout, "Iteration = %i, MaxDiffG = %12.5e\n", iter, MaxDiffG );


  while ( MaxDiffG>Tol ) { /* The error estimate is still too large */
    /* Replace the old values with the new values */
    tmp = new;
    new = old;
    old = tmp;


    /* For all processes except the last: send the last of its updated rows to the next
     * process and receive the first of the next process's updated rows in to the buffer
     * layer (asynchronously)                                                            */
    req_send10 = req_recv20 = MPI_REQUEST_NULL;
    if ( nid<noprocs-1 ) {
      MPI_Isend( &old[nrows][1],   ncols, MPI_DOUBLE, nid+1, 10, MPI_COMM_WORLD, &req_send10 );
      MPI_Irecv( &old[nrows+1][1], ncols, MPI_DOUBLE, nid+1, 20, MPI_COMM_WORLD, &req_recv20 );
    }

    /* For all processes except the first: send the first of its updated rows to the previous
     * process and receive the last of the previous process's updated rows in to the buffer
     * layer (asynchronously)                                                                 */
    req_send20 = req_recv10 = MPI_REQUEST_NULL;
    if ( nid>0 ) {
      MPI_Isend( &old[1][1], ncols, MPI_DOUBLE, nid-1, 20, MPI_COMM_WORLD, &req_send20 );
      MPI_Irecv( &old[0][1], ncols, MPI_DOUBLE, nid-1, 10, MPI_COMM_WORLD, &req_recv10 );
    }

    /* Update the solution values which don't require communication with other processes */
    MaxDiff = iteration( old, new, 2, nrows-1, 1, N-1 );


    /* Wait until solution values from the next process have been received... */
    if ( nid<noprocs-1 )
      MPI_Wait( &req_recv20, &status );
    /* ...before updating the solution values that require this information */
    diff = iteration( old, new, nrows, nrows, 1, N-1 );
    if ( diff>MaxDiff )
      MaxDiff = diff;


    /* Wait until solution values from the previous process have been received... */
    if ( nid>0 )
      MPI_Wait( &req_recv10, &status );
    /* ...before updating the solution values that require this information */
    diff = iteration( old, new, 1, 1, 1, N-1 );
    if ( diff>MaxDiff )
      MaxDiff = diff;
   

    /* Estimate the error for the complete update */
    MPI_Allreduce( &MaxDiff, &MaxDiffG, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );


    /* Write out the error estimate every 1000 iterations */
    iter++;
    if ( nid==0 && iter%1000==0 )
      fprintf( stdout, "Iteration = %i, MaxDiffG = %12.5e\n", iter, MaxDiffG );

 
    /* Make sure that the Isends are tidied up by MPI */
    if ( nid<noprocs-1 )
      MPI_Wait( &req_send10, &status );
    if ( nid>0 )
      MPI_Wait( &req_send20, &status );
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
  if ( nid==0 ) {
    rowmin--;
  }
  if ( nid==noprocs-1 ) {
    rowmax++;
  }
  colmin--;
  colmax++;

  /* Write the data to file */
  write_file( nid, new, rowmin, rowmax, colmin, colmax, first_row, first_col, noprocs );


  /* Free allocated memory */
  free( new );
  free( old );


  /* Finish up */
  MPI_Finalize();

  return (0);
}
