/*
 * Parallel implementation of Mandelbrot program.
 *
 * This program computes and displays the Mandelbrot set.
 * By default, it examines all points in the complex plane
 * that have both real and imaginary parts between -2 and 2.  
 * 
 * Code originally obtained from Web site for Wilkinson and
 * Allen's text on parallel programming:
 *   http://www.cs.uncc.edu/~abw/parallel/par_prog/
 * 
 * Reformatted and revised by B.Massingill.
 * Further reformatted by MEH.
 *   -- Cyclic decomposition implemented by allocating row by row.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

/* Define functions */
#define MAX(a,b) ( ((a)>(b)) ? (a):(b) )
#define MIN(a,b) ( ((a)<(b)) ? (a):(b) )

#include "mandelbrot_zoom-gui.h" /* contains setup(), interact(), zoom() */

/* Default values */
#define N           2           /* size of problem space (x, y from -N to N) */
#define NPIXELS     1000        /* size of display window in pixels */

/* Constants for message passing */
#define start_tag   1           /* ``start'' message (master to worker) */
#define data_tag    2           /* ``data''  message (worker to master) */


/* Structure definition for complex numbers */
typedef struct {
  double real, imag;
} complex;

/* Structure definition for the Xwindow */
typedef struct {
    Display *display;
    Window win;
    GC gc;
} XDATA;


/* Shorthand for some commonly-used types */
typedef unsigned int uint;
typedef unsigned long ulong;


/* Function declarations */
int master_pgm( int nworkers, uint width, uint height, double *xylim_orig, uint maxiter);
int worker_pgm( int nworkers, int nid, uint width, uint height, uint maxiter);



/* Draw a row of pixels in the display */
void drawLine( int row, int width, ulong *ilines, Display *display, Window window, GC gc )
{
  int col;

  for ( col=0; col<width; ++col ) {
    XSetForeground( display, gc, ilines[width*row+col] );
    XDrawPoint( display, window, gc, col, row );
  }
}



/* Main program */
int main( int argc, char **argv )
{
  int nid, noprocs;
  int returnval;
  uint maxiter;
  double xylim_orig[4];
  uint width=NPIXELS, height=NPIXELS;         /* dimensions of display window */
  double x0=0.0, y0=0.0;
  double size=2.0;


  /* Initialize for MPI */
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &nid );
  MPI_Comm_size( MPI_COMM_WORLD, &noprocs );
  if ( noprocs<2 ) {
    printf( "Number of processes must be at least 2\n" );
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }


  /* Initialise parameters */
  maxiter = 10000;
  xylim_orig[0] = x0 - size;
  xylim_orig[1] = x0 + size;
  xylim_orig[2] = y0 - size;
  xylim_orig[3] = y0 + size;


  /* Call master or worker code as appropriate */
  if ( nid==0 ) {
    returnval = master_pgm( noprocs-1, width, height, xylim_orig, maxiter);
    fprintf( stdout, "Master finished\n" );
  }
  else {
    returnval = worker_pgm( noprocs-1, nid, width, height, maxiter);
    fprintf( stdout, "Worker %d finished\n", nid );
  }


  /* Finish up */
  MPI_Finalize();


  return returnval;
}



/* Program for master process.
 *  -- Returns EXIT_SUCCESS or EXIT_FAILURE as appropriate. */
int master_pgm( int nworkers, uint width, uint height, double *xylim_orig, uint maxiter )
{
  XDATA windata;

  ulong min_colour, max_colour;
  int setup_return, stop_msg, i;
  double xylim[4];


  /* Initialise the current bounds for the region to be drawn */
  for ( i=0; i<4; ++i )
    xylim[i] = xylim_orig[i];


  /* Initialise for graphical display */
  setup_return = setup( width, height, &windata.display, &windata.win,
                        &windata.gc, &min_colour, &max_colour );


  /* If not successful, continue but don't display results */
  if ( setup_return!=EXIT_SUCCESS ) {
    fprintf( stderr, "Unable to initialize display, continuing\n" );
  }
  else {
    /* Trigger the workers to compute the pixel values */
    stop_msg = 0;
    MPI_Bcast( &stop_msg, 1, MPI_INT, 0, MPI_COMM_WORLD );

    /* Draw the initial region in the graphical display */
    draw( nworkers, width, height, xylim, maxiter, min_colour, max_colour,
          windata.display, windata.win, windata.gc, setup_return );

    /* Wait for interaction by the user and respond to it */
    interact( windata.display, windata.win, windata.gc, width, height,
              xylim, xylim_orig, nworkers, maxiter, min_colour, max_colour, setup_return );
  }


  /* Tell the workers to finish */
  stop_msg = 1;
  MPI_Bcast( &stop_msg, 1, MPI_INT, 0, MPI_COMM_WORLD );

  return EXIT_SUCCESS;
}



/* Draw the specified region in the graphical display */
void draw( int nworkers, uint width, uint height, double *xylim,
           uint maxiter, ulong min_colour, ulong max_colour,
           Display *display, Window window, GC gc, int setup_return )
{
  uint this_row, row_counter, col;
  double start_time, end_time;
  ulong *ilines, *row_data;
  MPI_Status status;

  double real_min = xylim[0];
  double real_max = xylim[1];
  double imag_min = xylim[2];
  double imag_max = xylim[3];


  /* Start timing */
  start_time = MPI_Wtime();


  /* Allocate memory */
  ilines   = calloc( width*height, sizeof(ulong) );
  row_data = calloc( width+1,      sizeof(ulong) );


  /* Broadcast values workers need to find colours and region boundary */
  MPI_Bcast( &min_colour, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
  MPI_Bcast( &max_colour, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
  MPI_Bcast( xylim, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD );


  /* The rows each worker should work are calculated by the worker */


  /* Draw the rows of pixels as they return to the master */
  for ( row_counter=0; row_counter<height; ++row_counter ) {
    MPI_Recv( row_data, width+1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, data_tag, MPI_COMM_WORLD, &status );

    /* Extract the row number for plotting */
    this_row = row_data[0];

    /* Put the data in to the array ilines */
    for ( col=0; col<width; ++col )
      ilines[width*this_row+col] = row_data[col+1];

    /* Draw the points for this row */
    drawLine( this_row, width, ilines, display, window, gc );
  }

  /* Be sure everything is written out */
  if ( setup_return==EXIT_SUCCESS ) {
    XFlush(display);
  }


  /* End timing */
  end_time = MPI_Wtime();


  /* Produce text output */
  fprintf( stdout, "Parallel program\n" );
  fprintf( stdout, "Number of worker processes = %d\n", nworkers );
  fprintf( stdout, "Centre = (%g, %g), Size = %g\n",
           (real_max+real_min)/2, (imag_max+imag_min)/2, (real_max-real_min)/2 );
  fprintf( stdout, "Maximum iterations = %d\n", maxiter );
  fprintf( stdout, "Execution time in seconds = %g\n\n", end_time-start_time );


  /* Free array memory */
  free( ilines );
  free( row_data );
}



/* Program for worker process.
 *  -- Returns EXIT_SUCCESS or EXIT_FAILURE as appropriate. */
int worker_pgm( int nworkers, int nid, uint width, uint height, uint maxiter )
{
  ulong min_colour, max_colour;
  double real_min, real_max, imag_min, imag_max;
  double xylim[4];
  uint col, row;

  ulong colour;
  double scale_real, scale_imag, scale_colour;
  uint count;
  double lengthsq, temp;
  complex z, c;
  ulong *row_data;

  int stop_msg;


  while (True) {
    /* Wait to receive a trigger message broadcast by the master before
     * attempting any computations */
    MPI_Bcast( &stop_msg, 1, MPI_INT, 0, MPI_COMM_WORLD );
    if ( stop_msg==1 ) {
      /* The user wants to exit the program */
      return EXIT_SUCCESS;
    }


    /* Receive broadcast values for computing colours and region boundary */
    MPI_Bcast( &min_colour, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( &max_colour, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( xylim, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    real_min = xylim[0];
    real_max = xylim[1];
    imag_min = xylim[2];
    imag_max = xylim[3];


    /* Compute factors to scale computational region to window */
    scale_real = (double) (real_max-real_min) / (double) width;
    scale_imag = (double) (imag_max-imag_min) / (double) height; 

    /* Compute factor for colour scaling */
    scale_colour = (double) (max_colour-min_colour) / (double) (maxiter-1);


    /* Allocate memory */
    row_data = calloc( width+1, sizeof(ulong) );


    /* Compute points and draw them ONE ROW AT A TIME
     *   -- Note that the loop counter implements cyclic decomposition */
    for ( row=nid-1; row<height; row+=nworkers ) {
      /* Store the row number so that the master can plot it */
      row_data[0] = row;

      /* Scale vertical display coordinates to actual region */
      c.imag = imag_min + ((double) (height-1-row) * scale_imag);
        /* height-1-row used so y-axis displays with larger values at top */

      for ( col=0; col<width; ++col ) {
        /* Scale horizontal display coordinates to actual region */
        c.real = real_min + ((double) col * scale_real);

        /* Compute z0,z1,... until divergence or maximum iterations */
        z.real = z.imag = 0.0;
        count = 0;
        do {
          temp   = z.real*z.real - z.imag*z.imag + c.real;
          z.imag = 2.0*z.real*z.imag + c.imag;
          z.real = temp;

          lengthsq = z.real*z.real + z.imag*z.imag;
          ++count;
        } while ( lengthsq<(N*N) && count<maxiter );

        /* Scale colour and store in array row_data */
        colour = (ulong) ((count-1) * scale_colour) + min_colour;
        row_data[col+1] = colour;
      }

      /* Send the data back to the master process for plotting */
      MPI_Send( row_data, width+1, MPI_UNSIGNED_LONG, 0, data_tag, MPI_COMM_WORLD );
    }


    /* Free array memory */
    free( row_data );
  }


  return EXIT_SUCCESS;
}
