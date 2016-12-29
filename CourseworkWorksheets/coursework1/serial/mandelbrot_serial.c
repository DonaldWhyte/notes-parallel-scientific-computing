/*
 * Serial implememntation of Mandelbrot program.
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
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

#include "mandelbrot-gui.h"     /* contains setup(), interact() */


/* Default values */
#define N           2           /* size of problem space (x, y from -N to N) */
#define NPIXELS     1000        /* size of display window in pixels */


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
int serial_pgm( uint width, uint height, double real_min, double real_max,
                double imag_min, double imag_max, uint maxiter);



/* Draw a row of pixels in the display */
void drawLine( int row, int width, ulong *ilines, XDATA *windata )
{
  int col;

  for ( col=0; col<width; ++col ) {
    XSetForeground( windata->display, windata->gc, ilines[width*row+col] );
    XDrawPoint( windata->display, windata->win, windata->gc, col, row );
  }
}



/* Main program */
int main( int argc, char **argv )
{
  int returnval;
  uint maxiter;
  double real_min=-N, real_max=N, imag_min=-N, imag_max=N;
  uint width=NPIXELS, height=NPIXELS;         /* dimensions of display window */
  double x0=0.0, y0=0.0;
  double size=2.0;


  /* Process command-line arguments */
  maxiter = 10000;
  real_min = x0 - size;
  real_max = x0 + size;
  imag_min = y0 - size;
  imag_max = y0 + size;


  /* Call serial code */
  returnval = serial_pgm( width, height, real_min, real_max, imag_min, imag_max, maxiter);
  fprintf( stdout, "Process finished\n" );


  return returnval;
}



/* Program for serial process.
 *  -- Returns EXIT_SUCCESS or EXIT_FAILURE as appropriate. */
int serial_pgm( uint width, uint height, double real_min, double real_max,
                double imag_min, double imag_max, uint maxiter ) {

  XDATA windata;

  ulong min_colour, max_colour;
  uint col, row;
  int setup_return;

  ulong colour;
  double scale_real, scale_imag, scale_colour;
  uint count;
  double lengthsq, temp;
  complex z, c;
  ulong *ilines;


  /* Initialise for graphical display */
  setup_return = setup( width, height, &windata.display, &windata.win,
                        &windata.gc, &min_colour, &max_colour );
  /* If not successful, continue but don't display results */
  if ( setup_return!=EXIT_SUCCESS ) {
    fprintf( stderr, "Unable to initialise display, continuing\n" );
  }


  /* Compute factors to scale computational region to window */
  scale_real = (double) (real_max-real_min) / (double) width;
  scale_imag = (double) (imag_max-imag_min) / (double) height; 

  /* Compute factor for colour scaling */
  scale_colour = (double) (max_colour-min_colour) / (double) (maxiter-1);


  /* Allocate memory */
  ilines = calloc( width*height, sizeof(ulong) );


  /* Compute points and draw them ONE ROW AT A TIME */
  for ( row=0; row<height; ++row ) {

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


      /* Scale colour and store */
      colour = (ulong) ((count-1) * scale_colour) + min_colour;

      ilines[width*row+col] = colour;
    }


    /* Draw the points for this row */
    drawLine( row, width, ilines, &windata );
  }


  /* Be sure everything is written out */
  if (setup_return==EXIT_SUCCESS) {
    XFlush(windata.display);
  }


  /* Produce text output */
  fprintf( stdout, "Sequential program\n" );
  fprintf( stdout, "Centre = (%g, %g), Size = %g\n",
           (real_max+real_min)/2, (imag_max+imag_min)/2, (real_max-real_min)/2 );
  fprintf( stdout, "Maximum iterations = %d\n", maxiter );


  /* Wait for user response, then exit program */
  if (setup_return==EXIT_SUCCESS) {
    interact( windata.display, &windata.win, windata.gc, width, height,
              real_min, real_max, imag_min, imag_max );
  }


  /* Free array memory */
  free(ilines);


  return EXIT_SUCCESS;
}
