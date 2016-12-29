/* GUI-related functions for Mandelbrot program (C version) */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>


/* Function declarations */
void draw( int nworkers, uint width, uint height, double *xylim,
           uint maxiter, ulong min_colour, ulong max_colour,
           Display *display, Window window, GC gc, int setup_return );
void zoom( int nworkers, uint width, uint height, double *xylim,
           uint maxiter, ulong min_colour, ulong max_colour,
           Display *display, Window window, GC gc, int setup_return );



/* Initialize for graphical display. 
 *   -- Width, height are dimensions of display, in pixels. */
int setup( uint width, uint height, Display **display, Window *win,
           GC *gc, ulong *min_colour, ulong *max_colour )
{
  /* Variables for graphical display */
  uint x=0, y=0;                      /* window position */
  uint border_width=4;                /* border width in pixels */
  uint disp_width, disp_height;       /* size of screen */
  uint screen;                        /* which screen */

  char *window_name="Mandelbrot Set", *disp_name=NULL;
  ulong valuemask=0;
  XGCValues values;

  ulong white, black;                 /* white, black pixel values */

  XEvent report;


  /* Connect to Xserver */
  if ( (*display = XOpenDisplay (disp_name)) == NULL ) {
    fprintf( stderr, "Cannot connect to X server %s\n", XDisplayName(disp_name) );
    return EXIT_FAILURE;
  }


  /* Initialise for graphical display  */
  screen = DefaultScreen( *display );
  disp_width  = DisplayWidth ( *display, screen );
  disp_height = DisplayHeight( *display, screen );
  *win = XCreateSimpleWindow( *display, RootWindow (*display, screen), x, y, width, height,
                border_width, BlackPixel (*display, screen), WhitePixel (*display, screen) );
  XStoreName( *display, *win, window_name );
  *gc = XCreateGC( *display, *win, valuemask, &values ); /* graphics context */
  white = WhitePixel( *display, screen );                /* colour value for white */
  black = BlackPixel( *display, screen );                /* colour value for black */
  XSetBackground( *display, *gc, white );
  XSetForeground( *display, *gc, black );
  XMapWindow( *display, *win );
  XSync( *display, False );


  /* Get min and max for range of colour values
   *   -- assumed to be defined by "white", "black" */
  *min_colour = (white > black) ? black : white;
  *max_colour = (white > black) ? white : black;


  /* Wait for keyboard input before starting program */
  fprintf( stdout, "Press any key (with cursor in display window) to start the program\n\n" );
  fflush( stdout );


  /*  Choose which events we want to handle   */
  XSelectInput( *display, *win, KeyPressMask );

  /* Wait for event */
  XNextEvent( *display, &report );


  return EXIT_SUCCESS;
}



/* Wait for user response before ending program.
 *  -- Also allows user to discover coordinates of points. */
void interact( Display *display, Window window, GC gc, uint width, uint height,
               double *xylim, double *xylim_orig, int nworkers, uint maxiter,
               ulong min_colour, ulong max_colour, int setup_return )
{
  double scale_real, scale_imag; 
  XEvent report;
  Window root_return, child_return;
  int root_x_return, root_y_return;
  int win_x_return, win_y_return;
  int stop_msg, i;
  uint mask_return;
  char text[10];
  KeySym key;

  double real_min = xylim[0];
  double real_max = xylim[1];
  double imag_min = xylim[2];
  double imag_max = xylim[3];


  fprintf( stdout, "To get the coordinates of a point in the display:\n" );
  fprintf( stdout, "     Use the mouse to click on the point\n" );
  fprintf( stdout, "To zoom in on a region of the display:\n" );
  fprintf( stdout, "     Press z on the keyboard\n" );
  fprintf( stdout, "     Use the mouse to click on one corner of the region you want to look at\n" );
  fprintf( stdout, "     Drag the mouse to the opposite corner and release the mouse button\n" );
  fprintf( stdout, "To return to the original view of the Mandelbrot set\n" );
  fprintf( stdout, "     Press r on the keyboard\n" );
  fprintf( stdout, "To end the program:\n" );
  fprintf( stdout, "     Press q on the keyboard\n\n" );
  fflush( stdout );

  /* Choose which events we want to handle */
  XSelectInput( display, window, ExposureMask | KeyPressMask | ButtonPressMask );


  /* Compute scaling factors (for processing mouse clicks) */
  scale_real = (double) (real_max-real_min) / (double) width;
  scale_imag = (double) (imag_max-imag_min) / (double) height;


  /* Event loop */
  stop_msg = 0;
  while ( stop_msg==0 ) { /* Monitor and respond to events until told to stop */
    XNextEvent( display, &report );

    switch ( report.type ) {

      case MappingNotify:

        /* Make sure that the keyboard mapping is up to date */
        XRefreshKeyboardMapping( &report.xmapping );
        break;

      case Expose:

        /* Redraw the data in the graphical display after being uncovered */
        if ( report.xexpose.count==0 ) {
          /* Tell the workers to reccompute the pixel values */
          MPI_Bcast( &stop_msg, 1, MPI_INT, 0, MPI_COMM_WORLD );

          /* Redraw the current region of the Mandelbrot set */
          draw( nworkers, width, height, xylim, maxiter,
                min_colour, max_colour, display, window, gc, setup_return );
        }
        break;

      case ButtonPress:

        /* Display the coordinates of the point the cursor is indicating when
         * the mouse button is pressed */
        XQueryPointer( display, window, &root_return, &child_return, &root_x_return,
                       &root_y_return, &win_x_return, &win_y_return, &mask_return );
        fprintf( stdout, "coordinates = (%g, %g)\n",
                 real_min + ((double) win_x_return * scale_real),
                 imag_min + ((double) (height-1-win_y_return) * scale_imag));
                   /* height-1-row used so y-axis displays with larger values at top */
        fflush( stdout );
        break;

      case KeyPress:

        /* Respond according to the key pressed */
        i = XLookupString( &report.xkey, text, 10, &key, 0 );
        if ( i==1 && text[0]=='q' )
          /* Pressing q stops the program and exits */
          stop_msg = 1;
        else if ( i==1 && text[0]=='z' ) {
          /* Pressing z allows the user to zoom in using the mouse */
          zoom( nworkers, width, height, xylim, maxiter,
                min_colour, max_colour, display, window, gc, setup_return );

          /* Tell the workers to recompute the pixel values */
          MPI_Bcast( &stop_msg, 1, MPI_INT, 0, MPI_COMM_WORLD );

          /* Clear the window then draw the new region of the Mandelbrot set */
          XClearWindow( display, window );
          draw( nworkers, width, height, xylim, maxiter,
                min_colour, max_colour, display, window, gc, setup_return );
        }
        else if ( i==1 && text[0]=='r' ) {
          /* Pressing r returns the display to its original region */
          for ( i=0; i<4; ++i )
            xylim[i] = xylim_orig[i];

          /* Tell the workers to recompute the pixel values */
          MPI_Bcast( &stop_msg, 1, MPI_INT, 0, MPI_COMM_WORLD );

          /* Clear the window then draw the original region of the Mandelbrot set */
          XClearWindow( display, window );
          draw( nworkers, width, height, xylim, maxiter,
                min_colour, max_colour, display, window, gc, setup_return );
        }

        /* Recompute scaling factors (for processing mouse clicks) */
        real_min = xylim[0];
        real_max = xylim[1];
        imag_min = xylim[2];
        imag_max = xylim[3];

        scale_real = (double) (real_max-real_min) / (double) width;
        scale_imag = (double) (imag_max-imag_min) / (double) height;
        break;
    }
  }


  /* Close the display */
  XFreeGC( display, gc );
  XDestroyWindow( display, window );
  XCloseDisplay( display );


  return;
}



/* Use the mouse to zoom in on a region of the graphical display */
void zoom( int nworkers, uint width, uint height, double *xylim,
           uint maxiter, ulong min_colour, ulong max_colour,
           Display *display, Window window, GC gc, int setup_return )
{
  XEvent zoom_report;
  double scale_real, scale_imag, sc;
  int stop_msg;
  Window root_return, child_return;
  int root_x_return, root_y_return;
  int win_x_return1, win_y_return1;
  int win_x_return2, win_y_return2;
  uint mask_return;
  int press=0;
  double xmin, ymin, temp;

  double real_min = xylim[0];
  double real_max = xylim[1];
  double imag_min = xylim[2];
  double imag_max = xylim[3];


  /* Compute scaling factors (for processing mouse clicks) */
  scale_real = (double) (real_max-real_min) / (double) width;
  scale_imag = (double) (imag_max-imag_min) / (double) height;

  /* Choose which events we want to handle */
  XSelectInput( display, window, ExposureMask | KeyPressMask | ButtonPressMask | ButtonReleaseMask );


  /* Event loop */
  stop_msg = 0;
  while ( stop_msg==0 ) { /* Monitor and respond to events until told to stop */
    XNextEvent( display, &zoom_report );

    switch ( zoom_report.type ) {

      case MappingNotify:

        /* Make sure that the keyboard mapping is up to date */
        XRefreshKeyboardMapping( &zoom_report.xmapping );
        break;

      case Expose:

        /* Redraw the data in the graphical display after being uncovered */
        if ( zoom_report.xexpose.count==0 ) {
          /* Tell the workers to recompute the pixel values */
          MPI_Bcast( &stop_msg, 1, MPI_INT, 0, MPI_COMM_WORLD );

          /* Redraw the current region of the Mandelbrot set */
          draw( nworkers, width, height, xylim, maxiter,
                min_colour, max_colour, display, window, gc, setup_return );
        }
        break;

      case ButtonPress:

        /* Extract the coordinates of the cursor position when the mouse button
         * is pressed */
        if ( press==0 ) {
          XQueryPointer( display, window, &root_return, &child_return, &root_x_return,
                         &root_y_return, &win_x_return1, &win_y_return1, &mask_return );
          press = 1;
        }
        break;

      case ButtonRelease:

        /* Extract the coordinates of the cursor position when the mouse button
         * is released, then tell the event loop to finish */
        XQueryPointer( display, window, &root_return, &child_return, &root_x_return,
                       &root_y_return, &win_x_return2, &win_y_return2, &mask_return );
        stop_msg = 1;
        break;
    }
  }


  /* Recompute the dimensions of the region being drawn */
  xmin = xylim[0];
  ymin = xylim[2];

  xylim[0] = xmin + scale_real * (double) win_x_return1;
  xylim[1] = xmin + scale_real * (double) win_x_return2;
  xylim[2] = ymin + scale_imag * (double) (height-1-win_y_return1);
  xylim[3] = ymin + scale_imag * (double) (height-1-win_y_return2);

  /* Reorder the limits so that it doesn't matter which vertices of the
   * rectangle the mouse button was pressed and released at */
  if ( xylim[0]>xylim[1] ) {
    temp     = xylim[0];
    xylim[0] = xylim[1];
    xylim[1] = temp;
  }
  if ( xylim[2]>xylim[3] ) {
    temp     = xylim[2];
    xylim[2] = xylim[3];
    xylim[3] = temp;
  }

  /* Ensure that the new region drawn is a square that contains the specified
   * rectangle */
  sc = MAX( xylim[1]-xylim[0], xylim[3]-xylim[2] );

  xylim[1] = xylim[0] + sc;
  xylim[3] = xylim[2] + sc;


  return;
}
