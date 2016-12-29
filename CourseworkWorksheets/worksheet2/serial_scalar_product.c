#include <stdlib.h>
#include <stdio.h>

int main( int argc, char **argv )
{
  int i, n;
  float *a, *b, sum=0.0;
  FILE *fp;

  // Read two vectors from file
  fp = fopen( "vecs_small.dat", "rt" );
  fscanf( fp, "%d", &n );

  a = (float *) calloc( n, sizeof(float) );
  b = (float *) calloc( n, sizeof(float) );

  for ( i=0; i<n; i++ )
    fscanf( fp, "%f %f", &a[i], &b[i] );
  fclose( fp );

  // Compute scalar product of vector
  for ( i=0; i<n; i++ )
    sum += a[i]*b[i];

  printf( "The scalar product is %f \n", sum );

  free( a ); free( b );

  return 0;
}
