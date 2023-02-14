#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

# define NPOINTS 2000
# define MAXITER 2000

struct complex{
  double real;
  double imag;
};

int main(void)
{
  int i, j, iter, numoutside = 0;
  int numThreads, currentThread;
  int neach, gridStart, gridEnd;
  double area, error, ztemp;
  struct complex z, c;

/* ****************************************************************************
 *     Outer loops run over npoints, initialise z=c
 *     Inner loop has the iteration z=z*z+c, and threshold test
 * ***************************************************************************/

#pragma omp parallel default(none), \
private(i, j, iter, numThreads, currentThread, gridStart, gridEnd, c, z, ztemp,\
neach), reduction(+:numoutside)
{
  /* get the current number of threads available */
  numThreads = omp_get_num_threads();

  #pragma omp single
  {
    printf("The number of threads being used is: %i\n", numThreads);
  }
  /* get which thread is current used */
  currentThread = omp_get_thread_num();

  /* calculate the number of points for each thread */
  neach = NPOINTS/numThreads;
  /* assign start and end point for each thread */
  gridStart = currentThread * neach;
  gridEnd = (currentThread+1) * neach;

  for (i=gridStart; i<gridEnd; i++)
  {
    for (j=0; j<NPOINTS; j++)
    {
      c.real = -2.0+2.5*(double)(i)/(double)(NPOINTS)+1.0e-7;
      c.imag = 1.125*(double)(j)/(double)(NPOINTS)+1.0e-7;
      z=c;
      for (iter=0; iter<MAXITER; iter++)
      {
	      ztemp=(z.real*z.real)-(z.imag*z.imag)+c.real;
	      z.imag=z.real*z.imag*2+c.imag;
	      z.real=ztemp;
	      if ((z.real*z.real+z.imag*z.imag)>4.0e0)
        {
	        numoutside++;
	        break;
	      }
      }
    }
  }
}

 /* Calculate area and error and output the results */

  area=2.0*2.5*1.125*(double)(NPOINTS*NPOINTS-numoutside)/(double)(NPOINTS*NPOINTS);
  error=area/(double)NPOINTS;
  printf("Area of Mandlebrot set = %12.8f +/- %12.8f\n",area,error);

  return 0;
}
