
/*
 * model-dependent routines for integrating geodesics,
 * including:
 *
 * stop_backward_integration
 * stepsize
 */

#include "decs.h"
#include "model.h"
#include "coordinates.h"

/* condition for stopping the backwards-in-lambda
   integration of the photon geodesic */

int stop_backward_integration(double X[NDIM], double Xhalf[NDIM], double Kcon[NDIM],
  double Xcam[NDIM])
{
  // Necessary geometric stop conditions
  if ((X[1] > log (1.1 * rmax_geo) && Kcon[1] < 0.) || // Stop either beyond rmax_geo
      X[1] < log (1.05 * Rh)) { // Or right near the horizon
    return (1);
  }

  // Additional stop condition for thin disks: the disk is opaque
#if THIN_DISK
  static int n_left = -1;
#pragma omp threadprivate(n_left)

  if (thindisk_region(X, Xhalf) && n_left < 0) { // Set timer when we reach disk
    n_left = 2;
    return 0;
  } else if (n_left < 0) { // Otherwise continue normally
    return 0;
  } else if (n_left > 0) { // Or decrement the timer if we need
    n_left--;
    return 0;
  } else {  // If timer is 0, stop and reset
    n_left = -1;
    return 1;
  }
#endif

  return (0);
}

double stepsize(double X[NDIM], double Kcon[NDIM])
{
  double eps = 0.01; // TODO take as parameter w/MAXNSTEP?

  double dl, dlx1, dlx2, dlx3;
  double idlx1,idlx2,idlx3 ;

  dlx1 = eps / (fabs(Kcon[1]) + SMALL*SMALL) ;
  //dlx2 = EPS * GSL_MIN(X[2], 1. - X[2]) / (fabs(Kcon[2]) + SMALL*SMALL) ;
  dlx2 = eps * MIN(X[2], 1. - X[2]) / (fabs(Kcon[2]) + SMALL*SMALL) ;
  dlx3 = eps / (fabs(Kcon[3]) + SMALL*SMALL) ;

  idlx1 = 1./(fabs(dlx1) + SMALL*SMALL) ;
  idlx2 = 1./(fabs(dlx2) + SMALL*SMALL) ;
  idlx3 = 1./(fabs(dlx3) + SMALL*SMALL) ;

  dl = 1. / (idlx1 + idlx2 + idlx3) ;

  //dl = MIN(MIN(1./dlx1, 1./idlx2), 1./dlx3);
  //printf("idlx = %e %e %e\n", 1./idlx1, 1./idlx2, 1./idlx3);

  //dl = MIN(dl, 0.5*DTd/Kcon[0]);

  return dl;
}

