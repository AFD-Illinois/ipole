
/* 

model-dependent routines for integrating geodesics, 
including:

stop_backward_integration
stepsize

*/

#include "decs.h"

/* condition for stopping the backwards-in-lambda
   integration of the photon geodesic */

//#define LRMAX (log(1.1*Rout))
#define LRMAX (log(rmax))
#define LRMIN (log(1.05*Rh))
int stop_backward_integration(double X[NDIM], double Kcon[NDIM], 
  double Xcam[NDIM])
{
    if ((X[1] > LRMAX && Kcon[1] < 0.) ||       /* out far */
        X[1] < LRMIN            /* in deep */
        )
        return (1);
    else
        return (0);             /* neither out far nor in deep */

}
#undef LRMIN
#undef LRMAX

//#define EPS     0.01

double stepsize(double X[NDIM], double Kcon[NDIM])
{
  double eps = 0.5;
  double small = 1.e-40; 

  double dl, dlx1, dlx2, dlx3;
  double idlx1,idlx2,idlx3 ;

  dlx1 = eps / (fabs(Kcon[1]) + small*small) ;
  //dlx2 = EPS * GSL_MIN(X[2], 1. - X[2]) / (fabs(Kcon[2]) + SMALL*SMALL) ;
  dlx2 = eps * MIN(X[2], 1. - X[2]) / (fabs(Kcon[2]) + small*small) ;
  dlx3 = eps / (fabs(Kcon[3]) + small*small) ;

  idlx1 = 1./(fabs(dlx1) + small*small) ;
  idlx2 = 1./(fabs(dlx2) + small*small) ;
  idlx3 = 1./(fabs(dlx3) + small*small) ;

  dl = 1. / (idlx1 + idlx2 + idlx3) ;

  //dl = MIN(MIN(1./dlx1, 1./idlx2), 1./dlx3);
  //printf("idlx = %e %e %e\n", 1./idlx1, 1./idlx2, 1./idlx3);

  //dl = MIN(dl, 0.5*DTd/Kcon[0]);

  return dl;
}
//#undef EPS

