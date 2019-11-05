#include "geodesics.h"

#include "decs.h"
#include "geometry.h"

/*
 * This is the main photon orbit integrator
 */
void push_photon (double X[NDIM], double Kcon[NDIM], double dl, double Xhalf[NDIM],
             double Kconhalf[NDIM])
{
  double lconn[NDIM][NDIM][NDIM];
  double dKcon[NDIM];
  double Xh[NDIM], Kconh[NDIM];
  int i, j, k;

  /* 2nd order: scheme
   take half-step and then evaluate derivatives
   at half-step. */

  /** half-step **/
  get_connection(X, lconn);

  /* advance K */
  for (k = 0; k < 4; k++)
    dKcon[k] = 0.;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      for (k = 0; k < 4; k++)
        dKcon[k] -= 0.5 * dl * lconn[k][i][j] * Kcon[i] * Kcon[j];
  for (k = 0; k < 4; k++)
    Kconh[k] = Kcon[k] + dKcon[k];

  /* advance X */
  for (i = 0; i < 4; i++)
    Xh[i] = X[i] + 0.5 * dl * Kcon[i];

  /********************/
  for (i = 0; i < 4; i++)
    Xhalf[i] = Xh[i];
  for (i = 0; i < 4; i++)
    Kconhalf[i] = Kconh[i];
  /********************/

  /** full step **/
  get_connection(Xh, lconn);

  /* advance K */
  for (k = 0; k < 4; k++)
    dKcon[k] = 0.;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      for (k = 0; k < 4; k++)
        dKcon[k] -= dl * lconn[k][i][j] * Kconh[i] * Kconh[j];
  for (k = 0; k < 4; k++)
    Kcon[k] += dKcon[k];

  /* advance X */
  for (k = 0; k < 4; k++)
    X[k] += dl * Kconh[k];

  /*
   gcov_func(X,gcov);
   lower(Kcon, gcov, Kcov);
   E0 = -Kcov[0];
   L0 = Kcov[3];
   fprintf(stdout,"along g: r=%g E0=%g L0=%g \n",exp(X[1]),E0,L0);
   */
  /* done! */
}

