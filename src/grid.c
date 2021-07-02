
#include "grid.h"
#include "coordinates.h"
#include "geometry.h"
#include <math.h>

int N1, N2, N3;

/********************************************************************
        Interpolation routines
 ********************************************************************/

/* return scalar in cgs units */
double interp_scalar(double X[NDIM], double ***var)
{
  double del[NDIM] = {0}, b1 = 0, b2 = 0, interp = 0;
  int i = 0, j = 0, k = 0, ip1 = 0, jp1 = 0, kp1 = 0;

  // zone and offset from X
  Xtoijk(X, &i, &j, &k, del);

  // since we read from data, adjust i,j,k for ghost zones
  i += 1;
  j += 1;
  k += 1;

  ip1 = i+1;
  jp1 = j+1;
  kp1 = k+1;

  b1 = 1.-del[1];
  b2 = 1.-del[2];

  // interpolate in x1 and x2
  interp = var[i][j][k]*b1*b2 +
    var[i][jp1][k]*b1*del[2] +
    var[ip1][j][k]*del[1]*b2 +
    var[ip1][jp1][k]*del[1]*del[2];

  // then interpolate in x3
  interp = (1.-del[3])*interp +
        del[3]*(var[i  ][j  ][kp1]*b1*b2 +
      var[i  ][jp1][kp1]*del[2]*b1 +
      var[ip1][j  ][kp1]*del[1]*b2 +
      var[ip1][jp1][kp1]*del[1]*del[2]);

  return interp;
}

/* return scalar interpolated in time */
double interp_scalar_time(double X[NDIM], double ***varA, double ***varB, double tfac)
{
  double vA = interp_scalar(X, varA);

#if SLOWLIGHT
  double vB = interp_scalar(X, varB);
  return tfac*vA + (1. - tfac)*vB;
#endif

  return vA;
}

/*
 *  returns geodesic coordinates associated with center of zone i,j,k
 */
void ijktoX(int i, int j, int k, double X[NDIM])
{
  // first do the naive thing
  X[1] = startx[1] + (i+0.5)*dx[1];
  X[2] = startx[2] + (j+0.5)*dx[2];
  X[3] = startx[3] + (k+0.5)*dx[3];

  // now transform to geodesic coordinates if necessary by first
  // converting to KS and then to destination coordinates (eKS).
  if (use_eKS_internal) {
      double xKS[4] = { 0 };
    if (metric == METRIC_MKS3) {
      double x0 = X[0];
      double x1 = X[1];
      double x2 = X[2];
      double x3 = X[3];

      double H0 = mks3H0;
      double MY1 = mks3MY1;
      double MY2 = mks3MY2;
      double MP0 = mks3MP0;

      xKS[0] = x0;
      xKS[1] = exp(x1) + mks3R0;
      xKS[2] = (M_PI*(1+1./tan((H0*M_PI)/2.)*tan(H0*M_PI*(-0.5+(MY1+(pow(2,MP0)*(-MY1+MY2))/pow(exp(x1)+R0,MP0))*(1-2*x2)+x2))))/2.;
      xKS[3] = x3;
    }

    double r, th;
    bl_coord(X, &r, &th);
    X[0] = xKS[0];
    X[1] = r;
    X[2] = th;
    X[3] = xKS[3];
  }
}

/*
 *  translates geodesic coordinates to a grid zone and returns offset
 *  for interpolation purposes. integer index corresponds to the zone
 *  center "below" the desired point and del[i] \in [0,1) returns the
 *  offset from that zone center.
 *
 *  0    0.5    1
 *  [     |     ]
 *  A  B  C DE  F
 *
 *  startx = 0.
 *  dx = 0.5
 *
 *  A -> (-1, 0.5)
 *  B -> ( 0, 0.0)
 *  C -> ( 0, 0.5)
 *  D -> ( 0, 0.9)
 *  E -> ( 1, 0.0)
 *  F -> ( 1, 0.5)
 *
 */
void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM])
{
  // unless we're reading from data, i,j,k are the normal expected thing
  double phi = 0;
  double XG[4] = {0};

  if (use_eKS_internal) {
    // the geodesics are evolved in eKS so invert through KS -> zone coordinates
    double r, th;
    bl_coord(X, &r, &th);
    double Xks[4] = { X[0], r, th, X[3] };
    if (metric == METRIC_MKS3) {
      double H0 = mks3H0, MY1 = mks3MY1, MY2 = mks3MY2, MP0 = mks3MP0;
      double KSx1 = Xks[1], KSx2 = Xks[2];
      XG[0] = Xks[0];
      XG[1] = log(Xks[1] - mks3R0);
      XG[2] = (-(H0*pow(KSx1,MP0)*M_PI) - pow(2.,1. + MP0)*H0*MY1*M_PI +
        2.*H0*pow(KSx1,MP0)*MY1*M_PI + pow(2.,1. + MP0)*H0*MY2*M_PI +
        2.*pow(KSx1,MP0)*atan(((-2.*KSx2 + M_PI)*tan((H0*M_PI)/2.))/M_PI))/(2.*
        H0*(-pow(KSx1,MP0) - pow(2.,1 + MP0)*MY1 + 2.*pow(KSx1,MP0)*MY1 +
          pow(2.,1. + MP0)*MY2)*M_PI);
      XG[3] = Xks[3];
    }
  } else {
    MULOOP XG[mu] = X[mu];
  }

  // the X[3] coordinate is allowed to vary so first map it to [0, cstopx[3])
  phi = fmod(XG[3], cstopx[3]);
  if(phi < 0.0) phi = cstopx[3]+phi;

  // get provisional zone index. see note above function for details. note we
  // shift to zone centers because that's where variables are most exact.
  *i = (int) ((XG[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
  *j = (int) ((XG[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
  *k = (int) ((phi  - startx[3]) / dx[3] - 0.5 + 1000) - 1000;

  // exotic coordinate systems sometime have issues. use this block to enforce
  // reasonable limits on *i,*j and *k. in the normal coordinate systems, this
  // block should never fire.
  if (*i < -1) *i = -1;
  if (*j < -1) *j = -1;
  if (*k < -1) *k = -1;
  if (*i >= N1) *i = N1-1;
  if (*j >= N2) *j = N2-1;
  if (*k >= N3) *k = N3-1;

  // now construct del
  del[1] = (XG[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
  del[2] = (XG[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
  del[3] = (phi - ((*k + 0.5) * dx[3] + startx[3])) / dx[3];

  // and enforce limits on del (for exotic coordinate systems)
  for (int i=0; i<4; ++i) {
    if (del[i] < 0.) del[i] = 0.;
    if (del[i] >= 1.) del[i] = 1.;
  }

}

int X_in_domain(double X[NDIM]) {
  // returns 1 if X is within the computational grid.
  // checks different sets of coordinates depending on
  // specified grid coordinates

  if (use_eKS_internal) {
    double XG[4] = { 0 };
    double r, th;
    bl_coord(X, &r, &th);
    double Xks[4] = { X[0], r, th, X[3] };

    if (metric == METRIC_MKS3) {
      // if METRIC_MKS3, ignore theta boundaries
      double H0 = mks3H0, MY1 = mks3MY1, MY2 = mks3MY2, MP0 = mks3MP0;
      double KSx1 = Xks[1], KSx2 = Xks[2];
      XG[0] = Xks[0];
      XG[1] = log(Xks[1] - mks3R0);
      XG[2] = (-(H0*pow(KSx1,MP0)*M_PI) - pow(2,1 + MP0)*H0*MY1*M_PI +
        2*H0*pow(KSx1,MP0)*MY1*M_PI + pow(2,1 + MP0)*H0*MY2*M_PI +
        2*pow(KSx1,MP0)*atan(((-2*KSx2 + M_PI)*tan((H0*M_PI)/2.))/M_PI))/(2.*
        H0*(-pow(KSx1,MP0) - pow(2,1 + MP0)*MY1 + 2*pow(KSx1,MP0)*MY1 +
          pow(2,1 + MP0)*MY2)*M_PI);
      XG[3] = Xks[3];

      if (XG[1] < startx[1] || XG[1] > stopx[1]) return 0;
    }

  } else {
    if(X[1] < startx[1] ||
       X[1] > stopx[1]  ||
       X[2] < startx[2] ||
       X[2] > stopx[2]) {
      return 0;
    }
  }

  return 1;
}

/*
 * return the gdet associated with zone coordinates for the zone at i,j,k
 */
double gdet_zone(int i, int j, int k)
{
  // get the X for the zone (in geodesic coordinates for bl_coord)
  // and in zone coordinates (for set_dxdX)
  double X[NDIM], Xzone[NDIM];
  ijktoX(i,j,k, X);
  Xzone[0] = 0.;
  Xzone[1] = startx[1] + (i+0.5)*dx[1];
  Xzone[2] = startx[2] + (j+0.5)*dx[2];
  Xzone[3] = startx[3] + (k+0.5)*dx[3];

  if (metric == METRIC_MINKOWSKI || metric == METRIC_EMINKOWSKI) {
    double gcov[NDIM][NDIM];
    gcov_func(Xzone, gcov);
    return gdet_func(gcov);
  }

  // then get gcov for the zone (in zone coordinates)
  double gcovKS[NDIM][NDIM], gcov[NDIM][NDIM];
  double r, th;
  double dxdX[NDIM][NDIM];
  MUNULOOP gcovKS[mu][nu] = 0.;
  MUNULOOP gcov[mu][nu] = 0.;
  bl_coord(X, &r, &th);
  gcov_ks(r, th, gcovKS);
  set_dxdX(Xzone, dxdX);
  MUNULOOP {
    for (int lam=0; lam<NDIM; ++lam) {
      for (int kap=0; kap<NDIM; ++kap) {
        gcov[mu][nu] += gcovKS[lam][kap]*dxdX[lam][mu]*dxdX[kap][nu];
      }
    }
  }

  return gdet_func(gcov);
}
