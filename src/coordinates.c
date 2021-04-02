// Coordinate-dependent (but not radiation model-dependent) functions for ipole

#include "coordinates.h"

#include "decs.h"
#include "geometry.h"

int METRIC_eKS;
int METRIC_MKS, METRIC_FMKS, METRIC_MKS3;
double a, hslope;
double poly_norm, poly_xt, poly_alpha, mks_smooth; // mmks
double mks3R0, mks3H0, mks3MY1, mks3MY2, mks3MP0; // mks3
double startx[NDIM], stopx[NDIM], dx[NDIM];
double R0, Rin, Rout, Rh;

/*
 * Despite the name, this returns r, th coordinates for a KS or BL
 * coordinate system (since they're equal), from a set of "modified"
 * native coordinates X
 */
void bl_coord(double X[NDIM], double *r, double *th)
{
  *r = exp(X[1]);

  if (METRIC_eKS) {
    *r = exp(X[1]);
    *th = M_PI * X[2];
  } else if (METRIC_MKS3) {
    *r = exp(X[1]) + mks3R0;
    *th = (M_PI
        * (1.
            + 1. / tan((mks3H0 * M_PI) / 2.)
                * tan(
                    mks3H0 * M_PI
                        * (-0.5
                            + (mks3MY1
                                + (pow(2., mks3MP0) * (-mks3MY1 + mks3MY2))
                                    / pow(exp(X[1]) + mks3R0, mks3MP0))
                                * (1. - 2. * X[2]) + X[2])))) / 2.;
  } else if (METRIC_FMKS) {
    double thG = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
    double y = 2 * X[2] - 1.;
    double thJ = poly_norm * y
        * (1. + pow(y / poly_xt, poly_alpha) / (poly_alpha + 1.)) + 0.5 * M_PI;
    *th = thG + exp(mks_smooth * (startx[1] - X[1])) * (thJ - thG);
  } else {
    *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
  }
}

void bl_to_ks(double X[NDIM], double ucon_bl[NDIM], double ucon_ks[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);

  double trans[NDIM][NDIM];
  MUNULOOP
    trans[mu][nu] = delta(mu, nu);

  trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
  trans[3][1] = a / (r * r - 2. * r + a * a);

  MULOOP
    ucon_ks[mu] = 0.;
  MUNULOOP
    ucon_ks[mu] += trans[mu][nu] * ucon_bl[nu];
}

void ks_to_bl(double X[NDIM], double ucon_ks[NDIM], double ucon_bl[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);

  double trans[NDIM][NDIM], rev_trans[NDIM][NDIM];
  MUNULOOP
    trans[mu][nu] = delta(mu, nu);

  trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
  trans[3][1] = a / (r * r - 2. * r + a * a);

  invert_matrix(trans, rev_trans);

  MULOOP
    ucon_bl[mu] = 0.;
  MUNULOOP
    ucon_bl[mu] += rev_trans[mu][nu] * ucon_ks[nu];
}

/*
 * returns g_{munu} at location specified by X
 */
void gcov_func(double X[NDIM], double gcov[NDIM][NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);

  // compute ks metric
  double Gcov_ks[NDIM][NDIM];
  gcov_ks(r, th, Gcov_ks);

  // convert from ks metric to mks/mmks
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);

  MUNULOOP
  {
    gcov[mu][nu] = 0;
    for (int lam = 0; lam < NDIM; ++lam) {
      for (int kap = 0; kap < NDIM; ++kap) {
        gcov[mu][nu] += Gcov_ks[lam][kap] * dxdX[lam][mu] * dxdX[kap][nu];
      }
    }
  }
}

// compute KS metric at point (r,th) in KS coordinates (cyclic in t, ph)
inline void gcov_ks(double r, double th, double gcov[NDIM][NDIM])
{
  double cth = cos(th);
  double sth = sin(th);

  double s2 = sth * sth;
  double rho2 = r * r + a * a * cth * cth;

  MUNULOOP gcov[mu][nu] = 0.;
  // Compute KS metric from KS coordinates (cyclic in t,phi)
  gcov[0][0] = -1. + 2. * r / rho2;
  gcov[0][1] = 2. * r / rho2;
  gcov[0][3] = -2. * a * r * s2 / rho2;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = 1. + 2. * r / rho2;
  gcov[1][3] = -a * s2 * (1. + 2. * r / rho2);

  gcov[2][2] = rho2;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2));
}

inline void gcov_bl(double r, double th, double gcov[NDIM][NDIM])
{
  double sth, cth, s2, a2, r2, DD, mu;
  sth = fabs(sin(th));
  s2 = sth * sth;
  cth = cos(th);
  a2 = a * a;
  r2 = r * r;
  DD = 1. - 2. / r + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;

  MUNULOOP gcov[mu][nu] = 0.;
  // Compute BL metric from BL coordinates
  gcov[0][0] = -(1. - 2. / (r * mu));
  gcov[0][3] = -2. * a * s2 / (r * mu);
  gcov[3][0] = gcov[0][3];
  gcov[1][1] = mu / DD;
  gcov[2][2] = r2 * mu;
  gcov[3][3] = r2 * sth * sth * (1. + a2 / r2 + 2. * a2 * s2 / (r2 * r * mu));

}

void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{
  // Jacobian with respect to KS basis where X is given in
  // non-KS basis
  MUNULOOP
    dxdX[mu][nu] = 0.;

  if (METRIC_eKS) { //  && metric == 0 // TODO eKS switch
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][2] = M_PI;
    dxdX[3][3] = 1.;

  } else if (METRIC_MKS3) {

    // mks3 ..
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][1] = -(pow(2., -1. + mks3MP0) * exp(X[1]) * mks3H0 * mks3MP0
        * (mks3MY1 - mks3MY2) * pow(M_PI, 2)
        * pow(exp(X[1]) + mks3R0, -1 - mks3MP0) * (-1 + 2 * X[2]) * 1.
        / tan((mks3H0 * M_PI) / 2.)
        * pow(
            1.
                / cos(
                    mks3H0 * M_PI
                        * (-0.5
                            + (mks3MY1
                                + (pow(2, mks3MP0) * (-mks3MY1 + mks3MY2))
                                    / pow(exp(X[1]) + mks3R0, mks3MP0))
                                * (1 - 2 * X[2]) + X[2])),
            2));
    dxdX[2][2] = (mks3H0 * pow(M_PI, 2)
        * (1
            - 2
                * (mks3MY1
                    + (pow(2, mks3MP0) * (-mks3MY1 + mks3MY2))
                        / pow(exp(X[1]) + mks3R0, mks3MP0))) * 1.
        / tan((mks3H0 * M_PI) / 2.)
        * pow(
            1.
                / cos(
                    mks3H0 * M_PI
                        * (-0.5
                            + (mks3MY1
                                + (pow(2, mks3MP0) * (-mks3MY1 + mks3MY2))
                                    / pow(exp(X[1]) + mks3R0, mks3MP0))
                                * (1 - 2 * X[2]) + X[2])),
            2)) / 2.;
    dxdX[3][3] = 1.;

  } else if (METRIC_FMKS) {

    // fmks
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][1] = -exp(mks_smooth * (startx[1] - X[1])) * mks_smooth
        * (
        M_PI / 2. -
        M_PI * X[2]
            + poly_norm * (2. * X[2] - 1.)
                * (1
                    + (pow((-1. + 2 * X[2]) / poly_xt, poly_alpha))
                        / (1 + poly_alpha))
            - 1. / 2. * (1. - hslope) * sin(2. * M_PI * X[2]));
    dxdX[2][2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2])
        + exp(mks_smooth * (startx[1] - X[1]))
            * (-M_PI
                + 2. * poly_norm
                    * (1.
                        + pow((2. * X[2] - 1.) / poly_xt, poly_alpha)
                            / (poly_alpha + 1.))
                + (2. * poly_alpha * poly_norm * (2. * X[2] - 1.)
                    * pow((2. * X[2] - 1.) / poly_xt, poly_alpha - 1.))
                    / ((1. + poly_alpha) * poly_xt)
                - (1. - hslope) * M_PI * cos(2. * M_PI * X[2]));
    dxdX[3][3] = 1.;

  } else {
    // mks
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][2] = M_PI - (hslope - 1.) * M_PI * cos(2. * M_PI * X[2]);
    dxdX[3][3] = 1.;
  }
}

void set_dXdx(double X[NDIM], double dXdx[NDIM][NDIM]) {
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);
  invert_matrix(dxdX, dXdx);
}

void vec_to_ks(double X[NDIM], double v_nat[NDIM], double v_ks[NDIM]) {
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);

  MULOOP v_ks[mu] = 0.;
  MUNULOOP v_ks[mu] += dxdX[mu][nu] * v_nat[nu];
}

void vec_from_ks(double X[NDIM], double v_ks[NDIM], double v_nat[NDIM]) {
  double dXdx[NDIM][NDIM];
  set_dXdx(X, dXdx);

  MULOOP v_nat[mu] = 0.;
  MUNULOOP v_nat[mu] += dXdx[mu][nu] * v_ks[nu];
}

// Root-find the camera location in native coordinates
double root_find(double X[NDIM])
{
  double th = X[2];
  double tha, thb, thc;

  double Xa[NDIM], Xb[NDIM], Xc[NDIM];
  Xa[1] = log(X[1]);
  Xa[3] = X[3];
  Xb[1] = Xa[1];
  Xb[3] = Xa[3];
  Xc[1] = Xa[1];
  Xc[3] = Xa[3];

  if (X[2] < M_PI / 2.) {
    Xa[2] = 0.;
    Xb[2] = 0.5 + SMALL;
  } else {
    Xa[2] = 0.5 - SMALL;
    Xb[2] = 1.;
  }

  double tol = 1.e-9;
  tha = theta_func(Xa);
  thb = theta_func(Xb);

  // check limits first
  if (fabs(tha-th) < tol) {
    return Xa[2];
  } else if (fabs(thb-th) < tol) {
    return Xb[2];
  }

  // bisect for a bit
  for (int i = 0; i < 1000; i++) {
    Xc[2] = 0.5 * (Xa[2] + Xb[2]);
    thc = theta_func(Xc);

    if ((thc - th) * (thb - th) < 0.)
      Xa[2] = Xc[2];
    else
      Xb[2] = Xc[2];
    
    double err = theta_func(Xc) - th;
    if (fabs(err) < tol)
      break;
  }

  return Xc[2];
}
