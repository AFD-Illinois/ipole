// Coordinate-dependent (but not radiation model-dependent) functions for ipole

#include "coordinates.h"

#include "decs.h"
#include "geometry.h"

int use_eKS_internal = 0;
int metric = -1;
double a, hslope; // mks
double poly_norm, poly_xt, poly_alpha, mks_smooth; // fmks
double mks3R0, mks3H0, mks3MY1, mks3MY2, mks3MP0; // mks3

// Coordinate parameters
double startx[NDIM], stopx[NDIM], dx[NDIM];
double cstartx[NDIM], cstopx[NDIM];
double R0, Rin, Rout, Rh;
// Tracing parameters we need independent of model
double rmax_geo = 100.;
double rmin_geo = 1.;

/*
 * Despite the name, this returns r, th coordinates for a KS or BL
 * coordinate system (since they're equal), from a set of "modified"
 * native coordinates X
 * If METRIC_MINKOWSKI is set, returns spherical coordinates instead
 */
void bl_coord(double X[NDIM], double *r, double *th)
{

  if (metric == METRIC_MINKOWSKI) {
    *r = X[1];
    *th = X[2];
    return;
  } else if (metric == METRIC_EMINKOWSKI) {
    *r = exp(X[1]);
    *th = X[2];
    return;
  }

  *r = exp(X[1]);

  if (use_eKS_internal) {
    *th = M_PI * X[2];
  } else {
    double y, thG, thJ;
    switch (metric) {
      case METRIC_EKS:
        *th = X[2];
        break;
      case METRIC_MKS:
        *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
        break;
      case METRIC_BHACMKS:
        *th = X[2] + (hslope / 2.) * sin(2. * X[2]);
        break;
      case METRIC_FMKS:
        thG = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
        y = 2 * X[2] - 1.;
        thJ = poly_norm * y
            * (1. + pow(y / poly_xt, poly_alpha) / (poly_alpha + 1.)) + 0.5 * M_PI;
        *th = thG + exp(mks_smooth * (startx[1] - X[1])) * (thJ - thG);
        break;
      case METRIC_MKS3:
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
        break;
    }
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

  if (metric == METRIC_MINKOWSKI) {
    MUNULOOP gcov[mu][nu] = 0;
    gcov[0][0] = -1;
    gcov[1][1] = 1;
    gcov[2][2] = r*r;
    gcov[3][3] = r*r*sin(th)*sin(th);
    return;
  } else if (metric == METRIC_EMINKOWSKI) {
    MUNULOOP gcov[mu][nu] = 0;
    gcov[0][0] = -1;
    gcov[1][1] = r*r;
    gcov[2][2] = r*r;
    gcov[3][3] = r*r*sin(th)*sin(th);
    return;
  }

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

//Variables used in EDGB metric
  double cost2 = cth*cth;
  double cost4 = cost2*cost2;
  //double a = 1;
  //double pow(a,2) = a * a;
  //double pow(a,4) = pow(a,2)*pow(a,2);
  double zeta = 0.01;
  MUNULOOP gcov[mu][nu] = 0.;
  // Compute KS metric from KS coordinates (cyclic in t,phi)
  
  gcov[0][0] = -1. + 2. *r/rho2;

  //EDGB metric
  /*-1 + 2/r - (2*pow(a,2)*cost2)/pow(r,3) + 2*pow(a,4)*cost4/pow(r,5) + zeta*(-1/15*(-400+96*r+66*pow(r,2)+130*pow(r,3)+5*pow(r,4))/pow(r,7)
    + (pow(a,2)*(pow(r,7)*(444696-562338*cost2)  + 8820000*(-1 + 3*cost2) - 19600*r*(-467 + 1251*cost2) - 63*pow(r,8)*(-3267 + 8926*cost2) 
    + 1050*pow(r,3)*(-1465 + 11997*cost2) - 2100*pow(r,2)*(-955 + 22577*cost2) - 6*pow(r,6)*(-59329 + 82437*cost2) + 15*pow(r,5)*(-52533 + 455029*cost2) 
    + 10*pow(r,4)*(-281221 + 1218513*cost2)))/(110250*pow(r,11)) + (pow(a,4)*(675*pow(r,12)*(-19717 + 67726*cost2) + 450*pow(r,11)*(43312 + 101589*cost2) 
    - 164640000*(-70 + 585*cost2 + 156*cost4) - 8232000*r*(3949 - 11058*cost2 + 6948*cost4) - 39200*pow(r,2)*(189191 - 824825*cost2 + 972045*cost4) 
    + 60*pow(r,9)*(717867 - 13885852*cost2 + 12733507*cost4) + 30*pow(r,10)*(209773 - 9090216*cost2 + 16888370*cost4) + 1400*pow(r,3)*(648009 - 14691730*cost2 
    + 26074500*cost4) + 420*pow(r,4)*(553219 - 32471380*cost2 + 222891320*cost4) - 14*pow(r,7)*(-11393603 + 38599350*cost2 + 359928985*cost4) 
    + 2*pow(r,8)*(59964679 - 491173980*cost2 + 452824800*cost4) - 28*pow(r,5)*(21852807 + 12094180*cost2 + 762315255*cost4) 
    - 14*pow(r,6)*(42004007 - 226439060*cost2 + 1041312310*cost4)))/(30870000*pow(r,15)));*/
  //Original value:
  //-1. + 2. * r / rho2;

  gcov[0][1] =2. * r / rho2 +0.0001;

  //EDGB metric
  /*2/r - (2*pow(a,2)*cost2)/pow(r,3) + (2*pow(a,4)*cost4)/pow(r,5) + zeta*(-1/30*(-800 + 912*r + 516*pow(r,2) + 470*pow(r,3) + 50*pow(r,4) + 15*pow(r,5))/ pow(r,7) 
    + (pow(a,2)*(-566034*pow(r,8) + 55125*pow(r,9) +  17640000*(-1 + 3*cost2) + 12600*pow(r,3)*(4006 + 1877*cost2) + 78400*r*(-779 + 2412*cost2) + 42*pow(r,7)*(-94339 
    + 107112*cost2) - 1400*pow(r,2)*(20431 + 132243*cost2) + 36*pow(r,6)*(-436917 + 491281*cost2) + 20*pow(r,4)*(-875941 + 2263053*cost2) + 20*pow(r,5)*(-1122937 
    + 2632446*cost2)))/(220500*pow(r,11)) - (pow(a,4)*(-80247150*pow(r,12) - 5788125*pow(r,13) + 450*pow(r,11)*(-1196839 + 812712*cost2) 
    + 329280000*(-70 + 585*cost2 + 156*cost4) + 49392000*r*(-1717 + 21664*cost2 + 9076*cost4) + 60*pow(r,10)*(-23601289 + 13406112*cost2 + 2187250*cost4) 
    + 78400*pow(r,2)*(-2206069 - 11318105*cost2 + 13657725*cost4) + 14000*pow(r,3)*(22540153 - 88517480*cost2 + 62290230*cost4) 
    + 30*pow(r,9)*(-145291221 + 30934768*cost2 + 146683252*cost4) + 280*pow(r,4)*(793805393 - 2014699860*cost2 + 507428040*cost4) 
    + 28*pow(r,6)*(-711534337 - 1340707620*cost2 + 4191169150*cost4) + 28*pow(r,5)*(2484408549 - 2498856340*cost2 + 4680981810*cost4) 
    + 14*pow(r,7)*(-1576790517 - 1991632680*cost2 + 4890821060*cost4) + 4*pow(r,8)*(-2741883352 - 2122491690*cost2 + 5145464715*cost4)))/(61740000*pow(r,15)));*/
  //Original value:
  //2. * r / rho2 +0.0001;

  gcov[0][3] = -2. * a * r * s2 / rho2+0.0001;

  //EDGB metric
  /*(2*a*(-1 + cost2))/r - (2*a*pow(a,2)*cost2*(-1 + cost2))/pow(r,3) + zeta*(-1/15*((-400 + 144*r + 90*pow(r,2) + 140*pow(r,3) 
    + 9*pow(r,4))*a*(-1 + cost2))/pow(r,7) - (pow(a,2)*a*(-1 + cost2)* (pow(r,4)*(2736210 - 4410530*cost2) + pow(r,5)*(766015 - 3620183*cost2) - 
    8820000*(-1 + 3*cost2) + 19600*r*(-467 + 1551*cost2) - 12*pow(r,6)*(26511 +  6310*cost2) + 750*pow(r,3)*(2051 + 8733*cost2) +  2100*pow(r,2)*(-955 + 21233*cost2) 
    + 3*pow(r,8)*(-63529 + 262520*cost2) + pow(r,7)*(-406611 + 563055*cost2)))/(110250*pow(r,11)));*/
  //Original value:
  //-2. * a * r * s2 / rho2+0.0001;

  gcov[1][0] = gcov[0][1];

  gcov[1][1] = 1. + 2. * r / rho2;

  //EDGB metric
  /*(2 + r)/r - (2*pow(a,2)*cost2)/pow(r,3) + (2*pow(a,4)*cost4)/pow(r,5) + zeta*(-1/15*(-400 + 816*r + 450*pow(r,2) + 340*pow(r,3) + 45*pow(r,4) + 15*pow(r,5))/pow(r,7) 
    + (pow(a,2)*(55125*pow(r,9) + 8820000*(-1 + 3*cost2) + 1050*pow(r,3)*(49537 + 10527*cost2) + 19600*r*(-3583 + 10899*cost2) + 21*pow(r,8)*(-36755 + 26778*cost2) 
    + 42*pow(r,7)*(-104927 + 120501*cost2) - 700*pow(r,2)*(43727 + 196755*cost2) + 6*pow(r,6)*(-2680831 + 3030123*cost2) + 10*pow(r,4)*(-1470661 + 3307593*cost2) 
    + 5*pow(r,5)*(-4334149 + 9164697*cost2)))/(110250*pow(r,11)) - (pow(a,4)*(-5788125*pow(r,13) + 225*pow(r,12)*(-415805 + 203178*cost2) 
    + 1350*pow(r,11)*(-384509 + 304767*cost2) + 164640000*(-70 + 585*cost2 + 156*cost4) + 8232000*r*(-14251 + 141042*cost2 + 47508*cost4) 
    + 30*pow(r,10)*(-46992805 + 17722008*cost2 + 21262870*cost4) + 39200*pow(r,2)*(-4601329 - 21811385*cost2 + 26343405*cost4) 
    + 30*pow(r,9)*(-143855487 + 3163064*cost2 + 172150266*cost4) + 1400*pow(r,3)*(226049539 - 899866530*cost2 + 648976800*cost4) 
    + 140*pow(r,4)*(1589270443 - 4126813860*cost2 + 1683530040*cost4) + 10*pow(r,8)*(-1084760405 - 947231472*cost2 + 2148750846*cost4) 
    + 28*pow(r,5)*(2462555742 - 2510950520*cost2 + 3918666555*cost4) + 14*pow(r,7)*(-1565396914 - 2030232030*cost2 + 4530892075*cost4) 
    + 14*pow(r,6)*(-1465072681 - 2454976180*cost2 + 7341025990*cost4)))/(30870000*pow(r,15)));*/
  //Original value:
  //1. + 2. * r / rho2;

  gcov[1][3] = -a * s2 * (1. + 2. * r / rho2); 

  //EDGB metric
  /*((2 + r)*a*(-1 + cost2))/r - (2*a*pow(a,2)*cost2*(-1 + cost2))/pow(r,3) + zeta*(-1/36750*((16660000 - 5350800*r + 4797450*pow(r,2) + 3526600*pow(r,3) 
    + 2965560*pow(r,4) + 918855*pow(r,5) + 187446*pow(r,6))*a*(-1 + cost2))/pow(r,7) + (a*pow(a,2)*(-1 + cost2)*(22085775*pow(r,9) 
    + 4571505*pow(r,10) + 49392000*(580 + 327*cost2) + 548800*r*(-23041 + 74715*cost2) + 6300*pow(r,3)*(-446807 + 973501*cost2) 
    + 126*pow(r,7)*(1064483 + 1485790*cost2) + 9800*pow(r,2)*(-1223527 + 1991748*cost2) + 12*pow(r,8)*(4458631 + 3456783*cost2) 
    + 280*pow(r,4)*(3773463 + 15733496*cost2) + 42*pow(r,6)*(6767669 + 23527525*cost2) + 56*pow(r,5)*(17855552 + 49207893*cost2)))/(3087000*pow(r,11)));*/
  //Original value:
  //-a * s2 * (1. + 2. * r / rho2);

  gcov[2][2] = rho2;

  //EDGB metric
  /*pow(r,2) + pow(a,2)*cost2 + ((8820000 - 6213200*r - 3416700*pow(r,2) - 1855650*pow(r,3) + 887110*pow(r,4) + 800733*pow(r,5) + 435540*pow(r,6) 
   + 187446*pow(r,7))*pow(a,2)*zeta*(1 - 3*cost2))/(110250*pow(r,8)) + (pow(a,4)*zeta*(45715050*pow(r,11)*(-1 + 3*cost2) + 5625*pow(r,10)*(-20749 + 58131*cost2) 
   + 493920000*(-70 + 585*cost2 + 156*cost4) + 24696000*r*(3049 - 10698*cost2 + 8868*cost4) + 117600*pow(r,2)*(280331 - 1711445*cost2 + 1596165*cost4) 
   + 180*pow(r,9)*(-1286466 - 846865*cost2 + 5819941*cost4) + 4200*pow(r,3)*(2362411 - 16650910*cost2 + 14489100*cost4) 
   - 1260*pow(r,4)*(-3173281 - 5026080*cost2 + 26477920*cost4) + 42*pow(r,8)*(-18071967 - 940590*cost2 + 54146980*cost4) 
   + 42*pow(r,6)*(-19116713 - 46592740*cost2 + 138130070*cost4) - 28*pow(r,5)*(11804979 - 261030540*cost2 + 235282135*cost4) 
   + 6*pow(r,7)*(-259078241 - 99440670*cost2 + 857000595*cost4)))/(92610000*pow(r,12));*/
  //Original value:
  //rho2;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];

  gcov[3][3] = s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2));

  //EDGB metric
  /*-(pow(r,2)*(-1 + cost2)) - (pow(a,2)*(2 + r - 2*cost2)*(-1 + cost2))/r - (2*pow(a,4)*cost2*pow((-1 + cost2),2))/pow(r,3) 
    + zeta*(((8820000 - 6213200*r - 3416700*pow(r,2) - 1855650*pow(r,3) + 887110*pow(r,4) + 800733*pow(r,5) + 435540*pow(r,6) 
    + 187446*pow(r,7))*pow(a,2)* (-1 + cost2)*(-1 + 3*cost2))/(110250*pow(r,8)) -  (pow(a,4)*(-1 + cost2)*(45715050*pow(r,11)*(-1 + 3*cost2) 
    + 5625*pow(r,10)*(-20749 + 58131*cost2) + 493920000*(-70 + 585*cost2 + 156*cost4) + 24696000*r*(3649 - 8958*cost2 + 6528*cost4) 
    + 352800*pow(r,2)*(84857 - 350495*cost2 + 320655*cost4) + 12600*pow(r,3)*(-82303 - 1443030*cost2 + 1592200*cost4) 
    + 180*pow(r,9)*(-411718 - 4345857*cost2 + 8444185*cost4) - 1260*pow(r,4)*(1578719 - 11450880*cost2 + 28150720*cost4) 
    + 42*pow(r,8)*(-1863327 - 67980150*cost2 + 104977900*cost4) + 28*pow(r,5)*(-14247879 - 109560360*cost2 + 137751665*cost4) 
    + 42*pow(r,6)*(30654807 - 316973820*cost2 + 358739630*cost4) + 6*pow(r,7)*(-25024421 - 1143700950*cost2 + 1667207055*cost4)))/(92610000*pow(r,12)));*/
  //Original value: 
  //s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2));
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
    dxdX[mu][nu] = delta(mu, nu);

  dxdX[1][1] = exp(X[1]); // Overridden by one case below

  if (use_eKS_internal) {
    dxdX[2][2] = M_PI;
  } else {
    switch (metric) {
      case METRIC_EKS:
        break;
      case METRIC_MKS:
        dxdX[2][2] = M_PI + (1 - hslope) * M_PI * cos(2. * M_PI * X[2]);
        break;
      case METRIC_BHACMKS:
        dxdX[2][2] = 1 + hslope * cos(2. * X[2]);
        break;
      case METRIC_FMKS:
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
        break;
      case METRIC_MKS3:
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
        break;
      case METRIC_MINKOWSKI:
        // Blank transform: just override L_11
        dxdX[1][1] = 1.;
        break;
      case METRIC_EMINKOWSKI:
        // keep radial transformation element!
        break;
    }
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

/*
 * Translate the input camera angles into a canonical Xcam in native coordinates
 */
void native_coord(double r, double th, double phi, double X[NDIM]) {
  if (metric == METRIC_MINKOWSKI) {
    X[0] = 1; X[1] = r; X[2] = th/180*M_PI; X[3] = phi/180*M_PI;
  } else if (metric == METRIC_EMINKOWSKI) {
    X[0] = 1; X[1] = log(r); X[2] = th/180*M_PI; X[3] = phi/180*M_PI;
  } else {
    double x[NDIM] = {0., r, th/180.*M_PI, phi/180.*M_PI};
    X[0] = 0.0;
    X[1] = log(r);
    X[2] = root_find(x);
    X[3] = phi/180.*M_PI;
  }
}

/**
 * Root-find the camera theta in native coordinates
 * 
 * TODO switch this to native_coord func above...
 */
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
    Xa[2] = cstartx[2];
    Xb[2] = (cstopx[2] - cstartx[2])/2 + SMALL;
  } else {
    Xa[2] = (cstopx[2] - cstartx[2])/2 - SMALL;
    Xb[2] = cstopx[2];
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
