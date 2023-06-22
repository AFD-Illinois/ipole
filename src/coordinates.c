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

  //double a = 1;
  //double pow(a,2) = a * a;
  //double pow(a,4) = pow(a,2)*pow(a,2);
  double zeta = 0.3;
  MUNULOOP gcov[mu][nu] = 0.;
  // Compute KS metric from KS coordinates (cyclic in t,phi)
  
  gcov[0][0] = -1. + 2./r - (2.*pow(a,2)*pow(cth,2))/pow(r,3) + 2.*pow(a,4)*pow(cth,4)/pow(r,5) + zeta*(-1./15.*(-400.+96.*r+66.*pow(r,2)+130.*pow(r,3)+5.*pow(r,4))/pow(r,7)
    + (pow(a,2)*(pow(r,7)*(444696.-562338.*pow(cth,2))  + 8820000.*(-1. + 3.*pow(cth,2)) - 19600.*r*(-467. + 1251.*pow(cth,2)) - 63.*pow(r,8)*(-3267. + 8926.*pow(cth,2)) 
    + 1050.*pow(r,3)*(-1465. + 11997.*pow(cth,2)) - 2100.*pow(r,2)*(-955. + 22577.*pow(cth,2)) - 6.*pow(r,6)*(-59329. + 82437.*pow(cth,2)) + 15.*pow(r,5)*(-52533. + 455029.*pow(cth,2)) 
    + 10.*pow(r,4)*(-281221. + 1218513.*pow(cth,2))))/(110250.*pow(r,11)) + (pow(a,4)*(675.*pow(r,12)*(-19717. + 67726.*pow(cth,2)) + 450.*pow(r,11)*(43312. + 101589.*pow(cth,2)) 
    - 164640000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) - 8232000.*r*(3949. - 11058.*pow(cth,2) + 6948.*pow(cth,4)) - 39200.*pow(r,2)*(189191. - 824825.*pow(cth,2) + 972045.*pow(cth,4)) 
    + 60.*pow(r,9)*(717867. - 13885852.*pow(cth,2) + 12733507.*pow(cth,4)) + 30.*pow(r,10)*(209773. - 9090216.*pow(cth,2) + 16888370.*pow(cth,4)) + 1400.*pow(r,3)*(648009. - 14691730.*pow(cth,2) 
    + 26074500.*pow(cth,4)) + 420.*pow(r,4)*(553219. - 32471380.*pow(cth,2) + 222891320.*pow(cth,4)) - 14.*pow(r,7)*(-11393603. + 38599350.*pow(cth,2) + 359928985.*pow(cth,4)) 
    + 2.*pow(r,8)*(59964679. - 491173980.*pow(cth,2) + 452824800.*pow(cth,4)) - 28.*pow(r,5)*(21852807. + 12094180.*pow(cth,2) + 762315255.*pow(cth,4)) 
    - 14.*pow(r,6)*(42004007. - 226439060.*pow(cth,2) + 1041312310.*pow(cth,4))))/(30870000.*pow(r,15)));

  //EDGB metric
  /*-1. + 2./r - (2.*pow(a,2)*pow(cth,2))/pow(r,3) + 2.*pow(a,4)*pow(cth,4)/pow(r,5) + zeta*(-1./15.*(-400.+96.*r+66.*pow(r,2)+130.*pow(r,3)+5.*pow(r,4))/pow(r,7)
    + (pow(a,2)*(pow(r,7)*(444696.-562338.*pow(cth,2))  + 8820000.*(-1. + 3.*pow(cth,2)) - 19600.*r*(-467. + 1251.*pow(cth,2)) - 63.*pow(r,8)*(-3267. + 8926.*pow(cth,2)) 
    + 1050.*pow(r,3)*(-1465. + 11997.*pow(cth,2)) - 2100.*pow(r,2)*(-955. + 22577.*pow(cth,2)) - 6.*pow(r,6)*(-59329. + 82437.*pow(cth,2)) + 15.*pow(r,5)*(-52533. + 455029.*pow(cth,2)) 
    + 10.*pow(r,4)*(-281221. + 1218513.*pow(cth,2))))/(110250.*pow(r,11)) + (pow(a,4)*(675.*pow(r,12)*(-19717. + 67726.*pow(cth,2)) + 450.*pow(r,11)*(43312. + 101589.*pow(cth,2)) 
    - 164640000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) - 8232000.*r*(3949. - 11058.*pow(cth,2) + 6948.*pow(cth,4)) - 39200.*pow(r,2)*(189191. - 824825.*pow(cth,2) + 972045.*pow(cth,4)) 
    + 60.*pow(r,9)*(717867. - 13885852.*pow(cth,2) + 12733507.*pow(cth,4)) + 30.*pow(r,10)*(209773. - 9090216.*pow(cth,2) + 16888370.*pow(cth,4)) + 1400.*pow(r,3)*(648009. - 14691730.*pow(cth,2) 
    + 26074500.*pow(cth,4)) + 420.*pow(r,4)*(553219. - 32471380.*pow(cth,2) + 222891320.*pow(cth,4)) - 14.*pow(r,7)*(-11393603. + 38599350.*pow(cth,2) + 359928985.*pow(cth,4)) 
    + 2.*pow(r,8)*(59964679. - 491173980.*pow(cth,2) + 452824800.*pow(cth,4)) - 28.*pow(r,5)*(21852807. + 12094180.*pow(cth,2) + 762315255.*pow(cth,4)) 
    - 14.*pow(r,6)*(42004007. - 226439060.*pow(cth,2) + 1041312310.*pow(cth,4))))/(30870000.*pow(r,15)));*/
  //Original value:
  //-1. + 2. * r / rho2;

  gcov[0][1] = 2./r - (2.*pow(a,2)*pow(cth,2))/pow(r,3) + (2.*pow(a,4)*pow(cth,4))/pow(r,5) + zeta*(-1./30.*(-800. + 912.*r + 516.*pow(r,2) + 470.*pow(r,3) + 50.*pow(r,4) + 15.*pow(r,5))/ pow(r,7) 
    + (pow(a,2)*(-566034.*pow(r,8) + 55125.*pow(r,9) +  17640000.*(-1. + 3.*pow(cth,2)) + 12600.*pow(r,3)*(4006. + 1877.*pow(cth,2)) + 78400.*r*(-779. + 2412.*pow(cth,2)) + 42.*pow(r,7)*(-94339. 
    + 107112.*pow(cth,2)) - 1400.*pow(r,2)*(20431. + 132243.*pow(cth,2)) + 36.*pow(r,6)*(-436917. + 491281.*pow(cth,2)) + 20.*pow(r,4)*(-875941. + 2263053.*pow(cth,2)) + 20.*pow(r,5)*(-1122937. 
    + 2632446.*pow(cth,2))))/(220500.*pow(r,11)) - (pow(a,4)*(-80247150.*pow(r,12) - 5788125.*pow(r,13) + 450.*pow(r,11)*(-1196839. + 812712.*pow(cth,2)) 
    + 329280000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) + 49392000.*r*(-1717. + 21664.*pow(cth,2) + 9076.*pow(cth,4)) + 60.*pow(r,10)*(-23601289. + 13406112.*pow(cth,2) + 2187250.*pow(cth,4)) 
    + 78400.*pow(r,2)*(-2206069. - 11318105.*pow(cth,2) + 13657725.*pow(cth,4)) + 14000.*pow(r,3)*(22540153. - 88517480.*pow(cth,2) + 62290230.*pow(cth,4)) 
    + 30.*pow(r,9)*(-145291221. + 30934768.*pow(cth,2) + 146683252.*pow(cth,4)) + 280.*pow(r,4)*(793805393. - 2014699860.*pow(cth,2) + 507428040.*pow(cth,4)) 
    + 28.*pow(r,6)*(-711534337. - 1340707620.*pow(cth,2) + 4191169150.*pow(cth,4)) + 28.*pow(r,5)*(2484408549. - 2498856340.*pow(cth,2) + 4680981810.*pow(cth,4)) 
    + 14.*pow(r,7)*(-1576790517. - 1991632680.*pow(cth,2) + 4890821060.*pow(cth,4)) + 4.*pow(r,8)*(-2741883352. - 2122491690.*pow(cth,2) + 5145464715.*pow(cth,4))))/(61740000.*pow(r,15)));

  //EDGB metric
  /*2./r - (2.*pow(a,2)*pow(cth,2))/pow(r,3) + (2.*pow(a,4)*pow(cth,4))/pow(r,5) + zeta*(-1./30.*(-800. + 912.*r + 516.*pow(r,2) + 470.*pow(r,3) + 50.*pow(r,4) + 15.*pow(r,5))/ pow(r,7) 
    + (pow(a,2)*(-566034.*pow(r,8) + 55125.*pow(r,9) +  17640000.*(-1. + 3.*pow(cth,2)) + 12600.*pow(r,3)*(4006. + 1877.*pow(cth,2)) + 78400.*r*(-779. + 2412.*pow(cth,2)) + 42.*pow(r,7)*(-94339. 
    + 107112.*pow(cth,2)) - 1400.*pow(r,2)*(20431. + 132243.*pow(cth,2)) + 36.*pow(r,6)*(-436917. + 491281.*pow(cth,2)) + 20.*pow(r,4)*(-875941. + 2263053.*pow(cth,2)) + 20.*pow(r,5)*(-1122937. 
    + 2632446.*pow(cth,2))))/(220500.*pow(r,11)) - (pow(a,4)*(-80247150.*pow(r,12) - 5788125.*pow(r,13) + 450.*pow(r,11)*(-1196839. + 812712.*pow(cth,2)) 
    + 329280000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) + 49392000.*r*(-1717. + 21664.*pow(cth,2) + 9076.*pow(cth,4)) + 60.*pow(r,10)*(-23601289. + 13406112.*pow(cth,2) + 2187250.*pow(cth,4)) 
    + 78400.*pow(r,2)*(-2206069. - 11318105.*pow(cth,2) + 13657725.*pow(cth,4)) + 14000.*pow(r,3)*(22540153. - 88517480.*pow(cth,2) + 62290230.*pow(cth,4)) 
    + 30.*pow(r,9)*(-145291221. + 30934768.*pow(cth,2) + 146683252.*pow(cth,4)) + 280.*pow(r,4)*(793805393. - 2014699860.*pow(cth,2) + 507428040.*pow(cth,4)) 
    + 28.*pow(r,6)*(-711534337. - 1340707620.*pow(cth,2) + 4191169150.*pow(cth,4)) + 28.*pow(r,5)*(2484408549. - 2498856340.*pow(cth,2) + 4680981810.*pow(cth,4)) 
    + 14.*pow(r,7)*(-1576790517. - 1991632680.*pow(cth,2) + 4890821060.*pow(cth,4)) + 4.*pow(r,8)*(-2741883352. - 2122491690.*pow(cth,2) + 5145464715.*pow(cth,4))))/(61740000.*pow(r,15)));*/
  //Original value:
  //2. * r / rho2 +0.0001;

  gcov[0][3] = (2.*a*(-1. + pow(cth,2)))/r - (2.*a*pow(a,2)*pow(cth,2)*(-1. + pow(cth,2)))/pow(r,3) + zeta*(-1./15.*((-400. + 144.*r + 90.*pow(r,2) + 140.*pow(r,3) 
    + 9.*pow(r,4))*a*(-1. + pow(cth,2)))/pow(r,7) - (pow(a,2)*a*(-1. + pow(cth,2))* (pow(r,4)*(2736210. - 4410530.*pow(cth,2)) + pow(r,5)*(766015. - 3620183.*pow(cth,2)) - 
    8820000.*(-1. + 3.*pow(cth,2)) + 19600.*r*(-467. + 1551.*pow(cth,2)) - 12.*pow(r,6)*(26511. +  6310.*pow(cth,2)) + 750.*pow(r,3)*(2051. + 8733.*pow(cth,2)) +  2100.*pow(r,2)*(-955. + 21233.*pow(cth,2)) 
    + 3.*pow(r,8)*(-63529. + 262520.*pow(cth,2)) + pow(r,7)*(-406611. + 563055.*pow(cth,2))))/(110250.*pow(r,11)));

  //EDGB metric
  /*(2.*a*(-1. + pow(cth,2)))/r - (2.*a*pow(a,2)*pow(cth,2)*(-1. + pow(cth,2)))/pow(r,3) + zeta*(-1./15.*((-400. + 144.*r + 90.*pow(r,2) + 140.*pow(r,3) 
    + 9.*pow(r,4))*a*(-1. + pow(cth,2)))/pow(r,7) - (pow(a,2)*a*(-1. + pow(cth,2))* (pow(r,4)*(2736210. - 4410530.*pow(cth,2)) + pow(r,5)*(766015. - 3620183.*pow(cth,2)) - 
    8820000.*(-1. + 3.*pow(cth,2)) + 19600.*r*(-467. + 1551.*pow(cth,2)) - 12.*pow(r,6)*(26511. +  6310.*pow(cth,2)) + 750.*pow(r,3)*(2051. + 8733.*pow(cth,2)) +  2100.*pow(r,2)*(-955. + 21233.*pow(cth,2)) 
    + 3.*pow(r,8)*(-63529. + 262520.*pow(cth,2)) + pow(r,7)*(-406611. + 563055.*pow(cth,2))))/(110250.*pow(r,11)));*/
  //Original value:
  //-2. * a * r * s2 / rho2;

  gcov[1][0] = gcov[0][1];

  gcov[1][1] = (2. + r)/r - (2.*pow(a,2)*pow(cth,2))/pow(r,3) + (2.*pow(a,4)*pow(cth,4))/pow(r,5) + zeta*(-1./15.*(-400. + 816.*r + 450.*pow(r,2) + 340.*pow(r,3) + 45.*pow(r,4) + 15.*pow(r,5))/pow(r,7) 
    + (pow(a,2)*(55125.*pow(r,9) + 8820000.*(-1. + 3.*pow(cth,2)) + 1050.*pow(r,3)*(49537. + 10527.*pow(cth,2)) + 19600.*r*(-3583. + 10899.*pow(cth,2)) + 21.*pow(r,8)*(-36755. + 26778.*pow(cth,2)) 
    + 42.*pow(r,7)*(-104927. + 120501.*pow(cth,2)) - 700.*pow(r,2)*(43727. + 196755.*pow(cth,2)) + 6.*pow(r,6)*(-2680831. + 3030123.*pow(cth,2)) + 10.*pow(r,4)*(-1470661. + 3307593.*pow(cth,2)) 
    + 5.*pow(r,5)*(-4334149. + 9164697.*pow(cth,2))))/(110250.*pow(r,11)) - (pow(a,4)*(-5788125.*pow(r,13) + 225.*pow(r,12)*(-415805. + 203178.*pow(cth,2)) 
    + 1350.*pow(r,11)*(-384509. + 304767.*pow(cth,2)) + 164640000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) + 8232000.*r*(-14251. + 141042.*pow(cth,2) + 47508.*pow(cth,4)) 
    + 30.*pow(r,10)*(-46992805. + 17722008.*pow(cth,2) + 21262870.*pow(cth,4)) + 39200.*pow(r,2)*(-4601329. - 21811385.*pow(cth,2) + 26343405.*pow(cth,4)) 
    + 30.*pow(r,9)*(-143855487. + 3163064.*pow(cth,2) + 172150266.*pow(cth,4)) + 1400.*pow(r,3)*(226049539. - 899866530.*pow(cth,2) + 648976800.*pow(cth,4)) 
    + 140.*pow(r,4)*(1589270443. - 4126813860.*pow(cth,2) + 1683530040.*pow(cth,4)) + 10.*pow(r,8)*(-1084760405. - 947231472.*pow(cth,2) + 2148750846.*pow(cth,4)) 
    + 28.*pow(r,5)*(2462555742. - 2510950520.*pow(cth,2) + 3918666555.*pow(cth,4)) + 14.*pow(r,7)*(-1565396914. - 2030232030.*pow(cth,2) + 4530892075.*pow(cth,4)) 
    + 14.*pow(r,6)*(-1465072681. - 2454976180.*pow(cth,2) + 7341025990.*pow(cth,4))))/(30870000.*pow(r,15)));

  //EDGB metric
  /*(2. + r)/r - (2.*pow(a,2)*pow(cth,2))/pow(r,3) + (2.*pow(a,4)*pow(cth,4))/pow(r,5) + zeta*(-1./15.*(-400. + 816.*r + 450.*pow(r,2) + 340.*pow(r,3) + 45.*pow(r,4) + 15.*pow(r,5))/pow(r,7) 
    + (pow(a,2)*(55125.*pow(r,9) + 8820000.*(-1. + 3.*pow(cth,2)) + 1050.*pow(r,3)*(49537. + 10527.*pow(cth,2)) + 19600.*r*(-3583. + 10899.*pow(cth,2)) + 21.*pow(r,8)*(-36755. + 26778.*pow(cth,2)) 
    + 42.*pow(r,7)*(-104927. + 120501.*pow(cth,2)) - 700.*pow(r,2)*(43727. + 196755.*pow(cth,2)) + 6.*pow(r,6)*(-2680831. + 3030123.*pow(cth,2)) + 10.*pow(r,4)*(-1470661. + 3307593.*pow(cth,2)) 
    + 5.*pow(r,5)*(-4334149. + 9164697.*pow(cth,2))))/(110250.*pow(r,11)) - (pow(a,4)*(-5788125.*pow(r,13) + 225.*pow(r,12)*(-415805. + 203178.*pow(cth,2)) 
    + 1350.*pow(r,11)*(-384509. + 304767.*pow(cth,2)) + 164640000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) + 8232000.*r*(-14251. + 141042.*pow(cth,2) + 47508.*pow(cth,4)) 
    + 30.*pow(r,10)*(-46992805. + 17722008.*pow(cth,2) + 21262870.*pow(cth,4)) + 39200.*pow(r,2)*(-4601329. - 21811385.*pow(cth,2) + 26343405.*pow(cth,4)) 
    + 30.*pow(r,9)*(-143855487. + 3163064.*pow(cth,2) + 172150266.*pow(cth,4)) + 1400.*pow(r,3)*(226049539. - 899866530.*pow(cth,2) + 648976800.*pow(cth,4)) 
    + 140.*pow(r,4)*(1589270443. - 4126813860.*pow(cth,2) + 1683530040.*pow(cth,4)) + 10.*pow(r,8)*(-1084760405. - 947231472.*pow(cth,2) + 2148750846.*pow(cth,4)) 
    + 28.*pow(r,5)*(2462555742. - 2510950520.*pow(cth,2) + 3918666555.*pow(cth,4)) + 14.*pow(r,7)*(-1565396914. - 2030232030.*pow(cth,2) + 4530892075.*pow(cth,4)) 
    + 14.*pow(r,6)*(-1465072681. - 2454976180.*pow(cth,2) + 7341025990.*pow(cth,4))))/(30870000.*pow(r,15)));*/
  //Original value:
  //1. + 2. * r / rho2;

  gcov[1][3] = ((2. + r)*a*(-1. + pow(cth,2)))/r - (2.*a*pow(a,2)*pow(cth,2)*(-1. + pow(cth,2)))/pow(r,3) + zeta*(-1./36750.*((16660000. - 5350800.*r + 4797450.*pow(r,2) + 3526600.*pow(r,3) 
    + 2965560.*pow(r,4) + 918855.*pow(r,5) + 187446.*pow(r,6))*a*(-1. + pow(cth,2)))/pow(r,7) + (a*pow(a,2)*(-1. + pow(cth,2))*(22085775.*pow(r,9) 
    + 4571505.*pow(r,10) + 49392000.*(580. + 327.*pow(cth,2)) + 548800.*r*(-23041. + 74715.*pow(cth,2)) + 6300.*pow(r,3)*(-446807. + 973501.*pow(cth,2)) 
    + 126.*pow(r,7)*(1064483. + 1485790.*pow(cth,2)) + 9800.*pow(r,2)*(-1223527. + 1991748.*pow(cth,2)) + 12.*pow(r,8)*(4458631. + 3456783.*pow(cth,2)) 
    + 280.*pow(r,4)*(3773463. + 15733496.*pow(cth,2)) + 42.*pow(r,6)*(6767669. + 23527525.*pow(cth,2)) + 56.*pow(r,5)*(17855552. + 49207893.*pow(cth,2))))/(3087000.*pow(r,11))); 

  //EDGB metric
  /*((2. + r)*a*(-1. + pow(cth,2)))/r - (2.*a*pow(a,2)*pow(cth,2)*(-1. + pow(cth,2)))/pow(r,3) + zeta*(-1./36750.*((16660000. - 5350800.*r + 4797450.*pow(r,2) + 3526600.*pow(r,3) 
    + 2965560.*pow(r,4) + 918855.*pow(r,5) + 187446.*pow(r,6))*a*(-1. + pow(cth,2)))/pow(r,7) + (a*pow(a,2)*(-1. + pow(cth,2))*(22085775.*pow(r,9) 
    + 4571505.*pow(r,10) + 49392000.*(580. + 327.*pow(cth,2)) + 548800.*r*(-23041. + 74715.*pow(cth,2)) + 6300.*pow(r,3)*(-446807. + 973501.*pow(cth,2)) 
    + 126.*pow(r,7)*(1064483. + 1485790.*pow(cth,2)) + 9800.*pow(r,2)*(-1223527. + 1991748.*pow(cth,2)) + 12.*pow(r,8)*(4458631. + 3456783.*pow(cth,2)) 
    + 280.*pow(r,4)*(3773463. + 15733496.*pow(cth,2)) + 42.*pow(r,6)*(6767669. + 23527525.*pow(cth,2)) + 56.*pow(r,5)*(17855552. + 49207893.*pow(cth,2))))/(3087000.*pow(r,11)));
*/
  //Original value:
  //-a * s2 * (1. + 2. * r / rho2);

  gcov[2][2] = pow(r,2) + pow(a,2)*pow(cth,2) + ((8820000. - 6213200.*r - 3416700.*pow(r,2) - 1855650.*pow(r,3) + 887110.*pow(r,4) + 800733.*pow(r,5) + 435540.*pow(r,6) 
   + 187446.*pow(r,7))*pow(a,2)*zeta*(1. - 3.*pow(cth,2)))/(110250.*pow(r,8)) + (pow(a,4)*zeta*(45715050.*pow(r,11)*(-1. + 3.*pow(cth,2)) + 5625.*pow(r,10)*(-20749. + 58131.*pow(cth,2)) 
   + 493920000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) + 24696000.*r*(3049. - 10698.*pow(cth,2) + 8868.*pow(cth,4)) + 117600.*pow(r,2)*(280331. - 1711445.*pow(cth,2) + 1596165.*pow(cth,4)) 
   + 180.*pow(r,9)*(-1286466. - 846865.*pow(cth,2) + 5819941.*pow(cth,4)) + 4200.*pow(r,3)*(2362411. - 16650910.*pow(cth,2) + 14489100.*pow(cth,4)) 
   - 1260.*pow(r,4)*(-3173281. - 5026080.*pow(cth,2) + 26477920.*pow(cth,4)) + 42.*pow(r,8)*(-18071967. - 940590.*pow(cth,2) + 54146980.*pow(cth,4)) 
   + 42.*pow(r,6)*(-19116713. - 46592740.*pow(cth,2) + 138130070.*pow(cth,4)) - 28.*pow(r,5)*(11804979. - 261030540.*pow(cth,2) + 235282135.*pow(cth,4)) 
   + 6.*pow(r,7)*(-259078241. - 99440670.*pow(cth,2) + 857000595.*pow(cth,4))))/(92610000.*pow(r,12));

  //EDGB metric
  /*pow(r,2) + pow(a,2)*pow(cth,2) + ((8820000. - 6213200.*r - 3416700.*pow(r,2) - 1855650.*pow(r,3) + 887110.*pow(r,4) + 800733.*pow(r,5) + 435540.*pow(r,6) 
   + 187446.*pow(r,7))*pow(a,2)*zeta*(1. - 3.*pow(cth,2)))/(110250.*pow(r,8)) + (pow(a,4)*zeta*(45715050.*pow(r,11)*(-1. + 3.*pow(cth,2)) + 5625.*pow(r,10)*(-20749. + 58131.*pow(cth,2)) 
   + 493920000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) + 24696000.*r*(3049. - 10698.*pow(cth,2) + 8868.*pow(cth,4)) + 117600.*pow(r,2)*(280331. - 1711445.*pow(cth,2) + 1596165.*pow(cth,4)) 
   + 180.*pow(r,9)*(-1286466. - 846865.*pow(cth,2) + 5819941.*pow(cth,4)) + 4200.*pow(r,3)*(2362411. - 16650910.*pow(cth,2) + 14489100.*pow(cth,4)) 
   - 1260.*pow(r,4)*(-3173281. - 5026080.*pow(cth,2) + 26477920.*pow(cth,4)) + 42.*pow(r,8)*(-18071967. - 940590.*pow(cth,2) + 54146980.*pow(cth,4)) 
   + 42.*pow(r,6)*(-19116713. - 46592740.*pow(cth,2) + 138130070.*pow(cth,4)) - 28.*pow(r,5)*(11804979. - 261030540.*pow(cth,2) + 235282135.*pow(cth,4)) 
   + 6.*pow(r,7)*(-259078241. - 99440670.*pow(cth,2) + 857000595.*pow(cth,4))))/(92610000.*pow(r,12));*/
  //Original value:
  //rho2;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];

  gcov[3][3] = -(pow(r,2)*(-1. + pow(cth,2))) - (pow(a,2)*(2. + r - 2.*pow(cth,2))*(-1. + pow(cth,2)))/r - (2.*pow(a,4)*pow(cth,2)*pow((-1. + pow(cth,2)),2))/pow(r,3) 
    + zeta*(((8820000. - 6213200.*r - 3416700.*pow(r,2) - 1855650.*pow(r,3) + 887110.*pow(r,4) + 800733.*pow(r,5) + 435540.*pow(r,6) 
    + 187446.*pow(r,7))*pow(a,2)* (-1. + pow(cth,2))*(-1. + 3.*pow(cth,2)))/(110250.*pow(r,8)) -  (pow(a,4)*(-1. + pow(cth,2))*(45715050.*pow(r,11)*(-1. + 3.*pow(cth,2)) 
    + 5625.*pow(r,10)*(-20749. + 58131.*pow(cth,2)) + 493920000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) + 24696000.*r*(3649. - 8958.*pow(cth,2) + 6528.*pow(cth,4)) 
    + 352800.*pow(r,2)*(84857. - 350495.*pow(cth,2) + 320655.*pow(cth,4)) + 12600.*pow(r,3)*(-82303. - 1443030.*pow(cth,2) + 1592200.*pow(cth,4)) 
    + 180.*pow(r,9)*(-411718. - 4345857.*pow(cth,2) + 8444185.*pow(cth,4)) - 1260.*pow(r,4)*(1578719. - 11450880.*pow(cth,2) + 28150720.*pow(cth,4)) 
    + 42.*pow(r,8)*(-1863327. - 67980150.*pow(cth,2) + 104977900.*pow(cth,4)) + 28.*pow(r,5)*(-14247879. - 109560360.*pow(cth,2) + 137751665.*pow(cth,4)) 
    + 42.*pow(r,6)*(30654807. - 316973820.*pow(cth,2) + 358739630.*pow(cth,4)) + 6.*pow(r,7)*(-25024421. - 1143700950.*pow(cth,2) + 1667207055.*pow(cth,4))))/(92610000.*pow(r,12)));

  //EDGB metric
  /*-(pow(r,2)*(-1. + pow(cth,2))) - (pow(a,2)*(2. + r - 2.*pow(cth,2))*(-1. + pow(cth,2)))/r - (2.*pow(a,4)*pow(cth,2)*pow((-1. + pow(cth,2)),2))/pow(r,3) 
    + zeta*(((8820000. - 6213200.*r - 3416700.*pow(r,2) - 1855650.*pow(r,3) + 887110.*pow(r,4) + 800733.*pow(r,5) + 435540.*pow(r,6) 
    + 187446.*pow(r,7))*pow(a,2)* (-1. + pow(cth,2))*(-1. + 3.*pow(cth,2)))/(110250.*pow(r,8)) -  (pow(a,4)*(-1. + pow(cth,2))*(45715050.*pow(r,11)*(-1. + 3.*pow(cth,2)) 
    + 5625.*pow(r,10)*(-20749. + 58131.*pow(cth,2)) + 493920000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) + 24696000.*r*(3649. - 8958.*pow(cth,2) + 6528.*pow(cth,4)) 
    + 352800.*pow(r,2)*(84857. - 350495.*pow(cth,2) + 320655.*pow(cth,4)) + 12600.*pow(r,3)*(-82303. - 1443030.*pow(cth,2) + 1592200.*pow(cth,4)) 
    + 180.*pow(r,9)*(-411718. - 4345857.*pow(cth,2) + 8444185.*pow(cth,4)) - 1260.*pow(r,4)*(1578719. - 11450880.*pow(cth,2) + 28150720.*pow(cth,4)) 
    + 42.*pow(r,8)*(-1863327. - 67980150.*pow(cth,2) + 104977900.*pow(cth,4)) + 28.*pow(r,5)*(-14247879. - 109560360.*pow(cth,2) + 137751665.*pow(cth,4)) 
    + 42.*pow(r,6)*(30654807. - 316973820.*pow(cth,2) + 358739630.*fourcost) + 6.*pow(r,7)*(-25024421. - 1143700950.*pow(cth,2) + 1667207055.*pow(cth,4))))/(92610000.*pow(r,12)));*/
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
