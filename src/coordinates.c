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

  //Added variables for EdGB/DCS
  double zeta = 0.15;
  double c2 = cth * cth;
  double c4 = c2*c2;
  double a2 = a*a;
  double a4 = a2*a2;
 
  

  MUNULOOP gcov[mu][nu] = 0.;
  // Compute KS metric from KS coordinates (cyclic in t,phi)
  
  gcov[0][0] = -1. + 2.*r/rho2
    + zeta*(-1./15.*(-400.+96.*r+66.*pow(r,2)+130.*pow(r,3)+5.*pow(r,4))/pow(r,7)
    + (a2*(pow(r,7)*(444696.-562338.*c2)  + 8820000.*(-1. + 3.*c2) - 19600.*r*(-467. + 1251.*c2) - 63.*pow(r,8)*(-3267. + 8926.*c2) 
    + 1050.*pow(r,3)*(-1465. + 11997.*c2) - 2100.*pow(r,2)*(-955. + 22577.*c2) - 6.*pow(r,6)*(-59329. + 82437.*c2) + 15.*pow(r,5)*(-52533. + 455029.*c2) 
    + 10.*pow(r,4)*(-281221. + 1218513.*c2)))/(110250.*pow(r,11)) + (a4*(675.*pow(r,12)*(-19717. + 67726.*c2) + 450.*pow(r,11)*(43312. + 101589.*c2) 
    - 164640000.*(-70. + 585.*c2 + 156.*c4) - 8232000.*r*(3949. - 11058.*c2 + 6948.*c4) - 39200.*pow(r,2)*(189191. - 824825.*c2 + 972045.*c4) 
    + 60.*pow(r,9)*(717867. - 13885852.*c2 + 12733507.*c4) + 30.*pow(r,10)*(209773. - 9090216.*c2 + 16888370.*c4) + 1400.*pow(r,3)*(648009. - 14691730.*c2 
    + 26074500.*c4) + 420.*pow(r,4)*(553219. - 32471380.*c2 + 222891320.*c4) - 14.*pow(r,7)*(-11393603. + 38599350.*c2 + 359928985.*c4) 
    + 2.*pow(r,8)*(59964679. - 491173980.*c2 + 452824800.*c4) - 28.*pow(r,5)*(21852807. + 12094180.*c2 + 762315255.*c4) 
    - 14.*pow(r,6)*(42004007. - 226439060.*c2 + 1041312310.*c4)))/(30870000.*pow(r,15)));

  //EdGB metric
  /*-1. + 2.*r/rho2
    + zeta*(-1./15.*(-400.+96.*r+66.*pow(r,2)+130.*pow(r,3)+5.*pow(r,4))/pow(r,7)
    + (a2*(pow(r,7)*(444696.-562338.*c2)  + 8820000.*(-1. + 3.*c2) - 19600.*r*(-467. + 1251.*c2) - 63.*pow(r,8)*(-3267. + 8926.*c2) 
    + 1050.*pow(r,3)*(-1465. + 11997.*c2) - 2100.*pow(r,2)*(-955. + 22577.*c2) - 6.*pow(r,6)*(-59329. + 82437.*c2) + 15.*pow(r,5)*(-52533. + 455029.*c2) 
    + 10.*pow(r,4)*(-281221. + 1218513.*c2)))/(110250.*pow(r,11)) + (a4*(675.*pow(r,12)*(-19717. + 67726.*c2) + 450.*pow(r,11)*(43312. + 101589.*c2) 
    - 164640000.*(-70. + 585.*c2 + 156.*c4) - 8232000.*r*(3949. - 11058.*c2 + 6948.*c4) - 39200.*pow(r,2)*(189191. - 824825.*c2 + 972045.*c4) 
    + 60.*pow(r,9)*(717867. - 13885852.*c2 + 12733507.*c4) + 30.*pow(r,10)*(209773. - 9090216.*c2 + 16888370.*c4) + 1400.*pow(r,3)*(648009. - 14691730.*c2 
    + 26074500.*c4) + 420.*pow(r,4)*(553219. - 32471380.*c2 + 222891320.*c4) - 14.*pow(r,7)*(-11393603. + 38599350.*c2 + 359928985.*c4) 
    + 2.*pow(r,8)*(59964679. - 491173980.*c2 + 452824800.*c4) - 28.*pow(r,5)*(21852807. + 12094180.*c2 + 762315255.*c4) 
    - 14.*pow(r,6)*(42004007. - 226439060.*c2 + 1041312310.*c4)))/(30870000.*pow(r,15)));*/
  //DCS Metric
  /*-1. + 2.*r/rho2
    + zeta*((a2*(-338688.*(-1. + 3.*c2) + 4221.*pow(r,6.)*(-1. + 3.*c2) + 4221.*pow(r,7.)*(-1. + 3.*c2) 
    - 420.*pow(r,2.)*(263. + 321.*c2) + 20.*pow(r,3.)*(-5428. + 2025.*c2) - 20.*pow(r,4.)*(1523. + 2781.*c2) + 168.*r*(-275. + 3471.*c2) 
    + 2.*pow(r,5.)*(-2482. + 6711.*c2)))/(37632.*pow(r,10.)) + (a4*(-436560.*pow(r,10.)*(-1. + 3.*c2) - 436560.*pow(r,11.)*(-1. + 3.*c2) 
    + pow(r,5.)*(-27367138. + 338395440.*c2 - 79842960.*c4) + pow(r,7.)*(5994779. + 9180120.*c2 - 23096700.*c4) 
    + pow(r,8.)*(1835988. + 14588424.*c2 - 20269020.*c4) + pow(r,9.)*(609036. + 6695880.*c2 - 14028580.*c4) 
    - 3048192.*(299. - 1440.*c2 + 720.*c4) + 1512.*r*(200377. - 262200.*c2 + 546900.*c4) 
    + 540.*pow(r,2.)*(209695. - 1641704.*c2 + 3850652.*c4) + 36.*pow(r,4.)*(-1309617. + 14023212.*c2 + 5328620.*c4) 
    + pow(r,6.)*(1493875. + 66961128.*c2 + 66934880.*c4) - 12.*pow(r,3.)*(-7413845. - 9229392.*c2 + 156818740.*c4)))/(13547520.*pow(r,14.)));*/
  //Original value:
  //-1. + 2. * r / rho2;

  gcov[0][1] = 2.*r/rho2
    + zeta*(-1./30.*(-800. + 912.*r + 516.*pow(r,2) + 470.*pow(r,3) + 50.*pow(r,4) + 15.*pow(r,5))/ pow(r,7) 
    + (a2*(-566034.*pow(r,8) + 55125.*pow(r,9) +  17640000.*(-1. + 3.*c2) + 12600.*pow(r,3)*(4006. + 1877.*c2) + 78400.*r*(-779. + 2412.*c2) + 42.*pow(r,7)*(-94339. 
    + 107112.*c2) - 1400.*pow(r,2)*(20431. + 132243.*c2) + 36.*pow(r,6)*(-436917. + 491281.*c2) + 20.*pow(r,4)*(-875941. + 2263053.*c2) + 20.*pow(r,5)*(-1122937. 
    + 2632446.*c2)))/(220500.*pow(r,11)) - (a4*(-80247150.*pow(r,12) - 5788125.*pow(r,13) + 450.*pow(r,11)*(-1196839. + 812712.*c2) 
    + 329280000.*(-70. + 585.*c2 + 156.*c4) + 49392000.*r*(-1717. + 21664.*c2 + 9076.*c4) + 60.*pow(r,10)*(-23601289. + 13406112.*c2 + 2187250.*c4) 
    + 78400.*pow(r,2)*(-2206069. - 11318105.*c2 + 13657725.*c4) + 14000.*pow(r,3)*(22540153. - 88517480.*c2 + 62290230.*c4) 
    + 30.*pow(r,9)*(-145291221. + 30934768.*c2 + 146683252.*c4) + 280.*pow(r,4)*(793805393. - 2014699860.*c2 + 507428040.*c4) 
    + 28.*pow(r,6)*(-711534337. - 1340707620.*c2 + 4191169150.*c4) + 28.*pow(r,5)*(2484408549. - 2498856340.*c2 + 4680981810.*c4) 
    + 14.*pow(r,7)*(-1576790517. - 1991632680.*c2 + 4890821060.*c4) + 4.*pow(r,8)*(-2741883352. - 2122491690.*c2 + 5145464715.*c4)))/(61740000.*pow(r,15)));

  //EdGB metric
  /*2.*r/rho2
    + zeta*(-1./30.*(-800. + 912.*r + 516.*pow(r,2) + 470.*pow(r,3) + 50.*pow(r,4) + 15.*pow(r,5))/ pow(r,7) 
    + (a2*(-566034.*pow(r,8) + 55125.*pow(r,9) +  17640000.*(-1. + 3.*c2) + 12600.*pow(r,3)*(4006. + 1877.*c2) + 78400.*r*(-779. + 2412.*c2) + 42.*pow(r,7)*(-94339. 
    + 107112.*c2) - 1400.*pow(r,2)*(20431. + 132243.*c2) + 36.*pow(r,6)*(-436917. + 491281.*c2) + 20.*pow(r,4)*(-875941. + 2263053.*c2) + 20.*pow(r,5)*(-1122937. 
    + 2632446.*c2)))/(220500.*pow(r,11)) - (a4*(-80247150.*pow(r,12) - 5788125.*pow(r,13) + 450.*pow(r,11)*(-1196839. + 812712.*c2) 
    + 329280000.*(-70. + 585.*c2 + 156.*c4) + 49392000.*r*(-1717. + 21664.*c2 + 9076.*c4) + 60.*pow(r,10)*(-23601289. + 13406112.*c2 + 2187250.*c4) 
    + 78400.*pow(r,2)*(-2206069. - 11318105.*c2 + 13657725.*c4) + 14000.*pow(r,3)*(22540153. - 88517480.*c2 + 62290230.*c4) 
    + 30.*pow(r,9)*(-145291221. + 30934768.*c2 + 146683252.*c4) + 280.*pow(r,4)*(793805393. - 2014699860.*c2 + 507428040.*c4) 
    + 28.*pow(r,6)*(-711534337. - 1340707620.*c2 + 4191169150.*c4) + 28.*pow(r,5)*(2484408549. - 2498856340.*c2 + 4680981810.*c4) 
    + 14.*pow(r,7)*(-1576790517. - 1991632680.*c2 + 4890821060.*c4) + 4.*pow(r,8)*(-2741883352. - 2122491690.*c2 + 5145464715.*c4)))/(61740000.*pow(r,15)));*/
  //DCS metric
  /*2.*r/rho2
    + zeta*((a2*(8442.*pow(r,7.) + pow(r,4.)*(504690. - 543240.*c2) + pow(r,6.)*(55419. - 50652.*c2) - 338688.*(-1. + 3.*c2) 
    - 2016.*pow(r,2.)*(-747. + 1291.*c2) - 168.*r*(-7789. + 20721.*c2) - 20.*pow(r,3.)*(-50957. + 73764.*c2) 
    - 2.*pow(r,5.)*(-94018. + 97629.*c2)))/(37632.*pow(r,10.)) + (a4*(-1746240.*pow(r,11.) + 240.*pow(r,10.)*(-48727. + 43656.*c2) 
    - 6096384.*(299. - 1440.*c2 + 720.*c4) + 8.*pow(r,9.)*(-4902855. + 3171874.*c2 + 332305.*c4) 
    - 3456.*pow(r,2.)*(5955596. - 12196010.*c2 + 3150155.*c4) - 3024.*r*(3416327. - 17156040.*c2 + 8162220.*c4) 
    + 96.*pow(r,3.)*(-80010145. + 69893298.*c2 + 34916840.*c4) + 8.*pow(r,7.)*(-55039417. - 5982330.*c2 + 69327885.*c4) 
    + pow(r,8.)*(-164813529. + 45059288.*c2 + 116931460.*c4) + 36.*pow(r,4.)*(-27676299. - 98275560.*c2 + 239844820.*c4)
     + 8.*pow(r,6.)*(-76467017. - 93284976.*c2 + 253740770.*c4) + 4.*pow(r,5.)*(-83689817. - 868048872.*c2 + 1345549080.*c4)))/(27095040.*pow(r,14.)));*/
  //Original value:
  //2. * r / rho2 +0.0001;

  gcov[0][3] = -2.*a*r*s2/rho2
    + zeta*(-1./15.*((-400. + 144.*r + 90.*pow(r,2) + 140.*pow(r,3) 
    + 9.*pow(r,4))*a*(-1. + c2))/pow(r,7) - (a2*a*(-1. + c2)* (pow(r,4)*(2736210. - 4410530.*c2) + pow(r,5)*(766015. - 3620183.*c2) - 
    8820000.*(-1. + 3.*c2) + 19600.*r*(-467. + 1551.*c2) - 12.*pow(r,6)*(26511. +  6310.*c2) + 750.*pow(r,3)*(2051. + 8733.*c2) +  2100.*pow(r,2)*(-955. + 21233.*c2) 
    + 3.*pow(r,8)*(-63529. + 262520.*c2) + pow(r,7)*(-406611. + 563055.*c2)))/(110250.*pow(r,11)));

  //EdGB metric
  /*-2.*a*r*s2/rho2
    + zeta*(-1./15.*((-400. + 144.*r + 90.*pow(r,2) + 140.*pow(r,3) 
    + 9.*pow(r,4))*a*(-1. + c2))/pow(r,7) - (a2*a*(-1. + c2)* (pow(r,4)*(2736210. - 4410530.*c2) + pow(r,5)*(766015. - 3620183.*c2) - 
    8820000.*(-1. + 3.*c2) + 19600.*r*(-467. + 1551.*c2) - 12.*pow(r,6)*(26511. +  6310.*c2) + 750.*pow(r,3)*(2051. + 8733.*c2) +  2100.*pow(r,2)*(-955. + 21233.*c2) 
    + 3.*pow(r,8)*(-63529. + 262520.*c2) + pow(r,7)*(-406611. + 563055.*c2)))/(110250.*pow(r,11)));*/
  //DCS metric
  /*-2.*a*r*sqsint/rho2
     +zeta*(-1./112.*((189. + 120.*r + 70.*pow(r,2.))*a*(-1. + c2))/pow(r,6.) + (pow(a,3)*(-1. + c2)*(-10160640.*(-1. + 3.*c2) 
     + 105828.*pow(r,7.)*(-1. + 5.*c2) + 30.*pow(r,5.)*(2247. + 95.*c2) + 55440.*r*(-25. + 453.*c2) + 15.*pow(r,6.)*(-271. + 24875.*c2) 
     + 360.*pow(r,2.)*(-9205. + 64257.*c2) + 16.*pow(r,4.)*(-40883. + 112735.*c2) + 8.*pow(r,3.)*(-382792. + 1712045.*c2)))/(1128960.*pow(r,10.)));*/
  //Original value:
  //-2. * a * r * s2 / rho2;

  gcov[1][0] = gcov[0][1];

  gcov[1][1] = 1. + 2.*r/rho2 
    + zeta*(-1./15.*(-400. + 816.*r + 450.*pow(r,2) + 340.*pow(r,3) + 45.*pow(r,4) + 15.*pow(r,5))/pow(r,7) 
    + (a2*(55125.*pow(r,9) + 8820000.*(-1. + 3.*c2) + 1050.*pow(r,3)*(49537. + 10527.*c2) + 19600.*r*(-3583. + 10899.*c2) + 21.*pow(r,8)*(-36755. + 26778.*c2) 
    + 42.*pow(r,7)*(-104927. + 120501.*c2) - 700.*pow(r,2)*(43727. + 196755.*c2) + 6.*pow(r,6)*(-2680831. + 3030123.*c2) + 10.*pow(r,4)*(-1470661. + 3307593.*c2) 
    + 5.*pow(r,5)*(-4334149. + 9164697.*c2)))/(110250.*pow(r,11)) - (a4*(-5788125.*pow(r,13) + 225.*pow(r,12)*(-415805. + 203178.*c2) 
    + 1350.*pow(r,11)*(-384509. + 304767.*c2) + 164640000.*(-70. + 585.*c2 + 156.*c4) + 8232000.*r*(-14251. + 141042.*c2 + 47508.*c4) 
    + 30.*pow(r,10)*(-46992805. + 17722008.*c2 + 21262870.*c4) + 39200.*pow(r,2)*(-4601329. - 21811385.*c2 + 26343405.*c4) 
    + 30.*pow(r,9)*(-143855487. + 3163064.*c2 + 172150266.*c4) + 1400.*pow(r,3)*(226049539. - 899866530.*c2 + 648976800.*c4) 
    + 140.*pow(r,4)*(1589270443. - 4126813860.*c2 + 1683530040.*c4) + 10.*pow(r,8)*(-1084760405. - 947231472.*c2 + 2148750846.*c4) 
    + 28.*pow(r,5)*(2462555742. - 2510950520.*c2 + 3918666555.*c4) + 14.*pow(r,7)*(-1565396914. - 2030232030.*c2 + 4530892075.*c4) 
    + 14.*pow(r,6)*(-1465072681. - 2454976180.*c2 + 7341025990.*c4)))/(30870000.*pow(r,15)));

  //EdGB metric
  /*1. + 2.*r/rho2 
    + zeta*(-1./15.*(-400. + 816.*r + 450.*pow(r,2) + 340.*pow(r,3) + 45.*pow(r,4) + 15.*pow(r,5))/pow(r,7) 
    + (a2*(55125.*pow(r,9) + 8820000.*(-1. + 3.*c2) + 1050.*pow(r,3)*(49537. + 10527.*c2) + 19600.*r*(-3583. + 10899.*c2) + 21.*pow(r,8)*(-36755. + 26778.*c2) 
    + 42.*pow(r,7)*(-104927. + 120501.*c2) - 700.*pow(r,2)*(43727. + 196755.*c2) + 6.*pow(r,6)*(-2680831. + 3030123.*c2) + 10.*pow(r,4)*(-1470661. + 3307593.*c2) 
    + 5.*pow(r,5)*(-4334149. + 9164697.*c2)))/(110250.*pow(r,11)) - (a4*(-5788125.*pow(r,13) + 225.*pow(r,12)*(-415805. + 203178.*c2) 
    + 1350.*pow(r,11)*(-384509. + 304767.*c2) + 164640000.*(-70. + 585.*c2 + 156.*c4) + 8232000.*r*(-14251. + 141042.*c2 + 47508.*c4) 
    + 30.*pow(r,10)*(-46992805. + 17722008.*c2 + 21262870.*c4) + 39200.*pow(r,2)*(-4601329. - 21811385.*c2 + 26343405.*c4) 
    + 30.*pow(r,9)*(-143855487. + 3163064.*c2 + 172150266.*c4) + 1400.*pow(r,3)*(226049539. - 899866530.*c2 + 648976800.*c4) 
    + 140.*pow(r,4)*(1589270443. - 4126813860.*c2 + 1683530040.*c4) + 10.*pow(r,8)*(-1084760405. - 947231472.*c2 + 2148750846.*c4) 
    + 28.*pow(r,5)*(2462555742. - 2510950520.*c2 + 3918666555.*c4) + 14.*pow(r,7)*(-1565396914. - 2030232030.*c2 + 4530892075.*c4) 
    + 14.*pow(r,6)*(-1465072681. - 2454976180.*c2 + 7341025990.*c4)))/(30870000.*pow(r,15)));*/
  //DCS metric
  /*1. + 2.*r/rho2  
    + zeta*(-1./37632.*(a2*(4221.*pow(r,7.)*(-5. + 3.*c2) + 338688.*(-1. + 3.*c2) + 3820.*pow(r,3.)*(-562. + 783.*c2) 
    + 21.*pow(r,6.)*(-5479. + 5427.*c2) + 168.*r*(-15853. + 44913.*c2) + 20.*pow(r,4.)*(-51992. + 51543.*c2) 
    + 84.*pow(r,2.)*(-37171. + 60363.*c2) + 6.*pow(r,5.)*(-63506. + 67323.*c2)))/pow(r,10.) + (a4*(436560.*pow(r,11.)*(-5. + 3.*c2) 
    + 240.*pow(r,10.)*(-50546. + 49113.*c2) - 3048192.*(299. - 1440.*c2 + 720.*c4) 
    + 4.*pow(r,9.)*(-9957969. + 4669778.*c2 + 4171755.*c4) - 1512.*r*(7033031. - 34574280.*c2 + 16871340.*c4) 
    + 5.*pow(r,7.)*(-89262023. - 11407752.*c2 + 115543956.*c4) + 72.*pow(r,4.)*(-13183341. - 56149386.*c2 + 117258100.*c4) 
    - 108.*pow(r,2.)*(191627547. - 398480840.*c2 + 120058220.*c4) + pow(r,8.)*(-166649517. + 30470864.*c2 + 137200480.*c4) 
    + 12.*pow(r,3.)*(-647495005. + 549916992.*c2 + 436153460.*c4) + pow(r,6.)*(-613230011. - 813240936.*c2 + 1962991280.*c4) 
    + pow(r,5.)*(-307392130. - 3810590928.*c2 + 5462039280.*c4)))/(13547520.*pow(r,14.)));*/
  //Original value:
  //1. + 2. * r / rho2;

  gcov[1][3] = -a *s2*(1. + 2.*r/rho2)
    + zeta*(-1./36750.*((16660000. - 5350800.*r + 4797450.*pow(r,2) + 3526600.*pow(r,3) 
    + 2965560.*pow(r,4) + 918855.*pow(r,5) + 187446.*pow(r,6))*a*(-1. + c2))/pow(r,7) + (a*a2*(-1. + c2)*(22085775.*pow(r,9) 
    + 4571505.*pow(r,10) + 49392000.*(580. + 327.*c2) + 548800.*r*(-23041. + 74715.*c2) + 6300.*pow(r,3)*(-446807. + 973501.*c2) 
    + 126.*pow(r,7)*(1064483. + 1485790.*c2) + 9800.*pow(r,2)*(-1223527. + 1991748.*c2) + 12.*pow(r,8)*(4458631. + 3456783.*c2) 
    + 280.*pow(r,4)*(3773463. + 15733496.*c2) + 42.*pow(r,6)*(6767669. + 23527525.*c2) + 56.*pow(r,5)*(17855552. + 49207893.*c2)))/(3087000.*pow(r,11))); 

  //EdGB metric
  /*-a *s2*(1. + 2.*r/rho2)
    + zeta*(-1./36750.*((16660000. - 5350800.*r + 4797450.*pow(r,2) + 3526600.*pow(r,3) 
    + 2965560.*pow(r,4) + 918855.*pow(r,5) + 187446.*pow(r,6))*a*(-1. + c2))/pow(r,7) + (a*a2*(-1. + c2)*(22085775.*pow(r,9) 
    + 4571505.*pow(r,10) + 49392000.*(580. + 327.*c2) + 548800.*r*(-23041. + 74715.*c2) + 6300.*pow(r,3)*(-446807. + 973501.*c2) 
    + 126.*pow(r,7)*(1064483. + 1485790.*c2) + 9800.*pow(r,2)*(-1223527. + 1991748.*c2) + 12.*pow(r,8)*(4458631. + 3456783.*c2) 
    + 280.*pow(r,4)*(3773463. + 15733496.*c2) + 42.*pow(r,6)*(6767669. + 23527525.*c2) + 56.*pow(r,5)*(17855552. + 49207893.*c2)))/(3087000.*pow(r,11)));*/
//DCS metric
/*-a *sqsint*(1. + 2.*r/rho2) 
  +zeta*(((571536. + 393540.*r + 219380.*pow(r,2.) + 64660.*pow(r,3.) + 19880.*pow(r,4.) + 4221.*pow(r,5.))*a*(-1. + c2))/(12544.*pow(r,6.)) 
  - (pow(a,3)*(-1. + c2)*(1570950.*pow(r,8.) + 327420.*pow(r,9.) - 60963840.*(-49. + 48.*c2) - 7560.*r*(-133087. + 6297.*c2) 
  + 1620.*pow(r,2.)*(47721. + 449461.*c2) + 25.*pow(r,6.)*(711795. + 481847.*c2) + 30.*pow(r,5.)*(1256960. + 2241169.*c2) 
  + pow(r,7.)*(4312974. + 2699780.*c2) + 120.*pow(r,3.)*(-870409. + 5949739.*c2) + 12.*pow(r,4.)*(379503. + 21372890.*c2)))/(3386880.*pow(r,10.)));*/
  //Original value:
  //-a * s2 * (1. + 2. * r / rho2);

  gcov[2][2] = rho2
   + ((8820000. - 6213200.*r - 3416700.*pow(r,2) - 1855650.*pow(r,3) + 887110.*pow(r,4) + 800733.*pow(r,5) + 435540.*pow(r,6) 
   + 187446.*pow(r,7))*pow(a,2)*zeta*(1. - 3.*pow(cth,2)))/(110250.*pow(r,8)) + (pow(a,4)*zeta*(45715050.*pow(r,11)*(-1. + 3.*pow(cth,2)) + 5625.*pow(r,10)*(-20749. + 58131.*pow(cth,2)) 
   + 493920000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) + 24696000.*r*(3049. - 10698.*pow(cth,2) + 8868.*pow(cth,4)) + 117600.*pow(r,2)*(280331. - 1711445.*pow(cth,2) + 1596165.*pow(cth,4)) 
   + 180.*pow(r,9)*(-1286466. - 846865.*pow(cth,2) + 5819941.*pow(cth,4)) + 4200.*pow(r,3)*(2362411. - 16650910.*pow(cth,2) + 14489100.*pow(cth,4)) 
   - 1260.*pow(r,4)*(-3173281. - 5026080.*pow(cth,2) + 26477920.*pow(cth,4)) + 42.*pow(r,8)*(-18071967. - 940590.*pow(cth,2) + 54146980.*pow(cth,4)) 
   + 42.*pow(r,6)*(-19116713. - 46592740.*pow(cth,2) + 138130070.*pow(cth,4)) - 28.*pow(r,5)*(11804979. - 261030540.*pow(cth,2) + 235282135.*pow(cth,4)) 
   + 6.*pow(r,7)*(-259078241. - 99440670.*pow(cth,2) + 857000595.*pow(cth,4))))/(92610000.*pow(r,12));

  //EdGB metric
  /*pow(r,2) + pow(a,2)*pow(cth,2) + ((8820000. - 6213200.*r - 3416700.*pow(r,2) - 1855650.*pow(r,3) + 887110.*pow(r,4) + 800733.*pow(r,5) + 435540.*pow(r,6) 
   + 187446.*pow(r,7))*pow(a,2)*zeta*(1. - 3.*pow(cth,2)))/(110250.*pow(r,8)) + (pow(a,4)*zeta*(45715050.*pow(r,11)*(-1. + 3.*pow(cth,2)) + 5625.*pow(r,10)*(-20749. + 58131.*pow(cth,2)) 
   + 493920000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) + 24696000.*r*(3049. - 10698.*pow(cth,2) + 8868.*pow(cth,4)) + 117600.*pow(r,2)*(280331. - 1711445.*pow(cth,2) + 1596165.*pow(cth,4)) 
   + 180.*pow(r,9)*(-1286466. - 846865.*pow(cth,2) + 5819941.*pow(cth,4)) + 4200.*pow(r,3)*(2362411. - 16650910.*pow(cth,2) + 14489100.*pow(cth,4)) 
   - 1260.*pow(r,4)*(-3173281. - 5026080.*pow(cth,2) + 26477920.*pow(cth,4)) + 42.*pow(r,8)*(-18071967. - 940590.*pow(cth,2) + 54146980.*pow(cth,4)) 
   + 42.*pow(r,6)*(-19116713. - 46592740.*pow(cth,2) + 138130070.*pow(cth,4)) - 28.*pow(r,5)*(11804979. - 261030540.*pow(cth,2) + 235282135.*pow(cth,4)) 
   + 6.*pow(r,7)*(-259078241. - 99440670.*pow(cth,2) + 857000595.*pow(cth,4))))/(92610000.*pow(r,12));*/
  //DCS metric
  /*rho2  
     -(zeta*(1080.*pow(r,4.)*(338688. + 80808.*r + 67380.*pow(r,2.) + 10360.*pow(r,3.) + 18908.*pow(r,4.) + 9940.*pow(r,5.) 
     + 4221.*pow(r,6.))*a2*(1. - 3.*c2) + a4*(3141900.*pow(r,9.)*(-1. + 3.*c2) + 1309680.*pow(r,10.)*(-1. + 3.*c2) 
     - 9144576.*(299. - 1440.*c2 + 720.*c4) - 4536.*r*(41543. - 947400.*c2 + 420780.*c4) 
     + 1620.*pow(r,2.)*(-93041. - 85688.*c2 + 1454300.*c4) + 12.*pow(r,8.)*(-535899. - 116700.*c2 + 2367475.*c4) 
     + 18.*pow(r,5.)*(-7227915. + 2913384.*c2 + 17254600.*c4) + 36.*pow(r,3.)*(4260605. - 30961368.*c2 + 27684860.*c4) 
     + 3.*pow(r,6.)*(-20890109. + 18841800.*c2 + 46796480.*c4) + pow(r,7.)*(-20549463. + 6081360.*c2 + 59683820.*c4) 
     + 4.*pow(r,4.)*(-33335811. - 120639012.*c2 + 227696380.*c4))))/(40642560.*pow(r,11.));*/
  //Original value:
  //rho2;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];

  gcov[3][3] =s2 * (rho2+a2*s2*(1. + 2. *r/rho2))
    + zeta*(((8820000. - 6213200.*r - 3416700.*pow(r,2) - 1855650.*pow(r,3) + 887110.*pow(r,4) + 800733.*pow(r,5) + 435540.*pow(r,6) 
    + 187446.*pow(r,7))*a2* (-1. + c2)*(-1. + 3.*c2))/(110250.*pow(r,8)) -  (a4*(-1. + c2)*(45715050.*pow(r,11)*(-1. + 3.*c2) 
    + 5625.*pow(r,10)*(-20749. + 58131.*c2) + 493920000.*(-70. + 585.*c2 + 156.*c4) + 24696000.*r*(3649. - 8958.*c2 + 6528.*c4) 
    + 352800.*pow(r,2)*(84857. - 350495.*c2 + 320655.*c4) + 12600.*pow(r,3)*(-82303. - 1443030.*c2 + 1592200.*c4) 
    + 180.*pow(r,9)*(-411718. - 4345857.*c2 + 8444185.*c4) - 1260.*pow(r,4)*(1578719. - 11450880.*c2 + 28150720.*c4) 
    + 42.*pow(r,8)*(-1863327. - 67980150.*c2 + 104977900.*c4) + 28.*pow(r,5)*(-14247879. - 109560360.*c2 + 137751665.*c4) 
    + 42.*pow(r,6)*(30654807. - 316973820.*c2 + 358739630.*c4) + 6.*pow(r,7)*(-25024421. - 1143700950.*c2 + 1667207055.*c4)))/(92610000.*pow(r,12)));

  //EdGB metric
  /*-(pow(r,2)*(-1. + pow(cth,2))) - (pow(a,2)*(2. + r - 2.*pow(cth,2))*(-1. + pow(cth,2)))/r - (2.*pow(a,4)*pow(cth,2)*pow((-1. + pow(cth,2)),2))/pow(r,3) 
    + zeta*(((8820000. - 6213200.*r - 3416700.*pow(r,2) - 1855650.*pow(r,3) + 887110.*pow(r,4) + 800733.*pow(r,5) + 435540.*pow(r,6) 
    + 187446.*pow(r,7))*pow(a,2)* (-1. + pow(cth,2))*(-1. + 3.*pow(cth,2)))/(110250.*pow(r,8)) -  (pow(a,4)*(-1. + pow(cth,2))*(45715050.*pow(r,11)*(-1. + 3.*pow(cth,2)) 
    + 5625.*pow(r,10)*(-20749. + 58131.*pow(cth,2)) + 493920000.*(-70. + 585.*pow(cth,2) + 156.*pow(cth,4)) + 24696000.*r*(3649. - 8958.*pow(cth,2) + 6528.*pow(cth,4)) 
    + 352800.*pow(r,2)*(84857. - 350495.*pow(cth,2) + 320655.*pow(cth,4)) + 12600.*pow(r,3)*(-82303. - 1443030.*pow(cth,2) + 1592200.*pow(cth,4)) 
    + 180.*pow(r,9)*(-411718. - 4345857.*pow(cth,2) + 8444185.*pow(cth,4)) - 1260.*pow(r,4)*(1578719. - 11450880.*pow(cth,2) + 28150720.*pow(cth,4)) 
    + 42.*pow(r,8)*(-1863327. - 67980150.*pow(cth,2) + 104977900.*pow(cth,4)) + 28.*pow(r,5)*(-14247879. - 109560360.*pow(cth,2) + 137751665.*pow(cth,4)) 
    + 42.*pow(r,6)*(30654807. - 316973820.*pow(cth,2) + 358739630.*pow(cth,4)) + 6.*pow(r,7)*(-25024421. - 1143700950.*pow(cth,2) + 1667207055.*pow(cth,4))))/(92610000.*pow(r,12)));*/
  //DCS metric
  /*sqsint * (rho2+a2*sqsint*(1. + 2. *r/rho2)) 
     +zeta*(-1./37632.*((338688. + 80808.*r + 67380.*pow(r,2.) + 10360.*pow(r,3.) + 18908.*pow(r,4.) + 9940.*pow(r,5.) 
     + 4221.*pow(r,6.))*a2*(-1. + c2)*(-1. + 3.*c2))/pow(r,7.) + (a4*(-1. + c2)*(3141900.*pow(r,9.)*(-1. + 3.*c2) 
     + 1309680.*pow(r,10.)*(-1. + 3.*c2) - 9144576.*(299. - 1440.*c2 + 720.*c4) 
     + 4536.*r*(119737. + 302280.*c2 + 63060.*c4) + 1620.*pow(r,2.)*(240495. - 1419832.*c2 + 2454908.*c4) 
     + 12.*pow(r,8.)*(-156009. - 1636260.*c2 + 3507145.*c4) + 18.*pow(r,5.)*(-4337355. - 8648856.*c2 + 25926280.*c4) 
     + 36.*pow(r,3.)*(10727645. - 56829528.*c2 + 47085980.*c4) + 3.*pow(r,6.)*(-6926429. - 37012920.*c2 + 88687520.*c4) 
     + pow(r,7.)*(-696903. - 73328880.*c2 + 119241500.*c4) + 4.*pow(r,4.)*(-9548811. - 215787012.*c2 + 299057380.*c4)))/(40642560.*pow(r,11.)));*/
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
