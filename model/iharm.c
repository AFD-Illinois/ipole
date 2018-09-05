#include "decs.h"
#include "hdf5_utils.h"

#define NVAR (10)
#define SLOW_LIGHT (0)

void interp_fourv(double X[NDIM], double ****fourv, double Fourv[NDIM]) ;
double interp_scalar(double X[NDIM], double ***var) ;
static double poly_norm, poly_xt, poly_alpha, mks_smooth;
static double game, gamp;

static double MBH, Mdotedd, tp_over_te, Thetae_unit;

static int RADIATION, ELECTRONS, DEREFINE_POLES;
#define TPTE (3.)

#define NSUP (3)
struct of_data {
  double t;
  double ****bcon;
  double ****bcov;
  double ****ucon;
  double ****ucov;
  double ****p;
  double ***ne;
  double ***thetae;
  double ***b;
};
static int nloaded = 0;

struct of_data dataA, dataB, dataC;
struct of_data *data[NSUP];
  
void load_iharm_data(int n, char *);

void parse_input(int argc, char *argv[], Params *params)
{

  // if params has been loaded, just read from it
  if ( params->loaded ) {
    thetacam = params->thetacam;
    freqcgs = params->freqcgs;
    MBH = params->MBH * MSUN;
    M_unit = params->M_unit;
    strcpy(fnam, params->dump);
    tp_over_te = params->tp_over_te;
    return;
  }


  if (argc != 7) {
    fprintf(stderr, "ERROR format is\n");
    fprintf(stderr, "  ipole theta[deg] freq[cgs] MBH[Msolar] M_unit[g] filename counterjet\n");
    exit(-1);
  }

  sscanf(argv[1], "%lf", &thetacam);
  sscanf(argv[2], "%lf", &freqcgs);
  sscanf(argv[3], "%lf", &MBH);
  MBH *= MSUN;
  sscanf(argv[4], "%lf", &M_unit);
  strcpy(fnam, argv[5]);
  sscanf(argv[6], "%d",  &counterjet);
  tp_over_te = TPTE;
}

void set_tinterp_ns(double X[NDIM], int *nA, int *nB)
{
  if (X[0] < data[1]->t) {
    *nA = 0; *nB = 1;
  } else {
    *nA = 1; *nB = 2;
  }
}

void update_data()
{
  #pragma omp single
  {
    #if SLOW_LIGHT
    // Get new filename
    int len = strlen(fnam);
    char buf[STRLEN];
    memmove(buf, fnam+len-11, 8);
    buf[8] = '\0';
    int fnum = atoi(buf);
    fnum += 1;
    char newfnam[STRLEN]; 
    memmove(newfnam, fnam, len-11);
    newfnam[len-11] = '\0';
    sprintf(buf, "%08d", fnum);
    strcat(newfnam+len-11, buf);
    sprintf(buf, ".h5");
    strcat(newfnam+len-8, buf);
    strcpy(fnam, newfnam);
    //exit(-1);

    // Reorder dataA, dataB, dataC in data[]
    if (nloaded % 3 == 0) {
      data[0] = &dataA;
      data[1] = &dataB;
      data[2] = &dataC;
    } else if (nloaded % 3 == 1) {
      data[0] = &dataB;
      data[1] = &dataC;
      data[2] = &dataA;
    } else if (nloaded % 3 == 2) {
      data[0] = &dataC;
      data[1] = &dataA;
      data[2] = &dataB;
    } else {
      printf("Fail! nloaded = %i nloaded mod 3 = %i\n", nloaded, nloaded % 3);
    }
    
    load_iharm_data(2, fnam);
    #else // FAST LIGHT
    if (nloaded % 3 == 0) {
      data[0] = &dataA;
      data[1] = &dataB;
      data[2] = &dataC;
    } else if (nloaded % 3 == 1) {
      data[0] = &dataB;
      data[1] = &dataC;
      data[2] = &dataA;
    } else if (nloaded % 3 == 2) {
      data[0] = &dataC;
      data[1] = &dataA;
      data[2] = &dataB;
    } else {
      printf("Fail! nloaded = %i nloaded mod 3 = %i\n", nloaded, nloaded % 3);
    }
    data[2]->t = data[1]->t + DTd;
    #endif 
    //data[2]->t = data[1]->t + DTd;
  } // omp single
}

void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{
  MUNULOOP dxdX[mu][nu] = 0.;

  if (DEREFINE_POLES) {

    // mmks
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][1] = -exp(mks_smooth*(startx[1]-X[1]))*mks_smooth*(
      M_PI/2. -
      M_PI*X[2] +
      poly_norm*(2.*X[2]-1.)*(1+(pow((-1.+2*X[2])/poly_xt,poly_alpha))/(1 + poly_alpha)) -
      1./2.*(1. - hslope)*sin(2.*M_PI*X[2])
      );
    dxdX[2][2] = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) +
      exp(mks_smooth*(startx[1]-X[1]))*(
        -M_PI +
        2.*poly_norm*(1. + pow((2.*X[2]-1.)/poly_xt,poly_alpha)/(poly_alpha+1.)) +
        (2.*poly_alpha*poly_norm*(2.*X[2]-1.)*pow((2.*X[2]-1.)/poly_xt,poly_alpha-1.))/((1.+poly_alpha)*poly_xt) -
        (1.-hslope)*M_PI*cos(2.*M_PI*X[2])
        );
    dxdX[3][3] = 1.;

  } else {

    // mks
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][2] = M_PI - (hslope - 1.)*M_PI*cos(2.*M_PI*X[2]);
    dxdX[3][3] = 1.;

  }

}

void gcov_func(double X[NDIM], double gcov[NDIM][NDIM])
{
  // returns g_{munu} at location specified by X
 
  MUNULOOP gcov[mu][nu] = 0.;
    
  double sth, cth, s2, rho2;
  double r, th;

  // despite the name, get equivalent values for
  // r, th for kerr coordinate system
  bl_coord(X, &r, &th);

  cth = cos(th);
  sth = sin(th);

  s2 = sth*sth;
  rho2 = r*r + a*a*cth*cth;

  // compute ks metric for ks coordinates (cyclic in t,phi)
  gcov[0][0] = -1. + 2.*r/rho2;
  gcov[0][1] = 2.*r/rho2;
  gcov[0][3] = -2.*a*r*s2/rho2;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = 1. + 2.*r/rho2;
  gcov[1][3] = -a*s2*(1. + 2.*r/rho2);
  
  gcov[2][2] = rho2;
  
  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2));

  // convert from ks metric to mks/mmks
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);

  double gcov_ks[NDIM][NDIM];
  MUNULOOP {
    gcov_ks[mu][nu] = gcov[mu][nu];
    gcov[mu][nu] = 0.;
  }

  MUNULOOP {
    for (int lam=0; lam<NDIM; ++lam) {
      for (int kap=0; kap<NDIM; ++kap) {
        gcov[mu][nu] += gcov_ks[lam][kap]*dxdX[lam][mu]*dxdX[kap][nu];
      }
    }
  }

}

void get_connection(double X[4], double lconn[4][4][4])
{
  get_connection_num(X, lconn);
}

void init_model(char *args[])
{
  void init_iharm_grid(char *);

  // set up initial ordering of data[]
  data[0] = &dataA;
  data[1] = &dataB;
  data[2] = &dataC;

  // set up grid for fluid data
  fprintf(stderr, "reading data header...\n");
  init_iharm_grid(fnam);
  fprintf(stderr, "success\n");

  // set all dimensional quantities from loaded parameters
  set_units(args[4]);

  // read fluid data
  fprintf(stderr, "reading data...\n");
  load_iharm_data(0, fnam);
  update_data();
  load_iharm_data(1, fnam);
  update_data();
  load_iharm_data(2, fnam);
  data[2]->t =10000.;
  fprintf(stderr, "success\n");

  // horizon radius
  Rh = 1 + sqrt(1. - a * a) ;
}

/*
 
  these supply basic model data to ipole

*/
void get_model_ucov(double X[NDIM], double Ucov[NDIM])
{
  double gcov[NDIM][NDIM];

  gcov_func(X, gcov);

  if(X[1] < startx[1] || 
     X[1] > stopx[1]  || 
     X[2] < startx[2] || 
     X[2] > stopx[2]) {
     
    // sensible default value
    Ucov[0] = -1./sqrt(-gcov[0][0]) ;
    Ucov[1] = 0. ;
    Ucov[2] = 0. ;
    Ucov[3] = 0. ;

    return ;
  }

  double Ucon[NDIM];
  get_model_ucon(X, Ucon);
  lower(Ucon, gcov, Ucov);

}

void get_model_ucon(double X[NDIM], double Ucon[NDIM])
{
  double gcov[NDIM][NDIM] ;
  double gcon[NDIM][NDIM] ;
  double tmp[NDIM] ;

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  if(X[1] < startx[1] || 
     X[1] > stopx[1]  || 
     X[2] < startx[2] || 
     X[2] > stopx[2]) {
      /* sensible default value */
      gcov_func(X, gcov) ;

    tmp[0] = -1./sqrt(-gcov[0][0]) ;
    tmp[1] = 0. ;
    tmp[2] = 0. ;
    tmp[3] = 0. ;

      gcon_func(gcov, gcon) ;
    Ucon[0] = 
      tmp[0]*gcon[0][0] +
      tmp[1]*gcon[0][1] +
      tmp[2]*gcon[0][2] +
      tmp[3]*gcon[0][3] ;
    Ucon[1] = 
      tmp[0]*gcon[1][0] +
      tmp[1]*gcon[1][1] +
      tmp[2]*gcon[1][2] +
      tmp[3]*gcon[1][3] ;
    Ucon[2] = 
      tmp[0]*gcon[2][0] +
      tmp[1]*gcon[2][1] +
      tmp[2]*gcon[2][2] +
      tmp[3]*gcon[2][3] ;
    Ucon[3] = 
      tmp[0]*gcon[3][0] +
      tmp[1]*gcon[3][1] +
      tmp[2]*gcon[3][2] +
      tmp[3]*gcon[3][3] ;
  
    return ;
  }
 
  // safer version to recover four velocity

  // interpolate primitive variables first
  double U1A, U2A, U3A, U1B, U2B, U3B, tfac;
  double Vcon[NDIM];
  int nA, nB;
  set_tinterp_ns(X, &nA, &nB);
  tfac = (X[0] - data[nA]->t)/(data[nB]->t - data[nA]->t);
  U1A = interp_scalar(X, data[nA]->p[U1]);
  U2A = interp_scalar(X, data[nA]->p[U2]);
  U3A = interp_scalar(X, data[nA]->p[U3]);
  U1B = interp_scalar(X, data[nB]->p[U1]);
  U2B = interp_scalar(X, data[nB]->p[U2]);
  U3B = interp_scalar(X, data[nB]->p[U3]);
  Vcon[1] = tfac*U1A + (1. - tfac)*U1B;
  Vcon[2] = tfac*U2A + (1. - tfac)*U2B;
  Vcon[3] = tfac*U3A + (1. - tfac)*U3B;

  // translate to four velocity
  double VdotV = 0.;
  for (int i = 1; i < NDIM; i++)
    for (int j = 1; j < NDIM; j++)
      VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
  double Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
  Ucon[0] = -Vfac * gcon[0][0];
  for (int i = 1; i < NDIM; i++)
    Ucon[i] = Vcon[i] - Vfac * gcon[0][i];

}

void get_model_bcov(double X[NDIM], double Bcov[NDIM])
{
  if(X[1] < startx[1] || 
     X[1] > stopx[1]  || 
     X[2] < startx[2] || 
     X[2] > stopx[2]) {

      Bcov[0] = 0. ;
      Bcov[1] = 0. ;
      Bcov[2] = 0. ;
      Bcov[3] = 0. ;

    return ;
  }

  double Bcon[NDIM];
  double gcov[NDIM][NDIM];

  get_model_bcon(X, Bcon);
  gcov_func(X, gcov);

  lower(Bcon, gcov, Bcov);

}

void get_model_bcon(double X[NDIM], double Bcon[NDIM])
{
  if(X[1] < startx[1] || 
     X[1] > stopx[1]  || 
     X[2] < startx[2] || 
     X[2] > stopx[2]) {

      Bcon[0] = 0. ;
      Bcon[1] = 0. ;
      Bcon[2] = 0. ;
      Bcon[3] = 0. ;

    return;
  }

  int nA, nB;
  double tfac;

  // interpolate primitive variables first
  double B1A, B2A, B3A, B1B, B2B, B3B, Bcon1, Bcon2, Bcon3;
  set_tinterp_ns(X, &nA, &nB);
  tfac = (X[0] - data[nA]->t)/(data[nB]->t - data[nA]->t);
  B1A = interp_scalar(X, data[nA]->p[B1]);
  B2A = interp_scalar(X, data[nA]->p[B2]);
  B3A = interp_scalar(X, data[nA]->p[B3]);
  B1B = interp_scalar(X, data[nB]->p[B1]);
  B2B = interp_scalar(X, data[nB]->p[B2]);
  B3B = interp_scalar(X, data[nB]->p[B3]);
  Bcon1 = tfac*B1A + (1. - tfac)*B1B;
  Bcon2 = tfac*B2A + (1. - tfac)*B2B;
  Bcon3 = tfac*B3A + (1. - tfac)*B3B;

  double Ucon[NDIM], Ucov[NDIM];
  double gcov[NDIM][NDIM];

  gcov_func(X, gcov);
  get_model_ucon(X, Ucon);
  lower(Ucon, gcov, Ucov);

  Bcon[0] = Bcon1*Ucov[1] + Bcon2*Ucov[2] + Bcon3*Ucov[3];
  Bcon[1] = (Bcon1 + Ucon[1] * Bcon[0]) / Ucon[0];
  Bcon[2] = (Bcon2 + Ucon[2] * Bcon[0]) / Ucon[0];
  Bcon[3] = (Bcon3 + Ucon[3] * Bcon[0]) / Ucon[0];

}

double get_model_thetae(double X[NDIM])
{
  if(X[1] < startx[1] || 
     X[1] > stopx[1]  || 
     X[2] < startx[2] || 
     X[2] > stopx[2]) {
      return(0.) ;
  }
  
  double thetaeA, thetaeB, tfac;
  int nA, nB;
  set_tinterp_ns(X, &nA, &nB);
  tfac = (X[0] - data[nA]->t)/(data[nB]->t - data[nA]->t);
  thetaeA = interp_scalar(X, data[nA]->thetae);
  thetaeB = interp_scalar(X, data[nB]->thetae);


  double thetae = tfac*thetaeA + (1. - tfac)*thetaeB;
  if (thetae < 0.) {
    printf("thetae negative!\n");
    printf("X[] = %g %g %g %g\n", X[0], X[1], X[2], X[3]);
    printf("t = %e %e %e\n", data[0]->t, data[1]->t, data[2]->t);
    printf("thetae = %e tfac = %e thetaeA = %e thetaeB = %e nA = %i nB = %i\n",
    thetae, tfac, thetaeA, thetaeB, nA, nB);
  }
  return tfac*thetaeA + (1. - tfac)*thetaeB;
}

//b field strength in Gauss
double get_model_b(double X[NDIM])
{

  if(X[1] < startx[1] || 
     X[1] > stopx[1]  || 
     X[2] < startx[2] || 
     X[2] > stopx[2]) {
      return(0.) ;
  }
  
  double bA, bB, tfac;
  int nA, nB;
  set_tinterp_ns(X, &nA, &nB);
  tfac = (X[0] - data[nA]->t)/(data[nB]->t - data[nA]->t);
  bA = interp_scalar(X, data[nA]->b);
  bB = interp_scalar(X, data[nB]->b);

  return tfac*bA + (1. - tfac)*bB;
}

double get_model_ne(double X[NDIM])
{
  if(X[1] < startx[1] || 
     X[1] > stopx[1]  || 
     X[2] < startx[2] || 
     X[2] > stopx[2]) {
      return(0.) ;
  }

  double neA, neB, tfac;
  int nA, nB;
  set_tinterp_ns(X, &nA, &nB);
  tfac = (X[0] - data[nA]->t)/(data[nB]->t - data[nA]->t);
  neA = interp_scalar(X, data[nA]->ne);
  neB = interp_scalar(X, data[nB]->ne);
  return tfac*neA + (1. - tfac)*neB;
}


/** HARM utilities **/

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]) ;

/********************************************************************

        Interpolation routines

 ********************************************************************/

/* return fluid four-vector in simulation units */
void interp_fourv(double X[NDIM], double ****fourv, double Fourv[NDIM]){
  double del[NDIM],b1,b2,b3,d1,d2,d3,d4;
  int i, j, k, ip1, jp1, kp1;

  /* find the current zone location and offsets del[0], del[1] */
  Xtoijk(X, &i, &j, &k, del);

  // since we read from data, adjust i,j,k for ghost zones
  i += 1;
  j += 1;
  k += 1;

  ip1 = i + 1;
  jp1 = j + 1;
  kp1 = k + 1;
  
  b1 = 1.-del[1];
  b2 = 1.-del[2];
  b3 = 1.-del[3];

  d1 = b1*b2;
  d3 = del[1] * b2;
  d2 = del[2] * b1;
  d4 = del[1] * del[2];


  /* Interpolate along x1,x2 first */
  Fourv[0] = d1*fourv[i][j][k][0] + d2*fourv[i][jp1][k][0] + d3*fourv[ip1][j][k][0] + d4*fourv[ip1][jp1][k][0];
  Fourv[1] = d1*fourv[i][j][k][1] + d2*fourv[i][jp1][k][1] + d3*fourv[ip1][j][k][1] + d4*fourv[ip1][jp1][k][1];
  Fourv[2] = d1*fourv[i][j][k][2] + d2*fourv[i][jp1][k][2] + d3*fourv[ip1][j][k][2] + d4*fourv[ip1][jp1][k][2];
  Fourv[3] = d1*fourv[i][j][k][3] + d2*fourv[i][jp1][k][3] + d3*fourv[ip1][j][k][3] + d4*fourv[ip1][jp1][k][3];

  /* Now interpolate above in x3 */
  Fourv[0] = b3*Fourv[0] + del[3]*(d1*fourv[i][j][kp1][0] + d2*fourv[i][jp1][kp1][0] + d3*fourv[ip1][j][kp1][0] + d4*fourv[ip1][jp1][kp1][0]);
  Fourv[1] = b3*Fourv[1] + del[3]*(d1*fourv[i][j][kp1][1] + d2*fourv[i][jp1][kp1][1] + d3*fourv[ip1][j][kp1][1] + d4*fourv[ip1][jp1][kp1][1]);
  Fourv[2] = b3*Fourv[2] + del[3]*(d1*fourv[i][j][kp1][2] + d2*fourv[i][jp1][kp1][2] + d3*fourv[ip1][j][kp1][2] + d4*fourv[ip1][jp1][kp1][2]);
  Fourv[3] = b3*Fourv[3] + del[3]*(d1*fourv[i][j][kp1][3] + d2*fourv[i][jp1][kp1][3] + d3*fourv[ip1][j][kp1][3] + d4*fourv[ip1][jp1][kp1][3]);
  //new

  //no interpolation of vectors at all
 
  //Fourv[0]=fourv[i][j][k][0];
  //Fourv[1]=fourv[i][j][k][1];
  //Fourv[2]=fourv[i][j][k][2];
  //Fourv[3]=fourv[i][j][k][3];
  
}

/* return scalar in cgs units */
double interp_scalar(double X[NDIM], double ***var)
{
  double del[NDIM],b1,b2,interp;
  int i, j, k, ip1, jp1, kp1;

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

/***********************************************************************************

          End interpolation routines

 ***********************************************************************************/


void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM])
{
  // unless we're reading from data, i,j,k are the normal expected thing
  double phi;

  /* Map X[3] into sim range, assume startx[3] = 0 */
  phi = fmod(X[3], stopx[3]);
  //fold it to be positive and find index
  if(phi < 0.0) phi = stopx[3]+phi;
  
  //give index of a zone - zone index is moved to the grid zone center/
  //to account for the fact that all variables are reconstrucuted at zone centers?
  *i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
  *j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
  *k = (int) ((phi  - startx[3]) / dx[3] - 0.5 + 1000) - 1000;  

  del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
  del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
  del[3] = (phi - ((*k + 0.5) * dx[3] + startx[3])) / dx[3];
}

//#define SINGSMALL (1.E-20)
/* return boyer-lindquist coordinate of point */
void bl_coord(double X[NDIM], double *r, double *th)
{
  *r = exp(X[1]);

  if (DEREFINE_POLES) {
    double thG = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
    double y = 2*X[2] - 1.;
    double thJ = poly_norm*y*(1. + pow(y/poly_xt,poly_alpha)/(poly_alpha+1.)) + 0.5*M_PI;
    *th = thG + exp(mks_smooth*(startx[1] - X[1]))*(thJ - thG);
  } else {
    *th = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
  }
}

void coord(int i, int j, int k, double *X)
{
  // returns zone-centered X from i,j,k
  X[0] = startx[0];
  X[1] = startx[1] + (i + 0.5) * dx[1];
  X[2] = startx[2] + (j + 0.5) * dx[2];
  X[3] = startx[3] + (k + 0.5) * dx[3];
}


void set_units(char *munitstr)
{
  L_unit = GNEWT * MBH / (CL * CL);
  T_unit = L_unit / CL;
  RHO_unit = M_unit / pow(L_unit, 3);
  U_unit = RHO_unit * CL * CL;
  B_unit = CL * sqrt(4.*M_PI*RHO_unit);
  Mdotedd=4.*M_PI*GNEWT*MBH*MP/CL/0.1/SIGMA_THOMSON;

  fprintf(stderr,"L,T,M units: %g [cm] %g [s] %g [g]\n",L_unit,T_unit,M_unit) ;
  fprintf(stderr,"rho,u,B units: %g [g cm^-3] %g [g cm^-1 s^-2] %g [G] \n",RHO_unit,U_unit,B_unit) ;
}

void init_physical_quantities(int n)
{
  int i, j, k;
  double bsq,sigma_m;

  // cover everything, even ghost zones
  for (i = 0; i < N1+2; i++) {
    for (j = 0; j < N2+2; j++) {
      for (k = 0; k < N3+2; k++) {
        data[n]->ne[i][j][k] = data[n]->p[KRHO][i][j][k] * RHO_unit/(MP+ME) ;

        bsq = data[n]->bcon[i][j][k][0] * data[n]->bcov[i][j][k][0] +
              data[n]->bcon[i][j][k][1] * data[n]->bcov[i][j][k][1] +
              data[n]->bcon[i][j][k][2] * data[n]->bcov[i][j][k][2] +
              data[n]->bcon[i][j][k][3] * data[n]->bcov[i][j][k][3] ;

        data[n]->b[i][j][k] = sqrt(bsq)*B_unit;
        sigma_m = bsq/data[n]->p[KRHO][i][j][k];

        // beta presciption
        //beta=p[UU][i][j][k]*(gam-1.)/0.5/bsq;
        //b2=pow(beta,2);
        //trat = trat_d * b2/(1. + b2) + trat_j /(1. + b2);
        //Thetae_unit = (gam - 1.) * (MP / ME) / trat;
        
        if (ELECTRONS) {
          data[n]->thetae[i][j][k] = data[n]->p[KEL][i][j][k]*pow(data[n]->p[KRHO][i][j][k],game-1.)*Thetae_unit;
        } else {
          data[n]->thetae[i][j][k] = Thetae_unit*data[n]->p[UU][i][j][k]/data[n]->p[KRHO][i][j][k];
        }
        data[n]->thetae[i][j][k] = MAX(data[n]->thetae[i][j][k], 1.e-3);
       
        //thetae[i][j][k] = (gam-1.)*MP/ME*p[UU][i][j][k]/p[KRHO][i][j][k];
        //printf("rho = %e thetae = %e\n", p[KRHO][i][j][k], thetae[i][j][k]);

        //data[n]->thetae[i][j][k] = Thetae_unit*data[n]->p[UU][i][j][k]/data[n]->p[KRHO][i][j][k]/4.;
        //printf("Thetae_unit = %e Thetae = %e\n", Thetae_unit, thetae[i][j][k]);
        
        //strongly magnetized = empty, no shiny spine
        if (sigma_m > 1.0) data[n]->ne[i][j][k]=0.0;
      }
    }
  }
}

// malloc utilities
void *malloc_rank1(int n1, int size)
{
  void *A;

  if ((A = malloc(n1*size)) == NULL) {
    fprintf(stderr,"malloc failure in malloc_rank1\n");
    exit(123);
  }

  return A;
}

double **malloc_rank2(int n1, int n2)
{

  double **A;
  double *space;
  int i;

  space = malloc_rank1(n1*n2, sizeof(double));
  A = malloc_rank1(n1, sizeof(double *));
  for(i = 0; i < n1; i++) A[i] = &(space[i*n2]);

  return A;
}


double ***malloc_rank3(int n1, int n2, int n3)
{

  double ***A;
  double *space;
  int i,j;

  space = malloc_rank1(n1*n2*n3, sizeof(double));
  A = malloc_rank1(n1, sizeof(double *));
  for(i = 0; i < n1; i++){
    A[i] = malloc_rank1(n2,sizeof(double *));
    for(j = 0; j < n2; j++){
      A[i][j] = &(space[n3*(j + n2*i)]);
    }
  }

  return A;
}

float **malloc_rank2_float(int n1, int n2)
{

  float **A;
  float *space;
  int i;

  space = malloc_rank1(n1*n2, sizeof(float));
  A = malloc_rank1(n1, sizeof(float *));
  for(i = 0; i < n1; i++) A[i] = &(space[i*n2]);

  return A;
}


float ***malloc_rank3_float(int n1, int n2, int n3)
{

  float ***A;
  float *space;
  int i,j;

  space = malloc_rank1(n1*n2*n3, sizeof(float));
  A = malloc_rank1(n1, sizeof(float *));
  for(i = 0; i < n1; i++){
    A[i] = malloc_rank1(n2,sizeof(float *));
    for(j = 0; j < n2; j++){
      A[i][j] = &(space[n3*(j + n2*i)]);
    }
  }

  return A;
}

float ****malloc_rank4_float(int n1, int n2, int n3, int n4)
{

  float ****A;
  float *space;
  int i,j,k;

  space = malloc_rank1(n1*n2*n3*n4, sizeof(float));
  A = malloc_rank1(n1, sizeof(float *));
  for(i=0;i<n1;i++){
    A[i] = malloc_rank1(n2,sizeof(float *));
    for(j=0;j<n2;j++){
      A[i][j] = malloc_rank1(n3,sizeof(float *));
      for(k=0;k<n3;k++){
        A[i][j][k] = &(space[n4*(k + n3*(j + n2*i))]);
      }
    }
  }

  return A;
}


double ****malloc_rank4(int n1, int n2, int n3, int n4)
{

  double ****A;
  double *space;
  int i,j,k;

  space = malloc_rank1(n1*n2*n3*n4, sizeof(double));
  A = malloc_rank1(n1, sizeof(double *));
  for(i=0;i<n1;i++){
    A[i] = malloc_rank1(n2,sizeof(double *));
    for(j=0;j<n2;j++){
      A[i][j] = malloc_rank1(n3,sizeof(double *));
      for(k=0;k<n3;k++){
        A[i][j][k] = &(space[n4*(k + n3*(j + n2*i))]);
      }
    }
  }

  return A;
}

double *****malloc_rank5(int n1, int n2, int n3, int n4, int n5)
{

  double *****A;
  double *space;
  int i,j,k,l;

  space = malloc_rank1(n1*n2*n3*n4*n5, sizeof(double));
  A = malloc_rank1(n1, sizeof(double *));
  for(i=0;i<n1;i++){
    A[i] = malloc_rank1(n2, sizeof(double *));
    for(j=0;j<n2;j++){
      A[i][j] = malloc_rank1(n3, sizeof(double *));
      for(k=0;k<n3;k++){
        A[i][j][k] = malloc_rank1(n4, sizeof(double *));
        for(l=0;l<n4;l++){
          A[i][j][k][l] = &(space[n5*(l + n4*(k + n3*(j + n2*i)))]);
        }
      }
    }
  }

  return A;
}

void init_storage(void)
{
  // one ghost zone on each side of the domain
  for (int n = 0; n < NSUP; n++) {
    data[n]->bcon = malloc_rank4(N1+2,N2+2,N3+2,NDIM);
    data[n]->bcov = malloc_rank4(N1+2,N2+2,N3+2,NDIM);
    data[n]->ucon = malloc_rank4(N1+2,N2+2,N3+2,NDIM);
    data[n]->ucov = malloc_rank4(N1+2,N2+2,N3+2,NDIM);
    data[n]->p = malloc_rank4(NVAR,N1+2,N2+2,N3+2);
    //p = (double ****)malloc_rank1(NVAR,sizeof(double *));
    //for(i = 0; i < NVAR; i++) p[i] = malloc_rank3(N1,N2,N3);
    data[n]->ne = malloc_rank3(N1+2,N2+2,N3+2);
    data[n]->thetae = malloc_rank3(N1+2,N2+2,N3+2);
    data[n]->b = malloc_rank3(N1+2,N2+2,N3+2);
  }
}

/* HDF5 v1.8 API */
#include <hdf5.h>
#include <hdf5_hl.h>

void init_iharm_grid(char *fname)
{
  // called at the beginning of the run and sets the static parameters
  // along with setting up the grid
 
  fprintf(stderr, "filename: %s\n", fname);

  printf("init grid\n");

  if ( hdf5_open(fname) < 0 ) {
    fprintf(stderr, "! unable to open file %s. exiting!\n", fname);
    exit(-2);
  }

  hdf5_read_single_val(&t0, "t", H5T_IEEE_F64LE);

  hdf5_set_directory("/header/");

  if ( hdf5_exists("has_electrons") )
    hdf5_read_single_val(&ELECTRONS, "has_electrons", H5T_STD_I32LE);
  if ( hdf5_exists("has_radiation") ) 
    hdf5_read_single_val(&RADIATION, "has_radiation", H5T_STD_I32LE);
  if ( hdf5_exists("has_derefine_poles") )
    hdf5_read_single_val(&DEREFINE_POLES, "has_derefine_poles", H5T_STD_I32LE);

  ELECTRONS = 0; // force TP_OVER_TE to overwrite bad electrons

  char metric[20];
  hid_t HDF5_STR_TYPE = hdf5_make_str_type(20);
  hdf5_read_single_val(&metric, "metric", HDF5_STR_TYPE);

  if ( strncmp(metric, "MMKS", 19) == 0 ) 
    DEREFINE_POLES = 1;

  hdf5_read_single_val(&N1, "n1", H5T_STD_I32LE);
  hdf5_read_single_val(&N2, "n2", H5T_STD_I32LE);
  hdf5_read_single_val(&N3, "n3", H5T_STD_I32LE);
  hdf5_read_single_val(&gam, "gam", H5T_IEEE_F64LE);

  if (ELECTRONS) {
    fprintf(stderr, "custom electron model loaded...\n");
    hdf5_read_single_val(&game, "gam_e", H5T_IEEE_F64LE);
    hdf5_read_single_val(&gamp, "gam_p", H5T_IEEE_F64LE);
    Thetae_unit = MP/ME;
  } else {
    Thetae_unit = MP/ME*(gam-1.)*1./(1. + tp_over_te);
  }

  if (RADIATION) {
    fprintf(stderr, "custom radiation field tracking information loaded...\n");
    hdf5_set_directory("/header/units/");
    hdf5_read_single_val(&M_unit, "M_unit", H5T_IEEE_F64LE);
    hdf5_read_single_val(&T_unit, "T_unit", H5T_IEEE_F64LE);
    hdf5_read_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
    hdf5_read_single_val(&Thetae_unit, "Thetae_unit", H5T_IEEE_F64LE);
    hdf5_read_single_val(&MBH, "Mbh", H5T_IEEE_F64LE);
    hdf5_read_single_val(&tp_over_te, "tp_over_te", H5T_IEEE_F64LE);
  }

  hdf5_set_directory("/header/geom/");
  hdf5_read_single_val(&startx[1], "startx1", H5T_IEEE_F64LE);
  hdf5_read_single_val(&startx[2], "startx2", H5T_IEEE_F64LE);
  hdf5_read_single_val(&startx[3], "startx3", H5T_IEEE_F64LE);
  hdf5_read_single_val(&dx[1], "dx1", H5T_IEEE_F64LE);
  hdf5_read_single_val(&dx[2], "dx2", H5T_IEEE_F64LE);
  hdf5_read_single_val(&dx[3], "dx3", H5T_IEEE_F64LE);

  hdf5_set_directory("/header/geom/mks/");
  if ( DEREFINE_POLES ) hdf5_set_directory("/header/geom/mmks/");

  hdf5_read_single_val(&a, "a", H5T_IEEE_F64LE);
  hdf5_read_single_val(&hslope, "hslope", H5T_IEEE_F64LE);
  hdf5_read_single_val(&Rin, "Rin", H5T_IEEE_F64LE);
  hdf5_read_single_val(&Rout, "Rout", H5T_IEEE_F64LE);

  if (DEREFINE_POLES) {
    fprintf(stderr, "custom refinement at poles loaded...\n");
    hdf5_read_single_val(&poly_xt, "poly_xt", H5T_IEEE_F64LE);
    hdf5_read_single_val(&poly_alpha, "poly_alpha", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mks_smooth, "mks_smooth", H5T_IEEE_F64LE);
    poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha + 1.)*1./pow(poly_xt, poly_alpha));
 }

  rmax = MIN(50., Rout);

  hdf5_set_directory("/");
  hdf5_read_single_val(&DTd, "dump_cadence", H5T_IEEE_F64LE);
  
  //startx[1] += 3*dx[1];
  //startx[2] += 3*dx[2];
  //      startx[3] += 3*dx[3];
  
  fprintf(stdout,"start: %g %g %g \n",startx[1],startx[2],startx[3]);
  //below is equivalent to the above
  /*
  startx[1] = log(Rin-R0);       
        startx[2] = th_cutout/M_PI ; 
  */

  //fprintf(stdout,"th_cutout: %g  %d x %d x %d\n",th_cutout,N1,N2,N3);
  

  //th_beg=th_cutout;
  //th_end=M_PI-th_cutout;
  //th_len = th_end-th_beg;

  // Ignore radiation interactions within one degree of polar axis
  th_beg = 0.0174;
  //th_end = 3.1241;

  stopx[0] = 1.;
  stopx[1] = startx[1]+N1*dx[1];
  stopx[2] = startx[2]+N2*dx[2];
  stopx[3] = startx[3]+N3*dx[3];

  fprintf(stdout,"stop: %g %g %g \n",stopx[1],stopx[2],stopx[3]);

  init_storage();

  hdf5_close();
}

void load_iharm_data(int n, char *fname)
{
  // loads relevant information from fluid dump file stored at fname
  // to the n'th copy of data (e.g., for slow light)

  printf("LOADING DATA\n");
  nloaded++;

  int i,j,k,l,m;
  double X[NDIM],UdotU,ufac,udotB;
  double gcov[NDIM][NDIM],gcon[NDIM][NDIM], g;
  double dMact, Ladv;
  double r,th;
 
  if ( hdf5_open(fname) < 0 ) {
    fprintf(stderr, "! unable to open file %s. Exiting!\n", fname);
    exit(-1);
  }

  hdf5_set_directory("/");

  int n_prims;
  hdf5_read_single_val(&n_prims, "/header/n_prim", H5T_STD_I32LE);

  // load into "center" of data
  hsize_t fdims[] = { N1, N2, N3, n_prims };
  hsize_t fstart[] = { 0, 0, 0, 0 };
  hsize_t fcount[] = { N1, N2, N3, 1 };
  hsize_t mdims[] = { N1+2, N2+2, N3+2, 1 };
  hsize_t mstart[] = { 1, 1, 1, 0 };

  fstart[3] = 0;
  hdf5_read_array(data[n]->p[KRHO][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
  fstart[3] = 1;
  hdf5_read_array(data[n]->p[UU][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
  fstart[3] = 2;
  hdf5_read_array(data[n]->p[U1][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
  fstart[3] = 3;
  hdf5_read_array(data[n]->p[U2][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
  fstart[3] = 4;
  hdf5_read_array(data[n]->p[U3][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
  fstart[3] = 5;
  hdf5_read_array(data[n]->p[B1][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
  fstart[3] = 6;
  hdf5_read_array(data[n]->p[B2][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
  fstart[3] = 7;
  hdf5_read_array(data[n]->p[B3][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE); 

  if (ELECTRONS) {
    fstart[3] = 8;
    hdf5_read_array(data[n]->p[KEL][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
    fstart[3] = 9;
    hdf5_read_array(data[n]->p[KTOT][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
  }

  hdf5_read_single_val(&(data[n]->t), "t", H5T_IEEE_F64LE);

  hdf5_close();

  X[0] = 0.;
  //X[3] = 0.;

  dMact = Ladv = 0.;

  // construct four-vectors over "real" zones
  for(i = 1; i < N1+1; i++){
    X[1] = startx[1] + ( i - 0.5)*dx[1];
    for(j = 1; j < N2+1; j++){
      X[2] = startx[2] + (j-0.5)*dx[2];
      gcov_func(X, gcov); // in system with cut off
      gcon_func(gcov, gcon);
      g = gdet_func(gcov);

      bl_coord(X, &r, &th);

      for(k = 1; k < N3+1; k++){
        UdotU = 0.;
        X[3] = startx[3] + (k-0.5)*dx[3];
        
        // the four-vector reconstruction should have gcov and gcon and gdet using the modified coordinates
        // interpolating the four vectors to the zone center !!!!
        for(l = 1; l < NDIM; l++) 
          for(m = 1; m < NDIM; m++) 
            UdotU += gcov[l][m]*data[n]->p[U1+l-1][i][j][k]*data[n]->p[U1+m-1][i][j][k];
        ufac = sqrt(-1./gcon[0][0]*(1 + fabs(UdotU)));
        data[n]->ucon[i][j][k][0] = -ufac*gcon[0][0];
        for(l = 1; l < NDIM; l++) 
          data[n]->ucon[i][j][k][l] = data[n]->p[U1+l-1][i][j][k] - ufac*gcon[0][l];
        lower(data[n]->ucon[i][j][k], gcov, data[n]->ucov[i][j][k]);

        // reconstruct the magnetic field three vectors
        udotB = 0.;
        
        for (l = 1; l < NDIM; l++) {
          udotB += data[n]->ucov[i][j][k][l]*data[n]->p[B1+l-1][i][j][k];
        }
        
        data[n]->bcon[i][j][k][0] = udotB;

        for (l = 1; l < NDIM; l++) {
          data[n]->bcon[i][j][k][l] = (data[n]->p[B1+l-1][i][j][k] + data[n]->ucon[i][j][k][l]*udotB)/data[n]->ucon[i][j][k][0];
        }

        lower(data[n]->bcon[i][j][k], gcov, data[n]->bcov[i][j][k]);

        if(i <= 21) { dMact += g * data[n]->p[KRHO][i][j][k] * data[n]->ucon[i][j][k][1]; }
        if(i >= 21 && i < 41 && 0) Ladv += g * data[n]->p[UU][i][j][k] * data[n]->ucon[i][j][k][1] * data[n]->ucov[i][j][k][0] ;
        if(i <= 21) Ladv += g * data[n]->p[UU][i][j][k] * data[n]->ucon[i][j][k][1] * data[n]->ucov[i][j][k][0] ;

      }
    }
  }

  // now copy primitives and four-vectors according to boundary conditions

  // radial -- just extend zones
  for (j=1; j<N2+1; ++j) {
    for (k=1; k<N3+1; ++k) {
      for (l=0; l<NDIM; ++l) {
        data[n]->bcon[0][j][k][l] = data[n]->bcon[1][j][k][l];
        data[n]->bcon[N1+1][j][k][l] = data[n]->bcon[N1][j][k][l];
        data[n]->bcov[0][j][k][l] = data[n]->bcov[1][j][k][l];
        data[n]->bcov[N1+1][j][k][l] = data[n]->bcov[N1][j][k][l];
        data[n]->ucon[0][j][k][l] = data[n]->ucon[1][j][k][l];
        data[n]->ucon[N1+1][j][k][l] = data[n]->ucon[N1][j][k][l];
        data[n]->ucov[0][j][k][l] = data[n]->ucov[1][j][k][l];
        data[n]->ucov[N1+1][j][k][l] = data[n]->ucov[N1][j][k][l];
      }
      for (l=0; l<NVAR; ++l) {
        data[n]->p[l][0][j][k] = data[n]->p[l][1][j][k];
        data[n]->p[l][N1+1][j][k] = data[n]->p[l][N1][j][k];
      }
    }
  }

  // elevation -- flip (this is a rotation by pi)
  for (i=0; i<N1+2; ++i) {
    for (k=1; k<N3+1; ++k) {
      if (N3%2 == 0) {
        int kflip = ( k + (N3/2) ) % N3;
        for (l=0; l<NDIM; ++l) {
          data[n]->bcon[i][0][k][l] = data[n]->bcon[i][1][kflip][l];
          data[n]->bcon[i][N2+1][k][l] = data[n]->bcon[i][N2][kflip][l];
          data[n]->bcov[i][0][k][l] = data[n]->bcov[i][1][kflip][l];
          data[n]->bcov[i][N2+1][k][l] = data[n]->bcov[i][N2][kflip][l];
          data[n]->ucon[i][0][k][l] = data[n]->ucon[i][1][kflip][l];
          data[n]->ucon[i][N2+1][k][l] = data[n]->ucon[i][N2][kflip][l];
          data[n]->ucov[i][0][k][l] = data[n]->ucov[i][1][kflip][l];
          data[n]->ucov[i][N2+1][k][l] = data[n]->ucov[i][N2][kflip][l];
        }
        for (l=0; l<NVAR; ++l) {
          data[n]->p[l][i][0][k] = data[n]->p[l][i][1][kflip];
          data[n]->p[l][i][N2+1][k] = data[n]->p[l][i][N2][kflip];
        }
      } else {
        int kflip1 = ( k + (N3/2) ) % N3;
        int kflip2 = ( k + (N3/2) + 1 ) % N3;
        for (l=0; l<NDIM; ++l) {
          data[n]->bcon[i][0][k][l]    = ( data[n]->bcon[i][1][kflip1][l] 
                                         + data[n]->bcon[i][1][kflip2][l] ) / 2.;
          data[n]->bcon[i][N2+1][k][l] = ( data[n]->bcon[i][N2][kflip1][l]
                                         + data[n]->bcon[i][N2][kflip2][l] ) / 2.;
          data[n]->bcov[i][0][k][l]    = ( data[n]->bcov[i][1][kflip1][l]
                                         + data[n]->bcov[i][1][kflip2][l] ) / 2.;
          data[n]->bcov[i][N2+1][k][l] = ( data[n]->bcov[i][N2][kflip1][l] 
                                         + data[n]->bcov[i][N2][kflip2][l] ) / 2.;
          data[n]->ucon[i][0][k][l]    = ( data[n]->ucon[i][1][kflip1][l]
                                         + data[n]->ucon[i][1][kflip2][l] ) / 2.;
          data[n]->ucon[i][N2+1][k][l] = ( data[n]->ucon[i][N2][kflip1][l]
                                         + data[n]->ucon[i][N2][kflip2][l] ) / 2.;
          data[n]->ucov[i][0][k][l]    = ( data[n]->ucov[i][1][kflip1][l] 
                                         + data[n]->ucov[i][1][kflip2][l] ) / 2.;
          data[n]->ucov[i][N2+1][k][l] = ( data[n]->ucov[i][N2][kflip1][l] 
                                         + data[n]->ucov[i][N2][kflip2][l] ) / 2.;
        }
        for (l=0; l<NVAR; ++l) {
          data[n]->p[l][i][0][k]    = ( data[n]->p[l][i][1][kflip1]
                                      + data[n]->p[l][i][1][kflip2] ) / 2.;
          data[n]->p[l][i][N2+1][k] = ( data[n]->p[l][i][N2][kflip1]
                                      + data[n]->p[l][i][N2][kflip2] ) / 2.;
        }
      }
    }
  }

  // azimuth -- periodic
  for (i=0; i<N1+2; ++i) {
    for (j=0; j<N2+2; ++j) {
      for (l=0; l<NDIM; ++l) {
        data[n]->bcon[i][j][0][l] = data[n]->bcon[i][j][N3][l];
        data[n]->bcon[i][j][N3+1][l] = data[n]->bcon[i][j][1][l];
        data[n]->bcov[i][j][0][l] = data[n]->bcov[i][j][N3][l];
        data[n]->bcov[i][j][N3+1][l] = data[n]->bcov[i][j][1][l];
        data[n]->ucon[i][j][0][l] = data[n]->ucon[i][j][N3][l];
        data[n]->ucon[i][j][N3+1][l] = data[n]->ucon[i][j][1][l];
        data[n]->ucov[i][j][0][l] = data[n]->ucov[i][j][N3][l];
        data[n]->ucov[i][j][N3+1][l] = data[n]->ucov[i][j][1][l];
      }
      for (l=0; l<NVAR; ++l) {
        data[n]->p[l][i][j][0] = data[n]->p[l][i][j][N3];
        data[n]->p[l][i][j][N3+1] = data[n]->p[l][i][j][1];
      }
    }
  }

  dMact *= dx[3]*dx[2] ;
  dMact /= 21. ;
  Ladv *= dx[3]*dx[2] ;
  Ladv /= 21. ;

  fprintf(stderr,"dMact: %g [code]\n",dMact) ;
  fprintf(stderr,"Ladv: %g [code]\n",Ladv) ;
  fprintf(stderr,"Mdot: %g [g/s] \n",-dMact*M_unit/T_unit) ;
  fprintf(stderr,"Mdot: %g [MSUN/YR] \n",-dMact*M_unit/T_unit/(MSUN / YEAR)) ;
  fprintf(stderr,"Mdot: %g [Mdotedd]\n",-dMact*M_unit/T_unit/Mdotedd) ;
  fprintf(stderr,"Mdotedd: %g [g/s]\n",Mdotedd) ;
  fprintf(stderr,"Mdotedd: %g [MSUN/YR]\n",Mdotedd/(MSUN/YEAR)) ;

  // now construct useful scalar quantities (over full (+ghost) zones of data)
  init_physical_quantities(n);
}

double root_find(double x[NDIM])
{
    double th = x[2];
    double thb, thc;
    double dtheta_func(double X[NDIM]), theta_func(double X[NDIM]);

    double Xa[NDIM], Xb[NDIM], Xc[NDIM];
    Xa[1] = log(x[1]);
    Xa[3] = x[3];
    Xb[1] = Xa[1];
    Xb[3] = Xa[3];
    Xc[1] = Xa[1];
    Xc[3] = Xa[3];

    if (x[2] < M_PI / 2.) {
      Xa[2] = 0. - SMALL;
      Xb[2] = 0.5 + SMALL;
    } else {
      Xa[2] = 0.5 - SMALL;
      Xb[2] = 1. + SMALL;
    }

    //tha = theta_func(Xa);
    thb = theta_func(Xb);

    /* bisect for a bit */
    double tol = 1.e-6;
    for (int i = 0; i < 100; i++) {
      Xc[2] = 0.5 * (Xa[2] + Xb[2]);
      thc = theta_func(Xc);

      if ((thc - th) * (thb - th) < 0.)
        Xa[2] = Xc[2];
      else
        Xb[2] = Xc[2];

      double err = theta_func(Xc) - th;
      if (fabs(err) < tol) break;
    }

    return (Xa[2]);
}

/*this does not depend on theta cut-outs there is no squizzing*/
double theta_func(double X[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);
  return th;
}

int radiating_region(double X[NDIM])
{
  if (X[1] < log(rmax) && X[2]>th_beg/M_PI && X[2]<(1.-th_beg/M_PI) ) {
    return 1;
  } else {
    return 0;
  }
}

