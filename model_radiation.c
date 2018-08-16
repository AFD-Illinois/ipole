
/**********************************************************/

/*** all you need to make a polarized radiative transfer***/
/***** used in ipole to evolve complex tensor N *******/
/***** along with standard evolution for I scalar *********/
/**********************************************************/
/**** written by Monika Moscibrodzka on 09 July 2014 ******/
/************ @ Eindhoven Airport *************************/
/************  last update: 9 May 2017   ******************/

/****************and then rewritten by C.Gammie ***********/

/**********************************************************/

#include "decs.h"
#include "defs.h"

/* the following definitions are used only locally */
#define S2     (1.41421356237310)	//sqrt(2)
#define S3     (1.73205080756888)	//sqrt(3)

/* declarations of local functions */
/* thermal plasma emissivity, absorptivity and Faraday conversion and rotation */

double g(double Xe);
double h(double Xe);
double Je(double Xe);
double jffunc(double Xe);
double I_I(double x);
double I_Q(double x);
double I_V(double x);
double besselk_asym(int n, double x);


/*************************SUPPORTING FUNCTIONS******************************/

/*invariant plasma emissivities/abs/rot in tetrad frame */
void jar_calc(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV)
{
  double nu, Thetae, Ne, B, theta, nusq;
  double x, Xe, omega0, nuc;
  double Bnuinv;
  double Ucov[NDIM];
  double Thetaer, wp2;

  Ne = get_model_ne(X);
  get_model_ucov(X, Ucov);
  if (isnan(Ucov[0])) {
    printf("UCOV[0] is nan! thread = %i\n", omp_get_thread_num());
    printf("X[] = %e %e %e %e\n", X[0],X[1],X[2],X[3]);
    printf("K[] = %e %e %e %e\n", Kcon[0],Kcon[1],Kcon[2],Kcon[3]);
    printf("Ne = %e\n", Ne);
  }
  theta = get_bk_angle(X, Kcon, Ucov);	/* angle between k & b  */

  if (theta <= 0. || theta >= M_PI) {	/* no emission/absorption along field  */

    *jI = 0.0;
    *jQ = 0.0;
    *jU = 0.0;
    *jV = 0.0;

    *aI = 0.0;
    *aQ = 0.0;
    *aU = 0.0;
    *aV = 0.0;

    *rQ = 0.0;
    *rU = 0.0;

    nu = get_fluid_nu(Kcon, Ucov);	/* freqcgs1;  freq in Hz */
    wp2 = 4. * M_PI * Ne * EE * EE / ME;
    B = get_model_b(X);		/* field in G                */
    omega0 = EE * B / ME / CL;
    Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
    Thetaer = 1. / Thetae;
    /* Faraday rotativities for thermal plasma */
    Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));
    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);

    /* invariant rotativities */
    *rV *= nu;

    if (isnan(*rV)  || *rV > 1.e100 || *rV < -1.e100) {
      printf("NAN RV theta! rV = %e nu = %e Ne = %e Thetae = %e\n", *rV, nu, Ne, Thetae);
    }

    return;

  } else {

    nu = get_fluid_nu(Kcon, Ucov);	/* freqcgs1;  freq in Hz */
    nusq = nu * nu;
    B = get_model_b(X);	/* field in G                */
    Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
    Thetaer = 1. / Thetae;

    omega0 = EE * B / ME / CL;
    wp2 = 4. * M_PI * Ne * EE * EE / ME;

    /* Faraday rotativities for thermal plasma */
    Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));

    /* Here I use approximate bessel functions to match rhoqv with grtrans */
    *rQ = 2. * M_PI * nu / 2. / CL * wp2 * omega0 * omega0 / pow(2 * M_PI * nu, 4) *
      jffunc(Xe) * (besselk_asym(1, Thetaer) / besselk_asym(2, Thetaer) +
          6. * Thetae) * sin(theta) * sin(theta);
    *rU = 0.0;
    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);

    /*synchrotron emissivity */
    nuc = 3.0 * EE * B * sin(theta) / 4.0 / M_PI / ME / CL * Thetae * Thetae + 1.0;
    x = nu / nuc;

    *jI = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_I(x);	// [g/s^2/cm = ergs/s/cm^3]
    *jQ = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_Q(x);
    *jU = 0.0;		// convention; depends on tetrad
    *jV = 2. * Ne * EE * EE * nu / tan(theta) / 3. / S3 / CL / Thetae / Thetae / Thetae * I_V(x);

    /* invariant emissivity */
    *jI /= nusq;
    *jQ /= nusq;
    *jU /= nusq;
    *jV /= nusq;

    /* invariant synchrotron absorptivity */
    Bnuinv = Bnu_inv(nu, Thetae);   /* Planck function */
    *aI = *jI / Bnuinv;
    *aQ = *jQ / Bnuinv;
    *aU = *jU / Bnuinv;
    *aV = *jV / Bnuinv;

    /* invariant rotativities */
    *rQ *= nu;
    *rU *= nu;
    *rV *= nu;

    if (isnan(*rV) || *rV > 1.e100 || *rV < -1.e100) {
      printf("NAN RV! rV = %e nu = %e Ne = %e Thetae = %e x = %e\n", *rV, nu, Ne, Thetae, x);
      printf("B = %e\n", B);
      for (int mu = 0; mu < NDIM; mu++) {
        printf("[%i] X[] = %e Kcon[] = %e Ucov[] = %e\n", mu,X[mu],Kcon[mu],Ucov[mu]);
      }
    }

  }

}

/*emissivity functions and functions used for Faraday conversion and rotation*/
/*from J. Dexter PhD thesis (checked with Leung harmony program, and Huang & Shcherbakov 2011*/
double g(double Xe)
{
  return 1. - 0.11 * log(1 + 0.035 * Xe);
}

double h(double Xe)
{
  return 2.011 * exp(-pow(Xe, 1.035) / 4.7) -
    cos(Xe * 0.5) * exp(-pow(Xe, 1.2) / 2.73) -
    0.011 * exp(-Xe / 47.2);
}

double Je(double Xe)
{
  return 0.43793091 * log(1. + 0.00185777 * pow(Xe, 1.50316886));
}

double jffunc(double Xe)
{
  double extraterm;
  extraterm =
    (0.011 * exp(-Xe / 47.2) -
     pow(2., -1. / 3.) / pow(3.,
       23. / 6.) * M_PI * 1e4 * pow(Xe + 1e-16,
         -8. / 3.)) *
    (0.5 + 0.5 * tanh((log(Xe) - log(120.)) / 0.1));
  return 2.011 * exp(-pow(Xe, 1.035) / 4.7) -
    cos(Xe * 0.5) * exp(-pow(Xe, 1.2) / 2.73) -
    0.011 * exp(-Xe / 47.2) + extraterm;
}

double I_I(double x)
{
  return 2.5651 * (1 + 1.92 * pow(x, -1. / 3.) +
      0.9977 * pow(x, -2. / 3.)) * exp(-1.8899 * pow(x,
        1. /
        3.));
}

double I_Q(double x)
{
  return 2.5651 * (1 + 0.93193 * pow(x, -1. / 3.) +
      0.499873 * pow(x, -2. / 3.)) * exp(-1.8899 * pow(x,
        1. /
        3.));
}

double I_V(double x)
{
  return (1.81348 / x + 3.42319 * pow(x, -2. / 3.) +
      0.0292545 * pow(x, -0.5) + 2.03773 * pow(x,
        -1. / 3.)) *
    exp(-1.8899 * pow(x, 1. / 3.));
}

double besselk_asym(int n, double x)
{

  if (n == 0)
    return -log(x / 2.) - 0.5772;

  if (n == 1)
    return 1. / x;

  if (n == 2)
    return 2. / x / x;

  fprintf(stderr,"this cannot happen\n");
  exit(1);
}

/* end of emissivity functions */

#undef S2
#undef S3

/*int radiating_region(double X[4])
  {
  if (X[1] < log(rmax) && X[2]>th_beg/M_PI && X[2]<(1.-th_beg/M_PI) ) {
  return 1;
  } else {
  return 0;
  }
  }*/

/* 

   thermal synchrotron emissivity 

   interpolates between Petrosian limit and
   classical thermal synchrotron limit

   good for Thetae >~ 1

*/

double jnu_synch(double nu, double Ne, double Thetae, double B, double theta)
{
  double K2,nuc,nus,x,f,j,sth ;

  //K2 = gsl_sf_bessel_Kn(2,1./Thetae) ;
  K2 = 2.*Thetae*Thetae ;

  nuc = EE*B/(2.*M_PI*ME*CL) ;
  sth = sin(theta) ;
  nus = (2./9.)*nuc*Thetae*Thetae*sth ;
  if(nu > 1.e12*nus) return(0.) ;
  x = nu/nus ;
  f = pow( pow(x,1./2.) + pow(2.,11./12.)*pow(x,1./6.), 2 ) ;
  j = (sqrt(2.)*M_PI*EE*EE*Ne*nus/(3.*CL*K2)) * f * exp(-pow(x,1./3.)) ;

  return(j) ;
}

/* get the invariant emissivity and opacity at a given position
   for a given wavevector */

void get_jkinv(double X[NDIM], double Kcon[NDIM], double *jnuinv,
    double *knuinv)
{
  double nu, theta, B, Thetae, Ne, Bnuinv;
  double Ucov[NDIM], Bcov[NDIM];
  double Ucon[NDIM], Bcon[NDIM];
  double Kcov[NDIM], gcov[NDIM][NDIM];

  /* get fluid parameters */
  Ne = get_model_ne(X);	/* check to see if we're outside fluid model */

  if (Ne == 0.) {
    *jnuinv = 0.;
    *knuinv = 0.;
    return;
  }

  /* get covariant four-velocity of fluid for use in get_bk_angle and get_fluid_nu */
  get_model_ucov(X, Ucov);
  get_model_bcov(X, Bcov);

  /*extra: print out stuff to test tetrads */
  get_model_ucon(X, Ucon);
  get_model_bcon(X, Bcon);

  gcov_func(X, gcov);
  lower(Kcon, gcov, Kcov);

  //theta = M_PI/2.;//get_bk_angle(X,Kcon,Ucov) ; /* angle between k & b */
  theta = get_bk_angle(X, Kcon, Ucov);	/* angle between k & b */
  if (theta <= 0. || theta >= M_PI) {	/* no emission along field */
    *jnuinv = 0.;
    *knuinv = 0.;
    return;
  }

  nu = get_fluid_nu(Kcon, Ucov);	 /* freq in Hz */

  B = get_model_b(X);		/* field in G */
  Thetae = get_model_thetae(X);	/* temp in e rest-mass units */

  /* assume emission is thermal */
  Bnuinv = Bnu_inv(nu, Thetae);
  *jnuinv = jnu_inv(nu, Thetae, Ne, B, theta);

  if (Bnuinv < SMALL)
    *knuinv = SMALL;
  else
    *knuinv = *jnuinv / Bnuinv;

  if (isnan(*jnuinv) || isnan(*knuinv)) {
    fprintf(stderr, "\nisnan get_jkinv\n");
    fprintf(stderr, ">> %g %g %g %g %g %g %g %g\n", *jnuinv, *knuinv,
        Ne, theta, nu, B, Thetae, Bnuinv);
  }

  if (counterjet == 1) { // Emission from X[2] > 0.5 only
    if (X[2] < 0.5) {
      *jnuinv = 0.;
    }
  } else if (counterjet == 2) { // Emission from X[2] < 0.5 only
    if (X[2] > 0.5) {
      *jnuinv = 0.;
    }
  }

  return;
}

