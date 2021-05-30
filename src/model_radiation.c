
/**********************************************************/
/*** all you need to make a polarized radiative transfer **/
/***** used in ipole to evolve complex tensor N ***********/
/***** along with standard evolution for I scalar *********/
/**********************************************************/
/**** written by Monika Moscibrodzka on 09 July 2014 ******/
/************ @ Eindhoven Airport *************************/
/************  last update: 9 May 2017   ******************/
/****************and then rewritten by C.Gammie ***********/
/**********************************************************/

#include "model_radiation.h"

#include "model.h"
#include "radiation.h"
#include "coordinates.h"
#include "geometry.h"
#include "debug_tools.h"
#include "par.h"
#include "decs.h"

// Symphony
#include "fits.h"
#include "params.h"

#include <omp.h>
#include <gsl/gsl_sf_bessel.h>

// Emission prescriptions, described below
#define E_THERMAL        1
#define E_KAPPA          2
#define E_POWERLAW       3
#define E_DEXTER_THERMAL 4
#define E_CUSTOM        10
// Rotation 
#define ROT_OLD         11
#define ROT_PIECEWISE   12
#define ROT_SHCHERBAKOV 13
// Debugging
#define E_UNPOL         15

// Declarations of local fitting functions, for Dexter fits and old rotativities
void dexter_j_fit_thermal(double Ne, double nu, double Thetae, double B, double theta,
                    double *jI, double *jQ, double *jU, double *jV);
void shcherbakov_rho_fit(double Ne, double nu, double Thetae, double B, double theta,
                    double *rQ, double *rU, double *rV);
void piecewise_rho_fit(double Ne, double nu, double Thetae, double B, double theta,
                    double *rQ, double *rU, double *rV);
void old_rho_fit(double Ne, double nu, double Thetae, double B, double theta,
                    double *rQ, double *rU, double *rV);

// Thermal plasma emissivity, absorptivity and Faraday conversion and rotation
double g(double Xe);
double h(double Xe);
double Je(double Xe);
double jffunc(double Xe);
double I_I(double x);
double I_Q(double x);
double I_V(double x);
double besselk_asym(int n, double x);

/* Get coeffs from a specific distribution */
void jar_calc_dist(int dist, double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV);

/**
 * Optionally load radiation model parameters
 */
static double model_kappa = 0;
static double powerlaw_gamma_cut = 1e10;
static double powerlaw_gamma_min = 1e2;
static double powerlaw_gamma_max = 1e5;
static double powerlaw_p = 3.25;
static double max_pol_frac_e = 0.99;
static double max_pol_frac_a = 0.99;

void try_set_radiation_parameter(const char *word, const char *value)
{
  set_by_word_val(word, value, "kappa", &model_kappa, TYPE_DBL);

  set_by_word_val(word, value, "powerlaw_gamma_cut", &powerlaw_gamma_cut, TYPE_DBL);
  set_by_word_val(word, value, "powerlaw_gamma_min", &powerlaw_gamma_min, TYPE_DBL);
  set_by_word_val(word, value, "powerlaw_gamma_max", &powerlaw_gamma_max, TYPE_DBL);
  set_by_word_val(word, value, "powerlaw_p", &powerlaw_p, TYPE_DBL);

  set_by_word_val(word, value, "max_pol_frac_e", &powerlaw_p, TYPE_DBL);
  set_by_word_val(word, value, "max_pol_frac_a", &powerlaw_p, TYPE_DBL);
}

/**
 * Wrapper to call different distributions at different places in the simulation domain
 * See jar_calc_dist for distributions
 * TODO put the general emission zero criteria here
 */
void jar_calc(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV, Params *params)
{
#if INTEGRATOR_TEST
  jar_calc_dist(E_CUSTOM, X, Kcon, jI, jQ, jU, jV, aI, aQ, aU, aV, rQ, rU, rV);
#else
  if (params->emission_type == E_UNPOL) {
    get_jkinv(X, Kcon, jI, aI, params);
    *jQ = 0.0; *jU = 0.0; *jV = 0.0;
    *aQ = 0.0; *aU = 0.0; *aV = 0.0;
    *rQ = 0; *rU = 0; *rV = 0;
  } else {
    jar_calc_dist(params->emission_type, X, Kcon, jI, jQ, jU, jV, aI, aQ, aU, aV, rQ, rU, rV);
  }
#endif

  // Zero emission to isolate the jet/counterjet portion
  if (params->isolate_counterjet == 1) { // Allow emission from X[2] > midplane only
    if (X[2] < (cstopx[2] - cstartx[2]) / 2) {
      *jI = *jQ = *jU = *jV = 0.;
    }
  } else if (params->isolate_counterjet == 2) { // from X[2] < midplane only
    if (X[2] > (cstopx[2] - cstartx[2]) / 2) {
      *jI = *jQ = *jU = *jV = 0.;
    }
  }
}

/**
 * Get the invariant plasma emissivities, absorptivities, rotativities in tetrad frame
 * This is a wrapper to the fitting functions in this file and in the emissivity/ dir
 *
 * The dist argument controls fitting functions/distribution:
 * 1. Thermal (Pandya+ 2016)
 * 2. Kappa (Pandya+ 2016)
 * 3. Power-law (Pandya+ 2016)
 * 4. Thermal (Dexter 2016)
 * Rotativity rhoQ always from Dexter 2016, rhoV from Shcherbakov 2008 (rhoU == 0)
 *
 * To be implemented?
 * 5. Power-law (Dexter 2016)
 * 6. Thermal (Revised based on Pandya+ 2016 to better match Leung 2011)
 * 
 * Testing distributions:
 * 10. Pass through model-determined values. Currently used for constant-coefficient testing.
 * 11. Emulate original ipole, with Dexter 2016 emissivities and rotativities based on Bessel fits
 * 12. Emulate the old ipole temporary fix, with Dexter 2016 emissivities and rotativities patched into the constant limit
 * 
 */
void jar_calc_dist(int dist, double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV)
{
  // Four-vectors needed for most calculations
  double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
  // Common parameters
  double nu = 0, Ne = 0, B = 0, theta = 0;
  // Other parameters are filled as needed
  // Thermal:
  double Thetae = 0;
  // Powerlaw:
  double powerlaw_p = 0, gamma_min = 0, gamma_max = 0, gamma_cut = 0;
  // Kappa:
  double kappa = 0, kappa_width = 0;
  // Symphony parameters struct, not to be confused with ipole Params
  struct parameters paramsM; int fit = 0;

  // Ignore everything if this isn't our job
  if (dist == E_CUSTOM) {
    get_model_jar(X, Kcon, jI, jQ, jU, jV, aI, aQ, aU, aV, rQ, rU, rV);
    return;
  }
  // Also ignore everything if we're not in an emitting region
  // (or there are actually no e- to emit)
  Ne = get_model_ne(X);
  if (Ne <= 0.) {
    *jI = 0.0; *jQ = 0.0; *jU = 0.0; *jV = 0.0;
    *aI = 0.0; *aQ = 0.0; *aU = 0.0; *aV = 0.0;
    *rQ = 0; *rU = 0; *rV = 0;
    return;
  }

  // If we must do work, grab the 4-vectors...
  get_model_fourv(X, Kcon, Ucon, Ucov, Bcon, Bcov);
#if DEBUG
  if (isnan(Ucov[0])) {
    void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]);
    int i,j,k;
    double del[4];
    Xtoijk(X, &i,&j,&k, del);
    fprintf(stderr, "UCOV[0] (%d,%d,%d) is nan! thread = %i\n", i,j,k, omp_get_thread_num());
    print_vector("Ucon", Ucon);
    print_vector("Ucov", Ucov);
    fprintf(stderr, "X[] = %e %e %e %e\n", X[0],X[1],X[2],X[3]);
    fprintf(stderr, "K[] = %e %e %e %e\n", Kcon[0],Kcon[1],Kcon[2],Kcon[3]);
    fprintf(stderr, "Ne = %e\n", Ne);
  }
#endif

  // ...then the parameters needed for all distributions...
  theta = get_bk_angle(X, Kcon, Ucov, Bcon, Bcov); // angle between k & b
  nu = get_fluid_nu(Kcon, Ucov); // freqcgs in Hz
  B = get_model_b(X); // field in G

  //...and the ones for the specific distribution we'll be evaluating.
  switch (dist) {
  case E_THERMAL:
    setConstParams(&paramsM);
    fit = paramsM.MAXWELL_JUETTNER;
    Thetae = get_model_thetae(X);
    break;
  case E_KAPPA:
    setConstParams(&paramsM);
    fit = paramsM.KAPPA_DIST;
    Thetae = get_model_thetae(X);
    kappa = model_kappa;
    kappa_width = (kappa - 3.) / kappa * Thetae;
    break;
  case E_POWERLAW:
    setConstParams(&paramsM);
    fit = paramsM.POWER_LAW;
    // NOTE WE REPLACE Ne!!
    get_model_powerlaw_vals(X, &powerlaw_p, &Ne, &gamma_min, &gamma_max, &gamma_cut);
    break;
  case E_DEXTER_THERMAL:
    Thetae = get_model_thetae(X);
  }

  // EMISSIVITIES
  // Avoid issues directly along field lines
  if (theta <= 0 || theta >= M_PI) {
    *jI = 0.0; *jQ = 0.0; *jU = 0.0; *jV = 0.0;
    *aI = 0.0; *aQ = 0.0; *aU = 0.0; *aV = 0.0;
  } else {
    // EMISSIVITIES
    if (dist == E_DEXTER_THERMAL || dist > 10) {
      // Use fits from Dexter when called for or emulating old behavior
      dexter_j_fit_thermal(Ne, nu, Thetae, B, theta, jI, jQ, jU, jV);
    } else {
      // Call into bundled Symphony code
      // Symphony uses an... interesting coordinate system. Correct it.
      // TODO fix Symphony, jV fit is bad
      *jI = j_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_I, Thetae, powerlaw_p, gamma_min, gamma_max, gamma_cut, kappa, kappa_width);
      *jQ = -j_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_Q, Thetae, powerlaw_p, gamma_min, gamma_max, gamma_cut, kappa, kappa_width);
      *jU = -j_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_U, Thetae, powerlaw_p, gamma_min, gamma_max, gamma_cut, kappa, kappa_width);
      *jV = j_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_V, Thetae, powerlaw_p, gamma_min, gamma_max, gamma_cut, kappa, kappa_width);
    }
    // Check basic relationships
    if (*jI < 0) {
      fprintf(stderr, "Negative total emissivity! Exiting!\n");
      exit(-1);
    }
    if (*jI * *jI < *jQ * *jQ + *jU * *jU + *jV * *jV) {
      // Transport does not like 100% polarization...
      double pol_frac_e = *jI / jP * max_pol_frac_e;
      *jQ *= pol_frac_e;
      *jU *= pol_frac_e;
      *jV *= pol_frac_e;
#if DEBUG
      fprintf(stderr, "Polarized emissivities too large:\n %g vs %g, corrected by %g\n",
      sqrt(*jQ * *jQ + *jU * *jU + *jV * *jV), *jI, pol_frac_e);
#endif
    }

    // Make invariant
    double nusq = nu*nu;
    *jI /= nusq;
    *jQ /= nusq;
    *jU /= nusq;
    *jV /= nusq;

    // ABSORPTIVITIES
    if (dist == E_THERMAL || dist == E_DEXTER_THERMAL || dist > 10) { // Thermal distributions
      // Get absorptivities via Kirchoff's law
      // Already invariant, guaranteed to respect aI > aP
      double Bnuinv = Bnu_inv(nu, Thetae); // Planck function
      *aI = *jI / Bnuinv;
      *aQ = *jQ / Bnuinv;
      *aU = *jU / Bnuinv;
      *aV = *jV / Bnuinv;
    } else {
      *aI = alpha_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_I, Thetae, powerlaw_p, gamma_min, gamma_max, gamma_cut, kappa, kappa_width);
      *aQ = -alpha_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_Q, Thetae, powerlaw_p, gamma_min, gamma_max, gamma_cut, kappa, kappa_width);
      *aU = -alpha_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_U, Thetae, powerlaw_p, gamma_min, gamma_max, gamma_cut, kappa, kappa_width);
      *aV = alpha_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_V, Thetae, powerlaw_p, gamma_min, gamma_max, gamma_cut, kappa, kappa_width);

      // Check basic relationships
      if (*aI < 0) {
        fprintf(stderr, "Negative total absorptivity! Exiting!\n");
        exit(-1);
      }
      if (*aI * *aI < *aQ * *aQ + *aU * *aU + *aV * *aV) {
        // Transport does not like 100% polarization...
        double pol_frac_a = *aI / aP * max_pol_frac_a;
        *aQ *= pol_frac_a;
        *aU *= pol_frac_a;
        *aV *= pol_frac_a;
#if DEBUG
        fprintf(stderr, "Polarized absorptivities too large:\n %g vs %g, corrected by %g\n",
        sqrt(*aQ * *aQ + *aU * *aU + *aV * *aV), *aI, pol_frac_a);
#endif
      }

      // Make invariant
      *aI *= nu;
      *aQ *= nu;
      *aU *= nu;
      *aV *= nu;
    }
  }

  // ROTATIVITIES
  // Fill them by defualt since we will need rV in any case
  if (dist == ROT_PIECEWISE) { // Old piecewise distribution, for compatibility
    piecewise_rho_fit(Ne, nu, Thetae, B, theta, rQ, rU, rV);
  } else if (dist == ROT_OLD) { // Old incorrect distribution, for compatibility
    old_rho_fit(Ne, nu, Thetae, B, theta, rQ, rU, rV);
  } else if (dist == E_DEXTER_THERMAL) { // Dexter rQ, Shcherbakov rV
    // TODO All-Dexter default option w/Taylor series, make this compat 
    shcherbakov_rho_fit(Ne, nu, Thetae, B, theta, rQ, rU, rV);
  } else { // TODO Fix Symphony, these are currently equal to ROT_OLD
    *rQ = rho_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_Q, Thetae, powerlaw_p, gamma_min, gamma_max, gamma_cut, kappa, kappa_width);
    *rU = rho_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_U, Thetae, powerlaw_p, gamma_min, gamma_max, gamma_cut, kappa, kappa_width);
    *rV = rho_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_V, Thetae, powerlaw_p, gamma_min, gamma_max, gamma_cut, kappa, kappa_width);
  }
  // Make invariant
  *rQ *= nu;
  *rU *= nu;
  *rV *= nu;
  
  // *Then* handle field lines, to leave rV intact
  if (theta <= 0 || theta >= M_PI) {
    *rQ = 0;
    *rU = 0;
  }

#if DEBUG
  // Spot check for NaN coefficients
  if (isnan(*rV) || *rV > 1.e100 || *rV < -1.e100) {
    fprintf(stderr, "\nNAN RV! rV = %e nu = %e Ne = %e Thetae = %e\n", *rV, nu, Ne, Thetae);
  }
  if (isnan(*jV) || *jV > 1.e100 || *jV < -1.e100) {
    fprintf(stderr, "\nNAN jV! jV = %e nu = %e Ne = %e Thetae = %e B = %e theta = %e\n", *jV, nu, Ne, Thetae, B, theta);
  }
#endif
}

/*
 * emissivity functions and functions used for Faraday conversion and rotation
 * from J. Dexter PhD thesis (checked with Leung harmony program, and Huang & Shcherbakov 2011
 * Also see Dexter 2016 Appendix A1
 */
void dexter_j_fit_thermal(double Ne, double nu, double Thetae, double B, double theta,
                    double *jI, double *jQ, double *jU, double *jV)
{
  // Synchrotron emissivity
  double nus = 3.0 * EE * B * sin(theta) / 4.0 / M_PI / ME / CL * Thetae * Thetae + 1.0;
  double x = nu / nus;

  *jI = Ne * EE * EE * nu / 2. / sqrt(3) / CL / Thetae / Thetae * I_I(x); // [g/s^2/cm = ergs/s/cm^3]
  *jQ = Ne * EE * EE * nu / 2. / sqrt(3) / CL / Thetae / Thetae * I_Q(x);
  *jU = 0.0;    // convention; depends on tetrad
  *jV = 2. * Ne * EE * EE * nu / tan(theta) / 3. / sqrt(3) / CL / Thetae / Thetae / Thetae * I_V(x);
}

void shcherbakov_rho_fit(double Ne, double nu, double Thetae, double B, double theta,
                    double *rQ, double *rU, double *rV)
{
  double Thetaer = 1. / Thetae;

  double omega0 = EE * B / ME / CL;
  double wp2 = 4. * M_PI * Ne * EE * EE / ME;

  // Faraday rotativities for thermal plasma
  double Xe = Thetae * sqrt(sqrt(2) * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));

  // These are the Dexter (2016) fit actually
  *rQ = 2. * M_PI * nu / 2. / CL * wp2 * omega0 * omega0 / pow(2 * M_PI * nu, 4) *
        jffunc(Xe) * (gsl_sf_bessel_Kn(1, Thetaer) / gsl_sf_bessel_Kn(2, Thetaer) +
                      6. * Thetae) * sin(theta) * sin(theta);
  *rU = 0.0;

  // Shcherbakov fit for rV.  Possibly questionable at very low frequency
  // Note the real bessel functions. Slow?
  *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
        gsl_sf_bessel_Kn(0, Thetaer) / (gsl_sf_bessel_Kn(2, Thetaer)+SMALL) * g(Xe) * cos(theta);
}

void piecewise_rho_fit(double Ne, double nu, double Thetae, double B, double theta,
                    double *rQ, double *rU, double *rV)
{
  double Thetaer = 1. / Thetae;

  double omega0 = EE * B / ME / CL;
  double wp2 = 4. * M_PI * Ne * EE * EE / ME;

  // Faraday rotativities for thermal plasma
  double Xe = Thetae * sqrt(sqrt(2) * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));

  // Approximate bessel functions to match rhoq,v with grtrans
  *rQ = 2. * M_PI * nu / 2. / CL * wp2 * omega0 * omega0 / pow(2 * M_PI * nu, 4) *
    jffunc(Xe) * (besselk_asym(1, Thetaer) / besselk_asym(2, Thetaer) +
        6. * Thetae) * sin(theta) * sin(theta);
  *rU = 0.0;
  // Switch between three different fits for rho_V
  if (Thetae > 3.0) {
    // High temperature: use approximations to bessel
    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
  } else if (0.2 < Thetae && Thetae <= 3.0) {
    // Mid temperature: use real bessel functions (TODO fit?)
    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      (gsl_sf_bessel_Kn(0, Thetaer) - Je(Xe)) / gsl_sf_bessel_Kn(2, Thetaer) * cos(theta);
  } else if (Thetae <= 0.2) {
    // Use the constant low-temperature limit
    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) * cos(theta);
  }
}

void old_rho_fit(double Ne, double nu, double Thetae, double B, double theta,
                    double *rQ, double *rU, double *rV)
{
  double Thetaer = 1. / Thetae;

  double omega0 = EE * B / ME / CL;
  double wp2 = 4. * M_PI * Ne * EE * EE / ME;

  // Faraday rotativities for thermal plasma
  double Xe = Thetae * sqrt(sqrt(2) * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));

  // Approximate bessel functions to match rhoq,v with grtrans
  *rQ = 2. * M_PI * nu / 2. / CL * wp2 * omega0 * omega0 / pow(2 * M_PI * nu, 4) *
    jffunc(Xe) * (besselk_asym(1, Thetaer) / besselk_asym(2, Thetaer) +
        6. * Thetae) * sin(theta) * sin(theta);
  *rU = 0.0;
  // Use approximations to Bessel fns for all space. Dangerous at low temp!
  *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
}

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
  double extraterm =
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

// UNPOLARIZED VERSIONS

/*
 * get the invariant emissivity and opacity at a given position for a given wavevector
 */
void get_jkinv(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv, Params *params)
{
  if (params->emission_type == E_CUSTOM) {
    get_model_jk(X, Kcon, jnuinv, knuinv);
  } else {
    double nu, theta, B, Thetae, Ne, Bnuinv;
    double Ucov[NDIM], Ucon[NDIM], Bcon[NDIM], Bcov[NDIM];

    /* get fluid parameters */
    Ne = get_model_ne(X);	/* check to see if we're outside fluid model */
    if (Ne == 0.) {
      *jnuinv = 0.;
      *knuinv = 0.;
      return;
    }

    /* get covariant four-velocity of fluid for use in get_bk_angle and get_fluid_nu */
    get_model_fourv(X, Kcon, Ucon, Ucov, Bcon, Bcov);
    theta = get_bk_angle(X, Kcon, Ucov, Bcon, Bcov);	/* angle between k & b */
    // No emission along field
    if (theta <= 0. || theta >= M_PI) {	/* no emission along field */
      *jnuinv = 0.;
      *knuinv = 0.;
      return;
    }

    // Only compute these if we must
    B = get_model_b(X);		/* field in G */
    Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
    nu = get_fluid_nu(Kcon, Ucov);	 /* freq in Hz */

    // TODO: could be cleaner here
    struct parameters paramsM;
    int fit;
    double kappa, kappa_width;

    if (params->emission_type == E_KAPPA) {

      setConstParams(&paramsM);
      fit = paramsM.KAPPA_DIST;

      kappa = model_kappa;
      kappa_width = (kappa - 3.) / kappa * Thetae;

      *jnuinv = j_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_I, Thetae, 
                        powerlaw_p, powerlaw_gamma_min, powerlaw_gamma_max, powerlaw_gamma_cut, 
                        kappa, kappa_width) / nu/nu;
      *knuinv = alpha_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_I, Thetae, 
                        powerlaw_p, powerlaw_gamma_min, powerlaw_gamma_max, powerlaw_gamma_cut, 
                        kappa, kappa_width) * nu;

    } else if (params->emission_type == E_POWERLAW) {

      setConstParams(&paramsM);
      fit = paramsM.POWER_LAW;

      get_model_powerlaw_vals(X, &powerlaw_p, &Ne, &powerlaw_gamma_min, &powerlaw_gamma_max, &powerlaw_gamma_cut);

      *jnuinv = j_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_I, Thetae,
                        powerlaw_p, powerlaw_gamma_min, powerlaw_gamma_max, powerlaw_gamma_cut,
                        kappa, kappa_width) / nu/nu;
      *knuinv = alpha_nu_fit(nu, B, Ne, theta, fit, paramsM.STOKES_I, Thetae,
                        powerlaw_p, powerlaw_gamma_min, powerlaw_gamma_max, powerlaw_gamma_cut,
                        kappa, kappa_width) * nu * exp(-nu / 5.e13); // nu_cutoff = 5.e13

    } else {

      /* assume emission is thermal */
      Bnuinv = Bnu_inv(nu, Thetae);
      *jnuinv = jnu_inv(nu, Thetae, Ne, B, theta);

      if (Bnuinv < SMALL) {
        *knuinv = SMALL;
      } else {
        *knuinv = *jnuinv / Bnuinv;
      }

    }

#if DEBUG
    if (isnan(*jnuinv) || isnan(*knuinv)) {
      fprintf(stderr, "\nisnan get_jkinv\n");
      fprintf(stderr, ">> %g %g %g %g %g %g %g %g\n", *jnuinv, *knuinv,
          Ne, theta, nu, B, Thetae, Bnuinv);
    }
#endif

    if (params->isolate_counterjet == 1) { // Emission from X[2] > midplane only
      if (X[2] < (cstopx[2] - cstartx[2]) / 2) {
        *jnuinv = 0.;
      }
    } else if (params->isolate_counterjet == 2) { // Emission from X[2] < midplane only
      if (X[2] > (cstopx[2] - cstartx[2]) / 2) {
        *jnuinv = 0.;
      }
    }
  }
}

/*
 * thermal synchrotron emissivity
 *
 * Interpolates between Petrosian limit and
 * classical thermal synchrotron limit
 * Good for Thetae >~ 1
 * See Leung+ 2011, restated Pandya+ 2016
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

