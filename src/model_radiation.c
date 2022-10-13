
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
#include "grid.h"
#include "debug_tools.h"
#include "par.h"
#include "decs.h"

// Symphony
#include "fits.h"
#include "params.h"
#include "bremss_fits.h"

#include <omp.h>
#include <gsl/gsl_sf_bessel.h>

// Emission prescriptions, described below
#define E_THERMAL        1
#define E_KAPPA          2
#define E_POWERLAW       3
#define E_DEXTER_THERMAL 4
#define E_LEUNG          5
#define E_CUSTOM        10
// Rotation
#define ROT_OLD         11
#define ROT_PIECEWISE   12
#define ROT_SHCHERBAKOV 13
// Debugging/internal
#define E_UNPOL         15

#define RAD_OVERFLOW 1e100

// Local functions for declaring different kappa/powerlaw distributions
void get_model_powerlaw_vals(double X[NDIM], double *p, double *n,
                             double *gamma_min, double *gamma_max, double *gamma_cut);
void get_model_kappa(double X[NDIM], double *kappa, double *kappa_width);

/* Get coeffs from a specific distribution */
void jar_calc_dist(int dist, int pol, double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV);

/**
 * Optionally load radiation model parameters
 */
static double model_kappa = 3.5;
static double powerlaw_gamma_cut = 1e10;
static double powerlaw_gamma_min = 1e2;
static double powerlaw_gamma_max = 1e5;
static double powerlaw_p = 3.25;
static double powerlaw_eta = 0.02;
static int variable_kappa = 0;
static double variable_kappa_min = 3.1;
static double variable_kappa_interp_start = 1e20;
static double variable_kappa_max = 7.0;
static double max_pol_frac_e = 0.99;
static double max_pol_frac_a = 0.99;
static int do_bremss = 0;
static int bremss_type = 2;

void try_set_radiation_parameter(const char *word, const char *value)
{
  set_by_word_val(word, value, "kappa", &model_kappa, TYPE_DBL);
  set_by_word_val(word, value, "variable_kappa", &variable_kappa, TYPE_INT);
  set_by_word_val(word, value, "variable_kappa_min", &variable_kappa_min, TYPE_DBL);
  set_by_word_val(word, value, "variable_kappa_interp_start", &variable_kappa_interp_start, TYPE_DBL);
  set_by_word_val(word, value, "variable_kappa_max", &variable_kappa_max, TYPE_DBL);

  set_by_word_val(word, value, "powerlaw_gamma_cut", &powerlaw_gamma_cut, TYPE_DBL);
  set_by_word_val(word, value, "powerlaw_gamma_min", &powerlaw_gamma_min, TYPE_DBL);
  set_by_word_val(word, value, "powerlaw_gamma_max", &powerlaw_gamma_max, TYPE_DBL);
  set_by_word_val(word, value, "powerlaw_p", &powerlaw_p, TYPE_DBL);
  set_by_word_val(word, value, "powerlaw_eta", &powerlaw_eta, TYPE_DBL);

  set_by_word_val(word, value, "bremss", &do_bremss, TYPE_INT);
  set_by_word_val(word, value, "bremss_type", &bremss_type, TYPE_INT);

  set_by_word_val(word, value, "max_pol_frac_e", &powerlaw_p, TYPE_DBL);
  set_by_word_val(word, value, "max_pol_frac_a", &powerlaw_p, TYPE_DBL);
}

/**
 * Get polarized emission, absorption, and rotation coefficients
 * 
 * This is a wrapper to jar_calc_dist, see implementation there
 * Also checks for NaN coefficients when built with DEBUG, for quicker debugging
 */
void jar_calc(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV, Params *params)
{
  jar_calc_dist(params->emission_type, 1, X, Kcon, jI, jQ, jU, jV, aI, aQ, aU, aV, rQ, rU, rV);

  // This wrapper can be used to call jar_calc_dist differently in e.g. funnel vs jet, or
  // depending on local fluid parameters, whatever
}

/**
 * Get the emission and absorption coefficients 
 * This is a wrapper to jar_calc_dist, see implementation there
 */
void get_jkinv(double X[NDIM], double Kcon[NDIM], double *jI, double *aI, Params *params)
{
  double jQ, jU, jV, aQ, aU, aV, rQ, rU, rV;
  jar_calc_dist(params->emission_type, 0, X, Kcon, jI, &jQ, &jU, &jV, aI, &aQ, &aU, &aV, &rQ, &rU, &rV);
}

/**
 * Get the invariant plasma emissivities, absorptivities, rotativities in tetrad frame
 * Calls the appropriate fitting functions from the code in symphony/,
 * and ensures that the results obey basic consistency
 *
 * The dist argument controls fitting functions/distribution:
 * 1. Thermal (Pandya+ 2016)
 * 2. Kappa (Pandya+ 2016, no rotativities)
 * 3. Power-law (Pandya+ 2016, rhoQ/rhoV Marszewski+ 2021)
 * 4. Thermal (Dexter 2016)
 * Thermal rhoQ taken from Dexter 2016, rhoV from Shcherbakov 2008 (rhoU == 0)
 *
 * To be implemented?
 * 5. Power-law (Dexter 2016)
 * 
 * Testing distributions:
 * 10. Pass through model-determined values -- used in most analytic models
 * 
 * 
 */
void jar_calc_dist(int dist, int pol, double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV)
{
  // Don't emit where there are no electrons
  // This was also used as shorthand in some models to cut off emission
  // Please use radiating_region instead for model applicability cutoffs,
  // or see integrate_emission for zeroing just jN
  double Ne = get_model_ne(X);
  if (Ne <= 0.) {
    *jI = 0.0; *jQ = 0.0; *jU = 0.0; *jV = 0.0;
    *aI = 0.0; *aQ = 0.0; *aU = 0.0; *aV = 0.0;
    *rQ = 0; *rU = 0; *rV = 0;
    return;
  }

  // Call through to the model if it's responsible for this job
  if (dist == E_CUSTOM) {
    get_model_jar(X, Kcon, jI, jQ, jU, jV, aI, aQ, aU, aV, rQ, rU, rV);
    return;
  }

  struct parameters paramsM;
  setConstParams(&paramsM);
  // Set the parameters common to all distributions...
  double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
  get_model_fourv(X, Kcon, Ucon, Ucov, Bcon, Bcov);
  double nu    = get_fluid_nu(Kcon, Ucov);
  double nusq = nu*nu;
  double theta = get_bk_angle(X, Kcon, Ucov, Bcon, Bcov);
  paramsM.electron_density = Ne;
  paramsM.nu               = nu;
  paramsM.observer_angle   = theta;
  paramsM.magnetic_field   = get_model_b(X);

  // ...and then the specific distribution and its parameters.
  switch (dist) {
  case E_KAPPA: // Kappa fits (Pandya + Marszewski)
    paramsM.distribution = paramsM.KAPPA_DIST;
    // Fall back to Dexter starting at kappa > kappa_interp_start, completely at kappa_max
    paramsM.dexter_fit = 1; // This only affects choice of thermal fallback
    paramsM.kappa_interp_begin = fmin(variable_kappa_interp_start, variable_kappa_max);
    paramsM.kappa_interp_end = variable_kappa_max;
    paramsM.theta_e = get_model_thetae(X);
    get_model_kappa(X, &(paramsM.kappa), &(paramsM.kappa_width));
    break;
  case E_POWERLAW: // Powerlaw fits (Pandya, no rotativities!)
    paramsM.distribution = paramsM.POWER_LAW;
    // NOTE WE REPLACE Ne!!
    get_model_powerlaw_vals(X, &(paramsM.power_law_p), &(paramsM.electron_density),
                            &(paramsM.gamma_min), &(paramsM.gamma_max), &(paramsM.gamma_cutoff));
    break;
  case E_THERMAL: // Pandya thermal fits
    paramsM.distribution = paramsM.MAXWELL_JUETTNER;
    paramsM.theta_e = get_model_thetae(X);
    break;
  case E_DEXTER_THERMAL: // Dexter thermal fits (default)
  default:
    paramsM.dexter_fit = 1;
    paramsM.distribution = paramsM.MAXWELL_JUETTNER;
    paramsM.theta_e = get_model_thetae(X);
    break;
  }

  // First, enforce no emission/absorption along field lines,
  // but allow Faraday rotation in polarized context
  // TODO this skips any rho_V NaN/other checks
  if (theta <= 0.0 || theta >= M_PI) {
    *jI = 0.0; *jQ = 0.0; *jU = 0.0; *jV = 0.0;
    *aI = 0.0; *aQ = 0.0; *aU = 0.0; *aV = 0.0;
    *rQ = 0.0; *rU = 0.0;
    if (pol && !(dist == E_UNPOL)) {
      *rV = rho_nu_fit(&paramsM, paramsM.STOKES_V) * nu;
    } else {
      *rV = 0.0;
    }
    return;
  }

  // Then, if performing unpolarized transport, calculate only what we need
  if (!pol || dist == E_UNPOL) {
    if (paramsM.distribution == paramsM.MAXWELL_JUETTNER) {
      paramsM.dexter_fit = 2; // Signal symphony fits to use Leung+
      *jI = j_nu_fit(&paramsM, paramsM.STOKES_I);
      if(do_bremss) *jI += bremss_I(&paramsM, bremss_type);
      *jI /= nusq; // Avoids loss of precision in small numbers
      double Bnuinv = Bnu_inv(nu, paramsM.theta_e); // Planck function
      if (Bnuinv > 0) {
        *aI = *jI / Bnuinv;
      } else {
        *aI = 0;
      }
    } else {
      paramsM.dexter_fit = 2; // Signal symphony fits to use Leung+ as fallback
      *jI = j_nu_fit(&paramsM, paramsM.STOKES_I) / nusq;
      *aI = alpha_nu_fit(&paramsM, paramsM.STOKES_I) * nu;
    }
  } else { // Finally, calculate all coefficients normally
    // EMISSIVITIES
    *jI = j_nu_fit(&paramsM, paramsM.STOKES_I);
    // Bremsstrahlung is computed only for thermal;
    // silently drop Bremss+other dists, as absorptivities will be wrong
    if(paramsM.distribution == paramsM.MAXWELL_JUETTNER && do_bremss)
      *jI += bremss_I(&paramsM, bremss_type);
    *jI /= nusq; // Avoids loss of precision in small numbers

    *jQ = -j_nu_fit(&paramsM, paramsM.STOKES_Q) / nusq;
    *jU = -j_nu_fit(&paramsM, paramsM.STOKES_U) / nusq;
    *jV = j_nu_fit(&paramsM, paramsM.STOKES_V) / nusq;
    // Check basic relationships
    double jP = sqrt(*jQ * *jQ + *jU * *jU + *jV * *jV);
    if (*jI  < jP/max_pol_frac_e) {
      // Transport does not like 100% polarization...
      double pol_frac_e = *jI / jP * max_pol_frac_e;
      *jQ *= pol_frac_e;
      *jU *= pol_frac_e;
      *jV *= pol_frac_e;
    }

    // ABSORPTIVITIES
    if (paramsM.distribution == paramsM.MAXWELL_JUETTNER) { // Thermal distributions
      // Get absorptivities via Kirchoff's law
      // Already invariant, guaranteed to respect aI > aP
      // Faster than calling Symphony code since we know jS, Bnu
      double Bnuinv = Bnu_inv(nu, paramsM.theta_e); // Planck function
      if (Bnuinv > 0) {
        *aI = *jI / Bnuinv;
        *aQ = *jQ / Bnuinv;
        *aU = *jU / Bnuinv;
        *aV = *jV / Bnuinv;
      } else {
        *aI = 0.;
        *aQ = 0.;
        *aU = 0.;
        *aV = 0.;
      }
    } else {
      *aI = alpha_nu_fit(&paramsM, paramsM.STOKES_I) * nu;
      // Note Bremss emission is available for thermal dists *only*
      *aQ = -alpha_nu_fit(&paramsM, paramsM.STOKES_Q) * nu;
      *aU = -alpha_nu_fit(&paramsM, paramsM.STOKES_U) * nu;
      *aV = alpha_nu_fit(&paramsM, paramsM.STOKES_V) * nu;

      // Check basic relationships
      double aP = sqrt(*aQ * *aQ + *aU * *aU + *aV * *aV);
      if (*aI < aP/max_pol_frac_a) {
        // Transport does not like 100% polarization...
        double pol_frac_a = *aI / aP * max_pol_frac_a;
        *aQ *= pol_frac_a;
        *aU *= pol_frac_a;
        *aV *= pol_frac_a;
      }
    }

    // ROTATIVITIES
    paramsM.dexter_fit = 0;  // Don't use the Dexter rhoV, as it's unstable at low temperature
    *rQ = rho_nu_fit(&paramsM, paramsM.STOKES_Q) * nu;
    *rU = rho_nu_fit(&paramsM, paramsM.STOKES_U) * nu;
    *rV = rho_nu_fit(&paramsM, paramsM.STOKES_V) * nu;
  }

#if DEBUG
  // Check for NaN coefficients
  if (isnan(*jI) || *jI > RAD_OVERFLOW || *jI < -RAD_OVERFLOW ||
      isnan(*jQ) || *jQ > RAD_OVERFLOW || *jQ < -RAD_OVERFLOW ||
      isnan(*jU) || *jU > RAD_OVERFLOW || *jU < -RAD_OVERFLOW ||
      isnan(*jV) || *jV > RAD_OVERFLOW || *jV < -RAD_OVERFLOW ||
      isnan(*aI) || *aI > RAD_OVERFLOW || *aI < -RAD_OVERFLOW ||
      isnan(*aQ) || *aQ > RAD_OVERFLOW || *aQ < -RAD_OVERFLOW ||
      isnan(*aU) || *aU > RAD_OVERFLOW || *aU < -RAD_OVERFLOW ||
      isnan(*aV) || *aV > RAD_OVERFLOW || *aV < -RAD_OVERFLOW ||
      isnan(*rQ) || *rQ > RAD_OVERFLOW || *rQ < -RAD_OVERFLOW ||
      isnan(*rU) || *rU > RAD_OVERFLOW || *rU < -RAD_OVERFLOW ||
      isnan(*rV) || *rV > RAD_OVERFLOW || *rV < -RAD_OVERFLOW) {
#pragma omp critical
    {
      fprintf(stderr, "\nNAN in emissivities!\n");
      fprintf(stderr, "j = %g %g %g %g alpha = %g %g %g %g rho = %g %g %g\n", *jI, *jQ, *jU, *jV, *aI, *aQ, *aU, *aV, *rQ, *rU, *rV);
      fprintf(stderr, "nu = %g Ne = %g Thetae = %g B = %g theta = %g kappa = %g kappa_width = %g\n",
              paramsM.nu, paramsM.electron_density, paramsM.theta_e, paramsM.magnetic_field, paramsM.observer_angle,
              paramsM.kappa, paramsM.kappa_width);
      // Powerlaw?
      exit(-1);
    }
  }
#endif
}

// SUPPORTING FUNCTIONS

/*
 * Get values for the powerlaw distribution:
 * power & cutoffs from parameters, plus number density of nonthermals
 */
void get_model_powerlaw_vals(double X[NDIM], double *p, double *n,
                             double *gamma_min, double *gamma_max, double *gamma_cut)
{
  *gamma_min = powerlaw_gamma_min;
  *gamma_max = powerlaw_gamma_max;
  *gamma_cut = powerlaw_gamma_cut;
  *p = powerlaw_p;

  double b = get_model_b(X);
  double u_nth = powerlaw_eta*b*b/2;
  // Number density of nonthermals
  *n = u_nth * (*p - 2)/(*p - 1) * 1/(ME * CL*CL * *gamma_min);
}

/*
 * Get values kappa, kappa_width for the kappa distribution emissivities
 * Optionally variable over the domain based on sigma and beta
 */
void get_model_kappa(double X[NDIM], double *kappa, double *kappa_width) {

  if (variable_kappa) {
    double sigma = get_model_sigma(X);
    double beta = get_model_beta(X);

    *kappa = 2.8 + 0.7/sqrt(sigma) + 3.7 * (1.0/pow(sigma, 0.19)) * tanh(23.4 * pow(sigma, 0.26) * beta);

    if (isnan(*kappa)) {
      double B = get_model_b(X);
      fprintf(stderr, "kappa, sigma, beta, B: %g %g %g %g \n", *kappa, sigma, beta, B);
    }
    // Lower limit for kappa.  Upper limit switches away from kappa fit entirely -> thermal
    if (*kappa < variable_kappa_min) *kappa = variable_kappa_min;

    //fprintf(stderr, "sigma, beta -> kappa %g %g -> %g\n", sigma, beta, *kappa);
  } else {
    *kappa = model_kappa;
  }

  double Thetae = get_model_thetae(X);
  *kappa_width = (*kappa - 3.) / *kappa * Thetae;
#if DEBUG
  if (isnan(*kappa_width) || isnan(*kappa)) {
    fprintf(stderr, "NaN kappa val! kappa, kappa_width, Thetae: %g %g %g\n", *kappa, *kappa_width, Thetae);
  }
#endif
}
