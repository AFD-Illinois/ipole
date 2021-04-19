#include "model.h"

#include "decs.h"
#include "hdf5_utils.h"

#include "coordinates.h"
#include "geometry.h"
#include "grid.h"
#include "model_radiation.h" // Only for outputting emissivities
#include "par.h"
#include "utils.h"

#include <string.h>
#include <assert.h>

// used by other files
double rmax_geo;
double L_unit;

// model parameters
static double MODEL_R_0 = 100.;
static double MODEL_TAU_0 = 1.e-5;
static double MODEL_THETAE_0 = 10.;
static double MODEL_BETA_0 = 20.;
static double MODEL_MBH = 4.1e6;
static double MODEL_TP_OVER_TE = 3;
static double MODEL_GAM = 13./9;

// TODO: duplicated code!
static double powerlaw_gamma_min = 1e2;
static double powerlaw_gamma_max = 1e5;
static double powerlaw_gamma_cut = 1e10;
static double powerlaw_p = 3.25;

// derived model parameters
double model_Ne_0 = 1.;
double model_B_0 = 1.;

void try_set_model_parameter(const char *word, const char *value)
{
  set_by_word_val(word, value, "Thetae0", &MODEL_THETAE_0, TYPE_DBL);
  set_by_word_val(word, value, "R0", &MODEL_R_0, TYPE_DBL);
  set_by_word_val(word, value, "tau0", &MODEL_TAU_0, TYPE_DBL);
  set_by_word_val(word, value, "beta0", &MODEL_BETA_0, TYPE_DBL);
  set_by_word_val(word, value, "gam", &MODEL_GAM, TYPE_DBL);
  set_by_word_val(word, value, "tp_over_te", &MODEL_TP_OVER_TE, TYPE_DBL);
  set_by_word_val(word, value, "MBH", &MODEL_MBH, TYPE_DBL);

  // TODO: figure out how to make consistent with model_radiation.c
  set_by_word_val(word, value, "powerlaw_gamma_min", &powerlaw_gamma_min, TYPE_DBL);
  set_by_word_val(word, value, "powerlaw_gamma_max", &powerlaw_gamma_max, TYPE_DBL);
  set_by_word_val(word, value, "powerlaw_gamma_cut", &powerlaw_gamma_cut, TYPE_DBL);
  set_by_word_val(word, value, "powerlaw_p", &powerlaw_p, TYPE_DBL);
}

// used in slow light with real data
void update_data_until(double *tA, double *tB, double tgt) { }
void update_data(double *tA, double *tB) { }
double get_dump_t(char *fnam, int dumpidx) { return 0.; }

void init_model(double *tA, double *tB)
{
  // set nice numbers here
  *tA = 0.;
  *tB = 1.;

  // set metric
  use_eKS_internal = 0;
  metric = METRIC_EMINKOWSKI;

  // set units
  double MBH = MODEL_MBH * MSUN;
  L_unit = GNEWT * MBH / (CL * CL);

  // derive model Ne (in cgs) 
  model_Ne_0 = MODEL_TAU_0 / SIGMA_THOMSON / MODEL_R_0 / L_unit;

  // derive model B (in gauss)
  double THETAE_UNIT = 1.;

  // since B = B(pressure), we need to specify the thermodynamics to
  // find pressure = pressure(Thetae)
  double gam = MODEL_GAM;
  double game = 4./3;
  double gamp = 5./3;

  // as implemented in the Illinois suite
  THETAE_UNIT = MP/ME * (game-1.) * (gamp-1.) / ( (gamp-1.) + (game-1)*MODEL_TP_OVER_TE );

  // as implemented in RAPTOR + kmonty
  THETAE_UNIT = MP/ME * (gam-1.) / (1. + MODEL_TP_OVER_TE);

  // now we can find B (again, in gauss)
  model_B_0 = CL * sqrt(8 * M_PI * (gam-1.) * (MP+ME) / MODEL_BETA_0) * sqrt( model_Ne_0 * MODEL_THETAE_0 ) / sqrt( THETAE_UNIT );

  fprintf(stderr, "Running with isothermal sphere model.\n");
  fprintf(stderr, "MBH, L_unit: %g [Msun], %g\n",MBH/MSUN, L_unit);
  fprintf(stderr, "Ne, Thetae, B: %g %g %g\n", model_Ne_0, MODEL_THETAE_0, model_B_0);
  fprintf(stderr, "Rout: %g\n", MODEL_R_0 * L_unit);

  // not really used in ipole other than for coordinates
  N1 = 300;
  N2 = 100;
  N3 = 1;

  // set up coordinates
  startx[0] = 0.;
  startx[1] = 0.;
  startx[2] = 0.;
  startx[3] = 0.;

  dx[0] = 0.1;
  dx[1] = (MODEL_R_0 - 0) / N1;
  dx[2] = M_PI / N2;
  dx[3] = 2. * M_PI / N3;
  
  stopx[0] = 1.;
  stopx[1] = startx[1]+N1*dx[1];
  stopx[2] = startx[2]+N2*dx[2];
  stopx[3] = startx[3]+N3*dx[3];

  rmax_geo = fmax(100., MODEL_R_0);

  fprintf(stderr, "Native coordinate start: %g %g %g stop: %g %g %g\n",
                  startx[1], startx[2], startx[3], stopx[1], stopx[2], stopx[3]);

}

/*
  these supply basic model data to ipole
*/

// Calculate Ucon,Ucov,Bcon,Bcov from primitives at location X using 
// interpolation (on the primitives). This has been all wrapped into
// a single function because some calculations require each other.
void get_model_fourv(double X[NDIM], double Kcon[NDIM],
                     double Ucon[NDIM], double Ucov[NDIM],
                     double Bcon[NDIM], double Bcov[NDIM])
{
  double r, h;
  bl_coord(X, &r, &h);

  double gcov[NDIM][NDIM];
  gcov_func(X, gcov);

  Ucon[0] = 1.;
  Ucon[1] = 0.;
  Ucon[2] = 0.;
  Ucon[3] = 0.;

  Bcon[0] = 0.;
  Bcon[1] = model_B_0 * cos(h) / r;
  Bcon[2] = - model_B_0 * sin(h) / (r + 1.e-8);
  Bcon[3] = 0.;

  lower(Ucon, gcov, Ucov);
  lower(Bcon, gcov, Bcov);
}

// used in diagnostics IO. not implemented
void get_model_primitives(double X[NDIM], double *p) { }

double get_model_thetae(double X[NDIM])
{
  if (radiating_region(X) == 0) {
    return 0;
  }

  return MODEL_THETAE_0;
}

//b field strength in Gauss
double get_model_b(double X[NDIM])
{
  if (radiating_region(X) == 0) {
    return 0.;
  }

  return model_B_0;
}

double get_model_ne(double X[NDIM])
{
  if (radiating_region(X) == 0) {
    return 0.;
  }

  return model_Ne_0;
}

void output_hdf5()
{
  hdf5_set_directory("/");

  // TODO maybe output the model parameters here

  hdf5_set_directory("/header/");
  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  hdf5_write_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);

  hdf5_set_directory("/");
}

int radiating_region(double X[NDIM])
{
  double r, h;
  bl_coord(X, &r, &h);
  return (r<MODEL_R_0) ? 1 : 0;
}

void get_model_powerlaw_vals(double X[NDIM], double *p, double *n,
                             double *gamma_min, double *gamma_max, double *gamma_cut)
{ 
  // TODO figure out how to make this all consistent
  //assert(1 == 0); 

  *gamma_min = powerlaw_gamma_min;
  *gamma_max = powerlaw_gamma_max;
  *gamma_cut = powerlaw_gamma_cut;
  *p = powerlaw_p;
}

// In case we want to mess with emissivities directly
void get_model_jar(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV) {return;}
void get_model_jk(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv) {return;}

