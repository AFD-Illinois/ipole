// Template for analytic emission or electron spatial distributions
// Fill get_model_jar and get_model_fourv and you're off!

#include "model.h"

#include "decs.h"
#include "coordinates.h"
#include "geometry.h"
#include "radiation.h"
#include "par.h"

// Globals we're in charge of
double M_unit;
double L_unit;
double T_unit;
double RHO_unit;
double U_unit;
double B_unit;
double Te_unit;

// Model parameters: public
double rmax_geo;
int counterjet = 0;
// Model parameters: private
static double MBH_solar = 4.e6;
static double MBH;
static int model;

// e.g. parameterization from GRRT paper
double A, alpha, height, l0, freqcgs;

// Forward declarations for non-public functions
void set_units();

void try_set_model_parameter(const char *word, const char *value)
{
  // Test the given word against our parameters' names,
  // and if it matches set the corresponding global
  set_by_word_val(word, value, "MBH", &MBH_solar, TYPE_DBL);
  set_by_word_val(word, value, "counterjet", &counterjet, TYPE_DBL);

  set_by_word_val(word, value, "model", &model, TYPE_INT);
  // Normal ipole pulls this, but we also need it for the GRRT problems
  // and this is easier than grabbing it from the 'params' struct
  set_by_word_val(word, value, "freqcgs", &freqcgs, TYPE_DBL);
}

void init_model(double *tA, double *tB)
{
  // Set all the geometry globals we need
  // TODO do this in geometry?  Deal with model/geom interface...
  use_eKS_internal = 1;
  metric = 0; // Doesn't matter due to above
  hslope = 1.0;

  if (model == 1) {
    A = 0;
    alpha = -3;
    height = 0;
    l0 = 0;
    a = 0.9;
  } else if (model == 2) {
    A = 0;
    alpha = -2;
    height = 0;
    l0 = 1;
    a = 0;
  } else if (model == 3) {
    A = 0;
    alpha = 0;
    height = 10./3;
    l0 = 1;
    a = 0.9;
  } else if (model == 4) {
    A = 1.e5;
    alpha = 0;
    height = 10./3;
    l0 = 1;
    a = 0.9;
  } else if (model == 5) {
    A = 1.e6;
    alpha = 0;
    height = 100./3;
    l0 = 1;
    a = 0.9;
  }

  // We already set stuff from parameters, so set_units here
  set_units();

  printf("Running analytic model %d:\nMBH: %g\na: %g\n", model, MBH, a);
  printf("A: %g\nalpha: %g\nh: %g\nl0: %g\n\n", A, alpha, height, l0);
}

void set_units()
{
  // Derive units we actually need
  MBH = MBH_solar * MSUN;
  L_unit = GNEWT * MBH / (CL * CL);
  T_unit = L_unit / CL;

  // Set all the geometry for coordinates.c
  // TODO function like initialize_coordinates, that makes sure these are all set.
  R0 = 0.;
  Rh = 1 + sqrt(1. - a * a);
  Rin = Rh;
  Rout = 1000.0;
  rmax_geo = MIN(1000.0, Rout);
  startx[0] = 0.0;
  startx[1] = log(Rin);
  startx[2] = 0.0;
  startx[3] = 0.0;
  stopx[0] = 0.0;
  stopx[1] = log(Rout);
  stopx[2] = 1.0;
  stopx[3] = 2*M_PI;
}

void output_hdf5()
{
  hdf5_set_directory("/header/");
  double zero = 0;
  hdf5_write_single_val(&zero, "t", H5T_IEEE_F64LE);
  hdf5_write_single_val(&a, "a", H5T_IEEE_F64LE);

  hdf5_write_single_val(&model, "model", H5T_STD_I32LE);
  hdf5_write_single_val(&A, "A", H5T_IEEE_F64LE);
  hdf5_write_single_val(&alpha, "alpha", H5T_IEEE_F64LE);
  hdf5_write_single_val(&height, "height", H5T_IEEE_F64LE);
  hdf5_write_single_val(&l0, "l0", H5T_IEEE_F64LE);

  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  hdf5_write_single_val(&zero, "M_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&T_unit, "T_unit", H5T_IEEE_F64LE);

  hdf5_set_directory("/");
}

//// INTERFACE: Functions called from elsewhere in ipole ////
void get_model_jar(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV)
{
  // Define a model here relating X,K -> j_S, alpha_S, rho_S
  // (and below relating X,K -> u,B 4-vectors)
  // ipole will do the rest
  double r, th;
  bl_coord(X, &r, &th);
  double n = (3.e-18) * exp(-1./2 * (pow(r/10, 2) + pow(height * cos(th), 2)));

  double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
  get_model_fourv(X, Ucon, Ucov, Bcon, Bcov);
  double nu = get_fluid_nu(Kcon, Ucov);

  *jI = n * pow(nu / freqcgs, -alpha) / pow(nu, 2);
  *aI = A * n * pow(nu / freqcgs, -(2.5 + alpha)) * nu;

  *jQ = 0;
  *jU = 0;
  *jV = 0;

  *aQ = 0;
  *aU = 0;
  *aV = 0;

  *rQ = 0;
  *rU = 0;
  *rV = 0;

  return;
}

void get_model_jk(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv)
{
  // Currently just takes jI, aI from _jar, but can be defined differently for comparison/raw unpolarized problems
  double jI, jQ, jU, jV, aI, aQ, aU, aV, rQ, rU, rV;
  get_model_jar(X, Kcon, &jI, &jQ, &jU, &jV, &aI, &aQ, &aU, &aV, &rQ, &rU, &rV);
  *jnuinv = jI;
  *knuinv = aI;
}

int radiating_region(double X[NDIM])
{
  // If you don't want conditionals in get_model_jar, 
  // you can control here where the coefficients are applied

  return 1;
}

void get_model_fourv(double X[NDIM], double Ucon[NDIM], double Ucov[NDIM],
                     double Bcon[NDIM], double Bcov[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);
  double R = r * sin(th);
  double l = (l0 / (1 + R)) * pow(R, 1 + 0.5);

  double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  // normal observer velocity for Ucon/Ucov
  Ucov[0] =
      - sqrt(-1. / (gcon[0][0] - 2. * gcon[0][3] * l
                  + gcon[3][3] * l * l));
  Ucov[1] = 0.;
  Ucov[2] = 0.;
  Ucov[3] = - l * Ucov[0];

  flip_index(Ucov, gcon, Ucon);
}

//// STUBS: Functions for normal models which we don't use ////
// This is only called for trace file output
void get_model_primitives(double X[NDIM], double *p) {return;}
// Define these to specify a fluid model: e- density/temperature for
// synchrotron radiation based on an energy distribution
double get_model_thetae(double X[NDIM]) {return 0;}
double get_model_b(double X[NDIM]) {return 0;}
double get_model_ne(double X[NDIM]) {return 0;}
void get_model_powerlaw_vals(double X[NDIM], double *p, double *n,
          double *gamma_min, double *gamma_max, double *gamma_cut) {return;}