// Generic model for analytic formulae for emissivities
// Stub out all the 

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
static double Mdot = 0.01;
static double MBH_solar = 10.;
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

  set_by_word_val(word, value, "Mdot", &Mdot, TYPE_DBL);

  set_by_word_val(word, value, "model", &model, TYPE_INT);
  set_by_word_val(word, value, "freqcgs", &freqcgs, TYPE_DBL);
}

void init_model(double *tA, double *tB)
{
  // Set all the geometry globals we need
  // TODO do this in geometry?  Deal with model/geom interface...
  use_eKS_internal = 1;
  metric = 0; // Doesn't matter due to above
  hslope = 1.0;
  // Needed for camera rootfinding
  // TODO standard geometry init that handles this...
  startx[2] = 0.;
  stopx[2] = 1;

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

  printf("Running analytic model %d:\nMBH: %g\nMdot: %g\na: %g\n", model, MBH, Mdot, a);
  printf("GRRT Parameters:\nA: %g\nalpha: %g\nh: %g\nl0: %g\n\n", A, alpha, height, l0);
}

void set_units()
{
  // Convert parameters to consistent CGS units
  MBH = MBH_solar * MSUN;
  // Note definition of Mdotedd w/ 10% efficiency
  double Mdotedd = 4. * M_PI * GNEWT * MBH * MP / 0.1 / CL / SIGMA_THOMSON;
  Mdot *= Mdotedd;

  // Derive everything else
  M_unit = Mdot;
  L_unit = GNEWT * MBH / (CL * CL);
  T_unit = L_unit / CL;
  RHO_unit = M_unit / pow(L_unit, 3.);
  U_unit = RHO_unit * CL * CL;
  B_unit = CL * sqrt(4. * M_PI * RHO_unit);

  // Set all the geometry
  // TODO function like initialize_coordinates, that makes sure these are all set.
  R0 = 0.;
  Rh = 1 + sqrt(1. - a * a);
  double z1 = 1. + pow(1. - a * a, 1. / 3.) * (pow(1. + a, 1. / 3.) + pow(1. - a, 1. / 3.));
  double z2 = sqrt(3. * a * a + z1 * z1);
  double r_isco = 3. + z2 - copysign(sqrt((3. - z1) * (3. + z1 + 2. * z2)), a);
  Rin = Rh;
  Rout = 100.0;
  rmax_geo = MIN(1000., Rout);
  startx[0] = 0.0;
  startx[1] = log(Rin);
  startx[2] = 0.0;
  startx[3] = 0.0;
  stopx[0] = 0.0;
  stopx[1] = log(Rout);
  stopx[2] = 1.0;
  stopx[3] = 2*M_PI;

  fprintf(stderr, "L,T,M units: %g [cm] %g [s] %g [g]\n", L_unit, T_unit, M_unit);
  fprintf(stderr, "rho,u,B units: %g [g cm^-3] %g [g cm^-1 s^-2] %g [G] \n", RHO_unit, U_unit, B_unit);
  fprintf(stderr, "Rh, Rin, Rout, r_isco: %g %g %g %g\n", Rh, Rin, Rout, r_isco);
}

void output_hdf5()
{
  hdf5_set_directory("/header/");
  double zero = 0;
  hdf5_write_single_val(&zero, "t", H5T_IEEE_F64LE);
  hdf5_write_single_val(&a, "a", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Mdot, "Mdot", H5T_IEEE_F64LE);

  hdf5_write_single_val(&model, "model", H5T_STD_I32LE);
  hdf5_write_single_val(&A, "A", H5T_IEEE_F64LE);
  hdf5_write_single_val(&alpha, "alpha", H5T_IEEE_F64LE);
  hdf5_write_single_val(&height, "height", H5T_IEEE_F64LE);
  hdf5_write_single_val(&l0, "l0", H5T_IEEE_F64LE);

  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  hdf5_write_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&M_unit, "M_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&T_unit, "T_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Te_unit, "Thetae_unit", H5T_IEEE_F64LE);

  hdf5_set_directory("/");
}

//// INTERFACE: Functions called from elsewhere in ipole ////
void get_model_jar(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV)
{
  // Define a model here relating X,K -> j_S, alpha_S, rho_S
  // ipole will do the rest
  double r, th;
  bl_coord(X, &r, &th);
  double n = (3.e-18) * exp(-1./2 * (pow(r/10, 2) + pow(height * cos(th), 2)));

  double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
  get_model_fourv(X, Ucon, Ucov, Bcon, Bcov);
  double nu = get_fluid_nu(Kcon, Ucov);

  *jI = n * pow(nu / freqcgs, -alpha) / pow(nu, 2);
  *aI = A * n * pow(nu / freqcgs, -(2.5 + alpha)) * nu;

  //TODO?
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

  return;
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