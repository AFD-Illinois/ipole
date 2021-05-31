#include "model.h"

#include "decs.h"
#include "coordinates.h"
#include "geometry.h"
#include "radiation.h"
#include "par.h"
#include "tetrads.h"

// Globals we're in charge of
double M_unit;
double L_unit;
double T_unit;
double U_unit;
double B_unit;
double RHO_unit;
// These actually set problem scale
double Te_unit = 1e11;
double Ne_unit = 5e6;

// Model parameters: public
double rmax_geo = 1000.0;
// Model parameters: private
static double rmin_geo = 0;
static double MBH_solar = 4.3e6;

static double MBH;

// e.g. parameterization from GRRT paper
double nth0, Te0, disk_h, pow_nth, pow_T;
double keplerian_factor, infall_factor;
double r_isco;

/**
 * This is a template for analytic problems, which consist of prescription:
 * X,K -> 4-vectors Ucon/cov, Bcon/cov
 * And either:
 * X,K -> emission coefficients jS, alphaS, rhoS
 * or:
 * X,K -> e- density and temperature ne, Thetae
 */

// Forward declarations for non-public functions
void set_units();

//// INITIALIZATION: Functions called from elsewhere in ipole ////

/**
 * This function is called for each word/value pair ipole encounters,
 * either from a parfile or the command line.
 * You can define new pairs here, skim what you need of ipole's defaults, or whatever
 * 
 * ipole will not warn on unspecified parameters. Have good defaults (set on declaration)
 */
void try_set_model_parameter(const char *word, const char *value)
{
  // Test the given word against our parameters' names,
  // and if it matches set the corresponding global
  set_by_word_val(word, value, "MBH", &MBH_solar, TYPE_DBL);
  set_by_word_val(word, value, "Ne_unit", &Ne_unit, TYPE_DBL);
  set_by_word_val(word, value, "Te_unit", &Te_unit, TYPE_DBL);
  // TODO NEED to move this into main parameters
  set_by_word_val(word, value, "rmax_geo", &rmax_geo, TYPE_DBL);
  set_by_word_val(word, value, "rmin_geo", &rmin_geo, TYPE_DBL);

  // Geometry parameters
  set_by_word_val(word, value, "a", &a, TYPE_DBL);
  // RIAF model parameters, from Odyssey.
  // Note nth0,Te0 are *different* as they are unitless!
  set_by_word_val(word, value, "nth0", &nth0, TYPE_DBL);
  set_by_word_val(word, value, "Te0", &Te0, TYPE_DBL);
  set_by_word_val(word, value, "disk_h", &disk_h, TYPE_DBL);
  set_by_word_val(word, value, "pow_nth", &pow_nth, TYPE_DBL);
  set_by_word_val(word, value, "pow_T", &pow_T, TYPE_DBL);
  // Inflow model parameters
  set_by_word_val(word, value, "keplerian_factor", &keplerian_factor, TYPE_DBL);
  set_by_word_val(word, value, "infall_factor", &infall_factor, TYPE_DBL);

  // Normal ipole pulls this, but we also need it for the GRRT problems
  // and this is easier than grabbing it from the 'params' struct
  //set_by_word_val(word, value, "freqcgs", &freqcgs, TYPE_DBL);
}

/**
 * Initialization takes boundary times, for slow light.  Most analytic models won't use them.
 */
void init_model(double *tA, double *tB)
{
  // Set all the geometry globals we need
  // TODO do this in geometry?  Deal with model/geom interface...
  use_eKS_internal = 0;
  metric = METRIC_MKS;
  hslope = 1.0;

  // We already set stuff from parameters, so set_units here
  set_units();

  fprintf(stderr, "Running RIAF model with a=%g, nth0=%g, Te0=%g, disk_h=%g, pow_nth=%g, pow_T=%g\n",
          a, nth0, Te0, disk_h, pow_nth, pow_T);
  fprintf(stderr, "Velocity model: Keplerian by %g, infall rate %g\n",
          keplerian_factor, infall_factor);
  // TODO B field when there are params
}

void set_units()
{
  // Derive units we actually need
  // We already have Te_unit as a parameter
  MBH = MBH_solar * MSUN;
  L_unit = GNEWT * MBH / (CL * CL);
  T_unit = L_unit / CL;
  
  RHO_unit = Ne_unit * (MP + ME);
  B_unit = CL * sqrt(4.*M_PI*RHO_unit);
  M_unit = RHO_unit * pow(L_unit, 3);

  // Set all the geometry for coordinates.c
  // TODO function like initialize_coordinates, that makes sure these are all set.
  R0 = 0.;
  Rh = 1 + sqrt(1. - a * a);
  // These are for tracing.  We only *emit* inside *_geo parameters
  Rin = Rh;
  Rout = rmax_geo;

  double z1 = 1. + pow(1. - a * a, 1. / 3.) * (pow(1. + a, 1. / 3.) + pow(1. - a, 1. / 3.));
  double z2 = sqrt(3. * a * a + z1 * z1);
  r_isco = 3. + z2 - copysign(sqrt((3. - z1) * (3. + z1 + 2. * z2)), a);
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

  // TODO output all riaf terms

  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  // TODO set M_unit correctly
  hdf5_write_single_val(&zero, "M_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&T_unit, "T_unit", H5T_IEEE_F64LE);

  hdf5_set_directory("/");
}

//// INTERFACE: Functions called from elsewhere in ipole ////
double get_model_ne(double X[NDIM])
{
  // Matter model defined in Gold et al 2020 section 3
  double r, th;
  bl_coord(X, &r, &th);
  double zc=r*cos(th);
  double rc=r*sin(th);
  return nth0 * exp(-zc*zc/2./rc/rc/disk_h/disk_h) * pow(r,pow_nth) * Ne_unit;
}

double get_model_thetae(double X[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);
  return Te0 * pow(r, pow_T) * Te_unit * KBOL / (ME*CL*CL); 
}

void get_model_fourv(double X[NDIM], double Kcon[NDIM], double Ucon[NDIM], double Ucov[NDIM],
                     double Bcon[NDIM], double Bcov[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);

  // Metrics: BL
  double bl_gcov[NDIM][NDIM], bl_gcon[NDIM][NDIM];
  gcov_bl(r, th, bl_gcov);
  gcon_func(bl_gcov, bl_gcon);
  // Native
  double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  // Get the 4-velocity
  double bl_Ucon[NDIM];
  double omegaK, omegaFF, omega;
  double K, ur, ut;
  if (r < Rh) {
    // Inside r_h, none
    double bl_Ucov[NDIM];
    bl_Ucov[0] = -1;
    bl_Ucov[1] = 0.;
    bl_Ucov[2] = 0.;
    bl_Ucov[3] = 0.;
    flip_index(bl_Ucov, bl_gcon, bl_Ucon);
  } else if (r < r_isco) {
    // Inside r_isco, freefall
    double omegaK_isco = 1. / (pow(r_isco, 3./2) + a);
  
    // Get conserved quantities at the ISCO...
    double bl_Ucon_isco[NDIM], bl_Ucov_isco[NDIM];
    bl_Ucon_isco[0] = 1.0;
    bl_Ucon_isco[1] = 0.0;
    bl_Ucon_isco[2] = 0.0;
    bl_Ucon_isco[3] = omegaK_isco;

    double bl_gcov_isco[NDIM][NDIM];
    gcov_bl(r_isco, th, bl_gcov_isco);

    normalize(bl_Ucon_isco, bl_gcov_isco);
    flip_index(bl_Ucon_isco, bl_gcov_isco, bl_Ucov_isco);
    double e = bl_Ucov_isco[0];
    double l = bl_Ucov_isco[3];

    // ...then set the infall velocity and find omega
    double bl_Ucon_tmp[NDIM], bl_Ucov_tmp[NDIM];
    double K_con = bl_gcon[0][0] * e * e + 2.0 * bl_gcon[0][3] * e * l + bl_gcon[3][3] * l * l;
    double urk_precut = -(1.0 + K_con) / bl_gcon[1][1];
    double urk = -sqrt(fmax(0.0, urk_precut));
    bl_Ucov_tmp[0] = e;
    bl_Ucov_tmp[1] = urk;
    bl_Ucov_tmp[2] = 0.0;
    bl_Ucov_tmp[3] = l;
    flip_index(bl_Ucov_tmp, bl_gcon, bl_Ucon_tmp);
    omegaK = bl_Ucon_tmp[3] / bl_Ucon_tmp[0];

    omegaFF = bl_gcon[0][3] / bl_gcon[0][0];
    // Compromise
    omega = omegaK + (1 - keplerian_factor)*(omegaFF - omegaK);

    // Then set the infall rate
    double urFF = -sqrt(fmax(0.0, -(1.0 + bl_gcon[0][0]) * bl_gcon[1][1]));
    ur = bl_Ucon_tmp[1] + infall_factor * (urFF - bl_Ucon_tmp[1]);

#if DEBUG
    if (fabs(ur) < 1e-10) {
      fprintf(stderr, "Bad ur: ur is %g\n", ur);
      fprintf(stderr, "Ucon BL: %g %g %g %g\n",
              bl_Ucon_tmp[0], bl_Ucon_tmp[1], bl_Ucon_tmp[2], bl_Ucon_tmp[3]);
      fprintf(stderr, "Ucov BL: %g %g %g %g\n",
              bl_Ucov_tmp[0], bl_Ucov_tmp[1], bl_Ucov_tmp[2], bl_Ucov_tmp[3]);
      fprintf(stderr, "urk was %g (%g pre-cut), e & l were %g %g\n", urk, urk_precut, e, l);
    }
#endif

    // Finally, get Ucon in BL coordinates
    K = bl_gcov[0][0] + 2*omega*bl_gcov[0][3] + omega*omega*bl_gcov[3][3];
    ut = sqrt(fmax(0.0, -(1. + ur*ur*bl_gcov[1][1]) / K));
    bl_Ucon[0] = ut;
    bl_Ucon[1] = ur;
    bl_Ucon[2] = 0.;
    bl_Ucon[3] = omega * ut;
  } else {
    // Outside r_isco, Keplerian
    omegaK = 1. / (pow(r, 3./2) + a);
    omegaFF = bl_gcon[0][3] / bl_gcon[0][0];

    // Compromise
    omega = omegaK + (1 - keplerian_factor)*(omegaFF - omegaK);
    // Set infall rate
    ur = infall_factor * -sqrt(fmax(0.0, -(1.0 + bl_gcon[0][0]) * bl_gcon[1][1]));

    // Get the normal observer velocity for Ucon/Ucov, in BL coordinates
    K = bl_gcov[0][0] + 2*omega*bl_gcov[0][3] + omega*omega*bl_gcov[3][3];
    ut = sqrt(fmax(0.0, -(1. + ur*ur*bl_gcov[1][1]) / K));
    bl_Ucon[0] = ut;
    bl_Ucon[1] = ur;
    bl_Ucon[2] = 0.;
    bl_Ucon[3] = omega * ut;
  }

  // Transform to KS coordinates,
  double ks_Ucon[NDIM];
  bl_to_ks(X, bl_Ucon, ks_Ucon);
  // then to our coordinates,
  vec_from_ks(X, ks_Ucon, Ucon);

  // and grab Ucov
  flip_index(Ucon, gcov, Ucov);


  // Check
#if DEBUG
  //if (r < r_isco) { fprintf(stderr, "ur = %g\n", Ucon[1]); }
  double bl_Ucov[NDIM];
  double dot_U = Ucon[0]*Ucov[0] + Ucon[1]*Ucov[1] + Ucon[2]*Ucov[2] + Ucon[3]*Ucov[3];
  double sum_U = Ucon[0]+Ucon[1]+Ucon[2]+Ucon[3];
  // Following condition gets handled better above
  // (r < r_isco && fabs(Ucon[1]) < 1e-10) ||
  if (get_fluid_nu(Kcon, Ucov) == 1. ||
      fabs(fabs(dot_U) - 1.) > 1e-10 || sum_U < 0.1) {
    flip_index(bl_Ucon, bl_gcov, bl_Ucov);
    fprintf(stderr, "RIAF model problem at r, th, phi = %g %g %g\n", r, th, X[3]);
    fprintf(stderr, "Omega K: %g FF: %g Final: %g K: %g ur: %g ut: %g\n",
            omegaK, omegaFF, omega, K, ur, ut);
    fprintf(stderr, "K1: %g K2: %g K3: %g\n", bl_gcov[0][0], 2*omega*bl_gcov[0][3], omega*omega*bl_gcov[3][3]);
    fprintf(stderr, "Ucon BL: %g %g %g %g\n", bl_Ucon[0], bl_Ucon[1], bl_Ucon[2], bl_Ucon[3]);
    fprintf(stderr, "Ucon KS: %g %g %g %g\n", ks_Ucon[0], ks_Ucon[1], ks_Ucon[2], ks_Ucon[3]);
    fprintf(stderr, "Ucon native: %g %g %g %g\n", Ucon[0], Ucon[1], Ucon[2], Ucon[3]);
    fprintf(stderr, "Ucov: %g %g %g %g\n", Ucov[0], Ucov[1], Ucov[2], Ucov[3]);
    fprintf(stderr, "Ubl.Ubl: %g\n", bl_Ucov[0]*bl_Ucon[0]+bl_Ucov[1]*bl_Ucon[1]+
                                    bl_Ucov[2]*bl_Ucon[2]+bl_Ucov[3]*bl_Ucon[3]);
    fprintf(stderr, "U.U: %g\n", Ucov[0]*Ucon[0]+Ucov[1]*Ucon[1]+Ucov[2]*Ucon[2]+Ucov[3]*Ucon[3]);
  }
#endif

  // Use pure toroidal field,
  // See Themis src/VRT2/src/AccretionFlows/mf_toroidal_beta.cpp/h
  double bl_Bcon[NDIM];
  bl_Bcon[0] = 0.0;
  bl_Bcon[1] = 0.0;
  bl_Bcon[2] = 0.0;
  bl_Bcon[3] = 1.0;

  // Transform to KS coordinates,
  double ks_Bcon[NDIM];
  bl_to_ks(X, bl_Bcon, ks_Bcon);
  // then to our coordinates,
  vec_from_ks(X, ks_Bcon, Bcon);
  normalize(Bcon, gcov);

  // Compute u.b and subtract it, normalize to get_model_b
  //project_out(Bcon, Ucon, gcov); ?
  double BdotU = 0;
  MULOOP BdotU += Bcon[mu] * Ucov[mu];
  MULOOP Bcon[mu] += BdotU * Ucon[mu];
  flip_index(Bcon, gcov, Bcov);
  double Bsq = 0;
  MULOOP Bsq += Bcon[mu] * Bcov[mu];
  double bmag = fmax(get_model_b(X), 1e-10);
  MULOOP Bcon[mu] *= bmag / sqrt(Bsq);

  flip_index(Bcon, gcov, Bcov);
}

/**
 * This problem defines no field in emission, but we want to control ipole's worst
 * tendencies when making tetrads.  This will return a correct value even for a
 * possible fluid/field model later, too.
 */
double get_model_b(double X[NDIM])
{
  // double Ucon[NDIM],Bcon[NDIM];
  // double Ucov[NDIM],Bcov[NDIM];
  // double Kcon[NDIM] = {0}; // TODO interface change if we ever need a real one here
  // get_model_fourv(X, Kcon, Ucon, Ucov, Bcon, Bcov);
  // return sqrt(Bcon[0]*Bcov[0] + Bcon[1]*Bcov[1] + Bcon[2]*Bcov[2] + Bcon[3]*Bcov[3]) * B_unit;

  double r, th;
  bl_coord(X, &r, &th);
  double nth = get_model_ne(X);
  // Get local field strength
  double eps = 0.1;
  double b = sqrt(8. * M_PI * eps * nth * MP * CL * CL / 6. / r);
  if (b == 0) b = 1.e-6;
  return b;
}

int radiating_region(double X[NDIM])
{
  // If you don't want conditionals in get_model_jar, 
  // you can control here where the coefficients are applied
  double r, th;
  bl_coord(X, &r, &th);
  return r > Rh + 0.1 && r > rmin_geo && r < Rout;
}

//// STUBS: Functions for normal models which we don't use ////
// Define these to specify a fluid model: e- density/temperature for
// synchrotron radiation based on an energy distribution
void get_model_powerlaw_vals(double X[NDIM], double *p, double *n,
          double *gamma_min, double *gamma_max, double *gamma_cut) {return;}
void update_data(double *tA, double *tB) {return;}
void update_data_until(double *tA, double *tB, double tgt) {return;}
// This is only called for trace file output, and doesn't really apply to analytic models
void get_model_primitives(double X[NDIM], double *p) {return;}
void get_model_jk(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv) {return;}
void get_model_jar(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV) {return;}
