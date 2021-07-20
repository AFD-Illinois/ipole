// Flat space. Fast light.
// Just an integrator/program stack test, nothing to see here.

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

// TODO get rid of these in ipole proper to get rid of them here
double rmax_geo = 1e30;
double model_dl;
// TODO this default needs to be kept in sync with maxnstep,
// to avoid difficulties
static int max_stokes = 10000;

// Model parameters: private
double jIc, jQc, jUc, jVc;
double aIc, aQc, aUc, aVc;
double rQc, rUc, rVc;

double *stokes_I, *stokes_Q, *stokes_U, *stokes_V;
double *lambda;

void try_set_model_parameter(const char *word, const char *value)
{
  // Test the given word against our parameters' names,
  // and if it matches set the corresponding global
  // TODO stepsize?
  set_by_word_val(word, value, "jI", &jIc, TYPE_DBL);
  set_by_word_val(word, value, "jQ", &jQc, TYPE_DBL);
  set_by_word_val(word, value, "jU", &jUc, TYPE_DBL);
  set_by_word_val(word, value, "jV", &jVc, TYPE_DBL);

  set_by_word_val(word, value, "aI", &aIc, TYPE_DBL);
  set_by_word_val(word, value, "aQ", &aQc, TYPE_DBL);
  set_by_word_val(word, value, "aU", &aUc, TYPE_DBL);
  set_by_word_val(word, value, "aV", &aVc, TYPE_DBL);

  set_by_word_val(word, value, "rQ", &rQc, TYPE_DBL);
  set_by_word_val(word, value, "rU", &rUc, TYPE_DBL);
  set_by_word_val(word, value, "rV", &rVc, TYPE_DBL);

  set_by_word_val(word, value, "dl", &model_dl, TYPE_DBL);
  // TODO precedence
  set_by_word_val(word, value, "maxnstep", &max_stokes, TYPE_INT);
  set_by_word_val(word, value, "max_stokes", &max_stokes, TYPE_INT);
}

void init_model(double *tA, double *tB)
{
  // Set all the geometry globals we need
  metric = METRIC_MINKOWSKI;
  use_eKS_internal = 0;
  Rh = 0.5; // Don't integrate past the coordinate origin, spherical coordinates don't like that

  // Unitless
  M_unit = 1.0;
  L_unit = 1.0;
  T_unit = 1.0;
  RHO_unit = 1.0;
  U_unit = 1.0;
  B_unit = 1.0;

  stokes_I = calloc(max_stokes, sizeof(double));
  stokes_Q = calloc(max_stokes, sizeof(double));
  stokes_U = calloc(max_stokes, sizeof(double));
  stokes_V = calloc(max_stokes, sizeof(double));
  lambda = calloc(max_stokes, sizeof(double));;

}

void output_hdf5()
{
  hdf5_set_directory("/header/");
  hdf5_make_directory("emissivities");
  hdf5_set_directory("/header/emissivities/");

  hdf5_write_single_val(&jIc, "jI", H5T_IEEE_F64LE);
  hdf5_write_single_val(&jQc, "jQ", H5T_IEEE_F64LE);
  hdf5_write_single_val(&jUc, "jU", H5T_IEEE_F64LE);
  hdf5_write_single_val(&jVc, "jV", H5T_IEEE_F64LE);

  hdf5_write_single_val(&aIc, "aI", H5T_IEEE_F64LE);
  hdf5_write_single_val(&aQc, "aQ", H5T_IEEE_F64LE);
  hdf5_write_single_val(&aUc, "aU", H5T_IEEE_F64LE);
  hdf5_write_single_val(&aVc, "aV", H5T_IEEE_F64LE);

  hdf5_write_single_val(&rQc, "rQ", H5T_IEEE_F64LE);
  hdf5_write_single_val(&rUc, "rU", H5T_IEEE_F64LE);
  hdf5_write_single_val(&rVc, "rV", H5T_IEEE_F64LE);

  hdf5_set_directory("/");
  hsize_t stokes_dims[1] = {max_stokes};
  hdf5_write_full_array(stokes_I, "I", 1, stokes_dims, H5T_IEEE_F64LE);
  hdf5_write_full_array(stokes_Q, "Q", 1, stokes_dims, H5T_IEEE_F64LE);
  hdf5_write_full_array(stokes_U, "U", 1, stokes_dims, H5T_IEEE_F64LE);
  hdf5_write_full_array(stokes_V, "V", 1, stokes_dims, H5T_IEEE_F64LE);
  hdf5_write_full_array(lambda, "lam", 1, stokes_dims, H5T_IEEE_F64LE);
}

//// PUBLIC INTERFACE: shortcut all the normal emissivity functions to return constants////

void record_stokes_parameters(double SI, double SQ, double SU, double SV, double lam)
{
  static int nstep = 0;
  stokes_I[nstep] = SI;
  stokes_Q[nstep] = SQ;
  stokes_U[nstep] = SU;
  stokes_V[nstep] = SV;
  lambda[nstep] = lam;
  ++nstep;
}

void get_model_jar(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV)
{
  *jI = jIc;
  *jQ = jQc;
  *jU = jUc;
  *jV = jVc;

  *aI = aIc;
  *aQ = aQc;
  *aU = aUc;
  *aV = aVc;

  *rQ = rQc;
  *rU = rUc;
  *rV = rVc;
}

void get_model_jk(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv)
{
  *jnuinv = jIc;
  *knuinv = aIc;
}


//// STUBS: Keep the rest of ipole off our back ////
int radiating_region(double X[NDIM])
{
    return 1;
}

// B does NOT affect polarization calculation, since we control the coefficients
// U, however, would, so we set it to 0
void get_model_fourv(double X[NDIM], double Kcon[NDIM],
                     double Ucon[NDIM], double Ucov[NDIM],
                     double Bcon[NDIM], double Bcov[NDIM])
{
  double gcov[NDIM][NDIM], gcon[NDIM][NDIM];

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  // Normal observer velocity for Ucon/Ucov and default
  // Bcon/Bcov to zero.

  Ucov[0] = -1./sqrt(-gcov[0][0]);
  Ucov[1] = 0.;
  Ucov[2] = 0.;
  Ucov[3] = 0.;

  Bcon[0] = 0.;
  Bcon[1] = 1.;
  Bcon[2] = 1.;
  Bcon[3] = 1.;

  MULOOP {Ucon[mu] = 0.; Bcon[mu] = 0.;}
  MUNULOOP {
    Ucon[mu] += Ucov[nu] * gcon[mu][nu];
    Bcov[mu] += Bcon[nu] * gcov[mu][nu];
  }
}

double get_model_thetae(double X[NDIM]) {return 0;}
double get_model_b(double X[NDIM]) {return 1;}
double get_model_ne(double X[NDIM]) {return 1;} // Otherwise we trigger the "empty space" emissivity
void get_model_primitives(double X[NDIM], double *p) {return;}
void update_data(double *tA, double *tB) {return;}
void update_data_until(double *tA, double *tB, double tgt) {return;}