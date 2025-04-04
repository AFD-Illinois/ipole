#include "model.h"

#include "decs.h"
#include "hdf5_utils.h"
#include "debug_tools.h"

#include "coordinates.h"
#include "geometry.h"
#include "grid.h"
#include "model_radiation.h"  // Only for outputting emissivities
#include "simcoords.h"  // For interpolating arbitrary grids
#include "par.h"
#include "utils.h"

#include "debug_tools.h"

#include <assert.h>
#include <string.h>

#define NVAR (10)
#define USE_FIXED_TPTE (0)
#define USE_MIXED_TPTE (1)
#define NSUP (3)

#define USE_GEODESIC_SIGMACUT (1)

#define FORMAT_HAMR_EKS (3)
#define FORMAT_IHARM_v1 (1)
#define FORMAT_KORAL_v2 (4)
int dumpfile_format = 0;

// UNITS
double M_unit;
double L_unit;
double T_unit;
double RHO_unit;
double U_unit;
double B_unit;
double Te_unit;

double target_mdot = 0.0;  // if this value > 0, use to renormalize M_unit &c.

// MOLECULAR WEIGHTS
static double Ne_factor = 1.;  // e.g., used for He with 2 protons+neutrons per 2 electrons
static double mu_i, mu_e, mu_tot;

// MODEL PARAMETERS: PUBLIC
double DTd;
double sigma_cut = 1.;
double beta_crit = 1.0;
double sigma_cut_high = -1.0;

// MODEL PARAMETERS: PRIVATE
static char fnam[STRLEN] = "dump.h5";

static double tp_over_te = 3.;
static double trat_small = 1.;
static double trat_large = 40.;
// Minimum number of dynamical times the cooling time must
// undershoot to be considered "small"
// lower values -> higher max T_e, higher values are restrictive
static double cooling_dynamical_times = 1.e-20;

static int dumpskip = 1;
static int dumpmin, dumpmax, dumpidx;
static double MBH_solar = 6.2e9;
static double MBH; // Set from previous
static double Mdot_dump;
static double MdotEdd_dump;
static double Ladv_dump;

static int reverse_field = 0;

double tf;

// MAYBES
//static double t0;

// ELECTRONS
//    0 : constant TP_OVER_TE
//    1 : use dump file model (kawazura?)
//    2 : use mixed TP_OVER_TE (beta model)
//    3 : use mixed TP_OVER_TE (beta model) with fluid temperature
//    9 : load Te (in Kelvin) from dump file (KORAL etc.)
// TODO the way this is selected is horrid.  Make it a parameter.
#define ELECTRONS_TFLUID (3)
static int RADIATION, ELECTRONS;
static double gam = 1.444444444444444, game = 1.333333333333333, gamp = 1.666666666666667;
static double Thetae_unit, Mdotedd;

// Ignore radiation interactions within one degree of polar axis
static double polar_cut = -1;
static double th_beg = 0.0174;
static int nloaded = 0;


static hdf5_blob fluid_header = { 0 };


// Debug KHARMA reader
#define DEBUG_READER (0)

struct of_data {
  double t;
  double ****p;
  double ***ne;
  double ***thetae;
  double ***b;
  double ***sigma;
  double ***beta;
};
static struct of_data dataA, dataB, dataC;
static struct of_data *data[NSUP];

// Definitions for functions not in model.h interface
void set_units();
void load_data(int n, char *, int dumpidx, int verbose);
void load_iharm_data(int n, char *, int dumpidx, int verbose);
void load_koral_data(int n, char *, int dumpidx, int verbose);
void load_hamr_data(int n, char *, int dumpidx, int verbose);
double get_dump_t(char *fnam, int dumpidx);
void init_grid(char *fnam, int dumpidx);
void init_iharm_grid(char *fnam, int dumpidx);
void init_koral_grid(char *fnam, int dumpidx);
void init_hamr_grid(char *fnam, int dumpidx);
void init_physical_quantities(int n, double rescale_factor);
void init_storage(void);

void try_set_model_parameter(const char *word, const char *value)
{
  // TODO remember to set defaults!

  // ipole no longer supports fixed-order command line arguments!
  // assume params is populated
  set_by_word_val(word, value, "MBH", &MBH_solar, TYPE_DBL);
  set_by_word_val(word, value, "M_unit", &M_unit, TYPE_DBL);
  set_by_word_val(word, value, "mdot", &target_mdot, TYPE_DBL);

  set_by_word_val(word, value, "dump", (void *)fnam, TYPE_STR);

  set_by_word_val(word, value, "tp_over_te", &tp_over_te, TYPE_DBL);
  set_by_word_val(word, value, "trat_small", &trat_small, TYPE_DBL);
  set_by_word_val(word, value, "trat_large", &trat_large, TYPE_DBL);
  set_by_word_val(word, value, "sigma_cut", &sigma_cut, TYPE_DBL);
  set_by_word_val(word, value, "sigma_cut_high", &sigma_cut_high, TYPE_DBL);
  set_by_word_val(word, value, "beta_crit", &beta_crit, TYPE_DBL);
  set_by_word_val(word, value, "cooling_dynamical_times", &cooling_dynamical_times, TYPE_DBL);

  set_by_word_val(word, value, "rmax_geo", &rmax_geo, TYPE_DBL);
  set_by_word_val(word, value, "rmin_geo", &rmin_geo, TYPE_DBL);

  set_by_word_val(word, value, "reverse_field", &reverse_field, TYPE_INT);
  // allow cutting out the spine
  set_by_word_val(word, value, "polar_cut_deg", &polar_cut, TYPE_DBL);

  // for slow light
  set_by_word_val(word, value, "dump_min", &dumpmin, TYPE_INT);
  set_by_word_val(word, value, "dump_max", &dumpmax, TYPE_INT);
  set_by_word_val(word, value, "dump_skip", &dumpskip, TYPE_INT);
  dumpidx = dumpmin;

  // override parameters based on input
  if (target_mdot > 0.) {
    M_unit = 1.;
  }
}

// Advance through dumps until we are closer to the next set
// of dumps corresponding to tA == tgt. Used when attempting
// to restart from slowlight restart file.
void update_data_until(double *tA, double *tB, double tgt)
{
  double tC = data[2]->t;

  while (tC < tgt) {
    dumpidx += dumpskip;
    tC = get_dump_t(fnam, dumpidx);
  }

  // reset dump index, just to be safe, then load on out ...
  dumpidx -= dumpskip;
  while (*tA < tgt) update_data(tA, tB);
}

// Use internal dumpidx variable to load the next "expected"
// dump into memory (used for slowlight mode). After calling
// this function, it is guaranteed that data are ordered:
//
//   data[0]->t < data[1]->t < data[2]->t
//
// This function uses dataA, dataB, dataC to save the actual
// data locations and then swaps which members live where.
void update_data(double *tA, double *tB)
{
  #if SLOW_LIGHT
  // Reorder dataA, dataB, dataC in data[]
  if (nloaded % 3 == 0) {
    data[0] = &dataB;
    data[1] = &dataC;
    data[2] = &dataA;
  } else if (nloaded % 3 == 1) {
    data[0] = &dataC;
    data[1] = &dataA;
    data[2] = &dataB;
  } else {
    data[0] = &dataA;
    data[1] = &dataB;
    data[2] = &dataC;
  }
  int nextdumpidx = dumpidx;
  dumpidx += dumpskip;
  if (nextdumpidx > dumpmax) {
    load_data(2, fnam, --nextdumpidx, 0);
    data[2]->t += 1.;
  } else {
    load_data(2, fnam, nextdumpidx, 0);
  }
  *tA = data[0]->t;
  *tB = data[1]->t;
  fprintf(stderr, "loaded data (dump %d) (%g < t < %g)\n", nextdumpidx, *tA, *tB);
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
}

double get_dump_t(char *fnam, int dumpidx)
{
  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);
  double t = -1.;

  hdf5_open(fname);
  if (dumpfile_format == FORMAT_IHARM_v1) {
    hdf5_set_directory("/");
    hdf5_read_single_val(&t, "/t", H5T_IEEE_F64LE);
  } else if (dumpfile_format == FORMAT_HAMR_EKS) {
    hdf5_read_attr_num(&t, "t", "", H5T_IEEE_F64LE);
  }
  hdf5_close();

  return t;
}

void get_dumpfile_type(char *fnam, int dumpidx)
{
  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);

  hdf5_open(fname);
  int hamr_attr_exists = hdf5_attr_exists("", "dscale");

  if (!hamr_attr_exists) {
    if (!hdf5_exists("header/version")) {
      // Converted BHAC dumps and very old iharm3d output do not include a version
      // BHAC output requires nothing special from ipole so we mark it "iharm_v1"
      dumpfile_format = FORMAT_IHARM_v1;
      fprintf(stderr, "bhac! (or old iharm)\n");
    } else {
      char harmversion[256];
      hdf5_read_single_val(harmversion, "header/version", hdf5_make_str_type(255));

      if ( strcmp(harmversion, "KORALv2") == 0 ) {
        dumpfile_format = FORMAT_KORAL_v2;
        fprintf(stderr, "koral!\n");
      } else {
        dumpfile_format = FORMAT_IHARM_v1;
        fprintf(stderr, "iharm!\n");
      }
    }
  } else {
    // note this will return -1 if the "header" group does not exist
    dumpfile_format = FORMAT_HAMR_EKS;
    fprintf(stderr, "hamr!\n");
  }

  hdf5_close();
}

void init_model(double *tA, double *tB)
{
  // set up initial ordering of data[]
  data[0] = &dataA;
  data[1] = &dataB;
  data[2] = &dataC;

  fprintf(stderr, "Determining dump file type... ");
  get_dumpfile_type(fnam, dumpmin);

  // set up grid for fluid data
  fprintf(stderr, "Reading data header...\n");
  init_grid(fnam, dumpmin);

  // set all dimensional quantities from loaded parameters
  set_units();

  // read fluid data
  fprintf(stderr, "Reading data...\n");
  load_data(0, fnam, dumpidx, 2);  
  // replaced dumpmin -> 2 because apparently that argument was just .. removed
  dumpidx += dumpskip;
  #if SLOW_LIGHT
  update_data(tA, tB);
  update_data(tA, tB);
  tf = get_dump_t(fnam, dumpmax) - 1.e-5;
  #else // FAST LIGHT
  data[2]->t = 10000.;
  #endif // SLOW_LIGHT

  // horizon radius
  Rh = 1 + sqrt(1. - a * a);

  // possibly cut around the pole
  if (polar_cut >= 0) {
    th_beg = 0.0174 * polar_cut;
  } else {
    if (dumpfile_format == FORMAT_HAMR_EKS) {
      th_beg = 0.0174 * 2;
    }
  }

  #if DEBUG_READER
    /* Set filename */
    char debug_fname[256];
    snprintf(debug_fname, sizeof(debug_fname), "debug_reader_iharm.h5");

    /* Create HDF5 file*/
    hdf5_create(debug_fname);

    /* Compute gcov, gcon */
    size_t total_elements = NDIM * NDIM * (N2 + 2) * (N1 + 2);
    double *gcov_global = malloc(total_elements * sizeof(double));
    double *gcon_global = malloc(total_elements * sizeof(double));
    // Use the arrays via indexing. For example, to access element [mu][nu][j][i]:
    #define IDX(mu, nu, j, i) (((mu) * NDIM * (N2+2) * (N1+2)) + ((nu) * (N2+2) * (N1+2)) + ((j) * (N1+2)) + (i))

#pragma omp parallel for collapse(2)
  for (int i = 0; i < N1+2; i++) {
    for (int j = 0; j < N2+2; j++) {

      double X[NDIM] = {0.};
      ijktoX(i, j, 0, X);
      double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
      gcov_func(X, gcov);
      gcon_func(gcov, gcon);

      for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
          gcov_global[IDX(mu, nu, j, i)] = gcov[mu][nu];
          gcon_global[IDX(mu, nu, j, i)] = gcon[mu][nu];
        }
      }

    }
  }

    /* Write gcov, gcon to file */
    hsize_t dims_grid[4] = { NDIM, NDIM, N2+2, N1+2 };
    hdf5_write_full_array(gcov_global, "gcov", 4, dims_grid, H5T_NATIVE_DOUBLE);
    hdf5_write_full_array(gcon_global, "gcon", 4, dims_grid, H5T_NATIVE_DOUBLE);

    /* Write physical quantities */
    hsize_t dims_phys[3] = { N1+2, N2+2, N3+2 };
    hdf5_write_full_array(data[0]->ne[0][0], "ne", 3, dims_phys, H5T_NATIVE_DOUBLE);
    hdf5_write_full_array(data[0]->thetae[0][0], "thetae", 3, dims_phys, H5T_NATIVE_DOUBLE);
    hdf5_write_full_array(data[0]->b[0][0], "b", 3, dims_phys, H5T_NATIVE_DOUBLE);
    hdf5_write_full_array(data[0]->sigma[0][0], "sigma", 3, dims_phys, H5T_NATIVE_DOUBLE);
    hdf5_write_full_array(data[0]->beta[0][0], "beta", 3, dims_phys, H5T_NATIVE_DOUBLE);

    /* Close HDF5 file */
    hdf5_close();

  #endif // DEBUG_READER
}

/*
  these supply basic model data to ipole
*/

// In slowlight mode, we perform linear interpolation in time. This function tells
// us how far we've progressed from data[nA]->t to data[nB]->t but "in reverse" as
// tinterp == 1 -> we're exactly on nA and tinterp == 0 -> we're exactly on nB.
double set_tinterp_ns(double X[NDIM], int *nA, int *nB)
{
  #if SLOW_LIGHT
  if (X[0] < data[1]->t) {
    *nA = 0; *nB = 1;
  } else {
    *nA = 1; *nB = 2;
  }
  double tinterp = 1. - ( X[0] - data[*nA]->t ) / ( data[*nB]->t - data[*nA]->t );
  if (tinterp < 0.) tinterp = 0.; //  in slow light, when we reset based on tB, sometimes we overshoot
  if (tinterp > 1.) tinterp = 1.; //  TODO, this should really only happen at r >> risco, but still...
  return tinterp;
  #else
  *nA = 0;
  *nB = 0;
  return 0.;
  #endif // SLOW_LIGHT
}

// Calculate Ucon,Ucov,Bcon,Bcov from primitives at location X using 
// interpolation (on the primitives). This has been all wrapped into
// a single function because some calculations require each other.
void get_model_fourv(double X[NDIM], double Kcon[NDIM],
                     double Ucon[NDIM], double Ucov[NDIM],
                     double Bcon[NDIM], double Bcov[NDIM])
{
  double gcov[NDIM][NDIM], gcon[NDIM][NDIM];

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  // If we're outside of the logical domain, default to
  // normal observer velocity for Ucon/Ucov and default
  // Bcon/Bcov to zero.
  if ( X_in_domain(X) == 0 ) {

    Ucov[0] = -1./sqrt(-gcon[0][0]);
    Ucov[1] = 0.;
    Ucov[2] = 0.;
    Ucov[3] = 0.;
    Ucon[0] = 0.;
    Ucon[1] = 0.;
    Ucon[2] = 0.;
    Ucon[3] = 0.;

    for (int mu=0; mu<NDIM; ++mu) {
      Ucon[0] += Ucov[mu] * gcon[0][mu];
      Ucon[1] += Ucov[mu] * gcon[1][mu];
      Ucon[2] += Ucov[mu] * gcon[2][mu];
      Ucon[3] += Ucov[mu] * gcon[3][mu];
      Bcon[mu] = 0.;
      Bcov[mu] = 0.;
    }

    return;
  }

  // Set Ucon and get Ucov by lowering

  // interpolate primitive variables first
  double Vcon[NDIM];
  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);
  Vcon[1] = interp_scalar_time(X, data[nA]->p[U1], data[nB]->p[U1], tfac);
  Vcon[2] = interp_scalar_time(X, data[nA]->p[U2], data[nB]->p[U2], tfac);
  Vcon[3] = interp_scalar_time(X, data[nA]->p[U3], data[nB]->p[U3], tfac);

  // translate to four velocity
  double VdotV = 0.;
  for (int i = 1; i < NDIM; i++)
    for (int j = 1; j < NDIM; j++)
      VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
  double Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
  Ucon[0] = -Vfac * gcon[0][0];
  for (int i = 1; i < NDIM; i++)
    Ucon[i] = Vcon[i] - Vfac * gcon[0][i];

  // lower (needed for Bcon)
  lower(Ucon, gcov, Ucov);

  // Now set Bcon and get Bcov by lowering

  // interpolate primitive variables first
  double Bcon1 = interp_scalar_time(X, data[nA]->p[B1], data[nB]->p[B1], tfac);
  double Bcon2 = interp_scalar_time(X, data[nA]->p[B2], data[nB]->p[B2], tfac);
  double Bcon3 = interp_scalar_time(X, data[nA]->p[B3], data[nB]->p[B3], tfac);

  // get Bcon
  Bcon[0] = Bcon1*Ucov[1] + Bcon2*Ucov[2] + Bcon3*Ucov[3];
  Bcon[1] = (Bcon1 + Ucon[1] * Bcon[0]) / Ucon[0];
  Bcon[2] = (Bcon2 + Ucon[2] * Bcon[0]) / Ucon[0];
  Bcon[3] = (Bcon3 + Ucon[3] * Bcon[0]) / Ucon[0];

  // lower
  lower(Bcon, gcov, Bcov);
}

// Get the primitive variables interpolated to a point X,
// And fill them in the next 8 array slots after p
// Not used for transport but useful for plotting along a geodesic later
void get_model_primitives(double X[NDIM], double *p)
{
  if ( X_in_domain(X) == 0 ) return;

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  for (int np=0; np<8; np++) {
    p[np] = interp_scalar_time(X, data[nA]->p[np], data[nA]->p[np], tfac);
  }
}

double get_model_thetae(double X[NDIM])
{
  if ( X_in_domain(X) == 0 ) return 0.;
  
  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);
  double thetae = interp_scalar_time(X, data[nA]->thetae, data[nB]->thetae, tfac);

#if DEBUG
  if (thetae < 0. || isnan(thetae)) {
    printf("thetae negative or NaN!\n");
    printf("X[] = %g %g %g %g\n", X[0], X[1], X[2], X[3]);
    printf("t = %e %e %e\n", data[0]->t, data[1]->t, data[2]->t);
    double thetaeA = interp_scalar(X, data[nA]->thetae);
    double thetaeB = interp_scalar(X, data[nB]->thetae);
    printf("thetaeA, thetaeB = %e %e", thetaeA, thetaeB);
    printf("thetae, tfac = %e %e\n", thetae, tfac);
  }
#endif

  return thetae;
}

//b field strength in Gauss
double get_model_b(double X[NDIM])
{
  // TODO how *should* we handle exiting the domain consistently?
  if ( X_in_domain(X) == 0 ) return 0.;
  
  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  return interp_scalar_time(X, data[nA]->b, data[nB]->b, tfac);
}

double get_model_sigma(double X[NDIM]) 
{
  if ( X_in_domain(X) == 0 ) return 0.;

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  return interp_scalar_time(X, data[nA]->sigma, data[nB]->sigma, tfac);
}

double get_model_beta(double X[NDIM]) 
{
  if ( X_in_domain(X) == 0 ) return 0.;  // TODO inf?

  double betaA, betaB, tfac;
  int nA, nB;
  tfac = set_tinterp_ns(X, &nA, &nB);
  betaA = interp_scalar(X, data[nA]->beta);
  betaB = interp_scalar(X, data[nB]->beta);
  return tfac*betaA + (1. - tfac)*betaB;
}

double get_sigma_smoothfac(double sigma)
{
  double sigma_above = sigma_cut;
  if (sigma_cut_high > 0) sigma_above = sigma_cut_high;
  if (sigma < sigma_cut) return 1;
  if (sigma >= sigma_above) return 0;
  double dsig = sigma_above - sigma_cut;
  return cos(M_PI / 2. / dsig * (sigma - sigma_cut));
}

double get_model_ne(double X[NDIM])
{
  if ( X_in_domain(X) == 0 ) return 0.;

  double sigma_smoothfac = 1;

#if USE_GEODESIC_SIGMACUT
  double sigma = get_model_sigma(X);
  if (sigma > sigma_cut) return 0.;
  sigma_smoothfac = get_sigma_smoothfac(sigma);
#endif

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  return interp_scalar_time(X, data[nA]->ne, data[nB]->ne, tfac) * sigma_smoothfac;
}

void set_units()
{
  MBH = MBH_solar * MSUN; // Convert to CGS
  L_unit = GNEWT * MBH / (CL * CL);
  T_unit = L_unit / CL;
  RHO_unit = M_unit / pow(L_unit, 3);
  U_unit = RHO_unit * CL * CL;
  B_unit = CL * sqrt(4.*M_PI*RHO_unit);
  Mdotedd=4.*M_PI*GNEWT*MBH*MP/CL/0.1/SIGMA_THOMSON;

  fprintf(stderr,"MBH: %g [Msun]\n",MBH/MSUN);
  fprintf(stderr,"L,T,M units: %g [cm] %g [s] %g [g]\n",L_unit,T_unit,M_unit) ;
  fprintf(stderr,"rho,u,B units: %g [g cm^-3] %g [g cm^-1 s^-2] %g [G] \n",RHO_unit,U_unit,B_unit) ;
}

void init_physical_quantities(int n, double rescale_factor)
{
#if DEBUG
  int ceilings = 0;
#endif

  rescale_factor = sqrt(rescale_factor);

  // cover everything, even ghost zones
#pragma omp parallel for collapse(3)
  for (int i = 0; i < N1+2; i++) {
    for (int j = 0; j < N2+2; j++) {
      for (int k = 0; k < N3+2; k++) {
        data[n]->ne[i][j][k] = data[n]->p[KRHO][i][j][k] * RHO_unit/(MP+ME) * Ne_factor;

        data[n]->b[i][j][k] *= rescale_factor;

        double bsq = data[n]->b[i][j][k] / B_unit;
        bsq = bsq*bsq;

        double sigma_m = bsq/data[n]->p[KRHO][i][j][k];
        double beta_m = data[n]->p[UU][i][j][k]*(gam-1.)/0.5/bsq;
#if DEBUG
        if(isnan(sigma_m)) {
          sigma_m = 0;
          fprintf(stderr, "Setting zero sigma!\n");
        }
        if(isnan(beta_m)) {
          beta_m = INFINITY;
          fprintf(stderr, "Setting INF beta!\n");
        }
#endif
        if (ELECTRONS == 1) {
          data[n]->thetae[i][j][k] = data[n]->p[KEL][i][j][k] * 
                                     pow(data[n]->p[KRHO][i][j][k],game-1.)*Thetae_unit;
        } else if (ELECTRONS == 2) {
          double betasq = beta_m*beta_m / beta_crit/beta_crit;
          double trat = trat_large * betasq/(1. + betasq) + trat_small /(1. + betasq);
          //Thetae_unit = (gam - 1.) * (MP / ME) / trat;
          // see, e.g., Eq. 8 of the EHT GRRT formula list
          double lcl_Thetae_u = (MP/ME) * (game-1.) * (gamp-1.) / ( (gamp-1.) + (game-1.)*trat );
          Thetae_unit = lcl_Thetae_u;
          data[n]->thetae[i][j][k] = lcl_Thetae_u*data[n]->p[UU][i][j][k]/data[n]->p[KRHO][i][j][k];
        } else if (ELECTRONS == 9) {
          // convert Kelvin -> Thetae
          data[n]->thetae[i][j][k] = data[n]->p[TFLK][i][j][k] * KBOL / ME / CL / CL;
        } else if (ELECTRONS == ELECTRONS_TFLUID) {
          double beta = data[n]->p[UU][i][j][k]*(gam-1.)/0.5/bsq;
          double betasq = beta*beta / beta_crit/beta_crit;
          double trat = trat_large * betasq/(1. + betasq) + trat_small /(1. + betasq);
          double dfactor = mu_tot / mu_e + mu_tot / mu_i * trat;
          data[n]->thetae[i][j][k] = data[n]->p[THF][i][j][k] / dfactor;
        } else {
          data[n]->thetae[i][j][k] = Thetae_unit*data[n]->p[UU][i][j][k]/data[n]->p[KRHO][i][j][k];
        }
#if DEBUG
        if(isnan(data[n]->thetae[i][j][k])) {
          data[n]->thetae[i][j][k] = 0.0;
          fprintf(stderr, "\nNaN Thetae!  Prims %g %g %g %g %g %g %g %g\n", data[n]->p[KRHO][i][j][k], data[n]->p[UU][i][j][k],
                  data[n]->p[U1][i][j][k], data[n]->p[U2][i][j][k], data[n]->p[U3][i][j][k], data[n]->p[B1][i][j][k],
                  data[n]->p[B2][i][j][k], data[n]->p[B3][i][j][k]);
          fprintf(stderr, "Setting floor temp!\n");
        }
#endif

        // Enforce a max on Thetae based on cooling time == dynamical time
        if (cooling_dynamical_times > 1e-20) {
          double X[NDIM];
          ijktoX(i, j, k, X);
          double r, th;
          bl_coord(X, &r, &th);
          // Calculate thetae_max based on matching the cooling time w/dynamical time
          // Makes sure to use b w/units, but r has already been rescaled
          double Thetae_max_dynamical =  1 / cooling_dynamical_times * 7.71232e46 / 2 / MBH * pow(data[n]->b[i][j][k], -2) * pow(r * sin(th), -1.5);
#if DEBUG
          if (Thetae_max_dynamical < data[n]->thetae[i][j][k]) {
            if (r > 2) fprintf(stderr, "Ceiling on temp! %g < %g, r, th %g %g\n", Thetae_max_dynamical, data[n]->thetae[i][j][k], r, th);
            ceilings++;
          }
#endif
          data[n]->thetae[i][j][k] = fmin(data[n]->thetae[i][j][k], Thetae_max_dynamical);
        }

        // Apply floor last in case the above is a very restrictive ceiling
        data[n]->thetae[i][j][k] = fmax(data[n]->thetae[i][j][k], 1.e-3);

        // Preserve sigma for cutting along geodesics, and for variable-kappa model
        data[n]->sigma[i][j][k] = fmax(sigma_m, SMALL);
        // Also record beta, for variable-kappa model
        data[n]->beta[i][j][k] = fmax(beta_m, SMALL);

        // Cut Ne (i.e. emission) based on sigma, if we're not doing so along each geodesic
        // Strongly magnetized = empty, no shiny spine
        if (sigma_m > sigma_cut && !USE_GEODESIC_SIGMACUT) {
          data[n]->b[i][j][k]=0.0;
          data[n]->ne[i][j][k]=0.0;
          data[n]->thetae[i][j][k]=0.0;
        }
      }
    }
  }
#if DEBUG
  fprintf(stderr, "TOTAL TEMPERATURE CEILING ZONES: %d of %d\n", ceilings, (N1+2)*(N2+2)*(N3+2));
#endif

}

void init_storage(void)
{ 
  // one ghost zone on each side of the domain
  for (int n = 0; n < NSUP; n++) {
    data[n]->p = malloc_rank4(NVAR,N1+2,N2+2,N3+2);
    data[n]->ne = malloc_rank3(N1+2,N2+2,N3+2);
    data[n]->thetae = malloc_rank3(N1+2,N2+2,N3+2);
    data[n]->b = malloc_rank3(N1+2,N2+2,N3+2);
    data[n]->sigma = malloc_rank3(N1+2,N2+2,N3+2);
    data[n]->beta = malloc_rank3(N1+2,N2+2,N3+2);
  }
}

void init_grid(char *fnam, int dumpidx)
{
  if (dumpfile_format == FORMAT_IHARM_v1) {
    init_iharm_grid(fnam, dumpidx);
  } else if (dumpfile_format == FORMAT_KORAL_v2) {
    init_koral_grid(fnam, dumpidx);
  } else if (dumpfile_format == FORMAT_HAMR_EKS) {
    init_hamr_grid(fnam, dumpidx);
  }
}

void init_hamr_grid(char *fnam, int dumpidx) 
{
  // no support for electron thermodynamics or radiation
  ELECTRONS = 0;
  RADIATION = 0;

  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);
  fprintf(stderr, "filename: %s\n", fname);
  hdf5_open(fname);

  // TODO: maybe reconstruct header blob?

  // note that MKS assumes x2 \in (0, 1), so we need to 
  // deal with remapping when we load the data
  metric = METRIC_MKS;
  hslope = 1.;
  hdf5_read_attr_num(&a, "a", "", H5T_IEEE_F64LE);

  hdf5_read_attr_num(&N1, "N1", "", H5T_STD_I32LE);
  hdf5_read_attr_num(&N2, "N2", "", H5T_STD_I32LE);
  hdf5_read_attr_num(&N3, "N3", "", H5T_STD_I32LE);

  hdf5_read_attr_num(&gam, "gam", "", H5T_IEEE_F64LE);
  game = 4./3; gamp = 5./3;

  Te_unit = Thetae_unit;
  // we can override which electron model to use here. print results if we're
  // overriding anything. ELECTRONS should only be nonzero if we need to make
  // use of extra variables (instead of just UU and RHO) for thetae
  if (!USE_FIXED_TPTE && !USE_MIXED_TPTE) {
    fprintf(stderr, "! no electron temperature model specified in model/iharm.c\n");
    exit(-3);
  } else if (USE_FIXED_TPTE && !USE_MIXED_TPTE) {
    ELECTRONS = 0; // force TP_OVER_TE to overwrite bad electrons
    fprintf(stderr, "using fixed tp_over_te ratio = %g\n", tp_over_te);
    //Thetae_unit = MP/ME*(gam-1.)*1./(1. + tp_over_te);
    // see, e.g., Eq. 8 of the EHT GRRT formula list.
    // this formula assumes game = 4./3 and gamp = 5./3
    Thetae_unit = 2./3. * MP/ME / (2. + tp_over_te);
  } else if (USE_MIXED_TPTE && !USE_FIXED_TPTE) {
    ELECTRONS = 2;
    fprintf(stderr, "using mixed tp_over_te with trat_small = %g, trat_large = %g, and beta_crit = %g\n",
      trat_small, trat_large, beta_crit);
    // Thetae_unit set per-zone below
  } else {
    fprintf(stderr, "! please change electron model in model/iharm.c\n");
    exit(-3);
  }

  // by this point, we're sure that Thetae_unit is what we want so we can set
  // Te_unit which is what ultimately get written to the dump files
  Te_unit = Thetae_unit;

  // TODO is there a way to determine this?
  DTd = 10.;

  // the grid is a bit confusing, but I've tried my best to reproduce it. the
  // info from the attributes does not match, and since we need the metric to
  // reconstruct values, it's important to get it right. so rather than doing 
  // it blindly, it might be better to ensure the output is correct.
   
  // get domain information from attributes
  hdf5_read_attr_arr(&(startx[1]), "startx", "");
  hdf5_read_attr_arr(&(dx[1]), "dx", "");
  hdf5_read_attr_num(&Rin, "Rin", "", H5T_IEEE_F64LE);
  hdf5_read_attr_num(&Rout, "Rout", "", H5T_IEEE_F64LE);

  // it appears that all of the npz coordinates are actually zone centers
  // so we need to translate startx[i] -= dx[i]/2, but this does not seem
  // to work nicely for phi, so we overwrite
  startx[1] -= dx[1]/2;
  startx[2] -= dx[2]/2;
  startx[3] = 0.;

  stopx[0] = 1.;
  stopx[1] = startx[1]+N1*dx[1];
  stopx[2] = startx[2]+N2*dx[2];
  stopx[3] = startx[3]+N3*dx[3];

  // now translate from hamr x2 \in (-1, 1) -> mks x2 \in (0, 1)
  startx[2] = (startx[2] + 1)/2.;
  stopx[2] = (stopx[2] + 1)/2.;
  dx[2] /= 2;

  // set limit for tracking geodesic emission
  rmax_geo = fmin(rmax_geo, Rout);
  rmin_geo = fmax(rmin_geo, Rin);

  // set boundary of coordinate system
  MULOOP cstartx[mu] = 0.;
  cstopx[0] = 0;
  cstopx[1] = log(Rout);
  cstopx[2] = 1.0;
  cstopx[3] = 2*M_PI;

  fprintf(stderr, "Native coordinate start: %g %g %g stop: %g %g %g\n",
                  cstartx[1], cstartx[2], cstartx[3], cstopx[1], cstopx[2], cstopx[3]);
  fprintf(stderr, "Grid start: %g %g %g stop: %g %g %g\n",
                  startx[1], startx[2], startx[3], stopx[1], stopx[2], stopx[3]);

  init_storage();
  hdf5_close();
}

void init_iharm_grid(char *fnam, int dumpidx)
{
  // called at the beginning of the run and sets the static parameters
  // along with setting up the grid
  
  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);
  fprintf(stderr, "filename: %s\n", fname);
  hdf5_open(fname);

  // get dump info to copy to ipole output
  //  hdf5_read_single_val(&t0, "t", H5T_IEEE_F64LE);
  fluid_header = hdf5_get_blob("/header");

  hdf5_set_directory("/header/");

  if ( hdf5_exists("has_electrons") ) {
    hdf5_read_single_val(&ELECTRONS, "has_electrons", H5T_STD_I32LE);
  } else {
    ELECTRONS = 0;
  }
  if ( hdf5_exists("has_radiation") ) {
    hdf5_read_single_val(&RADIATION, "has_radiation", H5T_STD_I32LE);
  } else {
    RADIATION = 0;
  }
  if ( hdf5_exists("has_derefine_poles") ) {
    printf("Dump includes flag 'has_derefine_poles' and is therefore non-standard and not well-defined");
    exit(-3);
  }

  if ( hdf5_exists("weights") ) {
    hdf5_set_directory("/header/weights/");
    hdf5_read_single_val(&mu_i, "mu_i", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mu_e, "mu_e", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mu_tot, "mu_tot", H5T_IEEE_F64LE);
    fprintf(stderr, "Loaded molecular weights (mu_i, mu_e, mu_tot): %g %g %g\n", mu_i, mu_e, mu_tot);
    Ne_factor = 1. / mu_e;
    hdf5_set_directory("/header/");
    ELECTRONS = ELECTRONS_TFLUID;
  }

  char metric_name[20];
  hid_t HDF5_STR_TYPE = hdf5_make_str_type(20);
  hdf5_read_single_val(&metric_name, "metric", HDF5_STR_TYPE);

  metric = 0;
  use_eKS_internal = 0;

  // Read the coordinate system and set the only value about it not recorded in the dump
  if ( strncmp(metric_name, "MKS", 19) == 0) {
    metric = METRIC_MKS;
    cstopx[2] = 1.0;
  } else if ( strncmp(metric_name, "BHAC_MKS", 19) == 0) {
    metric = METRIC_BHACMKS;
    cstopx[2] = M_PI;
  } else if ( strncmp(metric_name, "MMKS", 19) == 0 || strncmp(metric_name, "FMKS", 19) == 0) {
    // TODO Handle MMKS vs FMKS signifiers in dumps somehow...
    metric = METRIC_FMKS;
    cstopx[2] = 1.0;
  } else if ( strncmp(metric_name, "MKS3", 19) == 0 ) {
    use_eKS_internal = 1;
    metric = METRIC_MKS3;
    cstopx[2] = 1.0;
  } else if ( strncmp(metric_name, "EKS", 19) == 0 ) {
    metric = METRIC_EKS;
    cstopx[2] = M_PI;
  } else {
    fprintf(stderr, "File is in unknown metric %s.  Cannot continue.\n", metric_name);
    exit(-1);
  }

  hdf5_read_single_val(&N1, "n1", H5T_STD_I32LE);
  hdf5_read_single_val(&N2, "n2", H5T_STD_I32LE);
  hdf5_read_single_val(&N3, "n3", H5T_STD_I32LE);
  hdf5_read_single_val(&gam, "gam", H5T_IEEE_F64LE);

  if (hdf5_exists("gam_e")) { /* pyharm-converted files save gam_e, gam_p even when electrons are not run. Which overwrites the default values*/
    fprintf(stderr, "custom electron model loaded from dump file...\n");
    hdf5_read_single_val(&game, "gam_e", H5T_IEEE_F64LE);
    hdf5_read_single_val(&gamp, "gam_p", H5T_IEEE_F64LE);
  } // Else use default values set above
  Te_unit = Thetae_unit;

  // we can override which electron model to use here. print results if we're
  // overriding anything. ELECTRONS should only be nonzero if we need to make
  // use of extra variables (instead of just UU and RHO) for thetae
  if (!USE_FIXED_TPTE && !USE_MIXED_TPTE) {
    if (ELECTRONS != 1) {
      fprintf(stderr, "! no electron temperature model specified! Cannot continue\n");
      exit(-3);
    }
    ELECTRONS = 1;
    Thetae_unit = MP/ME;
  } else if (ELECTRONS == ELECTRONS_TFLUID) {
    fprintf(stderr, "Using Ressler/Athena electrons with mixed tp_over_te and\n");
    fprintf(stderr, "trat_small = %g, trat_large = %g, and beta_crit = %g\n", trat_small, trat_large, beta_crit);
  } else if (USE_FIXED_TPTE && !USE_MIXED_TPTE) {
    ELECTRONS = 0; // force TP_OVER_TE to overwrite bad electrons
    fprintf(stderr, "Using fixed tp_over_te ratio = %g\n", tp_over_te);
    //Thetae_unit = MP/ME*(gam-1.)*1./(1. + tp_over_te);
    // see, e.g., Eq. 8 of the EHT GRRT formula list. 
    // this formula assumes game = 4./3 and gamp = 5./3
    Thetae_unit = 2./3. * MP/ME / (2. + tp_over_te);
  } else if (USE_MIXED_TPTE && !USE_FIXED_TPTE) {
    ELECTRONS = 2;
    fprintf(stderr, "Using mixed tp_over_te with trat_small = %g, trat_large = %g, and beta_crit = %g\n", 
      trat_small, trat_large, beta_crit);
    // Thetae_unit set per-zone below
  } else {
    fprintf(stderr, "Unknown electron model %d! Cannot continue.\n", ELECTRONS);
    exit(-3);
  }
  fprintf(stderr, "sigma_cut = %g\n", sigma_cut);

  // by this point, we're sure that Thetae_unit is what we want so we can set
  // Te_unit which is what ultimately get written to the dump files
  Te_unit = Thetae_unit;

  if (RADIATION) {
    fprintf(stderr, "custom radiation field tracking information loaded...\n");
    fprintf(stderr, "!! warning, this branch is not tested!\n");
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

  switch (metric) {
    case METRIC_MKS:
      hdf5_set_directory("/header/geom/mks/");
      fprintf(stderr, "Using Modified Kerr-Schild coordinates MKS\n");
      break;
    case METRIC_BHACMKS:
      hdf5_set_directory("/header/geom/bhac_mks/");
      fprintf(stderr, "Using BHAC-style Modified Kerr-Schild coordinates BHAC_MKS\n");
      break;
    case METRIC_FMKS:
      hdf5_set_directory("/header/geom/mmks/"); // For compat. TODO autodetect 'fmks' here?
      fprintf(stderr, "Using Funky Modified Kerr-Schild coordinates FMKS\n");
      break;
    case METRIC_MKS3:
      hdf5_set_directory("/header/geom/mks3/");
      fprintf(stderr, "Using logarithmic KS coordinates internally\n");
      fprintf(stderr, "Converting from KORAL-style Modified Kerr-Schild coordinates MKS3\n");
      break;
    case METRIC_EKS:
      hdf5_set_directory("/header/geom/eks/");
      fprintf(stderr, "Using Kerr-Schild coordinates with exponential radial coordiante\n");
      break;
  }
  
  if ( metric == METRIC_MKS3 ) {
    hdf5_read_single_val(&a, "a", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mks3R0, "R0", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mks3H0, "H0", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mks3MY1, "MY1", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mks3MY2, "MY2", H5T_IEEE_F64LE);
    hdf5_read_single_val(&mks3MP0, "MP0", H5T_IEEE_F64LE);
    Rout = 100.; 
  } else if ( metric == METRIC_EKS ) {
    hdf5_read_single_val(&a, "a", H5T_IEEE_F64LE);
    hdf5_read_single_val(&Rin, "r_in", H5T_IEEE_F64LE);
    hdf5_read_single_val(&Rout, "r_out", H5T_IEEE_F64LE);
    fprintf(stderr, "eKS parameters a: %f Rin: %f Rout: %f\n", a, Rin, Rout);
    
  } else { // Some brand of MKS.  All have the same parameters
    hdf5_read_single_val(&a, "a", H5T_IEEE_F64LE);
    hdf5_read_single_val(&hslope, "hslope", H5T_IEEE_F64LE);
    if (hdf5_exists("Rin")) {
      hdf5_read_single_val(&Rin, "Rin", H5T_IEEE_F64LE);
      hdf5_read_single_val(&Rout, "Rout", H5T_IEEE_F64LE);
    } else {
      hdf5_read_single_val(&Rin, "r_in", H5T_IEEE_F64LE);
      hdf5_read_single_val(&Rout, "r_out", H5T_IEEE_F64LE);
    }
    fprintf(stderr, "MKS parameters a: %f hslope: %f Rin: %f Rout: %f\n", a, hslope, Rin, Rout);

    if (metric == METRIC_FMKS) {
      hdf5_read_single_val(&poly_xt, "poly_xt", H5T_IEEE_F64LE);
      hdf5_read_single_val(&poly_alpha, "poly_alpha", H5T_IEEE_F64LE);
      hdf5_read_single_val(&mks_smooth, "mks_smooth", H5T_IEEE_F64LE);
      poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha + 1.)*1./pow(poly_xt, poly_alpha));
      fprintf(stderr, "MKS parameters poly_xt: %f poly_alpha: %f mks_smooth: %f poly_norm: %f\n", poly_xt, poly_alpha, mks_smooth, poly_norm);
    }
  }

  // Don't emit beyond specified limit or coordinate limit
  rmax_geo = fmin(rmax_geo, Rout);
  rmin_geo = fmax(rmin_geo, Rin);

  hdf5_set_directory("/");
  hdf5_read_single_val(&DTd, "dump_cadence", H5T_IEEE_F64LE);

  // The rest of the 
  stopx[0] = 1.;
  stopx[1] = startx[1]+N1*dx[1];
  stopx[2] = startx[2]+N2*dx[2];
  stopx[3] = startx[3]+N3*dx[3];
  // Start & stop for specifically the *coordinate system*, for step sizes & so on
  cstartx[0] = 0;
  cstartx[1] = 0;
  cstartx[2] = 0;
  cstartx[3] = 0;
  cstopx[0] = 0;
  cstopx[1] = log(Rout);
  cstopx[3] = 2*M_PI;

  fprintf(stderr, "Native coordinate start: %g %g %g stop: %g %g %g\n",
                  cstartx[1], cstartx[2], cstartx[3], cstopx[1], cstopx[2], cstopx[3]);
  fprintf(stderr, "Grid start: %g %g %g stop: %g %g %g\n",
                  startx[1], startx[2], startx[3], stopx[1], stopx[2], stopx[3]);

  init_storage();
  hdf5_close();
}

void init_koral_grid(char *fnam, int dumpidx)
{
  // called at the beginning of the run and sets the static parameters
  // along with setting up the grid

  // assert(42==0);
  // this version of the code has not been validated to have the
  // right units and four-vector recovery. use at your own peril
  
  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);
  fprintf(stderr, "filename: %s\n", fname);

  // always load ks rh from the file to get proper extents
  // and potentially as a fallback for when we cannnot get
  // the inverse/reverse transformation eKS -> simcoords.
  load_simcoord_info_from_file(fname);

  // beecause of legacy global state choices, we can't open two 
  // files at once, so we have to call this after the above.
  hdf5_open(fname);

  // get dump info to copy to ipole output
  fluid_header = hdf5_get_blob("/header");

  // get time information for slow light
  // currently unsupported
  //hdf5_read_single_val(&DTd, "dump_cadence", H5T_IEEE_F64LE);

  // in general, we assume KORAL should use the simcoords feature. first
  // ensure that metric_out is KS and then set simcoords
  hdf5_set_directory("/header/");
  char metric_out[20], metric_run[20];
  hid_t HDF5_STR_TYPE = hdf5_make_str_type(20);
  hdf5_read_single_val(&metric_out, "metric_out", HDF5_STR_TYPE);
  hdf5_read_single_val(&metric_run, "metric_run", HDF5_STR_TYPE);

  if (strcmp(metric_out, "KS") != 0) {
    fprintf(stderr, "! expected koral metric_out==KS but got %s instead. quitting.\n", metric_out);
    exit(5);
  }

  // at this point, assume we will load data as eKS and trace as eKS. 
  use_simcoords = 1;
  metric = METRIC_MKS;
  hslope = 1.;

  // get simulation grid coordinates
  if (strcmp(metric_run, "MKS3") == 0) {
    simcoords = SIMCOORDS_KORAL_MKS3;
    hdf5_read_single_val(&a, "bhspin", H5T_IEEE_F64LE);
    hdf5_set_directory("/header/geom/mks3/");
    hdf5_read_single_val(&(mp_koral_mks3.r0), "mksr0", H5T_IEEE_F64LE);
    hdf5_read_single_val(&(mp_koral_mks3.h0), "mksh0", H5T_IEEE_F64LE);
    hdf5_read_single_val(&(mp_koral_mks3.my1), "mksmy1", H5T_IEEE_F64LE);
    hdf5_read_single_val(&(mp_koral_mks3.my2), "mksmy2", H5T_IEEE_F64LE);
    hdf5_read_single_val(&(mp_koral_mks3.mp0), "mksmp0", H5T_IEEE_F64LE);
    fprintf(stderr, "KORAL simulation was run with MKS3 coordinates.\n");
  } else if (strcmp(metric_run, "MKS2") == 0) {
    simcoords = SIMCOORDS_KORAL_MKS3;
    hdf5_read_single_val(&a, "bhspin", H5T_IEEE_F64LE);
    hdf5_set_directory("/header/geom/mks2/");
    hdf5_read_single_val(&(mp_koral_mks3.r0), "mksr0", H5T_IEEE_F64LE);
    hdf5_read_single_val(&(mp_koral_mks3.h0), "mksh0", H5T_IEEE_F64LE);
    mp_koral_mks3.my1 = 0;
    mp_koral_mks3.my2 = 0;
    mp_koral_mks3.mp0 = 0;
    fprintf(stderr, "KORAL simulation was run with MKS2 coordinates.\n");
  } else if (strcmp(metric_run, "JETCOORDS") == 0) {
    simcoords = SIMCOORDS_KORAL_JETCOORDS;
    hdf5_read_single_val(&a, "bhspin", H5T_IEEE_F64LE);
    hdf5_set_directory("/header/geom/jetcoords/");
    hdf5_read_single_val(&(mp_koral_jetcoords.alpha_1), "alpha1", H5T_IEEE_F64LE);
    hdf5_read_single_val(&(mp_koral_jetcoords.alpha_2), "alpha2", H5T_IEEE_F64LE);
    hdf5_read_single_val(&(mp_koral_jetcoords.cylindrify), "cylindrify", H5T_IEEE_F64LE);
    hdf5_set_directory("/header/geom/");
    fprintf(stderr, "KORAL simulation was run with JETCOORDS coordinates.\n");
  } else {
    fprintf(stderr, "! unknown koral metric_run (%s). quitting.\n", metric_run);
    exit(5);
  }

  // get grid information
  hdf5_set_directory("/header/");
  hdf5_read_single_val(&N1, "n1", H5T_STD_I32LE);
  hdf5_read_single_val(&N2, "n2", H5T_STD_I32LE);
  hdf5_read_single_val(&N3, "n3", H5T_STD_I32LE);
  hdf5_read_single_val(&gam, "gam", H5T_IEEE_F64LE);

  hdf5_set_directory("/header/geom/");
  hdf5_read_single_val(&startx[1], "startx1", H5T_IEEE_F64LE);
  hdf5_read_single_val(&startx[2], "startx2", H5T_IEEE_F64LE);
  hdf5_read_single_val(&startx[3], "startx3", H5T_IEEE_F64LE);
  hdf5_read_single_val(&dx[1], "dx1", H5T_IEEE_F64LE);
  hdf5_read_single_val(&dx[2], "dx2", H5T_IEEE_F64LE);
  hdf5_read_single_val(&dx[3], "dx3", H5T_IEEE_F64LE);
  hdf5_set_directory("/header/");

  // reset x3 grid. this makes it easier for simcoords + Xtoijk to handle the
  // translation, but it also means ipole and KORAL disagree about where x3=0
  if ( fabs(2.*M_PI - dx[3]*N3) >= dx[3] ) {
    fprintf(stderr, "! base koral domain extent in x3 is not 2 PI. quitting.\n");
    exit(5);
  } else {
    startx[3] = 0.;
  }

  stopx[0] = 1.;
  stopx[1] = startx[1]+N1*dx[1];
  stopx[2] = startx[2]+N2*dx[2];
  stopx[3] = startx[3]+N3*dx[3];
  MULOOP cstartx[mu] = startx[mu];
  MULOOP cstopx[mu] = stopx[mu];

  fprintf(stderr, "KORAL coordinates dx: %g %g %g\n", dx[1], dx[2], dx[3]);
  fprintf(stderr, "Native coordinate start: %g %g %g stop: %g %g %g\n",
                  cstartx[1], cstartx[2], cstartx[3], cstopx[1], cstopx[2], cstopx[3]);
  fprintf(stderr, "Grid start: %g %g %g stop: %g %g %g\n",
                  startx[1], startx[2], startx[3], stopx[1], stopx[2], stopx[3]);

  // unsupported features / undocumented elements
  hdf5_set_directory("/header/geom/mks3/");
  if (hdf5_exists("gam_e")) {
    fprintf(stderr, "! found gam_e in KORAL simulation. undocumented. quitting.\n");
    exit(5);
  }

  // maybe load radiation units from dump file
  RADIATION = 0;
  hdf5_set_directory("/header/");
  if ( hdf5_exists("has_radiation") ) {
    hdf5_read_single_val(&RADIATION, "has_radiation", H5T_STD_I32LE);
    if (RADIATION) {
      // Note set_units(...) get called AFTER this function returns
      fprintf(stderr, "koral dump file was from radiation run. loading units...\n");
      hdf5_set_directory("/header/units/");
      hdf5_read_single_val(&MBH_solar, "M_bh", H5T_IEEE_F64LE);
      hdf5_read_single_val(&RHO_unit, "M_unit", H5T_IEEE_F64LE);
      hdf5_read_single_val(&T_unit, "T_unit", H5T_IEEE_F64LE);
      hdf5_read_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
      hdf5_read_single_val(&U_unit, "U_unit", H5T_IEEE_F64LE);
      hdf5_read_single_val(&B_unit, "B_unit", H5T_IEEE_F64LE);
      M_unit = RHO_unit * L_unit*L_unit*L_unit;
    }
  }

  ELECTRONS = 0;
  hdf5_set_directory("/header/");
  if ( hdf5_exists("has_electrons") ) {
    hdf5_read_single_val(&ELECTRONS, "has_electrons", H5T_STD_I32LE);
    if (ELECTRONS) {
      fprintf(stderr, "koral dump has native electron temperature. forcing Thetae...\n");
      ELECTRONS = 9;
    } else {
      if (USE_MIXED_TPTE && !USE_FIXED_TPTE) {
        fprintf(stderr, "Using mixed tp_over_te with trat_small = %g, trat_large = %g, and beta_crit = %g\n", trat_small, trat_large, beta_crit);
        ELECTRONS = 2;
      } else {
        fprintf(stderr, "! koral unsupported without native electrons or mixed tp_over_te.\n");
        exit(6);
      }
    }
  }

  fprintf(stderr, "sigma_cut = %g\n", sigma_cut);
  
  fprintf(stderr, "generating simcoords grid... ");
  int interp_n1 = 1024;
  int interp_n2 = 1024;
  initialize_simgrid(interp_n1, interp_n2, 
                     startx[1], startx[1]+N1*dx[1], 
                     startx[2], startx[2]+N2*dx[2]);
  fprintf(stderr, "done!\n");

  init_storage();
  hdf5_close();
}

void output_hdf5()
{
  hdf5_set_directory("/");

  if (dumpfile_format == FORMAT_IHARM_v1) {
    hdf5_write_blob(fluid_header, "/fluid_header");
  }

  hdf5_write_single_val(&Mdot_dump, "Mdot", H5T_IEEE_F64LE);
  hdf5_write_single_val(&MdotEdd_dump, "MdotEdd", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Ladv_dump, "Ladv", H5T_IEEE_F64LE);

  hdf5_set_directory("/header/");
#if SLOW_LIGHT
  hdf5_write_single_val(&(data[1]->t), "t", H5T_IEEE_F64LE);
#else // FAST LIGHT
  hdf5_write_single_val(&(data[0]->t), "t", H5T_IEEE_F64LE);
#endif

  hdf5_write_single_val(&sigma_cut, "sigma_cut", H5T_IEEE_F64LE);
  hdf5_make_directory("electrons");
  hdf5_set_directory("/header/electrons/");
  if (ELECTRONS == 0) {
    hdf5_write_single_val(&tp_over_te, "tp_over_te", H5T_IEEE_F64LE);
  } else if (ELECTRONS == 2) {
    hdf5_write_single_val(&trat_small, "rlow", H5T_IEEE_F64LE);
    hdf5_write_single_val(&trat_large, "rhigh", H5T_IEEE_F64LE);
    hdf5_write_single_val(&beta_crit, "beta_crit", H5T_IEEE_F64LE);
  } else if (ELECTRONS == ELECTRONS_TFLUID) {
    hdf5_write_single_val(&mu_i, "mu_i", H5T_IEEE_F64LE);
    hdf5_write_single_val(&mu_e, "mu_e", H5T_IEEE_F64LE);
    hdf5_write_single_val(&mu_tot, "mu_tot", H5T_IEEE_F64LE);
  }
  hdf5_write_single_val(&ELECTRONS, "type", H5T_STD_I32LE);

  hdf5_set_directory("/header/");
  hdf5_write_single_val(&reverse_field,"field_config",H5T_STD_I32LE);
  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  hdf5_write_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&M_unit, "M_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&T_unit, "T_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Te_unit, "Thetae_unit", H5T_IEEE_F64LE);

  hdf5_set_directory("/");

  //fprintf(stderr, "Wrote model header\n");
}

void load_data(int n, char *fnam, int dumpidx, int verbose)
{
  if (dumpfile_format == FORMAT_IHARM_v1) {
    load_iharm_data(n, fnam, dumpidx, verbose);
  } else if (dumpfile_format == FORMAT_KORAL_v2) {
    load_koral_data(n, fnam, dumpidx, verbose);
  } else if (dumpfile_format == FORMAT_HAMR_EKS) {
    load_hamr_data(n, fnam, dumpidx, verbose);
  }
}

void populate_boundary_conditions(int n)
{
  // radial -- just extend zones
#pragma omp parallel for collapse(2)
  for (int j=1; j<N2+1; ++j) {
    for (int k=1; k<N3+1; ++k) {
      for (int l=0; l<NVAR; ++l) {
        data[n]->p[l][0][j][k] = data[n]->p[l][1][j][k];
        data[n]->p[l][N1+1][j][k] = data[n]->p[l][N1][j][k];
      }
      data[n]->b[0][j][k] = data[n]->b[1][j][k];
      data[n]->b[N1+1][j][k] = data[n]->b[N1][j][k];
    }
  }

  // elevation -- flip (this is a rotation by pi)
#pragma omp parallel for collapse(2)
  for (int i=0; i<N1+2; ++i) {
    for (int k=1; k<N3+1; ++k) {
      if (N3%2 == 0) {
        int kflip = ( (k - 1) + (N3/2) ) % N3 + 1;
        for (int l=0; l<NVAR; ++l) {
          data[n]->p[l][i][0][k] = data[n]->p[l][i][1][kflip];
          data[n]->p[l][i][N2+1][k] = data[n]->p[l][i][N2][kflip];
        }
        data[n]->b[i][0][k] = data[n]->b[i][1][kflip];
        data[n]->b[i][N2+1][k] = data[n]->b[i][N2][kflip];
      } else {
        int kflip1 = ( k + (N3/2) ) % N3;
        int kflip2 = ( k + (N3/2) + 1 ) % N3;
        for (int l=0; l<NVAR; ++l) {
          data[n]->p[l][i][0][k]    = ( data[n]->p[l][i][1][kflip1]
                                      + data[n]->p[l][i][1][kflip2] ) / 2.;
          data[n]->p[l][i][N2+1][k] = ( data[n]->p[l][i][N2][kflip1]
                                      + data[n]->p[l][i][N2][kflip2] ) / 2.;
        }
        data[n]->b[i][0][k]    = ( data[n]->b[i][1][kflip1]
                                 + data[n]->b[i][1][kflip2] ) / 2.;
        data[n]->b[i][N2+1][k] = ( data[n]->b[i][N2][kflip1]
                                 + data[n]->b[i][N2][kflip2] ) / 2.;
      }
    }
  }

  // azimuth -- periodic
#pragma omp parallel for collapse(2)
  for (int i=0; i<N1+2; ++i) {
    for (int j=0; j<N2+2; ++j) {
      for (int l=0; l<NVAR; ++l) {
        data[n]->p[l][i][j][0] = data[n]->p[l][i][j][N3];
        data[n]->p[l][i][j][N3+1] = data[n]->p[l][i][j][1];
      }
      data[n]->b[i][j][0] = data[n]->b[i][j][N3];
      data[n]->b[i][j][N3+1] = data[n]->b[i][j][1];
    }
  }
}


void remap_hamr(double *buffer, double ***memory, int n1, int n2, int n3, int ng)
{
  for (int i=0; i<n1; ++i) {
    for (int j=0; j<n2; ++j) {
      for (int k=0; k<n3; ++k) {
        memory[ng+i][ng+j][ng+k] = buffer[k+n3*j+n2*n3*i];
      } 
    }
  }
}

void load_hamr_data(int n, char *fnam, int dumpidx, int verbose)
{
  double dMact, Ladv;

  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);

  nloaded++;

  if ( hdf5_open(fname) < 0 ) {
    fprintf(stderr, "! unable to open file %s. Exiting!\n", fname);
    exit(-1);
  }

  // load 1d data into proper locations
  hsize_t fdims[] = { N1 * N2 * N3 };
  hsize_t fstart[] = { 0 };
  hsize_t fcount[] = { N1 * N2 * N3 };

  double *buffer = calloc(N1*N2*N3, sizeof(*buffer));

  hdf5_read_array(buffer, "RHO", 1, fdims, fstart, fcount, fdims, fstart, H5T_IEEE_F64LE);
  remap_hamr(buffer, data[n]->p[KRHO], N1, N2, N3, 1);

  hdf5_read_array(buffer, "UU", 1, fdims, fstart, fcount, fdims, fstart, H5T_IEEE_F64LE);
  remap_hamr(buffer, data[n]->p[UU], N1, N2, N3, 1);

  hdf5_read_array(buffer, "U1", 1, fdims, fstart, fcount, fdims, fstart, H5T_IEEE_F64LE);
  remap_hamr(buffer, data[n]->p[U1], N1, N2, N3, 1);

  hdf5_read_array(buffer, "U2", 1, fdims, fstart, fcount, fdims, fstart, H5T_IEEE_F64LE);
  remap_hamr(buffer, data[n]->p[U2], N1, N2, N3, 1);

  hdf5_read_array(buffer, "U3", 1, fdims, fstart, fcount, fdims, fstart, H5T_IEEE_F64LE);
  remap_hamr(buffer, data[n]->p[U3], N1, N2, N3, 1);

  hdf5_read_array(buffer, "B1", 1, fdims, fstart, fcount, fdims, fstart, H5T_IEEE_F64LE);
  remap_hamr(buffer, data[n]->p[B1], N1, N2, N3, 1);

  hdf5_read_array(buffer, "B2", 1, fdims, fstart, fcount, fdims, fstart, H5T_IEEE_F64LE);
  remap_hamr(buffer, data[n]->p[B2], N1, N2, N3, 1);

  hdf5_read_array(buffer, "B3", 1, fdims, fstart, fcount, fdims, fstart, H5T_IEEE_F64LE);
  remap_hamr(buffer, data[n]->p[B3], N1, N2, N3, 1);

  free(buffer);

  hdf5_read_attr_num(&(data[n]->t), "t", "", H5T_IEEE_F64LE);

  hdf5_close();

  dMact = Ladv = 0.;

  // construct four-vectors over "real" zones
#pragma omp parallel for collapse(2) reduction(+:dMact,Ladv)
  for(int i = 1; i < N1+1; i++) {
    for(int j = 1; j < N2+1; j++) {

      // we need to be tricky here. there are two coordinate systems that 
      // we have to translate between.
      double X[NDIM] = { 0. };
      double gcov_eKS[NDIM][NDIM], gcon_eKS[NDIM][NDIM];
      double gcov_hamr[NDIM][NDIM], gcon_hamr[NDIM][NDIM];
      double ucon[NDIM], ucov[NDIM];
      double bcon[NDIM], bcov[NDIM];
      double g, r, th, alpha;

      // this assumes axisymmetry in the coordinates
      ijktoX(i-1,j-1,0, X);
      bl_coord(X, &r, &th);

      // use gcov_eKS as temporary storage for gcov_KS during gcov_hamr calculation
      gcov_ks(r, th, gcov_eKS);
      //fprintf(stderr, "KS %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n", i,j,r,th,
      //  gcov_eKS[0][0],gcov_eKS[0][1],gcov_eKS[0][2],gcov_eKS[0][3],
      //  gcov_eKS[1][1],gcov_eKS[1][2],gcov_eKS[1][3],
      //  gcov_eKS[2][2],gcov_eKS[2][3],
      //  gcov_eKS[3][3]);

      double dxdX[NDIM][NDIM];
      MUNULOOP dxdX[mu][nu] = delta(mu, nu);
      dxdX[1][1] = r;
      dxdX[2][2] = M_PI/2.;
      
      MUNULOOP {
        gcov_hamr[mu][nu] = 0.;
        for (int lam=0; lam<NDIM; ++lam) {
          for (int kap=0; kap<NDIM; ++kap) {
            gcov_hamr[mu][nu] += gcov_eKS[lam][kap] * dxdX[lam][mu] * dxdX[kap][nu];
          }
        }
      }
      gcon_func(gcov_hamr, gcon_hamr);

      //fprintf(stderr, "HAMR %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n", i,j,r,th,
      //  gcov_hamr[0][0],gcov_hamr[0][1],gcov_hamr[0][2],gcov_hamr[0][3],
      //  gcov_hamr[1][1],gcov_hamr[1][2],gcov_hamr[1][3],
      //  gcov_hamr[2][2],gcov_hamr[2][3],
      //  gcov_hamr[3][3]);

      // now we can freely use gcov_eKS
      gcov_func(X, gcov_eKS);
      gcon_func(gcov_eKS, gcon_eKS);
      g = gdet_zone(i-1,j-1,0);
      alpha = sqrt(-1. / gcon_eKS[0][0]);

      for (int k=1; k<N3+1; ++k) {

        ijktoX(i-1,j-1,k,X);
        double UdotU = 0.;

        // construct ucon/ucov with everything hamr
        for(int l = 1; l < NDIM; l++) {
          for(int m = 1; m < NDIM; m++) {
            UdotU += gcov_hamr[l][m]*data[n]->p[U1+l-1][i][j][k]*data[n]->p[U1+m-1][i][j][k];
          }
        }
        double ufac = sqrt(-1./gcon_hamr[0][0]*(1 + fabs(UdotU)));
        ucon[0] = -ufac * gcon_hamr[0][0];
        for(int l = 1; l < NDIM; l++) {
          ucon[l] = data[n]->p[U1+l-1][i][j][k] - ufac*gcon_hamr[0][l];
        }
        flip_index(ucon, gcov_hamr, ucov);

        //fprintf(stderr, "ucon %d %d %d -> %g %g %g %g\n", i,j,k, ucon[0],ucon[1],ucon[2],ucon[3]);
        //fprintf(stderr, "ucov %d %d %d -> %g %g %g %g\n", i,j,k, ucov[0],ucov[1],ucov[2],ucov[3]);

        // reconstruct the magnetic field with everything hamr
        double udotB = 0.;
        for (int l = 1; l < NDIM; l++) {
          udotB += ucov[l]*data[n]->p[B1+l-1][i][j][k] / alpha;  // note alpha factor!
        }
        bcon[0] = udotB;
        for (int l = 1; l < NDIM; l++) {
          bcon[l] = (data[n]->p[B1+l-1][i][j][k]/alpha + ucon[l]*udotB)/ucon[0];  // note alpha!
        }
        flip_index(bcon, gcov_hamr, bcov);

        double bsq = 0.;
        for (int l=0; l<NDIM; ++l) bsq += bcon[l] * bcov[l];
        data[n]->b[i][j][k] = sqrt(bsq) * B_unit;
 
        //fprintf(stderr, "bcon %d %d %d -> %g %g %g %g\n", i,j,k, bcon[0],bcon[1],bcon[2],bcon[3]);
        //fprintf(stderr, "ucon %d %d %d -> %g %g %g %g\n", i,j,k, ucon[0],ucon[1],ucon[2],ucon[3]);
        //fprintf(stderr, "B %d %d %d -> %g %g %g\n", i,j,k, data[n]->p[B1][i][j][k],data[n]->p[B2][i][j][k],data[n]->p[B3][i][j][k]);
        //fprintf(stderr, "bsq = %g  %d %d %d\n", bsq, i, j, k);

        // translate from hamr -> eKS
        ucon[2] /= 2.;
        ucov[2] *= 2.;
        bcon[2] /= 2.;
        bcov[2] *= 2.;

        // resynthesize the primitives. note we have assumed that B^i = *F^{ti}
        data[n]->p[U1][i][j][k] = (gcon_eKS[0][1]*alpha*alpha + ucon[1]/ucon[0]) * ucon[0];
        data[n]->p[U2][i][j][k] = (gcon_eKS[0][2]*alpha*alpha + ucon[2]/ucon[0]) * ucon[0];
        data[n]->p[U3][i][j][k] = (gcon_eKS[0][3]*alpha*alpha + ucon[3]/ucon[0]) * ucon[0];
        data[n]->p[B1][i][j][k] = ucon[0] * bcon[1] - bcon[0] * ucon[1];
        data[n]->p[B2][i][j][k] = ucon[0] * bcon[2] - bcon[0] * ucon[2];
        data[n]->p[B3][i][j][k] = ucon[0] * bcon[3] - bcon[0] * ucon[3];

        // now do more checks. left in for posterity
        if (1 == 0) {
          fprintf(stderr, "init ucon %d %d %d %g %g %g %g\n", i,j,k,ucon[0],ucon[1],ucon[2],ucon[3]);
          double UdotU = 0.;
          for(int l = 1; l < NDIM; l++) {
            for(int m = 1; m < NDIM; m++) {
              UdotU += gcov_eKS[l][m]*data[n]->p[U1+l-1][i][j][k]*data[n]->p[U1+m-1][i][j][k];
            }
          }
          double ufac = sqrt(-1./gcon_eKS[0][0]*(1 + fabs(UdotU)));
          ucon[0] = -ufac * gcon_eKS[0][0];
          for(int l = 1; l < NDIM; l++) {
            ucon[l] = data[n]->p[U1+l-1][i][j][k] - ufac*gcon_eKS[0][l];
          }
          flip_index(ucon, gcov_eKS, ucov);
          fprintf(stderr, "end ucon %d %d %d %g %g %g %g\n", i,j,k,ucon[0],ucon[1],ucon[2],ucon[3]);

          fprintf(stderr, "init bcon %d %d %d %g %g %g %g\n", i,j,k,bcon[0],bcon[1],bcon[2],bcon[3]);
          double udotB = 0.;
          for (int l = 1; l < NDIM; l++) {
            udotB += ucov[l]*data[n]->p[B1+l-1][i][j][k];
          }
          bcon[0] = udotB;
          for (int l = 1; l < NDIM; l++) {
            bcon[l] = (data[n]->p[B1+l-1][i][j][k] + ucon[l]*udotB)/ucon[0];
          }
          fprintf(stderr, "end bcon %d %d %d %g %g %g %g\n", i,j,k,bcon[0],bcon[1],bcon[2],bcon[3]);

          bsq = 0.;
          MULOOP bsq += bcon[mu]*bcov[mu];
          fprintf(stderr, "end bsq %d %d %d %g\n", i,j,k, bsq);
        }

        // compute diagnostics about the dump. note g is in eKS, so ucon+ must be too!
        if(i <= 21) { dMact += g * data[n]->p[KRHO][i][j][k] * ucon[1]; }
        if(i >= 21 && i < 41 && 0) Ladv += g * data[n]->p[UU][i][j][k] * ucon[1] * ucov[0];
        if(i <= 21) Ladv += g * data[n]->p[UU][i][j][k] * ucon[1] * ucov[0] ;
      }
    }
  }

  // now copy primitives and four-vectors according to boundary conditions
  populate_boundary_conditions(n);

  dMact *= dx[3]*dx[2] ;
  dMact /= 21. ;
  Ladv *= dx[3]*dx[2] ;
  Ladv /= 21. ;

  Mdot_dump = -dMact*M_unit/T_unit;
  MdotEdd_dump = Mdotedd;
  Ladv_dump =  Ladv;

  double rescale_factor = 1.;

  if (target_mdot > 0) {
    fprintf(stderr, "Resetting M_unit to match target_mdot = %g ", target_mdot);

    double current_mdot = Mdot_dump/MdotEdd_dump;
    fprintf(stderr, "... is now %g\n", M_unit * fabs(target_mdot / current_mdot));
    rescale_factor = fabs(target_mdot / current_mdot);
    M_unit *= rescale_factor;

    set_units();
  }

  Mdot_dump = -dMact*M_unit/T_unit;
  MdotEdd_dump = Mdotedd;
  Ladv_dump =  Ladv;

  if (verbose == 2) {
    fprintf(stderr,"dMact: %g [code]\n",dMact);
    fprintf(stderr,"Ladv: %g [code]\n",Ladv_dump);
    fprintf(stderr,"Mdot: %g [g/s] \n",Mdot_dump);
    fprintf(stderr,"Mdot: %g [MSUN/YR] \n",Mdot_dump/(MSUN / YEAR));
    fprintf(stderr,"Mdot: %g [Mdotedd]\n",Mdot_dump/MdotEdd_dump);
    fprintf(stderr,"Mdotedd: %g [g/s]\n",MdotEdd_dump);
    fprintf(stderr,"Mdotedd: %g [MSUN/YR]\n",MdotEdd_dump/(MSUN/YEAR));
  } else if (verbose == 1) {
    fprintf(stderr,"Mdot: %g [MSUN/YR] \n",Mdot_dump/(MSUN / YEAR));
    fprintf(stderr,"Mdot: %g [Mdotedd]\n",Mdot_dump/MdotEdd_dump);
  }

  // now construct useful scalar quantities (over full (+ghost) zones of data)
  init_physical_quantities(n, rescale_factor);
}

void load_koral_data(int n, char *fnam, int dumpidx, int verbose)
{
  // loads relevant information from fluid dump file stored at fname
  // to the n'th copy of data (e.g., for slow light)

  double dMact, Ladv;

  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);

  nloaded++;

  if ( hdf5_open(fname) < 0 ) {
    fprintf(stderr, "! unable to open file %s. Exiting!\n", fname);
    exit(-1);
  }

  hdf5_set_directory("/");
  hdf5_read_single_val(&(data[n]->t), "t", H5T_IEEE_F64LE);

  hdf5_set_directory("/quants/");

  // load into "center" of data
  hsize_t fdims[] = { N1, N2, N3 };
  hsize_t fstart[] = { 0, 0, 0 };
  hsize_t fcount[] = { N1, N2, N3, 1 };
  hsize_t mdims[] = { N1+2, N2+2, N3+2 };
  hsize_t mstart[] = { 1, 1, 1 };

  hdf5_read_array(data[n]->p[KRHO][0][0], "rho", 3, fdims, fstart, fcount, 
                  mdims, mstart, H5T_IEEE_F64LE); 

  hdf5_read_array(data[n]->p[UU][0][0], "uint", 3, fdims, fstart, fcount, 
                  mdims, mstart, H5T_IEEE_F64LE); 

  hdf5_read_array(data[n]->p[U1][0][0], "U1", 3, fdims, fstart, fcount, 
                  mdims, mstart, H5T_IEEE_F64LE); 

  hdf5_read_array(data[n]->p[U2][0][0], "U2", 3, fdims, fstart, fcount, 
                  mdims, mstart, H5T_IEEE_F64LE); 

  hdf5_read_array(data[n]->p[U3][0][0], "U3", 3, fdims, fstart, fcount, 
                  mdims, mstart, H5T_IEEE_F64LE); 

  hdf5_read_array(data[n]->p[B1][0][0], "B1", 3, fdims, fstart, fcount, 
                  mdims, mstart, H5T_IEEE_F64LE); 

  hdf5_read_array(data[n]->p[B2][0][0], "B2", 3, fdims, fstart, fcount, 
                  mdims, mstart, H5T_IEEE_F64LE); 

  hdf5_read_array(data[n]->p[B3][0][0], "B3", 3, fdims, fstart, fcount, 
                  mdims, mstart, H5T_IEEE_F64LE); 

  if (ELECTRONS == 9) {
    hdf5_read_array(data[n]->p[TFLK][0][0], "te", 3, fdims, fstart, fcount, 
                    mdims, mstart, H5T_IEEE_F64LE); 
  }

  hdf5_close();

  dMact = Ladv = 0.;

  // construct four-vectors over "real" zones
#pragma omp parallel for collapse(2) reduction(+:dMact,Ladv)
  for(int i = 1; i < N1+1; i++) {
    for(int j = 1; j < N2+1; j++) {

      double X[NDIM] = { 0. };
      double gcov[NDIM][NDIM], gcon[NDIM][NDIM], gcov_KS[NDIM][NDIM], gcon_KS[NDIM][NDIM];
      double g, r, th;

      // this assumes axisymmetry in the coordinates
      ijktoX(i-1,j-1,0, X);
      gcov_func(X, gcov);
      gcon_func(gcov, gcon);
      g = gdet_zone(i-1,j-1,0);

      bl_coord(X, &r, &th);

      // the file is output in KS, so get KS metric
      gcov_ks(r, th, gcov_KS);
      gcon_func(gcov_KS, gcon_KS);

      for(int k = 1; k < N3+1; k++){

        ijktoX(i-1,j-1,k,X);
        double UdotU = 0.;
        
        for(int l = 1; l < NDIM; l++) 
          for(int m = 1; m < NDIM; m++) 
            UdotU += gcov_KS[l][m]*data[n]->p[U1+l-1][i][j][k]*data[n]->p[U1+m-1][i][j][k];
        double ufac = sqrt(-1./gcon_KS[0][0]*(1 + fabs(UdotU)));

        double ucon[NDIM] = { 0. };
        ucon[0] = -ufac * gcon_KS[0][0];

        for(int l = 1; l < NDIM; l++) 
          ucon[l] = data[n]->p[U1+l-1][i][j][k] - ufac*gcon_KS[0][l];

        double ucov[NDIM] = { 0. };
        flip_index(ucon, gcov_KS, ucov);

        // reconstruct the magnetic field three vectors
        double udotB = 0.;
        
        for (int l = 1; l < NDIM; l++) {
          udotB += ucov[l]*data[n]->p[B1+l-1][i][j][k];
        }
      
        double bcon[NDIM] = { 0. };
        double bcov[NDIM] = { 0. };

        bcon[0] = udotB;
        for (int l = 1; l < NDIM; l++) {
          bcon[l] = (data[n]->p[B1+l-1][i][j][k] + ucon[l]*udotB)/ucon[0];
        }
        flip_index(bcon, gcov_KS, bcov);

        double bsq = 0.;
        for (int l=0; l<NDIM; ++l) bsq += bcon[l] * bcov[l];
        data[n]->b[i][j][k] = sqrt(bsq) * B_unit;

        if(i <= 21) { dMact += g * data[n]->p[KRHO][i][j][k] * ucon[1]; }
        if(i >= 21 && i < 41 && 0) Ladv += g * data[n]->p[UU][i][j][k] * ucon[1] * ucov[0];
        if(i <= 21) Ladv += g * data[n]->p[UU][i][j][k] * ucon[1] * ucov[0];

        // trust ...
        //double udb1 = 0., udu1 = 0., bdb1 = 0.;
        //MULOOP { udb1 += ucon[mu]*bcov[mu]; udu1 += ucon[mu]*ucov[mu]; bdb1 += bcon[mu]*bcov[mu]; }

        // now translate from KS (outcoords) -> MKS:hslope=1
        ucon[1] /= r;
        ucon[2] /= M_PI;
        bcon[1] /= r;
        bcon[2] /= M_PI;

        // resynthesize the primitives. note we have assumed that B^i = *F^{ti}
        double alpha = sqrt(-1. / gcon[0][0]);
        data[n]->p[U1][i][j][k] = (gcon[0][1]*alpha*alpha + ucon[1]/ucon[0]) * ucon[0];
        data[n]->p[U2][i][j][k] = (gcon[0][2]*alpha*alpha + ucon[2]/ucon[0]) * ucon[0];
        data[n]->p[U3][i][j][k] = (gcon[0][3]*alpha*alpha + ucon[3]/ucon[0]) * ucon[0];
        data[n]->p[B1][i][j][k] = ucon[0] * bcon[1] - bcon[0] * ucon[1];
        data[n]->p[B2][i][j][k] = ucon[0] * bcon[2] - bcon[0] * ucon[2];
        data[n]->p[B3][i][j][k] = ucon[0] * bcon[3] - bcon[0] * ucon[3];

        // ... but verify
        //flip_index(bcon, gcov, bcov);
        //flip_index(ucon, gcov, ucov);
        //double udb2 = 0., udu2 = 0., bdb2 = 0.;
        //MULOOP { udb2 += ucon[mu]*bcov[mu]; udu2 += ucon[mu]*ucov[mu]; bdb2 += bcon[mu]*bcov[mu]; }
        //fprintf(stderr, "u.u %g %g   u.b %g %g   b.b %g %g\n", udu1,udu2, udb1,udb2, bdb1,bdb2);
      }
    }
  }

  // now copy primitives and four-vectors according to boundary conditions
  populate_boundary_conditions(n);

  dMact *= dx[3]*dx[2] ;
  dMact /= 21. ;
  Ladv *= dx[3]*dx[2] ;
  Ladv /= 21. ;

  Mdot_dump = -dMact*M_unit/T_unit;
  MdotEdd_dump = Mdotedd;
  Ladv_dump =  Ladv;

  double rescale_factor = 1.;

  if (target_mdot > 0) {
    fprintf(stderr, "Resetting M_unit to match target_mdot = %g ", target_mdot);

    double current_mdot = Mdot_dump/MdotEdd_dump;
    fprintf(stderr, "... is now %g\n", M_unit * fabs(target_mdot / current_mdot));
    rescale_factor = fabs(target_mdot / current_mdot);
    M_unit *= rescale_factor;

    set_units();
  }

  Mdot_dump = -dMact*M_unit/T_unit;
  MdotEdd_dump = Mdotedd;
  Ladv_dump =  Ladv;

  if (verbose == 2) {
    fprintf(stderr,"dMact: %g [code]\n",dMact);
    fprintf(stderr,"Ladv: %g [code]\n",Ladv_dump);
    fprintf(stderr,"Mdot: %g [g/s] \n",Mdot_dump);
    fprintf(stderr,"Mdot: %g [MSUN/YR] \n",Mdot_dump/(MSUN / YEAR));
    fprintf(stderr,"Mdot: %g [Mdotedd]\n",Mdot_dump/MdotEdd_dump);
    fprintf(stderr,"Mdotedd: %g [g/s]\n",MdotEdd_dump);
    fprintf(stderr,"Mdotedd: %g [MSUN/YR]\n",MdotEdd_dump/(MSUN/YEAR));
  } else if (verbose == 1) {
    fprintf(stderr,"Mdot: %g [MSUN/YR] \n",Mdot_dump/(MSUN / YEAR));
    fprintf(stderr,"Mdot: %g [Mdotedd]\n",Mdot_dump/MdotEdd_dump);
  }

  // now construct useful scalar quantities (over full (+ghost) zones of data)
  init_physical_quantities(n, rescale_factor);
}

// get dMact in the i'th radial zone (0 = 0 of the dump, so ignore ghost zones)
double get_code_dMact(int i, int n)
{
  i += 1;
  double dMact = 0;
#pragma omp parallel for collapse(1) reduction(+:dMact)
  for (int j=1; j<N2+1; ++j) {
    double X[NDIM] = { 0. };
    double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
    double g, r, th;

    // this assumes axisymmetry in the coordinates
    ijktoX(i-1,j-1,0, X);
    gcov_func(X, gcov);
    gcon_func(gcov, gcon);
    g = gdet_zone(i-1,j-1,0);
    bl_coord(X, &r, &th);

    for (int k=1; k<N3+1; ++k) {
      ijktoX(i-1,j-1,k,X);
      double UdotU = 0;

      for(int l = 1; l < NDIM; l++)
        for(int m = 1; m < NDIM; m++)
          UdotU += gcov[l][m]*data[n]->p[U1+l-1][i][j][k]*data[n]->p[U1+m-1][i][j][k];
      double ufac = sqrt(-1./gcon[0][0]*(1 + fabs(UdotU)));

      double ucon[NDIM] = { 0. };
      ucon[0] = -ufac * gcon[0][0];

      for(int l = 1; l < NDIM; l++)
        ucon[l] = data[n]->p[U1+l-1][i][j][k] - ufac*gcon[0][l];

      double ucov[NDIM] = { 0. };
      flip_index(ucon, gcov, ucov);

      dMact += g * data[n]->p[KRHO][i][j][k] * ucon[1];
    }

  }
  return dMact;
}

void load_iharm_data(int n, char *fnam, int dumpidx, int verbose)
{
  // loads relevant information from fluid dump file stored at fname
  // to the n'th copy of data (e.g., for slow light)

  double dMact, Ladv;

  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);

  nloaded++;

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

  if (ELECTRONS == 1) {
    fstart[3] = 8;
    hdf5_read_array(data[n]->p[KEL][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
    fstart[3] = 9;
    hdf5_read_array(data[n]->p[KTOT][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
  }

  //Reversing B Field
  if(reverse_field) {
    double multiplier = -1.0;
    for(int i=0;i<N1+2;i++){
      for(int j=0;j<N2+2;j++){
        for(int k=0;k<N3+2;k++){ 
          data[n]->p[B1][i][j][k] = multiplier*data[n]->p[B1][i][j][k];
          data[n]->p[B2][i][j][k] = multiplier*data[n]->p[B2][i][j][k];
          data[n]->p[B3][i][j][k] = multiplier*data[n]->p[B3][i][j][k];
        }
      }
    }
  }
  hdf5_read_single_val(&(data[n]->t), "t", H5T_IEEE_F64LE);

  if (ELECTRONS == ELECTRONS_TFLUID) {
    fstart[3] = 8;
    hdf5_read_array(data[n]->p[THF][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
  }

  hdf5_close();

  dMact = Ladv = 0.;

  // construct four-vectors over "real" zones
#pragma omp parallel for collapse(2) reduction(+:dMact,Ladv)
  for(int i = 1; i < N1+1; i++) {
    for(int j = 1; j < N2+1; j++) {

      double X[NDIM] = { 0. };
      double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
      double g, r, th;

      // this assumes axisymmetry in the coordinates
      ijktoX(i-1,j-1,0, X);
      gcov_func(X, gcov);
      gcon_func(gcov, gcon);
      g = gdet_zone(i-1,j-1,0);

      bl_coord(X, &r, &th);

      for(int k = 1; k < N3+1; k++){

        ijktoX(i-1,j-1,k,X);
        double UdotU = 0.;

        // the four-vector reconstruction should have gcov and gcon and gdet using the modified coordinates
        // interpolating the four vectors to the zone center !!!!
        for(int l = 1; l < NDIM; l++) 
          for(int m = 1; m < NDIM; m++) 
            UdotU += gcov[l][m]*data[n]->p[U1+l-1][i][j][k]*data[n]->p[U1+m-1][i][j][k];
        double ufac = sqrt(-1./gcon[0][0]*(1 + fabs(UdotU)));

        // TODO: move this above
        double ucon[NDIM] = { 0. };
        ucon[0] = -ufac * gcon[0][0];

        for(int l = 1; l < NDIM; l++) 
          ucon[l] = data[n]->p[U1+l-1][i][j][k] - ufac*gcon[0][l];

        double ucov[NDIM] = { 0. };
        flip_index(ucon, gcov, ucov);

        // reconstruct the magnetic field three vectors
        double udotB = 0.;

        for (int l = 1; l < NDIM; l++) {
          udotB += ucov[l]*data[n]->p[B1+l-1][i][j][k];
        }

        double bcon[NDIM] = { 0. };
        double bcov[NDIM] = { 0. };

        bcon[0] = udotB;
        for (int l = 1; l < NDIM; l++) {
          bcon[l] = (data[n]->p[B1+l-1][i][j][k] + ucon[l]*udotB)/ucon[0];
        }
        flip_index(bcon, gcov, bcov);

        double bsq = 0.;
        for (int l=0; l<NDIM; ++l) bsq += bcon[l] * bcov[l];
        data[n]->b[i][j][k] = sqrt(bsq) * B_unit;


        if(i <= 21) { dMact += g * data[n]->p[KRHO][i][j][k] * ucon[1]; }
        if(i >= 21 && i < 41 && 0) Ladv += g * data[n]->p[UU][i][j][k] * ucon[1] * ucov[0];
        if(i <= 21) Ladv += g * data[n]->p[UU][i][j][k] * ucon[1] * ucov[0];

      }
    }
  }

  // check if 21st zone (i.e., 20th without ghost zones) is beyond r_eh
  // otherwise recompute
  if (1==1) {
    double r_eh = 1. + sqrt(1. - a*a);
    int N2_by_2 = (int)(N2/2);

    double X[NDIM] = { 0. };
    ijktoX(20, N2_by_2, 0, X);

    double r, th;
    bl_coord(X, &r, &th);

    if (r < r_eh) {
      for (int i=1; i<N1+1; ++i) {
        ijktoX(i, N2_by_2, 0, X);
        bl_coord(X, &r, &th);
        if (r >= r_eh) {
          fprintf(stderr, "r_eh is beyond regular zones. recomputing at %g...\n", r);
          dMact = get_code_dMact(i, n) * 21;
          break;
        }
      }
    }
  }

  // now copy primitives and four-vectors according to boundary conditions
  populate_boundary_conditions(n);

  dMact *= dx[3]*dx[2] ;
  dMact /= 21. ;
  Ladv *= dx[3]*dx[2] ;
  Ladv /= 21. ;

  Mdot_dump = -dMact*M_unit/T_unit;
  MdotEdd_dump = Mdotedd;
  Ladv_dump =  Ladv;

  double rescale_factor = 1.;

  if (target_mdot > 0) {
    fprintf(stderr, "Resetting M_unit to match target_mdot = %g ", target_mdot);

    double current_mdot = Mdot_dump/MdotEdd_dump;
    fprintf(stderr, "... is now %g\n", M_unit * fabs(target_mdot / current_mdot));
    rescale_factor = fabs(target_mdot / current_mdot);
    M_unit *= rescale_factor;

    set_units();
  }

  Mdot_dump = -dMact*M_unit/T_unit;
  MdotEdd_dump = Mdotedd;
  Ladv_dump =  Ladv;

  if (verbose == 2) {
    fprintf(stderr,"dMact: %g [code]\n",dMact);
    fprintf(stderr,"Ladv: %g [code]\n",Ladv_dump);
    fprintf(stderr,"Mdot: %g [g/s] \n",Mdot_dump);
    fprintf(stderr,"Mdot: %g [MSUN/YR] \n",Mdot_dump/(MSUN / YEAR));
    fprintf(stderr,"Mdot: %g [Mdotedd]\n",Mdot_dump/MdotEdd_dump);
    fprintf(stderr,"Mdotedd: %g [g/s]\n",MdotEdd_dump);
    fprintf(stderr,"Mdotedd: %g [MSUN/YR]\n",MdotEdd_dump/(MSUN/YEAR));
  } else if (verbose == 1) {
    fprintf(stderr,"Mdot: %g [MSUN/YR] \n",Mdot_dump/(MSUN / YEAR));
    fprintf(stderr,"Mdot: %g [Mdotedd]\n",Mdot_dump/MdotEdd_dump);
  }

  // now construct useful scalar quantities (over full (+ghost) zones of data)
  init_physical_quantities(n, rescale_factor);
}


int radiating_region(double X[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);
  return (r > rmin_geo && r < rmax_geo && th > th_beg && th < (M_PI-th_beg));
}

// In case we want to mess with emissivities directly
void get_model_jar(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV) {return;}
void get_model_jk(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv) {return;}
