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

// used by other files but defined in model.c for some reason
double rmax_geo;
double L_unit;

// often we want to stop computing coefficients beyond some radius
double simulation_rout;

// local variables. custom to this model
static char fnam[STRLEN] = "dump.vtu";

static int DO_RENORMALIZE_SRP = 1;
static double MBH_solar = 4.e6;
static double Ne_max = 1.e6;
static double Thetae_max = 10.;
static double sigma = 0.5;
static double H0 = 0.01;
static double adiabatic_gamma = 5./3;

static double **midplane_Srp = NULL;
static double **midplane_v1 = NULL;
static double **midplane_v2 = NULL;
static double **midplane_v3 = NULL;

// defined at bottom
int radiating_region_rh(double r, double h); 

// convert from spherical KS to cylindrical t,rho,z,phi coordiantes
// ... some choices to be made here. using simplest option.
double get_coordinate_z(double r, double h)
{
  return r * cos(h);  
}
double get_coordinate_rho(double r, double h)
{
  return r * sin(h);
}

void try_set_model_parameter(const char *word, const char *value)
{
  set_by_word_val(word, value, "dump", (void *)fnam, TYPE_STR);

  set_by_word_val(word, value, "RENORMALIZE", &DO_RENORMALIZE_SRP, TYPE_INT);

  set_by_word_val(word, value, "MBH", &MBH_solar, TYPE_DBL);
  set_by_word_val(word, value, "a", &a, TYPE_DBL);

  set_by_word_val(word, value, "Ne_max", &Ne_max, TYPE_DBL);
  set_by_word_val(word, value, "Thetae_max", &Thetae_max, TYPE_DBL);
  set_by_word_val(word, value, "sigma", &sigma, TYPE_DBL);

  set_by_word_val(word, value, "H0", &H0, TYPE_DBL);
  set_by_word_val(word, value, "gamma", &adiabatic_gamma, TYPE_DBL);
}

void populate_boundary_conditions() 
{
  // copy over first and last radial zones
  for (int k=1; k<N3+1; ++k) {

    midplane_Srp[0][k] = midplane_Srp[1][k];
    midplane_Srp[N1-1][k] = midplane_Srp[N1-2][k];

    midplane_v1[0][k] = midplane_v1[1][k];
    midplane_v1[N1-1][k] = midplane_v1[N1-2][k];

    midplane_v2[0][k] = midplane_v2[1][k];
    midplane_v2[N1-1][k] = midplane_v2[N1-2][k];

    midplane_v3[0][k] = midplane_v3[1][k];
    midplane_v3[N1-1][k] = midplane_v3[N1-2][k];

  }

  // now copy over in azimuth
  for (int i=0; i<N1+2; ++i) {
    midplane_Srp[i][0] = midplane_Srp[i][1];
    midplane_Srp[i][N3-1] = midplane_Srp[i][N3-2];

    midplane_v1[i][0] = midplane_v1[i][1];
    midplane_v1[i][N3-1] = midplane_v1[i][N3-2];

    midplane_v2[i][0] = midplane_v2[i][1];
    midplane_v2[i][N3-1] = midplane_v2[i][N3-2];

    midplane_v3[i][0] = midplane_v3[i][1];
    midplane_v3[i][N3-1] = midplane_v3[i][N3-2];
  }
}

int get_N1_N3(const char *fnam, int *N1, int *N3)
{
  FILE * fp;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;

  fp = fopen(fnam, "r");
  if (fp == NULL) {
    fprintf(stderr, "! unable to open dump %s. quitting.\n", fnam);
    exit(3);
  }

  double Phi, tr, tPhi; 

  int total_zones = 0;
  int line_counter = 0;
  while ((read = getline(&line, &len, fp)) != -1) {
    line_counter += 1;

    // ignore comments and the first line
    if (line[0] == '#') continue;
    if (line_counter < 2) continue;
    
    // get total number of zones
    if (total_zones < 1) {
      sscanf(line, "%d", &total_zones);
      continue;
    }

    // load data
    sscanf(line, "%lg %lg ", &tr, &tPhi);
    
    if (line_counter == 3) {
      Phi = tPhi;
    } else {
      if (Phi != tPhi) {
        fclose(fp);
        *N1 = line_counter - 3;
        *N3 = total_zones / (*N1);
        return 0;
      }
    }
  }

  fclose(fp);
  return -1;
}

void init_model(double *tA, double *tB)
{
  // set nice numbers here
  *tA = 0.;
  *tB = 1.;

  // set metric
  use_eKS_internal = 0;
  metric = METRIC_MKS;
  hslope = 1.0;

  // get grid dimension
  int rv = get_N1_N3(fnam, &N1, &N3);
  N2 = 1;
  if (rv < 0) {
    fprintf(stderr, "! unable to get file dimensions.\n");
    exit(3);
  }
  fprintf(stderr, "file has dimensions: %d x %d\n", N1, N3);

  // initialize storage
  midplane_Srp = malloc_rank2(N1+2, N3+2);
  midplane_v1 = malloc_rank2(N1+2, N3+2);
  midplane_v2 = malloc_rank2(N1+2, N3+2);
  midplane_v3 = malloc_rank2(N1+2, N3+2);

  // load data
  int zoffset = 3;  // zones start on the zoffset'th line
  int line_counter = 0;

  FILE * fp;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;

  fp = fopen(fnam, "r");
  if (fp == NULL) {
    fprintf(stderr, "! unable to open dump %s. quitting.\n", fnam);
    exit(3);
  }

  int ti, tk;
  double tr, tPhi, trho, tflrho1, tflrho2, tv1, tv2, tv3, tpp, tlfac, tp;
  double max_Srp = 0;

  double *rslice = malloc(N1 * sizeof(*rslice));
  double *pslice = malloc(N3 * sizeof(*pslice));

  while ((read = getline(&line, &len, fp)) != -1) {
    line_counter += 1;

    // ignore comments and the first line
    if (line[0] == '#') continue;
    if (line_counter < 3) continue;

    // load data
    sscanf(line, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", 
           &tr, &tPhi, &trho, &tflrho1, &tflrho2,
           &tv1, &tv2, &tv3, &tpp, &tlfac, &tp);
    
    ti = (line_counter-zoffset) % N1;
    tk = (line_counter-zoffset) / N1;

    if (ti == 0) pslice[tk] = tPhi;
    if (tk == 0) rslice[ti] = tr;

    double Srp = trho / H0 / rslice[ti];
    if (max_Srp < Srp) max_Srp = Srp;

    midplane_Srp[ti+1][tk+1] = Srp;
    midplane_v1[ti+1][tk+1] = tv1;
    midplane_v2[ti+1][tk+1] = tv2;
    midplane_v3[ti+1][tk+1] = tv3;
  }

  fclose(fp);

  if (line) free(line);

  if (DO_RENORMALIZE_SRP) {
    for (int i=0; i<N1; ++i) {
      for (int k=0; k<N3; ++k) {
        midplane_Srp[i+1][k+1] /= max_Srp;
      }
    }
  }

  populate_boundary_conditions();

  // set units
  double MBH = MBH_solar * MSUN;
  L_unit = GNEWT * MBH / (CL * CL);
  
  // set up coordinates
  double dx1 = ( log(rslice[N1-1]) - log(rslice[0]) ) / ( N1 - 1 );
  simulation_rout = rslice[N1-1];
  rmax_geo = fmax(simulation_rout, 1000.);
  Rh = 1. + sqrt(1. - a*a);

  dx[0] = 0.1;
  dx[1] = dx1;
  dx[2] = 1. / N2;
  dx[3] = 2. * M_PI / N3;
 
  startx[0] = 0.;
  startx[1] = log(rslice[0]) - dx1/2.;
  startx[2] = 0.;
  startx[3] = 0.;
 
  stopx[0] = 1.;
  stopx[1] = startx[1]+N1*dx[1];
  stopx[2] = startx[2]+N2*dx[2];
  stopx[3] = startx[3]+N3*dx[3];

  MULOOP {
    cstartx[mu] = startx[mu];
    cstopx[mu] = stopx[mu];
  }

}

// called by ipole during ray tracing, should return
// Ucon, Ucov as regular and Bcon, Bcov. for safety,
// it makes sense to units of Bcon,Bcov such that we
// have Bcon.Bcov = bsq in Gauss^2
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

  double model_b = get_model_b(X);

  Bcon[0] = 0.;
  Bcon[1] = model_b * cos(h) / r;
  Bcon[2] = - model_b * sin(h) / (r + 1.e-8);
  Bcon[3] = 0.;

  lower(Ucon, gcov, Ucov);
  lower(Bcon, gcov, Bcov);
}

// dimensionless electron temperature
double get_model_thetae(double X[NDIM])
{
  double r, h;
  bl_coord(X, &r, &h);

  if (radiating_region_rh(r, h) == 0) {
    return 0;
  }

  double ne = get_model_ne(X);

  return sqrt(4. * M_PI * sigma * (MP+ME) * CL*CL * ne);
}

// magnetic field strength in Gauss
double get_model_b(double X[NDIM])
{
  double r, h;
  bl_coord(X, &r, &h);

  if (radiating_region_rh(r, h) == 0) {
    return 0.;
  }

  double z = get_coordinate_z(r, h);
  double zhsq = z*z / (H0*H0*r*r);

  // get Srp in midplane
  double Xp[NDIM] = { X[0], log(get_coordinate_rho(r, h)), X[2], X[3] };
  double tsrp = interp_scalar_2d(Xp, midplane_Srp);

  return Thetae_max * pow(tsrp, adiabatic_gamma - 1.) * (1. - zhsq);
}

// electron number density in cgs
double get_model_ne(double X[NDIM])
{
  double r, h;
  bl_coord(X, &r, &h);

  if (radiating_region_rh(r, h) == 0) {
    return 0.;
  }

  double z = get_coordinate_z(r, h);
  double zhsq = z*z / (H0*H0*r*r);

  // get Srp in midplane
  double Xp[NDIM] = { X[0], log(get_coordinate_rho(r, h)), X[2], X[3] };
  double tsrp = interp_scalar_2d(Xp, midplane_Srp);

  return Ne_max * tsrp * pow(1. - zhsq, 1. / (adiabatic_gamma - 1.));
}

// called when writing a new image file
void output_hdf5()
{
  hdf5_set_directory("/header/");
  
  hdf5_write_single_val(&Ne_max, "Ne_max", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Thetae_max, "Thetae_max", H5T_IEEE_F64LE);
  hdf5_write_single_val(&sigma, "sigma", H5T_IEEE_F64LE);

  hdf5_write_single_val(&H0, "H0", H5T_IEEE_F64LE);
  hdf5_write_single_val(&adiabatic_gamma, "gamma", H5T_IEEE_F64LE);

  hdf5_write_single_val(&N1, "N1", H5T_STD_I32LE);
  hdf5_write_single_val(&N3, "N3", H5T_STD_I32LE);

  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  hdf5_write_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
  
  hdf5_set_directory("/");
}

// used in ipolarray.c to turn off emission in "invalid" parts
// of the domain (e.g., if we shouldn't try to map from the KS
// coordinates into the grid)
int radiating_region(double X[NDIM])
{
  double r, h;
  bl_coord(X, &r, &h);
  return radiating_region_rh(r, h);
}

int radiating_region_rh(double r, double h)
{
  double z = get_coordinate_z(r, h);
  if (z*z >= H0*H0*r*r) return 0; 
  return ( r <= simulation_rout ) ? 1 : 0;
}

// used in slow light with real data
void update_data_until(double *tA, double *tB, double tgt) { }
void update_data(double *tA, double *tB) { }
double get_dump_t(char *fnam, int dumpidx) { return 0.; }

// not necessary to implement unless you know you want it
void get_model_primitives(double X[NDIM], double *p) { }
void get_model_powerlaw_vals(double X[NDIM], double *p, double *n, double *gamma_min, double *gamma_max, double *gamma_cut) { }

// in case we want to mess with emissivities directly
void get_model_jar(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV) {return;}
void get_model_jk(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv) {return;}

