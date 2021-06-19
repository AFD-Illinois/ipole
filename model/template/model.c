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

void try_set_model_parameter(const char *word, const char *value)
{
  // set_by_word_val(word, value, "parameter_name", &parameter_variable, TYPE_DBL);
}

void init_model(double *tA, double *tB)
{
  // set nice numbers here
  *tA = 0.;
  *tB = 1.;
  simulation_rout = 100.;

  // set metric
  use_eKS_internal = 0;
  metric = METRIC_MKS;

  // set units
  double MBH = 4.e6 * MSUN;
  L_unit = GNEWT * MBH / (CL * CL);
  
  // set up coordinate system. may be useful if you want to cut emisison outside
  // of some startx, stopx range using radiating_region(...)

  // not really used in ipole other than for coordinates
  N1 = 128;
  N2 = 128;
  N3 = 128;

  // set up coordinates
  startx[0] = 0.;
  startx[1] = 0.;
  startx[2] = 0.;
  startx[3] = 0.;

  dx[0] = 0.1;
  dx[1] = ( log(simulation_rout) - log(1) ) / N1;
  dx[2] = 1. / N2;
  dx[3] = 2. * M_PI / N3;
  
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

  Bcon[0] = 0.;
  Bcon[1] = 0.;
  Bcon[2] = 0.;
  Bcon[3] = 0.;

  lower(Ucon, gcov, Ucov);
  lower(Bcon, gcov, Bcov);
}

// dimensionless electron temperature
double get_model_thetae(double X[NDIM])
{
  if (radiating_region(X) == 0) {
    return 0;
  }

  return 0.;
}

// magnetic field strength in Gauss
double get_model_b(double X[NDIM])
{
  if (radiating_region(X) == 0) {
    return 0.;
  }

  return 0.;
}

// electron number density in cgs
double get_model_ne(double X[NDIM])
{
  if (radiating_region(X) == 0) {
    return 0.;
  }

  return 0.;
}

// called when writing a new image file
void output_hdf5()
{
  hdf5_set_directory("/");

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

