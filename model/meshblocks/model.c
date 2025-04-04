#include "model.h"

#include "decs.h"
#include "hdf5_utils.h"
#include "debug_tools.h"

#include "dict.h"
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
#include <ctype.h>

// macros
#define NSUP (3)  // how many files to load for slow light tracing
#define USE_GEODESIC_SIGMACUT (1)

#define P_RHO (0)
#define P_U1 (1)
#define P_U2 (2)
#define P_U3 (3)
#define P_UU (4)
#define P_B1 (5)
#define P_B2 (6)
#define P_B3 (7)

// used by other files
double L_unit;

/////////////////////
// used internally //
/////////////////////

// used internally but for scaling / intermediates
double MBH = 0.;
double MBH_solar = 0.;
double Mdotedd = 0.;
double M_unit = 0.;
double T_unit = 0.;
double B_unit = 0.;
double U_unit = 0.;
double RHO_unit = 0.;

// future-looking file format support
#define FORMAT_NONE (0)
#define FORMAT_ATHENAK_1P1 (1)
int dumpfile_format = FORMAT_NONE;

// thermodynamics
static double adiabatic_gamma = 13./9;
static double beta_crit = 1.0;
static double trat_small = 1.0;
static double trat_large = 40.0;
static double game = 4./3;
static double gamp = 5./3;

// geodesic parameters
static double polar_cut = -1.0;
double sigma_cut = 1.0;
static double sigma_cut_high = -1.0;

double tf = -1.;

// file/model data
dict *model_params = NULL;
// TODO: come up with way to save this

// for loading files
static int nloaded = 0;
static char fnam[STRLEN] = "dump.h5";
static int dumpidx;
static int dumpmin;
static int dumpmax;
static int dumpskip;
static double DTd = 0.;
static double last_time_loaded = 0.;

// geometry for the mesh. assumes static
double *meshblock_geometry;
int *meshblock_levels;
int last_mb = -1;
#pragma omp threadprivate(last_mb)

// data structures for the mesh
int n_variables = 0;
size_t n_meshblocks = 0;
size_t mb_nx1 = 0;
size_t mb_nx2 = 0;
size_t mb_nx3 = 0;
struct of_data {
  double t;
  double *****p;  // mbi, nvar, mb_nx1, mb_nx2, mb_nx3
  double ****ne;  // mbi, mb_nx1, mb_nx2, mb_nx3
  double ****thetae;
  double ****b;
  double ****sigma;
  double ****beta;
};
static struct of_data dataA, dataB, dataC;
static struct of_data *data[NSUP];

void try_set_model_parameter(const char *word, const char *value)
{
  set_by_word_val(word, value, "dump", (void *)fnam, TYPE_STR);

  // slow light
  set_by_word_val(word, value, "dump_min", &dumpmin, TYPE_INT);
  set_by_word_val(word, value, "dump_max", &dumpmax, TYPE_INT);
  set_by_word_val(word, value, "dump_skip", &dumpskip, TYPE_INT);
  dumpidx = dumpmin;

  // geodesic parameters
  set_by_word_val(word, value, "rmax_geo", &rmax_geo, TYPE_DBL);
  set_by_word_val(word, value, "rmin_geo", &rmin_geo, TYPE_DBL);
  set_by_word_val(word, value, "sigma_cut", &sigma_cut, TYPE_DBL);
  set_by_word_val(word, value, "sigma_cut_high", &sigma_cut_high, TYPE_DBL);
  set_by_word_val(word, value, "polar_cut_deg", &polar_cut, TYPE_DBL);

  // units
  set_by_word_val(word, value, "MBH", &MBH_solar, TYPE_DBL);
  set_by_word_val(word, value, "M_unit", &M_unit, TYPE_DBL);

  // thermodynamics
  set_by_word_val(word, value, "trat_small", &trat_small, TYPE_DBL);
  set_by_word_val(word, value, "trat_large", &trat_large, TYPE_DBL);
  set_by_word_val(word, value, "beta_crit", &beta_crit, TYPE_DBL);

  // for slow light
  set_by_word_val(word, value, "dump_min", &dumpmin, TYPE_INT);
  set_by_word_val(word, value, "dump_max", &dumpmax, TYPE_INT);
  set_by_word_val(word, value, "dump_skip", &dumpskip, TYPE_INT);
}

// advance through dumps until we are closer to the next set
// of dumps corresponding to tA == tgt. Used when attempting
// to restart from slowlight restart file.
void update_data_until(double *tA, double *tB, double tgt) {
  double tC = data[2]->t;
  while (tC < tgt) {
    dumpidx += dumpskip;
    tC = get_dump_t(fnam, dumpidx);
  }
  // reset dump index, just to be safe, then load on out ...
  dumpidx -= dumpskip;
  while (*tA < tgt) update_data(tA, tB);
}

// use internal dumpidx variable to load the next "expected"
// dump into memory (used for slowlight mode). After calling
// this function, it is guaranteed that data are ordered:
//
//   data[0]->t < data[1]->t < data[2]->t
//
// this function uses dataA, dataB, dataC to save the actual
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
    load_athenak_data(2, fnam, --nextdumpidx, 0);
    data[2]->t += 1.;
  } else {
    load_athenak_data(2, fnam, nextdumpidx, 0);
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

double get_athenak_dump_time(char *fname)
{
  double time = -1.0;

  char *res = NULL;
  char *line = NULL;

  FILE *fp;
  fp = fopen(fname, "rb");

  fseek(fp, 0, SEEK_SET);
  if (read_line(fp, &line, "! unable to read file header.\n") == -1) return -1;
  if (strcmp(line, "Athena binary output version=1.1\n") != 0) {
    fprintf(stderr, "! got bad file header: %s\n", line);
    return -1;
  }

  if (read_line(fp, &line, "! unable to read header line.\n") == -1) return -1;
  get_line_token(line, "=", 1, &res);
  int preheader_size = atoi(res);
  for (int i=0; i<preheader_size-1; ++i) {
    char *key = NULL;
    char *value = NULL;
    read_line(fp, &line, "! unable to read preheader line\n");
    get_line_token(line, "=", 1, &value);
    get_line_token(line, "=", 0, &key);
    strim(key);
    strim(value);
    if (strcmp(key, "time") == 0) {
      time = atof(value);
    }
  }

  fclose(fp);

  return time;
}

double get_dump_t(char *fnam, int dumpidx) 
{
  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);

  return get_athenak_dump_time(fname);
}

// right now this only supports "athenak_v1p1" but could in
// principle be extended to support other meshblock formats
void set_dumpfile_type(char *fnam, int index)
{
  char fname[256];
  snprintf(fname, 255, fnam, index);

  dumpfile_format = FORMAT_ATHENAK_1P1;
  fprintf(stderr, "set to athenak (v1.1)\n");
}

void strim(char *s)
{
  char *p = s;
  int l = strlen(p);
  while(isspace(p[l-1])) p[--l] = 0;
  while(*p && isspace(*p)) ++p, --l;
  memmove(s, p, l+1);
}

void get_line_token(char *line, char *delimiter, int index, char **result)
{
  char *save = NULL;
  *result = strtok_r(line, delimiter, &save);
  if (index-- == 0) return;
  while (*result != NULL) {
    *result = strtok_r(NULL, delimiter, &save);
    if (index-- == 0) return;
  }
}

int read_line(FILE *fp, char **line, char *message)
{
  int read;
  size_t line_length;
  if ((read = getline(line, &line_length, fp)) == -1) {
    fprintf(stderr, message);
    return -1;
  }
  return 0;
}

int read_int(FILE *fp, size_t bytes)
{
  // Warning: this is not safe across different sizeof(int)
  char values[bytes];
  if (fread(values, 1, bytes, fp)) ;
  int ival = 0;
  if (bytes == 4) {
    memcpy(&ival, values, sizeof(int));
  } else {
    fprintf(stderr, "unsupported integer size\n");
    exit(9);
  }
  return ival;
}

double read_double(FILE *fp, size_t bytes)
{
  char values[bytes];
  if (fread(values, 1, bytes, fp)) ;
  double dval = 0.;
  if (bytes == 8) {
    memcpy(&dval, values, sizeof(double));
  } else if (bytes == 4) {
    float fval;
    memcpy(&fval, values, sizeof(float));
    dval = fval;
  }
  return dval;
}

int is_within_meshblock(size_t mb, double x1, double x2, double x3)
{
  if (x1 < meshblock_geometry[6*mb+0]) return 0;
  if (x1 >= meshblock_geometry[6*mb+1]) return 0;
  if (x2 < meshblock_geometry[6*mb+2]) return 0;
  if (x2 >= meshblock_geometry[6*mb+3]) return 0;
  if (x3 < meshblock_geometry[6*mb+4]) return 0;
  if (x3 >= meshblock_geometry[6*mb+5]) return 0;
  return 1;
}

// assumes x1, x2, x3 in the same coordinates as the logical meshblock grid
int get_meshblock(double x1, double x2, double x3)
{
  if (is_within_meshblock(last_mb, x1, x2, x3)) {
    return last_mb;
  }

  for (int mb=0; mb<n_meshblocks; ++mb) {
    if (is_within_meshblock(mb, x1, x2, x3)) {
      last_mb = mb;
      return mb;
    }
  }

  return -1;
}

void set_units()
{
  MBH = MBH_solar * MSUN;
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

void set_meshblock_indices(double X[NDIM], int *mbn, int *i, int *j, int *k, double *di, double *dj, double *dk)
{
  double x, y, z;
  sks_to_cks(exp(X[1]), X[2], X[3], &x, &y, &z);
  int mb = get_meshblock(x, y, z);
  *mbn = mb;

  double xmin = meshblock_geometry[6*mb + 0];
  double xmax = meshblock_geometry[6*mb + 1];
  double ymin = meshblock_geometry[6*mb + 2];
  double ymax = meshblock_geometry[6*mb + 3];
  double zmin = meshblock_geometry[6*mb + 4];
  double zmax = meshblock_geometry[6*mb + 5];

  double dx = (xmax - xmin) / mb_nx1;
  double dy = (ymax - ymin) / mb_nx2;
  double dz = (zmax - zmin) / mb_nx3;
  
  // shift to zone centers because that's where variables are most exact.
  *i = (int)((x - xmin) / dx - 0.5 + 1000) - 1000;
  *j = (int)((y - ymin) / dy - 0.5 + 1000) - 1000;
  *k = (int)((z - zmin) / dz - 0.5 + 1000) - 1000;

  // now construct del
  *di = (x - ((*i + 0.5) * dx + xmin)) / dx;
  *dj = (y - ((*j + 0.5) * dy + ymin)) / dy;
  *dk = (z - ((*k + 0.5) * dz + zmin)) / dz;
}

double mb_interp_scalar(int mb, double x, double y, double z, double ****var)
{
  double xmin = meshblock_geometry[6*mb + 0];
  double xmax = meshblock_geometry[6*mb + 1];
  double ymin = meshblock_geometry[6*mb + 2];
  double ymax = meshblock_geometry[6*mb + 3];
  double zmin = meshblock_geometry[6*mb + 4];
  double zmax = meshblock_geometry[6*mb + 5];

  double dx = (xmax - xmin) / mb_nx1;
  double dy = (ymax - ymin) / mb_nx2;
  double dz = (zmax - zmin) / mb_nx3;

  // shift to zone centers because that's where variables are most exact.
  int i = (int)((x - xmin) / dx - 0.5 + 1000) - 1000;
  int j = (int)((y - ymin) / dy - 0.5 + 1000) - 1000;
  int k = (int)((z - zmin) / dz - 0.5 + 1000) - 1000;

  // now construct del
  double del1 = (x - ((i + 0.5) * dx + xmin)) / dx;
  double del2 = (y - ((j + 0.5) * dy + ymin)) / dy;
  double del3 = (z - ((k + 0.5) * dz + zmin)) / dz;

  // since we read from data, adjust i,j,k for ghost zones
  i += 1;
  j += 1;
  k += 1;

  int ip1 = i+1;
  int jp1 = j+1;
  int kp1 = k+1;

  double b1 = 1.-del1;
  double b2 = 1.-del2;

  // interpolate in x1 and x2
  double interp = var[mb][i  ][j  ][k]*b1*b2 +
                  var[mb][i  ][jp1][k]*b1*del2 +
                  var[mb][ip1][j  ][k]*del1*b2 +
                  var[mb][ip1][jp1][k]*del1*del2;

  // then interpolate in x3
  interp = (1.-del3)*interp +
  del3*(var[mb][i  ][j  ][kp1]*b1*b2 +
        var[mb][i  ][jp1][kp1]*b1*del2 +
        var[mb][ip1][j  ][kp1]*del1*b2 +
        var[mb][ip1][jp1][kp1]*del1*del2);

  return interp;
}

double mb_interp_scalar_time(int mb, double x, double y, double z, double ****varA, double ****varB, double tfac)
{
  double vA = mb_interp_scalar(mb, x, y, z, varA);

#if SLOWLIGHT
  double vB = mb_interp_scalar(mb, x, y, z, varB);
  return tfac*vA + (1. - tfac)*vB;
#endif

  return vA;
}

void read_set_athenak_header(FILE *fp, dict *model_params)
{
  char *line = NULL;
  char *res = NULL;

  if (read_line(fp, &line, "! unable to get header size.\n") == -1) return;
  get_line_token(line, "=", 1, &res);
  size_t header_offset = atoi(res);

  char header[header_offset+10];
  memset(header, 0, header_offset+10);
  if (fread(header, 1, header_offset, fp)) ;

  char group[128];
  char keyname[256];

  char *save = NULL;
  char *p = strtok_r(header, "\n", &save);

  do {  
    strim(p);
    if (p[0] == '#') continue;
    char *res = NULL;
    get_line_token(p, "#", 0, &res);
    strim(res);
    if (res[0] == '<') {
      strncpy(group, res+1, strlen(res)-2);
      group[strlen(res)-2] = 0;
    } else {
      char *key = NULL;
      char *value = NULL;
      get_line_token(res, "=", 1, &value);
      get_line_token(res, "=", 0, &key);
      if (strlen(key) < 1 || strlen(value) < 1) continue;
      strim(key);
      strim(value);
      snprintf(keyname, 255, "%s/%s", group, key);
      dict_add(model_params, keyname, value);
    }
  } while((p = strtok_r(NULL, "\n", &save)) != NULL);
}

void init_physical_quantities(int n, double rescale_factor)
{
  rescale_factor = sqrt(rescale_factor);

#pragma omp parallel for collapse(4)
  for (int mb=0; mb<n_meshblocks; ++mb) {
    for (int i=0; i<mb_nx1+2; ++i) {
      for (int j=0; j<mb_nx2+2; ++j) {
        for (int k=0; k<mb_nx3+2; ++k) {
          data[n]->ne[mb][i][j][k] = data[n]->p[P_RHO][mb][i][j][k] * RHO_unit/(MP+ME);
          data[n]->b[mb][i][j][k] *= rescale_factor;

          double bsq = data[n]->b[mb][i][j][k] / B_unit;
          bsq *= bsq;

          double sigma_m = bsq / data[n]->p[P_RHO][mb][i][j][k];
          double beta_m = data[n]->p[P_UU][mb][i][j][k] * (adiabatic_gamma-1.) / (0.5 * bsq);

          // right now we only support rhigh prescription
          double betasq = beta_m*beta_m / beta_crit/beta_crit;
          double trat = trat_large * betasq/(1. + betasq) + trat_small / (1. + betasq);
          double Thetae_conversion = (MP/ME) * (game-1.) * (gamp-1.) / ( (gamp-1.) + (game-1.)*trat );
          data[n]->thetae[mb][i][j][k] = Thetae_conversion * data[n]->p[P_UU][mb][i][j][k] / data[n]->p[P_RHO][mb][i][j][k];

          // apply per-zone floors/ceilings
          data[n]->thetae[mb][i][j][k] = fmax(data[n]->thetae[mb][i][j][k], 1.e-3);
          data[n]->sigma[mb][i][j][k] = fmax(sigma_m, SMALL);
          data[n]->beta[mb][i][j][k] = fmax(beta_m, SMALL);
          
          // apply per-zone sigma_cut if desired
          if (sigma_m > sigma_cut && !USE_GEODESIC_SIGMACUT) {
            data[n]->b[mb][i][j][k] = 0.;
            data[n]->ne[mb][i][j][k] = 0.;
            data[n]->thetae[mb][i][j][k] = 0.;
          }
        }
      }
    }
  }
}

void populate_boundary_conditions_ordered(int n)
{
  // get whenever we're pulling from zones at or above our refinement
  // level to ensure that we're not pulling from ghost zones
  for (int mb=0; mb<n_meshblocks; ++mb) {

    double xmin = meshblock_geometry[mb*6 + 0];
    double xmax = meshblock_geometry[mb*6 + 1];
    double ymin = meshblock_geometry[mb*6 + 2];
    double ymax = meshblock_geometry[mb*6 + 3];
    double zmin = meshblock_geometry[mb*6 + 4];
    double zmax = meshblock_geometry[mb*6 + 5];

    double dx = (xmax - xmin) / mb_nx1;
    double dy = (ymax - ymin) / mb_nx2;
    double dz = (zmax - zmin) / mb_nx3;

    int level = meshblock_levels[mb];

    for (int i=0; i<mb_nx1+2; ++i) {
      for (int j=0; j<mb_nx2+2; ++j) {
        for (int k=0; k<mb_nx3+2; ++k) {
          if (i==0 || i==mb_nx1+1 || j==0 || j==mb_nx2+1 || k==0 || k==mb_nx3+1) {

            double x = xmin + dx*(i-0.5);
            double y = ymin + dy*(j-0.5);
            double z = zmin + dz*(k-0.5);

            int omb = get_meshblock(x, y, z);
            if (omb < 0) continue;
            if (meshblock_levels[omb] < level) continue;

            data[n]->b[mb][i][j][k] = mb_interp_scalar(omb, x, y, z, data[n]->b);

            for (int l=0; l<n_variables; ++l) {
              data[n]->p[l][mb][i][j][k] = mb_interp_scalar(omb, x, y, z, data[n]->p[l]);
            }

          }
        }
      }
    }
  }

  // now set only when we need to read from a larger zone
  for (int mb=0; mb<n_meshblocks; ++mb) {

    double xmin = meshblock_geometry[mb*6 + 0];
    double xmax = meshblock_geometry[mb*6 + 1];
    double ymin = meshblock_geometry[mb*6 + 2];
    double ymax = meshblock_geometry[mb*6 + 3];
    double zmin = meshblock_geometry[mb*6 + 4];
    double zmax = meshblock_geometry[mb*6 + 5];

    double dx = (xmax - xmin) / mb_nx1;
    double dy = (ymax - ymin) / mb_nx2;
    double dz = (zmax - zmin) / mb_nx3;

    int level = meshblock_levels[mb];

    for (int i=0; i<mb_nx1+2; ++i) {
      for (int j=0; j<mb_nx2+2; ++j) {
        for (int k=0; k<mb_nx3+2; ++k) {
          if (i==0 || i==mb_nx1+1 || j==0 || j==mb_nx2+1 || k==0 || k==mb_nx3+1) {

            double x = xmin + dx*(i-0.5);
            double y = ymin + dy*(j-0.5);
            double z = zmin + dz*(k-0.5);

            int omb = get_meshblock(x, y, z);
            if (omb < 0) continue;
            if (meshblock_levels[omb] >= level) continue;

            data[n]->b[mb][i][j][k] = mb_interp_scalar(omb, x, y, z, data[n]->b);

            for (int l=0; l<n_variables; ++l) {
              data[n]->p[l][mb][i][j][k] = mb_interp_scalar(omb, x, y, z, data[n]->p[l]);
            }

          }
        }
      }
    }
  }
}

void populate_boundary_conditions(int n)
{
  // we only need to set the primitives and b, which is derived from them, since
  // everything else will be correctly handled in the init_physical_quantities(...)
  // function

  // simplest version just directly copies
  if (1==1) {
    
    for (int mb=0; mb<n_meshblocks; ++mb) {
      // six faces
      for (int j=1; j<=mb_nx2; ++j) {
        for (int k=1; k<=mb_nx3; ++k) {
          for (int l=0; l<n_variables; ++l) {
            data[n]->p[l][mb][0][j][k] = data[n]->p[l][mb][1][j][k];
            data[n]->p[l][mb][mb_nx1 + 1][j][k] = data[n]->p[l][mb][mb_nx1][j][k];
          }
          data[n]->b[mb][0][j][k] = data[n]->b[mb][1][j][k];
          data[n]->b[mb][mb_nx1 + 1][j][k] = data[n]->b[mb][mb_nx1][j][k];
        }
      }

      for (int i=1; i<=mb_nx1; ++i) {
        for (int k=1; k<=mb_nx3; ++k) {
          for (int l=0; l<n_variables; ++l) {
            data[n]->p[l][mb][i][0][k] = data[n]->p[l][mb][i][1][k];
            data[n]->p[l][mb][i][mb_nx2 + 1][k] = data[n]->p[l][mb][i][mb_nx2][k];
          }
          data[n]->b[mb][i][0][k] = data[n]->b[mb][i][1][k];
          data[n]->b[mb][i][mb_nx2 + 1][k] = data[n]->b[mb][i][mb_nx2][k];
        }
      }

      for (int i=1; i<=mb_nx1; ++i) {
        for (int j=1; j<=mb_nx2; ++j) {
          for (int l=0; l<n_variables; ++l) {
            data[n]->p[l][mb][i][j][0] = data[n]->p[l][mb][i][j][1];
            data[n]->p[l][mb][i][j][mb_nx3 + 1] = data[n]->p[l][mb][i][j][mb_nx3];
          }
          data[n]->b[mb][i][j][0] = data[n]->b[mb][i][j][1];
          data[n]->b[mb][i][j][mb_nx3 + 1] = data[n]->b[mb][i][j][mb_nx3];
        }
      }

      // edges
      for (int s=1; s<=mb_nx1; ++s) {
        for (int l=0; l<n_variables; ++l) {
          data[n]->p[l][mb][s][0][0] = (data[n]->p[l][mb][s][1][0] + data[n]->p[l][mb][s][0][1]) / 2.;
          data[n]->p[l][mb][mb_nx1 + 1][s][0] = (data[n]->p[l][mb][mb_nx1][s][0] + data[n]->p[l][mb][mb_nx1 + 1][s][1]) / 2.;
          data[n]->p[l][mb][s][mb_nx2 + 1][0] = (data[n]->p[l][mb][s][mb_nx2][0] + data[n]->p[l][mb][s][mb_nx2 + 1][1]) / 2.;
          data[n]->p[l][mb][s][0][mb_nx3 + 1] = (data[n]->p[l][mb][s][0][mb_nx3] + data[n]->p[l][mb][s][1][mb_nx3 + 1]) / 2.;
        }
        data[n]->b[mb][s][0][0] = (data[n]->b[mb][s][1][0] + data[n]->b[mb][s][0][1]) / 2.;
        data[n]->b[mb][mb_nx1 + 1][s][0] = (data[n]->b[mb][mb_nx1][s][0] + data[n]->b[mb][mb_nx1 + 1][s][1]) / 2.;
        data[n]->b[mb][s][mb_nx2 + 1][0] = (data[n]->b[mb][s][mb_nx2][0] + data[n]->b[mb][s][mb_nx2 + 1][1]) / 2.;
        data[n]->b[mb][s][0][mb_nx3 + 1] = (data[n]->b[mb][s][0][mb_nx3] + data[n]->b[mb][s][1][mb_nx3 + 1]) / 2.;
      }
      for (int s = 1; s <= mb_nx2; ++s) {
        for (int l = 0; l < n_variables; ++l) {
          data[n]->p[l][mb][0][s][0] = (data[n]->p[l][mb][1][s][0] + data[n]->p[l][mb][0][s][1]) / 2.;
          data[n]->p[l][mb][mb_nx1 + 1][s][0] = (data[n]->p[l][mb][mb_nx1][s][0] + data[n]->p[l][mb][mb_nx1 + 1][s][1]) / 2.;
          data[n]->p[l][mb][0][s][mb_nx3 + 1] = (data[n]->p[l][mb][1][s][mb_nx3 + 1] + data[n]->p[l][mb][0][s][mb_nx3]) / 2.;
          data[n]->p[l][mb][mb_nx1 + 1][s][mb_nx3 + 1] = (data[n]->p[l][mb][mb_nx1][s][mb_nx3 + 1] + data[n]->p[l][mb][mb_nx1 + 1][s][mb_nx3]) / 2.;
        }
        data[n]->b[mb][0][s][0] = (data[n]->b[mb][1][s][0] + data[n]->b[mb][0][s][1]) / 2.;
        data[n]->b[mb][mb_nx1 + 1][s][0] = (data[n]->b[mb][mb_nx1][s][0] + data[n]->b[mb][mb_nx1 + 1][s][1]) / 2.;
        data[n]->b[mb][0][s][mb_nx3 + 1] = (data[n]->b[mb][1][s][mb_nx3 + 1] + data[n]->b[mb][0][s][mb_nx3]) / 2.;
        data[n]->b[mb][mb_nx1 + 1][s][mb_nx3 + 1] = (data[n]->b[mb][mb_nx1][s][mb_nx3 + 1] + data[n]->b[mb][mb_nx1 + 1][s][mb_nx3]) / 2.;
      }
      for (int s = 1; s <= mb_nx3; ++s) {
        for (int l = 0; l < n_variables; ++l) {
          data[n]->p[l][mb][0][0][s] = (data[n]->p[l][mb][1][0][s] + data[n]->p[l][mb][0][1][s]) / 2.;
          data[n]->p[l][mb][mb_nx1 + 1][0][s] = (data[n]->p[l][mb][mb_nx1][0][s] + data[n]->p[l][mb][mb_nx1 + 1][1][s]) / 2.;
          data[n]->p[l][mb][0][mb_nx2 + 1][s] = (data[n]->p[l][mb][1][mb_nx2 + 1][s] + data[n]->p[l][mb][0][mb_nx2][s]) / 2.;
          data[n]->p[l][mb][mb_nx1 + 1][mb_nx2 + 1][s] = (data[n]->p[l][mb][mb_nx1][mb_nx2 + 1][s] + data[n]->p[l][mb][mb_nx1 + 1][mb_nx2][s]) / 2.;
        }
        data[n]->b[mb][0][0][s] = (data[n]->b[mb][1][0][s] + data[n]->b[mb][0][1][s]) / 2.;
        data[n]->b[mb][mb_nx1 + 1][0][s] = (data[n]->b[mb][mb_nx1][0][s] + data[n]->b[mb][mb_nx1 + 1][1][s]) / 2.;
        data[n]->b[mb][0][mb_nx2 + 1][s] = (data[n]->b[mb][1][mb_nx2 + 1][s] + data[n]->b[mb][0][mb_nx2][s]) / 2.;
        data[n]->b[mb][mb_nx1 + 1][mb_nx2 + 1][s] = (data[n]->b[mb][mb_nx1][mb_nx2 + 1][s] + data[n]->b[mb][mb_nx1 + 1][mb_nx2][s]) / 2.;
      }
      
      // corners
      for (int l=0; l < n_variables; ++l) {
        data[n]->p[l][mb][0][0][0] = (data[n]->p[l][mb][1][0][0] + data[n]->p[l][mb][0][1][0] + data[n]->p[l][mb][0][0][1]) / 3.;
        data[n]->p[l][mb][mb_nx1 + 1][0][0] = (data[n]->p[l][mb][mb_nx1][0][0] + data[n]->p[l][mb][mb_nx1 + 1][1][0] + data[n]->p[l][mb][mb_nx1 + 1][0][1]) / 3.;
        data[n]->p[l][mb][0][mb_nx2 + 1][0] = (data[n]->p[l][mb][1][mb_nx2 + 1][0] + data[n]->p[l][mb][0][mb_nx2][0] + data[n]->p[l][mb][0][mb_nx2 + 1][1]) / 3.;
        data[n]->p[l][mb][mb_nx1 + 1][mb_nx2 + 1][0] = (data[n]->p[l][mb][mb_nx1][mb_nx2 + 1][0] + data[n]->p[l][mb][mb_nx1 + 1][mb_nx2][0] + data[n]->p[l][mb][mb_nx1 + 1][mb_nx2 + 1][1]) / 3.;
        data[n]->p[l][mb][0][0][mb_nx3 + 1] = (data[n]->p[l][mb][1][0][mb_nx3 + 1] + data[n]->p[l][mb][0][1][mb_nx3 + 1] + data[n]->p[l][mb][0][0][mb_nx3]) / 3.;
        data[n]->p[l][mb][mb_nx1 + 1][0][mb_nx3 + 1] = (data[n]->p[l][mb][mb_nx1][0][mb_nx3 + 1] + data[n]->p[l][mb][mb_nx1 + 1][1][mb_nx3 + 1] + data[n]->p[l][mb][mb_nx1 + 1][0][mb_nx3]) / 3.;
        data[n]->p[l][mb][0][mb_nx2 + 1][mb_nx3 + 1] = (data[n]->p[l][mb][1][mb_nx2 + 1][mb_nx3 + 1] + data[n]->p[l][mb][0][mb_nx2][mb_nx3 + 1] + data[n]->p[l][mb][0][mb_nx2 + 1][mb_nx3]) / 3.;
        data[n]->p[l][mb][mb_nx1 + 1][mb_nx2 + 1][mb_nx3 + 1] = (data[n]->p[l][mb][mb_nx1][mb_nx2 + 1][mb_nx3 + 1] + data[n]->p[l][mb][mb_nx1 + 1][mb_nx2][mb_nx3 + 1] + data[n]->p[l][mb][mb_nx1 + 1][mb_nx2 + 1][mb_nx3]) / 3.;
      }
      data[n]->b[mb][0][0][0] = (data[n]->b[mb][1][0][0] + data[n]->b[mb][0][1][0] + data[n]->b[mb][0][0][1]) / 3.;
      data[n]->b[mb][mb_nx1 + 1][0][0] = (data[n]->b[mb][mb_nx1][0][0] + data[n]->b[mb][mb_nx1 + 1][1][0] + data[n]->b[mb][mb_nx1 + 1][0][1]) / 3.;
      data[n]->b[mb][0][mb_nx2 + 1][0] = (data[n]->b[mb][1][mb_nx2 + 1][0] + data[n]->b[mb][0][mb_nx2][0] + data[n]->b[mb][0][mb_nx2 + 1][1]) / 3.;
      data[n]->b[mb][mb_nx1 + 1][mb_nx2 + 1][0] = (data[n]->b[mb][mb_nx1][mb_nx2 + 1][0] + data[n]->b[mb][mb_nx1 + 1][mb_nx2][0] + data[n]->b[mb][mb_nx1 + 1][mb_nx2 + 1][1]) / 3.;
      data[n]->b[mb][0][0][mb_nx3 + 1] = (data[n]->b[mb][1][0][mb_nx3 + 1] + data[n]->b[mb][0][1][mb_nx3 + 1] + data[n]->b[mb][0][0][mb_nx3]) / 3.;
      data[n]->b[mb][mb_nx1 + 1][0][mb_nx3 + 1] = (data[n]->b[mb][mb_nx1][0][mb_nx3 + 1] + data[n]->b[mb][mb_nx1 + 1][1][mb_nx3 + 1] + data[n]->b[mb][mb_nx1 + 1][0][mb_nx3]) / 3.;
      data[n]->b[mb][0][mb_nx2 + 1][mb_nx3 + 1] = (data[n]->b[mb][1][mb_nx2 + 1][mb_nx3 + 1] + data[n]->b[mb][0][mb_nx2][mb_nx3 + 1] + data[n]->b[mb][0][mb_nx2 + 1][mb_nx3]) / 3.;
      data[n]->b[mb][mb_nx1 + 1][mb_nx2 + 1][mb_nx3 + 1] = (data[n]->b[mb][mb_nx1][mb_nx2 + 1][mb_nx3 + 1] + data[n]->b[mb][mb_nx1 + 1][mb_nx2][mb_nx3 + 1] + data[n]->b[mb][mb_nx1 + 1][mb_nx2 + 1][mb_nx3]) / 3.;
    }
  }

  // this function overwrites what happened above by getting prim/b values from
  // neighboring zones in neighbor.level >= self.level zones first (safe), then
  // from the remaining neighbor pairs, which should (by then) be populated, as
  // one of two different-level neighbors will always have the higher level.
  populate_boundary_conditions_ordered(n);
}

size_t get_athenak_datastart(char *fname)
{
  char *res = NULL;
  char *line = NULL;

  FILE *fp;
  fp = fopen(fname, "rb");

  fseek(fp, 0, SEEK_SET);
  if (read_line(fp, &line, "! unable to read file header.\n") == -1) return -1;
  if (strcmp(line, "Athena binary output version=1.1\n") != 0) {
    fprintf(stderr, "! got bad file header: %s\n", line);
    return -1;
  }

  if (read_line(fp, &line, "! unable to read header line.\n") == -1) return -1;
  get_line_token(line, "=", 1, &res);
  int preheader_size = atoi(res);
  for (int i=0; i<preheader_size-1; ++i) {
    read_line(fp, &line, "! unable to read preheader line\n");
  }

  if (read_line(fp, &line, "! unable to read number of variables.\n") == -1) return -1;
  if (read_line(fp, &line, "! unable to read variable names.\n") == -1) return -1;

  // deal with header
  read_set_athenak_header(fp, NULL);

  size_t datastart = ftell(fp);
  fclose(fp);

  return datastart;
}

size_t process_athenak_header(char *fnam, int index)
{
  char fname[256];
  snprintf(fname, 255, fnam, index);
  
  char *res = NULL;
  char *line = NULL;

  FILE *fp;
  fp = fopen(fname, "rb");

  fseek(fp, 0, SEEK_END);
  size_t filesize = ftell(fp);

  fseek(fp, 0, SEEK_SET);
  if (read_line(fp, &line, "! unable to read file header.\n") == -1) return -1;
  if (strcmp(line, "Athena binary output version=1.1\n") != 0) {
    fprintf(stderr, "! got bad file header: %s\n", line);
    return -1;
  }

  if (read_line(fp, &line, "! unable to read header line.\n") == -1) return -1;
  get_line_token(line, "=", 1, &res);
  int preheader_size = atoi(res);
  for (int i=0; i<preheader_size-1; ++i) {
    char *key = NULL;
    char *value = NULL;
    read_line(fp, &line, "! unable to read preheader line\n");
    get_line_token(line, "=", 1, &value);
    get_line_token(line, "=", 0, &key);
    strim(key);
    strim(value);
    dict_add(model_params, key, value);
  }

  if (read_line(fp, &line, "! unable to read number of variables.\n") == -1) return -1;
  get_line_token(line, "=", 1, &res);
  n_variables = atoi(res);
  if (read_line(fp, &line, "! unable to read variable names.\n") == -1) return -1;
  // TODO process variable names in "line"?
  
  // deal with header
  read_set_athenak_header(fp, model_params);
  adiabatic_gamma = atof(dict_get(model_params, "gamma", "1.4444444444444"));
 
  size_t datastart = ftell(fp);
  fclose(fp);

  // compute meshblock size
  size_t locsize = atoi(dict_get(model_params, "size of location", "4"));
  size_t varsize = atoi(dict_get(model_params, "size of variable", "4"));
  mb_nx1 = atoi(dict_get(model_params, "meshblock/nx1", "64"));
  mb_nx2 = atoi(dict_get(model_params, "meshblock/nx2", "64"));
  mb_nx3 = atoi(dict_get(model_params, "meshblock/nx3", "32"));

  size_t meshblock_bytes = 24 + 16 + 6*locsize + mb_nx1*mb_nx2*mb_nx3*n_variables*varsize;
  n_meshblocks = (filesize - datastart) / meshblock_bytes; 

  // initialize storage
  for (int n = 0; n < NSUP; n++) {
    data[n]->p = malloc_rank5(n_variables,n_meshblocks,mb_nx1+2,mb_nx2+2,mb_nx3+2);
    data[n]->ne = malloc_rank4(n_meshblocks,mb_nx1+2,mb_nx2+2,mb_nx3+2);
    data[n]->thetae = malloc_rank4(n_meshblocks,mb_nx1+2,mb_nx2+2,mb_nx3+2);
    data[n]->b = malloc_rank4(n_meshblocks,mb_nx1+2,mb_nx2+2,mb_nx3+2);
    data[n]->sigma = malloc_rank4(n_meshblocks,mb_nx1+2,mb_nx2+2,mb_nx3+2);
    data[n]->beta = malloc_rank4(n_meshblocks,mb_nx1+2,mb_nx2+2,mb_nx3+2);
  }
  meshblock_geometry = malloc(n_meshblocks * 6 * sizeof(double));
  meshblock_levels = malloc(n_meshblocks * sizeof(int));

  if (line) {
    free(line);
  }

  // force geometry
  a = atof(dict_get(model_params, "coord/a", "0"));
  Rh = 1. + sqrt(1. - a*a);

  use_eKS_internal = 0;  // since we're already using eKS, don't force translation
  metric = METRIC_EKS;

  double x1min = fabs(atof(dict_get(model_params, "mesh/x1min", "0")));
  double x1max = fabs(atof(dict_get(model_params, "mesh/x1max", "0")));
  double x2min = fabs(atof(dict_get(model_params, "mesh/x2min", "0")));
  double x2max = fabs(atof(dict_get(model_params, "mesh/x2max", "0")));
  double x3min = fabs(atof(dict_get(model_params, "mesh/x3min", "0")));
  double x3max = fabs(atof(dict_get(model_params, "mesh/x3max", "0")));
    
  x1max = x1max > x1min ? x1max : x1min;
  x2max = x2max > x2min ? x2max : x2min;
  x3max = x3max > x3min ? x3max : x3min;
  double rmax = x1max > x2max ? x1max : x2max;
  rmax = rmax > x3max ? rmax : x3max;

  startx[0] = 0.0;
  startx[1] = log(0.9*Rh);
  startx[2] = 0.0;
  startx[3] = 0.0;
  stopx[0] = 0.0;
  stopx[1] = log(rmax/0.9);
  stopx[2] = M_PI;
  stopx[3] = 2*M_PI;

  MULOOP cstartx[mu] = startx[mu];
  MULOOP cstopx[mu] = stopx[mu];

  fprintf(stderr, "done!\n");

  fprintf(stderr, "rmin_geo, rmax_geo: %g %g\n", rmin_geo, rmax_geo);
  fprintf(stderr, "Native coordinate start: %g %g %g stop: %g %g %g\n",
            cstartx[1], cstartx[2], cstartx[3], cstopx[1], cstopx[2], cstopx[3]);
  fprintf(stderr, "Grid start: %g %g %g stop: %g %g %g\n",
            startx[1], startx[2], startx[3], stopx[1], stopx[2], stopx[3]);

  set_units();

  return datastart;
}

// returns b in code units
double convert_ub_prim_cks_eks(double x, double y, double z, double Uprim[3], double Bprim[3])
{ 
  double gcov_cks[NDIM][NDIM];
  double gcon00 = get_g_cks(x, y, z, gcov_cks);
  double alpha = 1. / sqrt(-gcon00);

  // don't trust the pole
  double SMALLISH = 1.e-15;
  if (fabs(x) < SMALLISH) x = SMALLISH;
  if (fabs(y) < SMALLISH) y = SMALLISH;
  
  // construct CKS ucon and ucov
  double ucon_cks[4] = {0.};
  double ucov_cks[4] = {0.};

  double UdotU = 0.;
  for (int l=1; l<NDIM; ++l) {
    for (int m=1; m<NDIM; ++m) {
      UdotU += gcov_cks[l][m] * Uprim[l-1] * Uprim[m-1];
    }
  }
  double gamma = sqrt(1. + UdotU);
 
  ucon_cks[0] = gamma / alpha;
  ucon_cks[1] = Uprim[0] - gamma * alpha * gcov_cks[0][1];
  ucon_cks[2] = Uprim[1] - gamma * alpha * gcov_cks[0][2];
  ucon_cks[3] = Uprim[2] - gamma * alpha * gcov_cks[0][3];
  flip_index(ucon_cks, gcov_cks, ucov_cks);

  if (1==0) {
    double udu = 0.;
    for (int i=0; i<4; ++i) udu += ucon_cks[i] * ucov_cks[i];
    fprintf(stderr, "u.u = %g\n", udu);
  }
 
  // construct CKS bcon
  double bcon_cks[4] = {0.};
  bcon_cks[0] = Bprim[0]*ucov_cks[1] + Bprim[1]*ucov_cks[2] + Bprim[2]*ucov_cks[3];
  bcon_cks[1] = (Bprim[0] + ucon_cks[1]*bcon_cks[0]) / ucon_cks[0];
  bcon_cks[2] = (Bprim[1] + ucon_cks[2]*bcon_cks[0]) / ucon_cks[0];
  bcon_cks[3] = (Bprim[2] + ucon_cks[3]*bcon_cks[0]) / ucon_cks[0];

  if (1==0) {
    double bdu = 0.;
    for (int i=0; i<4; ++i) bdu += bcon_cks[i] * ucov_cks[i];
    fprintf(stderr, "b.u = %g\n", bdu);
  }

  // get eks geometry
  double X_eks[NDIM] = {0.};
  double gcov_eks[NDIM][NDIM], gcon_eks[NDIM][NDIM];
  cks_to_sks(x, y, z, X_eks);
  X_eks[1] = log(X_eks[1]);
  gcov_func(X_eks, gcov_eks);
  gcon_func(gcov_eks, gcon_eks);

  // convert to eKS four-vectors
  double ucon_eks[4] = {0.};
  double bcon_eks[4] = {0.};
  double bcov_eks[4] = {0.};

  double R = sqrt(x*x + y*y + z*z);
  double r = sqrt(R*R - a*a + sqrt((R*R-a*a)*(R*R-a*a) + 4.*a*a*z*z))/sqrt(2.);
  double sqrt_term = 2.*r*r - R*R + a*a;

  ucon_eks[0] = ucon_cks[0];
  ucon_eks[1] = ucon_cks[1] * (x*r)/sqrt_term +
                ucon_cks[2] * (y*r)/sqrt_term +
                ucon_cks[3] * z/r * (r*r+a*a)/sqrt_term;
  ucon_eks[2] = ucon_cks[1] * (x*z)/(r * sqrt_term * sqrt(1.-z*z/r/r) + SMALLISH) +
                ucon_cks[2] * (y*z)/(r * sqrt_term * sqrt(1.-z*z/r/r) + SMALLISH) +
                ucon_cks[3] * ( (z*z)*(r*r+a*a)/(r*r*r * sqrt_term * sqrt(1.0-z*z/r/r) + SMALLISH) 
                       - 1.0/(r*sqrt(1.0-z*z/r/r) + SMALLISH) );
  ucon_eks[3] = ucon_cks[1] * (-y/(x*x+y*y+SMALLISH) + a*r*x/((r*r+a*a)*sqrt_term)) +
                ucon_cks[2] * (x/(x*x+y*y+SMALLISH) + a*r*y/((r*r+a*a)*sqrt_term)) +
                ucon_cks[3] * (a*z/r/sqrt_term);
  ucon_eks[1] /= exp(X_eks[1]);

  bcon_eks[0] = bcon_cks[0];
  bcon_eks[1] = bcon_cks[1] * (x*r)/sqrt_term +
                bcon_cks[2] * (y*r)/sqrt_term +
                bcon_cks[3] * z/r * (r*r+a*a)/sqrt_term;
  bcon_eks[2] = bcon_cks[1] * (x*z)/(r * sqrt_term * sqrt(1.0-z*z/r/r) + SMALLISH) +
                bcon_cks[2] * (y*z)/(r * sqrt_term * sqrt(1.0-z*z/r/r) + SMALLISH) +
                bcon_cks[3] * ( (z*z)*(r*r+a*a)/(r*r*r * sqrt_term * sqrt(1.0-z*z/r/r) + SMALLISH) 
                       - 1.0/(r*sqrt(1.0-z*z/r/r) + SMALLISH) );
  bcon_eks[3] = bcon_cks[1] * (-y/(x*x+y*y+SMALLISH) + a*r*x/((r*r+a*a)*sqrt_term)) + \
                bcon_cks[2] * (x/(x*x+y*y+SMALLISH) + a*r*y/((r*r+a*a)*sqrt_term)) + \
                bcon_cks[3] * (a*z/r/sqrt_term);
  bcon_eks[1] /= exp(X_eks[1]);
  flip_index(bcon_eks, gcov_eks, bcov_eks);

  double bsq = 0.;
  MULOOP bsq += bcon_eks[mu] * bcov_eks[mu];

  if (1==0) {
    double ucov_eks[4] = {0};
    flip_index(ucon_eks, gcov_eks, ucov_eks);
    double udu = 0.;
    double bdu = 0.;
    double udu_cks = 0.;
    for (int i=0; i<4; ++i) {
      udu += ucon_eks[i] * ucov_eks[i];
      udu_cks += ucon_cks[i] * ucov_cks[i];
      bdu += bcon_eks[i] * ucov_eks[i];
    }
    fprintf(stderr, "u.u = %g  b.u = %g  u.u cks = %g\n", udu, bdu, udu_cks);
  }

  if (1==0) {
    double bcov_cks[4] = {0.};
    flip_index(bcon_cks, gcov_cks, bcov_cks);
    double bsq_cks = 0.;
    MULOOP bsq_cks += bcon_cks[mu] * bcov_cks[mu];
    fprintf(stderr, "bsq cks vs eks %g %g\n", bsq_cks, bsq);
  }

  // set primitives
  Uprim[0] = ucon_eks[1] - ucon_eks[0] * gcon_eks[0][1] / gcon_eks[0][0];
  Uprim[1] = ucon_eks[2] - ucon_eks[0] * gcon_eks[0][2] / gcon_eks[0][0];
  Uprim[2] = ucon_eks[3] - ucon_eks[0] * gcon_eks[0][3] / gcon_eks[0][0];
  Bprim[0] = bcon_eks[1] * ucon_eks[0] - bcon_eks[0] * ucon_eks[1];
  Bprim[1] = bcon_eks[2] * ucon_eks[0] - bcon_eks[0] * ucon_eks[2];
  Bprim[2] = bcon_eks[3] * ucon_eks[0] - bcon_eks[0] * ucon_eks[3];

  return bsq;
}

void load_athenak_data(int n, char *fnam, int index, size_t datastart)
{
  char fname[256];
  snprintf(fname, 255, fnam, index);
  fprintf(stderr, "Loading file %s... ", fname);

  data[n]->t = get_athenak_dump_time(fname);
  DTd = data[n]->t - last_time_loaded;
  last_time_loaded = data[n]->t;

  nloaded++;

  if (datastart == 0) {
    datastart = get_athenak_datastart(fname);
  }

  FILE *fp;
  fp = fopen(fname, "rb");
  fseek(fp, 0, SEEK_END);
  long int filesize = ftell(fp);
  fseek(fp, datastart, SEEK_SET);

  size_t locsize = atoi(dict_get(model_params, "size of location", "4"));
  size_t varsize = atoi(dict_get(model_params, "size of variable", "4"));

  size_t nmb = 0;

  // get more info on meshblocks and such
  while (ftell(fp) < filesize) {

    int si = read_int(fp, 4);
    int ei = read_int(fp, 4);
    int sj = read_int(fp, 4);
    int ej = read_int(fp, 4);
    int sk = read_int(fp, 4);
    int ek = read_int(fp, 4);
    int nx1_out = ei-si + 1;
    int nx2_out = ej-sj + 1;
    int nx3_out = ek-sk + 1;

    if (nx1_out != mb_nx1 || nx2_out != mb_nx2 || nx3_out != mb_nx3) {
      fprintf(stderr, "mesh block sizes do not agree: %d %d %d %ld %ld %ld\n",
              nx1_out, nx2_out, nx3_out, mb_nx1, mb_nx2, mb_nx3);
      exit(7);
    }
   
    // logical location (not used)
    int idum = read_int(fp, 4);
    idum = read_int(fp, 4);
    idum = read_int(fp, 4);
    idum = read_int(fp, 4);
    meshblock_levels[nmb] = idum;

    double xmin = read_double(fp, locsize);
    double xmax = read_double(fp, locsize);
    double ymin = read_double(fp, locsize);
    double ymax = read_double(fp, locsize);
    double zmin = read_double(fp, locsize);
    double zmax = read_double(fp, locsize);

    meshblock_geometry[nmb*6 + 0] = xmin;
    meshblock_geometry[nmb*6 + 1] = xmax;
    meshblock_geometry[nmb*6 + 2] = ymin;
    meshblock_geometry[nmb*6 + 3] = ymax;
    meshblock_geometry[nmb*6 + 4] = zmin;
    meshblock_geometry[nmb*6 + 5] = zmax;
 
    for (int v=0; v<n_variables; ++v) {
      for (int k=0; k<mb_nx3; ++k) {
        for (int j=0; j<mb_nx2; ++j) {
          for (int i=0; i<mb_nx1; ++i) {
            data[n]->p[v][nmb][i+1][j+1][k+1] = read_double(fp, varsize);
          }
        }
      }
    }

    nmb++;
  }

  fclose(fp);

  // construct four-vector quantities and convert from CKS -> eKS
#pragma omp parallel for
  for (int mb=0; mb<n_meshblocks; ++mb) {
    
    double xmin = meshblock_geometry[mb*6 + 0];
    double xmax = meshblock_geometry[mb*6 + 1];
    double ymin = meshblock_geometry[mb*6 + 2];
    double ymax = meshblock_geometry[mb*6 + 3];
    double zmin = meshblock_geometry[mb*6 + 4];
    double zmax = meshblock_geometry[mb*6 + 5];

    double dx = (xmax - xmin) / mb_nx1;
    double dy = (ymax - ymin) / mb_nx2;
    double dz = (zmax - zmin) / mb_nx3;

    for (int i=0; i<mb_nx1; ++i) {
      for (int j=0; j<mb_nx2; ++j) {
        for (int k=0; k<mb_nx3; ++k) {

          double x = xmin + dx*(i+0.5);
          double y = ymin + dy*(j+0.5);
          double z = zmin + dz*(k+0.5);

          double Uprim[3] = {data[n]->p[P_U1][mb][i+1][j+1][k+1], data[n]->p[P_U2][mb][i+1][j+1][k+1], data[n]->p[P_U3][mb][i+1][j+1][k+1]};
          double Bprim[3] = {data[n]->p[P_B1][mb][i+1][j+1][k+1], data[n]->p[P_B2][mb][i+1][j+1][k+1], data[n]->p[P_B3][mb][i+1][j+1][k+1]};

          double bsq = convert_ub_prim_cks_eks(x, y, z, Uprim, Bprim);
          data[n]->b[mb][i+1][j+1][k+1] = sqrt(bsq) * B_unit;

          data[n]->p[P_U1][mb][i+1][j+1][k+1] = Uprim[0];
          data[n]->p[P_U2][mb][i+1][j+1][k+1] = Uprim[1];
          data[n]->p[P_U3][mb][i+1][j+1][k+1] = Uprim[2];
          data[n]->p[P_B1][mb][i+1][j+1][k+1] = Bprim[0];
          data[n]->p[P_B2][mb][i+1][j+1][k+1] = Bprim[1];
          data[n]->p[P_B3][mb][i+1][j+1][k+1] = Bprim[2];
        }
      }
    }
  }

  // populate boundary conditions
  populate_boundary_conditions(n);

  // initialize scalars
  double rescale_factor = 1.0;  // can be used to rescale according to accretion rate
  init_physical_quantities(n, rescale_factor);

  fprintf(stderr, "done!\n");
}

// this function is run after reading the parameter file but
// before setting up anything else, so it is responsible for
// data structures and, e.g., setting geometry.
void init_model(double *tA, double *tB)
{
  // set data in order in memory
  data[0] = &dataA;
  data[1] = &dataB;
  data[2] = &dataC;

  model_params = dict_new();

  fprintf(stderr, "Determining dump file type... ");
  set_dumpfile_type(fnam, dumpmin);

  fprintf(stderr, "Reading file header... ");
  // this method ensures that the grid has been initialized
  size_t datastart = process_athenak_header(fnam, dumpmin);

  // read fluid data
  fprintf(stderr, "Reading data...\n");
  load_athenak_data(0, fnam, dumpmin, datastart);
  dumpidx += dumpskip;
#if SLOW_LIGHT
  update_data(tA, tB);
  update_data(tA, tB);
  tf = get_dump_t(fnam, dumpmax) - 1.e-5;
#else // FAST LIGHT
  data[2]->t = 10000.;
#endif // SLOW_LIGHT
  DTd = 10.;

  fprintf(stderr, "\n");
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
  double x, y, z;
  sks_to_cks(exp(X[1]), X[2], X[3], &x, &y, &z);
  int mb = get_meshblock(x, y, z);

  double gcov[NDIM][NDIM], gcon[NDIM][NDIM];

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  if (mb < 0 || radiating_region(X) == 0) {
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
  double Vcon[NDIM] = {0.};
  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);
  Vcon[1] = mb_interp_scalar_time(mb, x, y, z, data[nA]->p[P_U1], data[nB]->p[P_U1], tfac);
  Vcon[2] = mb_interp_scalar_time(mb, x, y, z, data[nA]->p[P_U2], data[nB]->p[P_U2], tfac);
  Vcon[3] = mb_interp_scalar_time(mb, x, y, z, data[nA]->p[P_U3], data[nB]->p[P_U3], tfac);

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
  double Bcon1 = mb_interp_scalar_time(mb, x, y, z, data[nA]->p[P_B1], data[nB]->p[P_B1], tfac);
  double Bcon2 = mb_interp_scalar_time(mb, x, y, z, data[nA]->p[P_B2], data[nB]->p[P_B2], tfac);
  double Bcon3 = mb_interp_scalar_time(mb, x, y, z, data[nA]->p[P_B3], data[nB]->p[P_B3], tfac);

  // get Bcon
  Bcon[0] = Bcon1*Ucov[1] + Bcon2*Ucov[2] + Bcon3*Ucov[3];
  Bcon[1] = (Bcon1 + Ucon[1] * Bcon[0]) / Ucon[0];
  Bcon[2] = (Bcon2 + Ucon[2] * Bcon[0]) / Ucon[0];
  Bcon[3] = (Bcon3 + Ucon[3] * Bcon[0]) / Ucon[0];

  // lower
  lower(Bcon, gcov, Bcov); 
}

// get primitives at X (including time). not used to compute
// coefficients so we don't run checks. useful for debugging
void get_model_primitives(double X[NDIM], double *p)
{
  double x, y, z;
  sks_to_cks(exp(X[1]), X[2], X[3], &x, &y, &z);
  int mb = get_meshblock(x, y, z);

  if (mb < 0 || radiating_region(X) == 0) {
    return;
  }

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  for (int l=0; l<n_variables; ++l) {
    p[l] = mb_interp_scalar_time(mb, x, y, z, data[nA]->p[l], data[nB]->p[l], tfac);
  }
}

double get_model_thetae(double X[NDIM])
{
  double x, y, z;
  sks_to_cks(exp(X[1]), X[2], X[3], &x, &y, &z);
  int mb = get_meshblock(x, y, z);

  if (mb < 0 || radiating_region(X) == 0) {
    return 0.;
  }

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  double interpv = mb_interp_scalar_time(mb, x, y, z, data[nA]->thetae, data[nB]->thetae, tfac);
  return interpv > 0 ? interpv : 0;
}

//b field strength in Gauss
double get_model_b(double X[NDIM])
{
  double x, y, z;
  sks_to_cks(exp(X[1]), X[2], X[3], &x, &y, &z);
  int mb = get_meshblock(x, y, z);

  if (mb < 0 || radiating_region(X) == 0) {
    return 0.;
  }

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  double interpv = mb_interp_scalar_time(mb, x, y, z, data[nA]->b, data[nB]->b, tfac);
  return interpv > 0 ? interpv : 0;
}

double get_model_sigma(double X[NDIM])
{
  double x, y, z;
  sks_to_cks(exp(X[1]), X[2], X[3], &x, &y, &z);
  int mb = get_meshblock(x, y, z);

  if (mb < 0 || radiating_region(X) == 0) {
    return 0.;
  }

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  double interpv = mb_interp_scalar_time(mb, x, y, z, data[nA]->sigma, data[nB]->sigma, tfac);
  return interpv > 0 ? interpv : 0;
}

double get_model_beta(double X[NDIM])
{
  double x, y, z;
  sks_to_cks(exp(X[1]), X[2], X[3], &x, &y, &z);
  int mb = get_meshblock(x, y, z);

  if (mb < 0 || radiating_region(X) == 0) {
    return 0.;
  }

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  double interpv = mb_interp_scalar_time(mb, x, y, z, data[nA]->beta, data[nB]->beta, tfac);
  return interpv > 0 ? interpv : 0;
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
  double x, y, z;
  sks_to_cks(exp(X[1]), X[2], X[3], &x, &y, &z);
  int mb = get_meshblock(x, y, z);

  if (mb < 0 || radiating_region(X) == 0) {
    return 0.;
  }

  double sigma_smoothfac = 1;

#if USE_GEODESIC_SIGMACUT
  double sigma = get_model_sigma(X);
  if (sigma > sigma_cut) return 0.;
  sigma_smoothfac = get_sigma_smoothfac(sigma);
#endif

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  double interpv = mb_interp_scalar_time(mb, x, y, z, data[nA]->ne, data[nB]->ne, tfac) * sigma_smoothfac;
  return interpv > 0 ? interpv : 0;
}

void output_hdf5()
{
  hdf5_set_directory("/");
  
  hdf5_set_directory("/header/");
  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  hdf5_write_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&M_unit, "M_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&T_unit, "T_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&B_unit, "B_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&RHO_unit, "RHO_unit", H5T_IEEE_F64LE);

  hdf5_set_directory("/header/");
  hdf5_make_directory("electrons");
  hdf5_set_directory("/header/electrons/");
  hdf5_write_single_val(&trat_small, "rlow", H5T_IEEE_F64LE);
  hdf5_write_single_val(&trat_large, "rhigh", H5T_IEEE_F64LE);
  hdf5_write_single_val(&beta_crit, "beta_crit", H5T_IEEE_F64LE);

  hdf5_set_directory("/");
}

// This should return 1 when we're in a region of space that is supported by the emission model
int radiating_region(double X[NDIM])
{
  double r, h;
  bl_coord(X, &r, &h);

  if (r < rmax_geo) return 1;

  return 0; 
}

// In case we want to mess with emissivities directly
void get_model_jar(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV) {return;}
void get_model_jk(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv) {return;}

