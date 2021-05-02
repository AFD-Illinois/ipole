#include "simcoords.h"
#include "decs.h"
#include "coordinates.h"
#include "model_params.h"
#include "koral_coords.h"
#include "debug_tools.h"
#include "hdf5_utils.h"

#include <math.h>
#include <assert.h>

// switches
int use_simcoords = 0;
int simcoords = 0;

// global metric defintiions
metric_params_koral_mks3 mp_koral_mks3 = { 0 };
metric_params_koral_jetcoords mp_koral_jetcoords = { 0 };

// invisible "members"
static double *simcoords_x1 = NULL;
static double *simcoords_x2 = NULL;
static size_t sc_n1 = 1;
static size_t sc_n2 = 1;
static double er0 = 0;
static double h0 = 0;
static double der = 1.e-5;
static double dh = 1.e-5;
static double x1i_oob = 1.e10;

static double *ks_r = NULL;
static double *ks_h = NULL;
static size_t n1 = 1;
static size_t n2 = 1;
static double startx1, startx2, dx1, dx2, maxks_r;

// coordinate function forward declarations
static void rev_KORAL_MKS3(double *xKS, double *xMKS);
static void rev_saved_rh(double *eKS, double *gridcoord);

#define ij2oned(i,j) ((size_t)(j+sc_n2*(i)))

// interface functions
void load_ks_rh_from_file(const char *fname)
{
  if ( hdf5_open(fname) < 0 ) {
    fprintf(stderr, "! unable to open file %s. exiting!\n", fname);
    exit(-1);
  }

  size_t n3;
  hdf5_set_directory("/header/");
  hdf5_read_single_val(&n1, "n1", H5T_STD_I32LE);
  hdf5_read_single_val(&n2, "n2", H5T_STD_I32LE);
  hdf5_read_single_val(&n3, "n3", H5T_STD_I32LE);

  hdf5_set_directory("/header/geom/");
  hdf5_read_single_val(&startx1, "startx1", H5T_IEEE_F64LE);
  hdf5_read_single_val(&startx2, "startx2", H5T_IEEE_F64LE);
  hdf5_read_single_val(&dx1, "dx1", H5T_IEEE_F64LE);
  hdf5_read_single_val(&dx2, "dx2", H5T_IEEE_F64LE);

  ks_r = calloc(n1*n2, sizeof(*ks_r));
  ks_h = calloc(n1*n2, sizeof(*ks_h));

  hdf5_set_directory("/grid_out/");

  hsize_t fdims[] = { n1, n2, n3 };
  hsize_t fstart[] = { 0, 0, 0 };
  hsize_t fcount[] = { n1, n2, 1 };
  hsize_t mdims[] = { n1, n2, 1 };
  hsize_t mstart[] = { 0, 0, 0 };

  hdf5_read_array(ks_r, "r", 3, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);
  hdf5_read_array(ks_h, "th", 3, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE);

  for (int i=0; i<n1; ++i) {
    for (int j=0; j<n2; ++j) {
      if (ks_r[j+n2*i] > maxks_r) {
        maxks_r = ks_r[j+n2*i];
      }
    }
  }

  hdf5_close();
}

void initialize_simgrid(size_t n1, size_t n2, double x1i, double x1f, double x2i, double x2f) 
{
  if (ks_r == NULL) {
    fprintf(stderr, "! must call load_ks_rh_from_file(...) before initialize_simgrid(...)\n");
    exit(-7);
  }

  if (use_simcoords == 0) return;

  sc_n1 = n1;
  sc_n2 = n2;

  simcoords_x1 = calloc(sc_n1*sc_n2, sizeof(*simcoords_x1));
  simcoords_x2 = calloc(sc_n1*sc_n2, sizeof(*simcoords_x2));

  // note we've made the assumption that x1i,x2i gives the left edge of the grid and
  // x1f,x2f gives the right edge. this means the range is (n1+1)*dx1,(n2+1)*dx2, so
  // that we cover the full domain. if we have an oob error when trying to determine
  // the ii,jj for interpolating, we return x1f+1, so x1f+1 must not be a valid grid
  // coordinate in the fluid model!
  double Rin = 1.05 * (1. + sqrt(1. - a*a));
  double Rout = maxks_r;
  double h_min = 0.;
  double h_max = M_PI;

  fprintf(stderr, "Rout %g\n\n", Rout);

  // set limit for tracking geodesic emission
  rmax_geo = fmin(rmax_geo, Rout);
  rmin_geo = fmax(rmin_geo, Rin);

  x1i_oob = x1f + 1.;
  er0 = log(Rin);
  der = (log(Rout) - log(Rin)) / sc_n1;
  h0 = h_min;
  dh = (h_max - h_min) / sc_n2;

  // coordinate system is MKS, so set reasonable values here
  cstartx[0] = 0;
  cstartx[1] = log(Rin);
  cstartx[2] = 0;
  cstartx[3] = 0;
  cstopx[0] = 0;
  cstopx[1] = log(Rout);
  cstopx[2] = 1.0;
  cstopx[3] = 2*M_PI;

#pragma omp parallel for schedule(dynamic,2) collapse(2) shared(simcoords_x1,simcoords_x2)
  for (size_t i=0; i<n1; ++i) {
    for (size_t j=0; j<n2; ++j) {
  
      double eKS[NDIM] = { 0 };
      double gridcoord[NDIM] = { 0 };

      eKS[1] = exp(er0 + der*i);
      eKS[2] = h0 + dh*j;

      if (simcoords == SIMCOORDS_KORAL_MKS3) {
        rev_KORAL_MKS3(eKS, gridcoord);
      } else if (simcoords == SIMCOORDS_KORAL_JETCOORDS) {
        rev_saved_rh(eKS, gridcoord);
      } else {
        assert(1==0);  // unknown simulation coordinates
      }

      simcoords_x1[ij2oned(i,j)] = gridcoord[1];
      simcoords_x2[ij2oned(i,j)] = gridcoord[2];
    }
  }
}

void finalize_simgrid()
{
  free(simcoords_x2);
  free(simcoords_x1);

  free(ks_h);
  free(ks_r);
}

void eks_to_simcoord(double eKS[NDIM], double gridcoord[NDIM])
{ 
  double i = (eKS[1] - er0)  / der;
  size_t ii = (size_t) i;
  double di = i - ii;

  double j = (eKS[2] - h0) / dh;
  size_t jj = (size_t) j;
  double dj = j - jj;

  if (ii < 0 || ii >= sc_n1 || jj < 0 || jj >= sc_n2) {
    gridcoord[1] = x1i_oob;
    return;
  }

  gridcoord[0] = eKS[0];

  gridcoord[1] = simcoords_x1[ ij2oned(ii,jj) ] * (1.-di)*(1.-dj)
               + simcoords_x1[ ij2oned(ii,jj+1) ] * (1.-di)*dj
               + simcoords_x1[ ij2oned(ii+1,jj) ] * di*(1.-dj)
               + simcoords_x1[ ij2oned(ii+1,jj+1) ] * di*dj;

  gridcoord[2] = simcoords_x2[ ij2oned(ii,jj) ] * (1.-di)*(1.-dj)
               + simcoords_x2[ ij2oned(ii,jj+1) ] * (1.-di)*dj
               + simcoords_x2[ ij2oned(ii+1,jj) ] * di*(1.-dj)
               + simcoords_x2[ ij2oned(ii+1,jj+1) ] * di*dj;

  gridcoord[3] = eKS[3];
}


//////
// function definitions for populating the interpolation grid
//////

static void rev_saved_rh(double *eKS, double *gridcoord)
{
  int gridi = -1;
  double griddx1 = 0.;
  double r = exp(eKS[1]);

  // here we assume x1 is independent of x2 and that grid zones 
  // are not incredibly skew
  for (int i=0; i<n1-1; ++i) {
    if (ks_r[n2*i] <= log(r) && log(r) < ks_r[n2*(i+1)]) {
      gridi = i;
      griddx1 = ( r - exp(ks_r[n2*i]) ) / ( exp(ks_r[n2*(i+1)]) - exp(ks_r[n2*i]) );
      break;
    }
  }

  int gridj = -1;
  double griddx2 = 0.;
  double h = eKS[2];
  for (int j=0; j<n2-1; ++j) {
    if (ks_h[j+n2*gridi] <= h && h < ks_h[(j+1)+n2*gridi]) {
      gridj = j;
      double a = ks_h[j+n2*gridi];
      double b = ks_h[1+j+n2*gridi];
      double c = ks_h[j+n2*(1+gridi)];
      double d = ks_h[1+j+n2*(1+gridi)];
      griddx2  = (h - a + a*griddx1 - c*griddx1);
      griddx2 /= (b - a + a*griddx1 - b*griddx1 + d*griddx1 - c*griddx1);
      break;
    }
  }

  gridcoord[0] = eKS[0];
  gridcoord[1] = startx1 + (gridi+griddx1+0.5) * dx1;
  gridcoord[2] = startx2 + (gridj+griddx2+0.5) * dx2;
  gridcoord[3] = eKS[3];
}

static void rev_KORAL_MKS3(double *xKS, double *xMKS)
{
  double KSx0=xKS[0];
  double KSx1=xKS[1];
  double KSx2=xKS[2];
  double KSx3=xKS[3];

  double R0 = mp_koral_mks3.r0;
  double H0 = mp_koral_mks3.h0;
  double MY1 = mp_koral_mks3.my1;
  double MY2 = mp_koral_mks3.my2;
  double MP0 = mp_koral_mks3.mp0;

  xMKS[0] = KSx0;

  xMKS[1] = log(KSx1 - R0);

  xMKS[2]  = -H0*M_PI * pow(KSx1, MP0) - H0*MY1*M_PI * pow(2., 1.+MP0);
  xMKS[2] += 2.*H0*MY1*M_PI * pow(KSx1, MP0) + H0*MY2*M_PI*pow(2., 1.+MP0);
  xMKS[2] += 2.*pow(KSx1, MP0)*atan(((-2.*KSx2 + M_PI)*tan((H0*M_PI)/2.))/M_PI);
  xMKS[2] /= 2.*M_PI*H0 * (-pow(KSx1, MP0) - pow(2., 1.+MP0)*MY1 + 2*pow(KSx1, MP0)*MY1 + pow(2., 1.+MP0)*MY2);

  xMKS[3] = KSx3;
}

