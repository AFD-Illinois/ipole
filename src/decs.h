#ifndef DECS_H
#define DECS_H

#include "constants.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NDIM	4

#define xstr(s) str(s)
#define str(s) #s

#define STRLEN (2048)

#define NIMG (4+1) // Stokes vector and faraday depth
#define imgindex(n,i,j) (((n)*NX+(i))*NY+(j))

#define SIGMA_CUT (1.)

#define sign(x) (((x) < 0) ? -1 : ((x) > 0))
#define delta(x, y) ( ((x) == (y)) ? 1 : 0 )

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

#define DLOOP  for(k=0;k<NDIM;k++)for(l=0;l<NDIM;l++)
#define IMLOOP for(int i=0;i<NDIM;i++) for(int j=0;j<NDIM;j++)
#define MULOOP for(int mu=0;mu<NDIM;mu++)                                        
#define MUNULOOP for(int mu=0;mu<NDIM;mu++) for(int nu=0;nu<NDIM;nu++)

#define P4VEC(S,X) fprintf(stderr, "%s: %g %g %g %g\n", S, X[0], X[1], X[2], X[3]);

#define MAXNSTEP 50000

// QU convention sets whether QU/EVPA is measured as
//   - 0  "East of North"  or  "observer convention"
//   - 1  "North of West"
#define QU_CONVENTION 0

struct of_traj {
  double dl;
  double X[NDIM];
  double Kcon[NDIM];
  double Xhalf[NDIM];
  double Kconhalf[NDIM];
} traj[MAXNSTEP];
#pragma omp threadprivate(traj)

// used for slow light
struct of_image {
  int nstep;
  double intensity;
  double tau;
  double tauF;
  double complex N_coord[NDIM][NDIM];
};

/** model-independent subroutines **/
/* core routines */
void init_XK (int i, int j, double Xcam[4], double fovx, double fovy,
              double X[4], double Kcon[4], double rotcam, double xoff,
              double yoff);

double approximate_solve (double Ii, double ji, double ki, double jf, double kf,
                   double dl, double *tau);
void get_jkinv (double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv);
int stop_forward_integration (double X[NDIM], double Kcon[NDIM], double Xcam[NDIM]);
int stop_backward_integration (double X[NDIM], double Kcon[NDIM],
                           double Xcam[NDIM]);

/* tetrad related */
void coordinate_to_tetrad (double Ecov[NDIM][NDIM], double K[NDIM],
                      double K_tetrad[NDIM]);
void tetrad_to_coordinate (double Ecov[NDIM][NDIM], double K_tetrad[NDIM],
                      double K[NDIM]);
void make_camera_tetrad (double X[NDIM], double Econ[NDIM][NDIM],
                    double Ecov[NDIM][NDIM]);
void make_plasma_tetrad (double Ucon[NDIM], double Kcon[NDIM], double Bcon[NDIM],
                    double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM],
                    double Ecov[NDIM][NDIM]);

/* radiation */
double Bnu_inv (double nu, double Thetae);
double jnu_inv (double nu, double Thetae, double Ne, double B, double theta);
double get_fluid_nu (double Kcon[NDIM], double Ucov[NDIM]);
double get_bk_angle (double X[NDIM], double Kcon[NDIM], double Ucov[NDIM],
              double Bcon[NDIM], double Bcov[NDIM]);

/* emissivity */
double jnu_synch (double nu, double Ne, double Thetae, double B, double theta);

void init_N(double Xi[NDIM],double Kconi[NDIM],double complex Ncon[NDIM][NDIM]);
void evolve_N(double Xi[NDIM],double Kconi[NDIM],
    double Xf[NDIM],double Kconf[NDIM],
    double Xhalf[NDIM],double Kconhalf[NDIM],
    double dlam,
    double complex N_coord[NDIM][NDIM],
    double *tauF);
void project_N(double X[NDIM],double Kcon[NDIM],
    double complex Ncon[NDIM][NDIM],
    double *Stokes_I, double *Stokes_Q,double *Stokes_U,double *Stokes_V, double rotcam);

#endif // DECS_H
