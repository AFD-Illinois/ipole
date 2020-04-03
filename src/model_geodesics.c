
/*
 * model-dependent routines for integrating geodesics,
 * including:
 *
 * stop_backward_integration
 * stepsize
 */

#include "model_geodesics.h"

#include "geodesics.h"
#include "coordinates.h"
#include "geometry.h"
#include "tetrads.h"
#include "decs.h"
#include "model.h"
#include "coordinates.h"

/*
 * Trace a single geodesic
 * Takes a starting location and frequency in *units of electron mass frequency*,
 * as well as a (long enough!) array of trajectory location structs of_traj
 */
int trace_geodesic(double X[NDIM], double Kcon[NDIM], struct of_traj *traj) {

  double Xhalf[NDIM], Kconhalf[NDIM];

  // Setup geodesic
  MULOOP {
    Xhalf[mu] = X[mu];
    Kconhalf[mu] = Kcon[mu];
  }
  int nstep = 0;

  // Integrate backwards
  while (!stop_backward_integration(X, Xhalf, Kcon)) {
    // Set the current position and wavevector
    MULOOP {
      traj[nstep].X[mu] = X[mu];
      traj[nstep].Kcon[mu] = Kcon[mu];
      traj[nstep].Xhalf[mu] = Xhalf[mu];
      traj[nstep].Kconhalf[mu] = Kconhalf[mu];
    }

    /* This stepsize function can be troublesome inside of R = 2M,
       and should be used cautiously in this region. */
    double dl = stepsize(X, Kcon);

    /* Geodesics in ipole are integrated using
     * dx^\mu/d\lambda = k^\mu
     * The positions x^mu are in simulation units, since different
     * coordinates sometimes have different units (e.g. x^r, x^\theta)

     * The convention we have adopted is that:
     * E = -u^\mu k_\mu,
     * which is always photon energy measured by an observer with four-velocity u^\mu,
     * is in units of *electron rest-mass energy*.
     * This implies that dl is *not* in cgs units, but in weird hybrid units.
     * This line sets dl to be in cgs units.
     *
     * Set the *next* step's size, so that when integrating back it is
     * the *last* step's size
     */
    traj[nstep+1].dl = dl * L_unit * HPL / (ME * CL * CL);

    /* move photon one step backwards, the procecure updates X
       and Kcon full step and returns also values in the middle */
    push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
    nstep++;

    if (nstep > MAXNSTEP - 2) {
      fprintf(stderr, "MAXNSTEP exceeded!\n");
      break;
    }
  }

  // Don't record final step because it violated "stop" condition
  nstep--;

  //fprintf(stderr, "Geodesic nsteps: %d\n", nstep);

  return nstep;
}

/*
 * Initialize a geodesic from the camera
 * This takes the parameters struct directly since most of them are
 * camera parameters anyway
 */
void init_XK(int i, int j, int nx, int ny, double Xcam[NDIM],
             Params params, double fovx, double fovy,
             double X[NDIM], double Kcon[NDIM])
{
  double Econ[NDIM][NDIM];
  double Ecov[NDIM][NDIM];
  double Kcon_tetrad[NDIM];

  make_camera_tetrad(Xcam, Econ, Ecov);

  // Construct outgoing wavevectors
  // xoff: allow arbitrary offset for e.g. ML training imgs
  // +0.5: project geodesics from px centers
  // -0.01: Prevent nasty artefact at 0 spin/phicam
  // xoff/yoff are separated to keep consistent behavior between refinement levels
  double dxoff = (i+0.5-0.01)/nx - 0.5 + params.xoff/params.nx;
  double dyoff = (j+0.5)/ny - 0.5 + params.yoff/params.ny;
  Kcon_tetrad[0] = 0.;
  Kcon_tetrad[1] = (dxoff*cos(params.rotcam) - dyoff*sin(params.rotcam)) * fovx;
  Kcon_tetrad[2] = (dxoff*sin(params.rotcam) + dyoff*cos(params.rotcam)) * fovy;
  Kcon_tetrad[3] = 1.;

  /* normalize */
  null_normalize(Kcon_tetrad, 1.);

  /* translate into coordinate frame */
  tetrad_to_coordinate(Econ, Kcon_tetrad, Kcon);

  /* set position */
  MULOOP X[mu] = Xcam[mu];
}


/* condition for stopping the backwards-in-lambda
   integration of the photon geodesic */

int stop_backward_integration(double X[NDIM], double Xhalf[NDIM], double Kcon[NDIM])
{
  // Necessary geometric stop conditions
  if ((X[1] > log (1.1 * rmax_geo) && Kcon[1] < 0.) || // Stop either beyond rmax_geo
      X[1] < log (1.05 * Rh)) { // Or right near the horizon
    return (1);
  }

  // Additional stop condition for thin disks: the disk is opaque
#if THIN_DISK
  static int n_left = -1;
#pragma omp threadprivate(n_left)

  if (thindisk_region(X, Xhalf) && n_left < 0) { // Set timer when we reach disk
    n_left = 2;
    return 0;
  } else if (n_left < 0) { // Otherwise continue normally
    return 0;
  } else if (n_left > 0) { // Or decrement the timer if we need
    n_left--;
    return 0;
  } else {  // If timer is 0, stop and reset
    n_left = -1;
    return 1;
  }
#endif

  return (0);
}

double stepsize(double X[NDIM], double Kcon[NDIM])
{
  double eps = 0.01; // TODO take as parameter w/MAXNSTEP?

  double dl, dlx1, dlx2, dlx3;
  double idlx1,idlx2,idlx3 ;

  dlx1 = eps / (fabs(Kcon[1]) + SMALL*SMALL) ;
  //dlx2 = EPS * GSL_MIN(X[2], 1. - X[2]) / (fabs(Kcon[2]) + SMALL*SMALL) ;
  dlx2 = eps * MIN(X[2], stopx[2] - X[2]) / (fabs(Kcon[2]) + SMALL*SMALL) ;
  dlx3 = eps / (fabs(Kcon[3]) + SMALL*SMALL) ;

  idlx1 = 1./(fabs(dlx1) + SMALL*SMALL) ;
  idlx2 = 1./(fabs(dlx2) + SMALL*SMALL) ;
  idlx3 = 1./(fabs(dlx3) + SMALL*SMALL) ;

  dl = 1. / (idlx1 + idlx2 + idlx3) ;

  //dl = MIN(MIN(1./dlx1, 1./idlx2), 1./dlx3);
  //printf("idlx = %e %e %e\n", 1./idlx1, 1./idlx2, 1./idlx3);

  //dl = MIN(dl, 0.5*DTd/Kcon[0]);

  return dl;
}

