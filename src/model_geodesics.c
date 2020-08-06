
/*
 * model-dependent routines for integrating geodesics,
 * including:
 *
 * stop_backward_integration
 * stepsize
 */

#include "model_geodesics.h"

#include "coordinates.h"
#include "decs.h"
#include "geodesics.h"
#include "geometry.h"
#include "model.h"
#include "model_tetrads.h"
#include "tetrads.h"

#include "debug_tools.h"

// TODO pick one or runtime
#define STEP_STRICT_MIN 0

/**
 * Trace a single geodesic
 * Takes a starting location and frequency in *units of electron mass frequency*,
 * as well as an array of trajectory location structs of_traj, size MAXNSTEP
 * 
 * TODO this routine definitely doesn't initialize traj[0].
 * 
 * Note convention changes for integrating backwards:
 * * Xhalf & Kconhalf trail X & Kcon in this function, and thereafter lead them
 * * In this function, dl is the length of the step *to* point N;
 *   afterward it is *from* point N onward
 */
int trace_geodesic(double X[NDIM], double Kcon[NDIM], struct of_traj *traj, double eps, int step_max) {

  double Xhalf[NDIM], Kconhalf[NDIM];

  // Save the first step
  // Note "half" values *trail* when integrating away from camera, so they don't
  // make sense to record for the starting location
  MULOOP {
    traj[0].X[mu] = X[mu];
    traj[0].Kcon[mu] = Kcon[mu];
  }
  // Initialize Xhalf for the first while loop check though
  MULOOP Xhalf[mu] = X[mu];

  int nstep = 0;

  // Integrate backwards
  while ( (!stop_backward_integration(X, Xhalf, Kcon)) && (nstep < step_max - 1) ) {
    nstep++;

    /* This stepsize function can be troublesome inside of R = 2M,
       and should be used cautiously in this region. */
    double dl = stepsize(X, Kcon, eps);

    /* Geodesics in ipole are integrated using
     * dx^\mu/d\lambda = k^\mu
     * The positions x^mu are in simulation units, since different
     * coordinates sometimes have different units (e.g. x^r, x^\theta)
     *
     * The convention we have adopted is that:
     * E = -u^\mu k_\mu,
     * which is always photon energy measured by an observer with four-velocity u^\mu,
     * is in units of *electron rest-mass energy*.
     * This implies that dl is *not* in cgs units, but in weird hybrid units.
     * This line sets dl to be in cgs units.
     */

#if INTEGRATOR_TEST
    traj[nstep].dl = dl;
#else
    traj[nstep].dl = dl * L_unit * HPL / (ME * CL * CL);
#endif

    // To print each point (TODO option?)
    // print_vector("X", X);
    // print_vector("Kcon", Kcon);
    // fprintf(stderr, "dl: %f\n", dl);
    // double r, th;
    // bl_coord(X, &r, &th);
    // fprintf(stderr, "BL r, th: %g, %g\n", r, th);
    // double Gcov[NDIM][NDIM];
    // gcov_func(X, Gcov);
    // print_matrix("Gcov", Gcov);
    // getchar();

    /* move photon one step backwards, the procecure updates X
       and Kcon full step and returns also values in the middle */
    push_photon(X, Kcon, -dl, Xhalf, Kconhalf);

    // Set the new position and wavevector
    MULOOP {
      traj[nstep].X[mu] = X[mu];
      traj[nstep].Kcon[mu] = Kcon[mu];
      traj[nstep].Xhalf[mu] = Xhalf[mu];
      traj[nstep].Kconhalf[mu] = Kconhalf[mu];
    }
  }

#if DEBUG
  if (nstep >= step_max-1) {
    fprintf(stderr, "\nMaxNStep exceeded!  Inaccuracies are likely!\n");
  }
#endif

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
  // xoff/yoff are separated to keep consistent behavior between refinement levels
  double dxoff = (i+0.5-0.01)/nx - 0.5 + params.xoff/nx;
  double dyoff = (j+0.5)/ny - 0.5 + params.yoff/ny;
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
  double r, th;
  bl_coord(X, &r, &th);
  if ((r > (1.1 * rmax_geo) && Kcon[1] < 0.) || // Stop either beyond rmax_geo
      r < (1.05 * Rh)) { // Or right near the horizon
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


/**
 * This is slightly modified from ipole pre-polyphony.
 * It is now a strict minimum, which may reduce some extra steps near the pole,
 * but runs the risk of overstepping there.
 * TODO test this clears Angelo's issue & maybe the old thin-disk issue
 */
double stepsize(double X[NDIM], double Kcon[NDIM], double eps)
{
  double dl;
  double dlx1 = eps / (fabs(Kcon[1]) + SMALL*SMALL);
  double dlx2 = eps * fmin(fabs(X[2]), fabs(stopx[2] - X[2])) / (fabs(Kcon[2]) + SMALL*SMALL);
  double dlx3 = eps / (fabs(Kcon[3]) + SMALL*SMALL);

  if (STEP_STRICT_MIN) {
    dl = fmin(fmin(dlx1, dlx2), dlx3);
  } else {
    double idlx1 = 1./(fabs(dlx1) + SMALL*SMALL) ;
    double idlx2 = 1./(fabs(dlx2) + SMALL*SMALL) ;
    double idlx3 = 1./(fabs(dlx3) + SMALL*SMALL) ;

    dl = 1. / (idlx1 + idlx2 + idlx3) ;
  }

  return dl;
}

