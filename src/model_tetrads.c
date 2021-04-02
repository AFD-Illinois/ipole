/*

 model-dependent functions related to creation and manipulation of tetrads

 */

#include "model_tetrads.h"

#include "decs.h"
#include "coordinates.h"
#include "geometry.h"
#include "tetrads.h"

#include "debug_tools.h"

/** tetrad making routines **/

/* 
 * econ/ecov index key:
 * Econ[k][l]
 * k: index attached to tetrad basis
 * index down
 * l: index attached to coordinate basis
 * index up
 * Ecov switches both indices
 */

/* 
 * make orthonormal basis for plasma frame.
 * e^0 along U
 * e^2 along b
 * e^3 along spatial part of K
 */
void make_plasma_tetrad(double Ucon[NDIM], double Kcon[NDIM], double Bcon[NDIM],
                        double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM],
                        double Ecov[NDIM][NDIM])
{
  int k, l;

  // start w/ time component parallel to U
  set_Econ_from_trial(Econ[0], 0, Ucon);
  normalize(Econ[0], Gcov);

  /*** done w/ basis vector 0 ***/

  /* now use the trial vector in basis vector 3 */
  /* cast a suspicious eye on the trial vector... */

  set_Econ_from_trial(Econ[3], 3, Kcon);

  /* project out econ0 */
  project_out(Econ[3], Econ[0], Gcov);
  normalize(Econ[3], Gcov);

  /*** done w/ basis vector 3 ***/

  // repeat for x2 unit basis vector
  set_Econ_from_trial(Econ[2], 2, Bcon);

  // project out econ0,3
  project_out(Econ[2], Econ[0], Gcov);
  project_out(Econ[2], Econ[3], Gcov);
  normalize(Econ[2], Gcov);

  /*** done w/ basis vector 2 ***/

  // whatever is left is econ1
  for (k = 0; k < 4; k++) /* trial vector */
    Econ[1][k] = 1.;
  // project out econ[0-2]
  project_out(Econ[1], Econ[0], Gcov);
  project_out(Econ[1], Econ[2], Gcov);
  project_out(Econ[1], Econ[3], Gcov);
  normalize(Econ[1], Gcov);

  /* check handedness */
  double dot = check_handedness(Econ, Gcov);

  // less restrictive condition on geometry for eKS coordinates which are
  // used when the exotic is expected.
  if ((fabs(fabs(dot) - 1.) > 1.e-10 && METRIC_eKS == 0)
      || (fabs(fabs(dot) - 1.) > 1.e-7 && METRIC_eKS == 1)) {
    fprintf(stderr, "that's odd: %g\n", fabs(dot) - 1.);
    fprintf(stderr, "Ucon[] = %e %e %e %e\n", Ucon[0], Ucon[1], Ucon[2],
            Ucon[3]);
    fprintf(stderr, "Kcon[] = %e %e %e %e\n", Kcon[0], Kcon[1], Kcon[2],
            Kcon[3]);
    fprintf(stderr, "Bcon[] = %e %e %e %e\n", Bcon[0], Bcon[1], Bcon[2],
            Bcon[3]);

    double ucov[4];
    flip_index(Ucon, Gcov, ucov);
    double udotu = 0., udotb = 0.;
    MULOOP
    {
      udotu += Ucon[mu] * ucov[mu];
      udotb += Bcon[mu] * ucov[mu];
    }
    fprintf(stderr, "u.u = %g  u.b = %g\n", udotu, udotb);

    int i, j, k;
    double del[NDIM];
    double X[NDIM];
    void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]);
    Xtoijk(X, &i, &j, &k, del);
    fprintf(stderr, "X[]: %g %g %g %g  (%d %d %d)\n", X[0], X[1], X[2], X[3], i,
            j, k);
    exit(-2);
  }

  /* we expect dot = 1. for right-handed system.
   If not, flip Econ[1] to make system right-handed. */
  if (dot < 0.) {
    for (k = 0; k < 4; k++)
      Econ[1][k] *= -1.;
  }

  /*** done w/ basis vector 3 ***/

  /* now make covariant version */
  for (k = 0; k < 4; k++) {

    /* lower coordinate basis index */
    flip_index(Econ[k], Gcov, Ecov[k]);
  }

  /* then raise tetrad basis index */
  for (l = 0; l < 4; l++) {
    Ecov[0][l] *= -1.;
  }

  // For paranoia could run check_ortho here

  return;
}

/*
 * Make orthonormal basis for camera frame.
 *
 * e^0 along Ucam
 * e^3 outward (!) along radius vector
 * e^2 toward north pole of coordinate system ("y" in the image plane)
 * e^1 in the remaining direction ("x" in the image plane)
 *
 * This combination measures the final Stokes parameters correctly (IEEE/IAS).
 * These values are then translated if a different convention is to be output.
 *
 * Points the camera so that the angular momentum k_phi at FOV center is 0
 */

void make_camera_tetrad(double X[NDIM], double Econ[NDIM][NDIM],
                        double Ecov[NDIM][NDIM])
{
  double Gcov[NDIM][NDIM], Gcon[NDIM][NDIM];
  gcov_func(X, Gcov);
  gcon_func(Gcov, Gcon);
  double Ucam[NDIM], Kcon[NDIM], trial[NDIM];

#ifdef CAMERA_CENTER_ZAMO
  // center the camera according to impact parameter, i.e., make it
  // so that Kcontetrad = ( 1, 0, 0, 1 ) corresponds to an outgoing
  // wavevector with zero angular momentum / zero impact parameter.

  // use normal observer velocity. this forces (Gcov.Econ[0])[3] = 0.
  trial[0] = -1.;
  trial[1] = 0.;
  trial[2] = 0.;
  trial[3] = 0.;
  flip_index(trial, Gcon, Ucam);

  // set Kcon (becomes Econ[3][mu]) outward directed with central 
  // pixel k_phi = 0. this ensures that a photon with zero impact 
  // parameter will be in the center of the field of view.
  trial[0] = 1.;
  trial[1] = 1.;
  trial[2] = 0.;
  trial[3] = 0.;
  flip_index(trial, Gcon, Kcon);

  // set the y camera direction to be parallel to the projected
  // spin axis of the black hole (on the image plane defined to
  // be normal to the Kcon vector above).
  trial[0] = 0.;
  trial[1] = 0.;
  trial[2] = 1.;
  trial[3] = 0.;

  make_plasma_tetrad(Ucam, Kcon, trial, Gcov, Econ, Ecov);

#else
  // old centering method

  Ucam[0] = 1.;
  Ucam[1] = 0.;
  Ucam[2] = 0.;
  Ucam[3] = 0.;

  trial[0] = 1.;
  trial[1] = 1.;
  trial[2] = 0.;
  trial[3] = 0.;
  flip_index(trial, Gcon, Kcon);

  trial[0] = 0.;
  trial[1] = 0.;
  trial[2] = 1.;
  trial[3] = 0.;

  make_plasma_tetrad(Ucam, Kcon, trial, Gcov, Econ, Ecov);

#endif

}

