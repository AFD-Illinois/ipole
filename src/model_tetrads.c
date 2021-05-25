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
 *
 * make orthonormal basis for plasma frame.
 * e^0 along U
 * e^2 along b
 * e^3 along spatial part of K
 * 
 * Returns flag for whether the tetrad is suspicious.
 * Ideally ipole should crash on these errors but there are a lot of corner cases...
 */
int make_plasma_tetrad(double ucon[NDIM], double Kcon[NDIM], double bcon[NDIM],
                        double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM],
                        double Ecov[NDIM][NDIM])
{
  // start w/ time component parallel to U
  set_Econ_from_trial(Econ[0], 0, ucon);
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
  set_Econ_from_trial(Econ[2], 2, bcon);

  // project out econ0,3
  project_out(Econ[2], Econ[0], Gcov);
  project_out(Econ[2], Econ[3], Gcov);
  normalize(Econ[2], Gcov);

  /*** done w/ basis vector 2 ***/

  // whatever is left is econ1
  for (int k = 0; k < 4; k++) /* trial vector */
    Econ[1][k] = 1.;
  // project out econ[0-2]
  project_out(Econ[1], Econ[0], Gcov);
  project_out(Econ[1], Econ[2], Gcov);
  project_out(Econ[1], Econ[3], Gcov);
  normalize(Econ[1], Gcov);

  int oddflag = 0;

  /* check handedness */
  double dot;
  if (check_handedness(Econ, Gcov, &dot)) {
    oddflag |= 16;
  }

  // less restrictive condition on geometry for eKS coordinates which are
  // used when the exotic is expected.
  if ((fabs(fabs(dot) - 1.) > 1.e-10 && use_eKS_internal == 0)
      || (fabs(fabs(dot) - 1.) > 1.e-7 && use_eKS_internal == 1)) {
    oddflag |= 1;
  }

  /* we expect dot = 1. for right-handed system.
   If not, flip Econ[1] to make system right-handed. */
  if (dot < 0.) {
    for (int k = 0; k < 4; k++)
      Econ[1][k] *= -1.;
  }

  /*** done w/ basis vector 3 ***/

  /* now make covariant version */
  for (int k = 0; k < 4; k++) {
    /* lower coordinate basis index */
    flip_index(Econ[k], Gcov, Ecov[k]);
  }

  /* then raise tetrad basis index */
  for (int l = 0; l < 4; l++) {
    Ecov[0][l] *= -1.;
  }

  // For paranoia could run check_ortho here

  return oddflag;
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
 * Points the camera so that the angular momentum k_{th,phi} at FOV center is 0
 */

int make_camera_tetrad(double X[NDIM], double Econ[NDIM][NDIM],
                        double Ecov[NDIM][NDIM])
{
  double Gcov[NDIM][NDIM], Gcon[NDIM][NDIM];
  gcov_func(X, Gcov);
  gcon_func(Gcov, Gcon);
  double Ucam[NDIM], Kcon[NDIM], trial[NDIM];

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

  int sing = make_plasma_tetrad(Ucam, Kcon, trial, Gcov, Econ, Ecov);
#if DEBUG
  if(sing) {
    fprintf(stderr, "\nError making Camera tetrad, something is wrong!\n");
    fprintf(stderr, "Used the following vectors:\n");
    print_vector("X", X);
    print_vector("Ucam", Ucam);
    print_vector("Kcon", Kcon);
    print_vector("trial", trial);
    print_matrix("gcov", Gcov);
    print_matrix("gcon", Gcon);
  }
#endif
  return sing;
}

/*
 * Make orthonormal basis for camera frame -- old implementation
 *
 * e^0 along Ucam
 * e^3 outward (!) along radius vector
 * e^2 toward north pole of coordinate system ("y" in the image plane)
 * e^1 in the remaining direction ("x" in the image plane)
 *
 * This combination measures the final Stokes parameters correctly (IEEE/IAS).
 * These values are then translated if a different convention is to be output.
 *
 * Points the camera so that the *contravariant wavevector* k^{th,phi} = 0
 */

int make_camera_tetrad_old(double X[NDIM], double Econ[NDIM][NDIM],
                        double Ecov[NDIM][NDIM])
{
  double Gcov[NDIM][NDIM], Gcon[NDIM][NDIM];
  gcov_func(X, Gcov);
  gcon_func(Gcov, Gcon);
  double Ucam[NDIM], Kcon[NDIM], trial[NDIM];

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

  int sing = make_plasma_tetrad(Ucam, Kcon, trial, Gcov, Econ, Ecov);
#if DEBUG
  if(sing) {
    fprintf(stderr, "\nError making Camera tetrad, something is wrong!\n");
    fprintf(stderr, "Used the following vectors:\n");
    print_vector("X", X);
    print_vector("Ucam", Ucam);
    print_vector("Kcon", Kcon);
    print_vector("trial", trial);
    print_matrix("gcov", Gcov);
    print_matrix("gcon", Gcon);
  }
#endif
  return sing;
}