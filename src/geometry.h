/*
 * geometry.h
 *
 *  Created on: Sep 9, 2019
 *      Author: bprather
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "decs.h"

#include "coordinates.h"
#include "matrix.h"
#include "debug_tools.h"

// #include <math.h>
// #include <stdio.h>
// #include <stdlib.h>

/* assumes gcov has been set first; returns sqrt{|g|} */
static inline double gdet_func(double gcov[][NDIM])
{
  int permute[NDIM];
  double gcovtmp[NDIM][NDIM];
  double gdet;

  for (int i = 0; i < NDIM; i++) {
    for (int j = 0; j < NDIM; j++) {
      gcovtmp[i][j] = gcov[i][j];
    }
  }
  if (LU_decompose(gcovtmp, permute) != 0) {
#if DEBUG
    fprintf(stderr, "\nSingluar matrix in gdet_func()!\n");
    print_matrix("gcov", gcov);
#endif
    return -1;
  }

  gdet = 1.;
  for (int i = 0; i < NDIM; i++)
    gdet *= gcovtmp[i][i];

  return (sqrt(fabs(gdet)));
}

/* invert gcov to get gcon */
static inline int gcon_func(double gcov[][NDIM], double gcon[][NDIM])
{
  int sing = invert_matrix(gcov, gcon);
#if DEBUG
  if (sing) {
    fprintf(stderr, "\nSingluar matrix in gcon_func()!\n");
    print_matrix("gcov", gcov);
  }
#endif
  return sing;
}

static inline void flip_index(double ucon[NDIM], double Gcov[NDIM][NDIM], double ucov[NDIM])
{

    ucov[0] = Gcov[0][0] * ucon[0]
  + Gcov[0][1] * ucon[1]
  + Gcov[0][2] * ucon[2]
  + Gcov[0][3] * ucon[3];
    ucov[1] = Gcov[1][0] * ucon[0]
  + Gcov[1][1] * ucon[1]
  + Gcov[1][2] * ucon[2]
  + Gcov[1][3] * ucon[3];
    ucov[2] = Gcov[2][0] * ucon[0]
  + Gcov[2][1] * ucon[1]
  + Gcov[2][2] * ucon[2]
  + Gcov[2][3] * ucon[3];
    ucov[3] = Gcov[3][0] * ucon[0]
  + Gcov[3][1] * ucon[1]
  + Gcov[3][2] * ucon[2]
  + Gcov[3][3] * ucon[3];

    return;
}
static inline void lower(double ucon[NDIM], double Gcov[NDIM][NDIM], double ucov[NDIM])
{ flip_index(ucon, Gcov, ucov); }

/* generic connection routine, using numerical derivatives */
/* Sets the spatial discretization in numerical derivatives : */
#define DEL 1.e-7

static inline void get_connection(double X[NDIM], double conn[NDIM][NDIM][NDIM])
{
  double tmp[NDIM][NDIM][NDIM];
  double Xh[NDIM], Xl[NDIM];
  double gcon[NDIM][NDIM];
  double gcov[NDIM][NDIM];
  double gh[NDIM][NDIM];
  double gl[NDIM][NDIM];

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  for (int k = 0; k < NDIM; k++) {
    for (int l = 0; l < NDIM; l++) {
      Xh[l] = X[l];
      Xl[l] = X[l];
    }

    Xh[k] += DEL;
    Xl[k] -= DEL;
    gcov_func(Xh, gh);
    gcov_func(Xl, gl);

    for (int i = 0; i < NDIM; i++) {
      for (int j = 0; j < NDIM; j++) {
        conn[i][j][k] = (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
      }
    }
  }

  /* now rearrange to find \Gamma_{ijk} */
  for (int i = 0; i < NDIM; i++)
    for (int j = 0; j < NDIM; j++)
      for (int k = 0; k < NDIM; k++)
        tmp[i][j][k] = 0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);

  /* finally, raise index */
  for (int i = 0; i < NDIM; i++) {
    for (int j = 0; j < NDIM; j++) {
      for (int k = 0; k < NDIM; k++) {
        conn[i][j][k] = 0.;
        for (int l = 0; l < NDIM; l++) {
          conn[i][j][k] += gcon[i][l] * tmp[l][j][k];
        }
      }
    }
  }

  /* done! */

}
#undef DEL

/*
 * Normalize input vector so that |v . v| = 1
 * Overwrites input
 */
static inline void normalize(double vcon[NDIM], double Gcov[NDIM][NDIM])
{
  double norm = 0.;

  for (int k = 0; k < 4; k++)
    for (int l = 0; l < 4; l++)
      norm += vcon[k] * vcon[l] * Gcov[k][l];

  norm = sqrt(fabs(norm));
  for (int k = 0; k < 4; k++)
    vcon[k] /= norm;

  return;
}

/* normalize null vector in a tetrad frame */
static inline void null_normalize(double Kcon[NDIM], double fnorm)
{
  double inorm =
    sqrt(Kcon[1] * Kcon[1] + Kcon[2] * Kcon[2] + Kcon[3] * Kcon[3]);

  Kcon[0] = fnorm;
  Kcon[1] *= fnorm / inorm;
  Kcon[2] *= fnorm / inorm;
  Kcon[3] *= fnorm / inorm;
}

#endif /* GEOMETRY_H */
