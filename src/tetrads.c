/*
 * all model-independent functions related to creation and
 * manipulation of tetrads
 */

#include "decs.h"
#include "geometry.h"

static inline double determinant3(double a00, double a01, double a02,
                                  double a10, double a11, double a12,
                                  double a20, double a21, double a22)
{
  return a00 * (a11 * a22 - a12 * a21)
       - a01 * (a10 * a22 - a12 * a20)
       + a02 * (a10 * a21 - a11 * a20);
}

static inline double determinant4(double A[NDIM][NDIM])
{
  return A[0][0] * determinant3(A[1][1], A[1][2], A[1][3],
                                A[2][1], A[2][2], A[2][3],
                                A[3][1], A[3][2], A[3][3])
       - A[0][1] * determinant3(A[1][0], A[1][2], A[1][3],
                                A[2][0], A[2][2], A[2][3],
                                A[3][0], A[3][2], A[3][3])
       + A[0][2] * determinant3(A[1][0], A[1][1], A[1][3],
                                A[2][0], A[2][1], A[2][3],
                                A[3][0], A[3][1], A[3][3])
       - A[0][3] * determinant3(A[1][0], A[1][1], A[1][2],
                                A[2][0], A[2][1], A[2][2],
                                A[3][0], A[3][1], A[3][2]);
}

/* input and vectors are contravariant (index up) */
void coordinate_to_tetrad(double Ecov[NDIM][NDIM], double K[NDIM],
                          double K_tetrad[NDIM])
{
  int k;

  for (k = 0; k < 4; k++) {
    K_tetrad[k] = Ecov[k][0] * K[0] + Ecov[k][1] * K[1] + Ecov[k][2] * K[2]
        + Ecov[k][3] * K[3];
  }
}

/* input and vectors are contravariant (index up) */
void tetrad_to_coordinate(double Econ[NDIM][NDIM], double K_tetrad[NDIM],
                          double K[NDIM])
{
  int l;

  for (l = 0; l < 4; l++) {
    K[l] = Econ[0][l] * K_tetrad[0] + Econ[1][l] * K_tetrad[1]
        + Econ[2][l] * K_tetrad[2] + Econ[3][l] * K_tetrad[3];
  }

  return;
}

/* 
 * Copy the trial vector into a tetrad basis vector,
 * checking to see if it is null, and if it is null
 * setting to some default value
 */
void set_Econ_from_trial(double Econ[4], int defdir, double trial[4])
{
  double norm = 0.;
  int k;

  for (k = 1; k < 4; k++)
    norm += fabs(trial[k]);
  for (k = 0; k < 4; k++) /* trial vector */
    if (norm <= SMALL) /* bad trial vector; default to defdir */
      Econ[k] = delta(k, defdir);
    else
      Econ[k] = trial[k];

  return;
}

/*
 * Check the handedness of a tetrad basis.
 * Basis is assumed to be in form e^\mu_{(a)} = Econ[a][mu]
 * levi_(ijkl) e0^i e1^j e2^k e3^l will be +1 if spatial
 * components are right-handed, -1 if left-handed.
 * experience suggests that roundoff produces errors of
 * order 10^{-12} in the result.
 */
int check_handedness(double Econ[NDIM][NDIM], double Gcov[NDIM][NDIM], double *dot)
{
  double g = gdet_func(Gcov);
  if(g < 0.) {
#if DEBUG
    fprintf(stderr, "\nEncountered singular gcov checking handedness!\n");
#endif
    return 1;
  }

  /* The Levi-Civita contraction is the determinant of Econ. */
  *dot = g * determinant4(Econ);

  return 0;
}

/*
 * project out vconb from vcona
 * both arguments are index up (contravariant)
 * covariant metric is third argument.
 * overwrite the first argument on return
 */
void project_out(double vcona[NDIM], double vconb[NDIM], double Gcov[4][4])
{

  double adotb, vconb_sq;
  int k, l;

  vconb_sq = 0.;
  for (k = 0; k < 4; k++)
    for (l = 0; l < 4; l++)
      vconb_sq += vconb[k] * vconb[l] * Gcov[k][l];

  adotb = 0.;
  for (k = 0; k < 4; k++)
    for (l = 0; l < 4; l++)
      adotb += vcona[k] * vconb[l] * Gcov[k][l];

  for (k = 0; k < 4; k++)
    vcona[k] -= vconb[k] * adotb / vconb_sq;

  return;
}
