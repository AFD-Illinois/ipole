#include "geometry.h"

#include "decs.h"
#include "coordinates.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* 
 find determinant of covariant metric gcov
 find gcon, as numerical inverse of gcov
 find connection, as numerical derivative of gcov
 These routines are taken mostly out of HARM.
 CFG 21 July 06
 MM 11 July 17
 */

/* assumes gcov has been set first; returns sqrt{|g|} */
double gdet_func(double gcov[][NDIM])
{

  int i, j;
  int permute[NDIM];
  double gcovtmp[NDIM][NDIM];
  double gdet;
  int LU_decompose(double A[][NDIM], int permute[]);

  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++) {
      gcovtmp[i][j] = gcov[i][j];
    }
  if (LU_decompose(gcovtmp, permute) != 0) {
    fprintf(stderr, "gdet_func(): singular matrix encountered! \n");
    exit(1);
  }
  gdet = 1.;

  for (i = 0; i < NDIM; i++)
    gdet *= gcovtmp[i][i];

  return (sqrt(fabs(gdet)));
}

/* invert gcov to get gcon */
int gcon_func(double gcov[][NDIM], double gcon[][NDIM])
{
  int invert_matrix(double Am[][NDIM], double Aminv[][NDIM]);

  int sing = invert_matrix(gcov, gcon);
  if (sing) {
    for (int mu = 0; mu < NDIM; mu++) {
      for (int nu = 0; nu < NDIM; nu++) {
        printf("gcov[%i][%i] = %e\n", mu, nu, gcov[mu][nu]);
      }
    }
  }

  return sing;

  /* done! */
}

inline void flip_index(double ucon[NDIM], double Gcov[NDIM][NDIM], double ucov[NDIM])
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

/* generic connection routine, using numerical derivatives */
/* Sets the spatial discretization in numerical derivatives : */
#define DEL 1.e-7

void get_connection(double X[NDIM], double conn[NDIM][NDIM][NDIM])
{
  int i, j, k, l;
  double tmp[NDIM][NDIM][NDIM];
  double Xh[NDIM], Xl[NDIM];
  double gcon[NDIM][NDIM];
  double gcov[NDIM][NDIM];
  double gh[NDIM][NDIM];
  double gl[NDIM][NDIM];

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  for (k = 0; k < NDIM; k++) {
    for (l = 0; l < NDIM; l++)
      Xh[l] = X[l];
    for (l = 0; l < NDIM; l++)
      Xl[l] = X[l];
    Xh[k] += DEL;
    Xl[k] -= DEL;
    gcov_func(Xh, gh);
    gcov_func(Xl, gl);

    for (i = 0; i < NDIM; i++) {
      for (j = 0; j < NDIM; j++) {
        conn[i][j][k] = (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
      }
    }
  }

  /* now rearrange to find \Gamma_{ijk} */
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      for (k = 0; k < NDIM; k++)
        tmp[i][j][k] = 0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);

  /* finally, raise index */
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      for (k = 0; k < NDIM; k++) {
        conn[i][j][k] = 0.;
        for (l = 0; l < NDIM; l++)
          conn[i][j][k] += gcon[i][l] * tmp[l][j][k];
      }

  /* done! */

}
#undef DEL

/*
 * Normalize input vector so that |v . v| = 1
 * Overwrites input
 */
void normalize(double vcon[NDIM], double Gcov[NDIM][NDIM])
{
  int k, l;
  double norm;

  norm = 0.;
  for (k = 0; k < 4; k++)
    for (l = 0; l < 4; l++)
      norm += vcon[k] * vcon[l] * Gcov[k][l];

  norm = sqrt(fabs(norm));
  for (k = 0; k < 4; k++)
    vcon[k] /= norm;

  return;
}

/* normalize null vector in a tetrad frame */
void null_normalize(double Kcon[NDIM], double fnorm)
{
  double inorm;

  inorm =
    sqrt(Kcon[1] * Kcon[1] + Kcon[2] * Kcon[2] + Kcon[3] * Kcon[3]);

  Kcon[0] = fnorm;
  Kcon[1] *= fnorm / inorm;
  Kcon[2] *= fnorm / inorm;
  Kcon[3] *= fnorm / inorm;
}

/* the completely antisymmetric symbol; not a tensor
 in the coordinate basis */
int levi_civita(int i, int j, int k, int l)
{
  int index[NDIM], n, do_sort, n_perm, val, n_swap;
  if (i == j || i == k || i == l || j == k || j == l || k == l) {
    return 0;
  } else {
    index[0] = i;
    index[1] = j;
    index[2] = k;
    index[3] = l;
    do_sort = 1;
    n_perm = 0;
    while (do_sort) {
      n_swap = 0;
      for (n = 0; n < NDIM - 1; n++) {
        if (index[n] > index[n + 1]) {
          n_perm++;
          n_swap++;
          val = index[n];
          index[n] = index[n + 1];
          index[n + 1] = val;
        }
      }
      do_sort = n_swap;
    }
    return (n_perm % 2) ? -1 : 1;
  }
}

// Coordinate functions not dependent on system
double theta_func(double X[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);
  return th;
}


/*
 invert_matrix():

 Uses LU decomposition and back substitution to invert a matrix
 A[][] and assigns the inverse to Ainv[][].  This routine does not
 destroy the original matrix A[][].

 Returns (1) if a singular matrix is found,  (0) otherwise.
 */
int invert_matrix(double Am[][NDIM], double Aminv[][NDIM])
{

  int i, j;
  int n = NDIM;
  int permute[NDIM];
  double dxm[NDIM], Amtmp[NDIM][NDIM];
  int LU_decompose(double A[][NDIM], int permute[]);
  void LU_substitution(double A[][NDIM], double B[], int permute[]);

  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++) {
      Amtmp[i][j] = Am[i][j];
    }

  // Get the LU matrix:
  if (LU_decompose(Amtmp, permute) != 0) {
    fprintf(stderr, "invert_matrix(): singular matrix encountered! \n");
    return (1);
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      dxm[j] = 0.;
    }
    dxm[i] = 1.;

    /* Solve the linear system for the i^th column of the inverse matrix: :  */
    LU_substitution(Amtmp, dxm, permute);

    for (j = 0; j < n; j++) {
      Aminv[j][i] = dxm[j];
    }

  }

  return (0);
}

/*
 LU_decompose():

 Performs a LU decomposition of the matrix A using Crout's method
 with partial implicit pivoting.  The exact LU decomposition of the
 matrix can be reconstructed from the resultant row-permuted form via
 the integer array permute[]

 The algorithm closely follows ludcmp.c of "Numerical Recipes
 in C" by Press et al. 1992.

 This will be used to solve the linear system  A.x = B

 Returns (1) if a singular matrix is found,  (0) otherwise.
*/
int LU_decompose(double A[][NDIM], int permute[])
{

  const double absmin = 1.e-30; /* Value used instead of 0 for singular matrices */

  //static double row_norm[NDIM];
  double row_norm[NDIM];
  double absmax, maxtemp;

  int i, j, k, max_row;
  int n = NDIM;

  max_row = 0;

  /* Find the maximum elements per row so that we can pretend later
   we have unit-normalized each equation: */

  for (i = 0; i < n; i++) {
    absmax = 0.;

    for (j = 0; j < n; j++) {

      maxtemp = fabs(A[i][j]);

      if (maxtemp > absmax) {
        absmax = maxtemp;
      }
    }

    /* Make sure that there is at least one non-zero element in this row: */
    if (absmax == 0.) {
      fprintf(stderr, "LU_decompose(): row-wise singular matrix!\n");
      return (1);
    }

    row_norm[i] = 1. / absmax; /* Set the row's normalization factor. */
  }

  /* The following the calculates the matrix composed of the sum
   of the lower (L) tridagonal matrix and the upper (U) tridagonal
   matrix that, when multiplied, form the original maxtrix.
   This is what we call the LU decomposition of the maxtrix.
   It does this by a recursive procedure, starting from the
   upper-left, proceding down the column, and then to the next
   column to the right.  The decomposition can be done in place
   since element {i,j} require only those elements with {<=i,<=j}
   which have already been computed.
   See pg. 43-46 of "Num. Rec." for a more thorough description.
   */

  /* For each of the columns, starting from the left ... */
  for (j = 0; j < n; j++) {

    /* For each of the rows starting from the top.... */

    /* Calculate the Upper part of the matrix:  i < j :   */
    for (i = 0; i < j; i++) {
      for (k = 0; k < i; k++) {
        A[i][j] -= A[i][k] * A[k][j];
      }
    }

    absmax = 0.0;

    /* Calculate the Lower part of the matrix:  i <= j :   */

    for (i = j; i < n; i++) {

      for (k = 0; k < j; k++) {
        A[i][j] -= A[i][k] * A[k][j];
      }

      /* Find the maximum element in the column given the implicit
       unit-normalization (represented by row_norm[i]) of each row:
       */
      maxtemp = fabs(A[i][j]) * row_norm[i];

      if (maxtemp >= absmax) {
        absmax = maxtemp;
        max_row = i;
      }

    }

    /* Swap the row with the largest element (of column j) with row_j.  absmax
     This is the partial pivoting procedure that ensures we don't divide
     by 0 (or a small number) when we solve the linear system.
     Also, since the procedure starts from left-right/top-bottom,
     the pivot values are chosen from a pool involving all the elements
     of column_j  in rows beneath row_j.  This ensures that
     a row  is not permuted twice, which would mess things up.
     */
    if (max_row != j) {

      /* Don't swap if it will send a 0 to the last diagonal position.
       Note that the last column cannot pivot with any other row,
       so this is the last chance to ensure that the last two
       columns have non-zero diagonal elements.
       */

      if ((j == (n - 2)) && (A[j][j + 1] == 0.)) {
        max_row = j;
      } else {
        for (k = 0; k < n; k++) {

          maxtemp = A[j][k];
          A[j][k] = A[max_row][k];
          A[max_row][k] = maxtemp;

        }

        /* Don't forget to swap the normalization factors, too...
         but we don't need the jth element any longer since we
         only look at rows beneath j from here on out.
         */
        row_norm[max_row] = row_norm[j];
      }
    }

    /* Set the permutation record s.t. the j^th element equals the
     index of the row swapped with the j^th row.  Note that since
     this is being done in successive columns, the permutation
     vector records the successive permutations and therefore
     index of permute[] also indexes the chronology of the
     permutations.  E.g. permute[2] = {2,1} is an identity
     permutation, which cannot happen here though.
     */

    permute[j] = max_row;

    if (A[j][j] == 0.) {
      A[j][j] = absmin;
    }

    /* Normalize the columns of the Lower tridiagonal part by their respective
     diagonal element.  This is not done in the Upper part because the
     Lower part's diagonal elements were set to 1, which can be done w/o
     any loss of generality.
     */
    if (j != (n - 1)) {
      maxtemp = 1. / A[j][j];

      for (i = (j + 1); i < n; i++) {
        A[i][j] *= maxtemp;
      }
    }

  }

  return (0);

  /* End of LU_decompose() */

}

/*
 LU_substitution():

 Performs the forward (w/ the Lower) and backward (w/ the Upper)
 substitutions using the LU-decomposed matrix A[][] of the original
 matrix A' of the linear equation:  A'.x = B.  Upon entry, A[][]
 is the LU matrix, B[] is the source vector, and permute[] is the
 array containing order of permutations taken to the rows of the LU
 matrix.  See LU_decompose() for further details.

 Upon exit, B[] contains the solution x[], A[][] is left unchanged.
*/
void LU_substitution(double A[][NDIM], double B[], int permute[])
{
  int i, j;
  int n = NDIM;
  double tmpvar;

  /* Perform the forward substitution using the LU matrix.
   */
  for (i = 0; i < n; i++) {

    /* Before doing the substitution, we must first permute the
     B vector to match the permutation of the LU matrix.
     Since only the rows above the currrent one matter for
     this row, we can permute one at a time.
     */
    tmpvar = B[permute[i]];
    B[permute[i]] = B[i];
    for (j = (i - 1); j >= 0; j--) {
      tmpvar -= A[i][j] * B[j];
    }
    B[i] = tmpvar;
  }

  /* Perform the backward substitution using the LU matrix.
   */
  for (i = (n - 1); i >= 0; i--) {
    for (j = (i + 1); j < n; j++) {
      B[i] -= A[i][j] * B[j];
    }
    B[i] /= A[i][i];
  }

  /* End of LU_substitution() */
}
