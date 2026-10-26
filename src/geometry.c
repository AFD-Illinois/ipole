#include "geometry.h"

#include "decs.h"

#include "coordinates.h"
#include "debug_tools.h"

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

int invert_matrix(double Am[NDIM][NDIM], double Aminv[NDIM][NDIM]);
int LU_decompose(double A[NDIM][NDIM], int permute[NDIM]);
void LU_substitution(double A[NDIM][NDIM], double B[NDIM], int permute[NDIM]);

/* assumes gcov has been set first; returns sqrt{|g|} */
double gdet_func(double gcov[NDIM][NDIM])
{
  int i, j;
  int permute[NDIM];
  double gcovtmp[NDIM][NDIM];
  double gdet;

  for (i = 0; i < NDIM; i++) {
    for (j = 0; j < NDIM; j++) {
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
  for (i = 0; i < NDIM; i++)
    gdet *= gcovtmp[i][i];

  return (sqrt(fabs(gdet)));
}

/* invert gcov to get gcon */
int gcon_func(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM])
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
  if(gcon_func(gcov, gcon)) {
    fprintf(stderr, "Encountered singluar gcov when getting connection!\n");
    print_vector("X", X);
    print_matrix("gcov", gcov);
  }

  /* Every metric implemented by gcov_func is stationary and axisymmetric:
   * derivatives with respect to t (0) and phi (3) vanish exactly.  Only
   * radial and polar finite differences are needed. */
  for (i = 0; i < NDIM; ++i)
    for (j = 0; j < NDIM; ++j) {
      conn[i][j][0] = 0.;
      conn[i][j][3] = 0.;
    }

  for (k = 1; k <= 2; k++) {
    for (l = 0; l < NDIM; l++)
      Xh[l] = X[l];
    for (l = 0; l < NDIM; l++)
      Xl[l] = X[l];
    Xh[k] += DEL;
    Xl[k] -= DEL;
    gcov_func(Xh, gh);
    gcov_func(Xl, gl);

    for (i = 0; i < NDIM; i++) {
      for (j = i; j < NDIM; j++) {
        conn[i][j][k] = (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
        conn[j][i][k] = conn[i][j][k];
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

/* Small forward-mode dual number used to differentiate the FMKS metric in
 * the two nontrivial coordinate directions without repeated metric calls. */
typedef struct { double v, d1, d2; } dual2;

static inline dual2 d2_const(double v) { return (dual2){v, 0., 0.}; }
static inline dual2 d2_add(dual2 a, dual2 b)
{ return (dual2){a.v+b.v, a.d1+b.d1, a.d2+b.d2}; }
static inline dual2 d2_sub(dual2 a, dual2 b)
{ return (dual2){a.v-b.v, a.d1-b.d1, a.d2-b.d2}; }
static inline dual2 d2_mul(dual2 a, dual2 b)
{ return (dual2){a.v*b.v, a.d1*b.v+a.v*b.d1, a.d2*b.v+a.v*b.d2}; }
static inline dual2 d2_scale(dual2 a, double b)
{ return (dual2){a.v*b, a.d1*b, a.d2*b}; }
static inline dual2 d2_div(dual2 a, dual2 b)
{
  double inv = 1. / b.v, inv2 = inv * inv;
  return (dual2){a.v*inv, (a.d1*b.v-a.v*b.d1)*inv2,
                  (a.d2*b.v-a.v*b.d2)*inv2};
}
static inline dual2 d2_exp(dual2 a)
{
  double v = exp(a.v);
  return (dual2){v, v*a.d1, v*a.d2};
}
static inline dual2 d2_sin(dual2 a)
{
  double s = sin(a.v), c = cos(a.v);
  return (dual2){s, c*a.d1, c*a.d2};
}
static inline dual2 d2_cos(dual2 a)
{
  double c = cos(a.v), s = -sin(a.v);
  return (dual2){c, s*a.d1, s*a.d2};
}
static inline dual2 d2_pow(dual2 a, double p)
{
  double v = pow(a.v, p), slope = p*pow(a.v, p-1.);
  return (dual2){v, slope*a.d1, slope*a.d2};
}

static void fmks_metric_and_derivs(double X[NDIM], double gcov[NDIM][NDIM],
                                   double gcon[NDIM][NDIM],
                                   double dg[NDIM][NDIM][NDIM])
{
  const dual2 x1 = {X[1], 1., 0.};
  const dual2 x2 = {X[2], 0., 1.};
  const dual2 one = d2_const(1.), half = d2_const(0.5);
  dual2 y = d2_sub(d2_scale(x2, 2.), one);
  dual2 yp = d2_div(y, d2_const(poly_xt));
  dual2 yalpha = d2_pow(yp, poly_alpha);
  dual2 yalpha_m1 = d2_pow(yp, poly_alpha - 1.);

  dual2 thg = d2_add(d2_scale(x2, M_PI),
      d2_scale(d2_sin(d2_scale(x2, 2.*M_PI)), (1.-hslope)/2.));
  dual2 thj = d2_add(d2_scale(d2_mul(y, d2_add(one,
      d2_scale(yalpha, 1./(poly_alpha+1.)))), poly_norm), d2_const(M_PI/2.));
  dual2 smooth = d2_exp(d2_scale(d2_sub(d2_const(startx[1]), x1), mks_smooth));
  dual2 th = d2_add(thg, d2_mul(smooth, d2_sub(thj, thg)));

  dual2 j21 = d2_scale(d2_mul(smooth, d2_sub(thj, thg)), -mks_smooth);
  dual2 c2p = d2_cos(d2_scale(x2, 2.*M_PI));
  dual2 j22inside = d2_add(d2_const(-M_PI),
      d2_scale(d2_add(one, d2_scale(yalpha, 1./(poly_alpha+1.))),
               2.*poly_norm));
  j22inside = d2_add(j22inside,
      d2_scale(d2_mul(y, yalpha_m1),
               2.*poly_alpha*poly_norm/((poly_alpha+1.)*poly_xt)));
  j22inside = d2_sub(j22inside,
      d2_scale(c2p, (1.-hslope)*M_PI));
  dual2 j22 = d2_add(d2_add(d2_const(M_PI),
                     d2_scale(c2p, (1.-hslope)*M_PI)),
                     d2_mul(smooth, j22inside));

  dual2 r = d2_exp(x1);
  dual2 sth = d2_sin(th), cth = d2_cos(th);
  dual2 s2 = d2_mul(sth, sth);
  dual2 rho2 = d2_add(d2_mul(r, r),
                      d2_scale(d2_mul(cth, cth), a*a));
  dual2 f = d2_div(d2_scale(r, 2.), rho2);

  dual2 g00 = d2_add(d2_const(-1.), f);
  dual2 g01 = d2_mul(f, r);
  dual2 g03 = d2_scale(d2_mul(s2, f), -a);
  dual2 g11 = d2_add(d2_mul(d2_add(one, f), d2_mul(r, r)),
                     d2_mul(rho2, d2_mul(j21, j21)));
  dual2 g12 = d2_mul(rho2, d2_mul(j21, j22));
  dual2 g13 = d2_scale(d2_mul(s2, d2_mul(d2_add(one, f), r)), -a);
  dual2 g22 = d2_mul(rho2, d2_mul(j22, j22));
  dual2 g33 = d2_mul(s2, d2_add(rho2,
      d2_scale(d2_mul(s2, d2_add(one, f)), a*a)));

  /* The inverse Kerr-Schild metric has a closed form.  Transform it back
   * with the inverse FMKS Jacobian instead of running a generic LU inverse. */
  double gkscon[NDIM][NDIM] = {{0.}};
  const double rv = r.v, rho2v = rho2.v, s2v = s2.v, fv = f.v;
  gkscon[0][0] = -(1. + fv);
  gkscon[0][1] = gkscon[1][0] = fv;
  gkscon[1][1] = (rv*rv - 2.*rv + a*a) / rho2v;
  gkscon[1][3] = gkscon[3][1] = a / rho2v;
  gkscon[2][2] = 1. / rho2v;
  gkscon[3][3] = 1. / (rho2v * s2v);
  double invJ[NDIM][NDIM] = {{0.}};
  invJ[0][0] = invJ[3][3] = 1.;
  invJ[1][1] = 1. / rv;
  invJ[2][1] = -j21.v / (rv * j22.v);
  invJ[2][2] = 1. / j22.v;
  for (int i = 0; i < NDIM; ++i)
    for (int j = 0; j < NDIM; ++j) {
      double value = 0.;
      for (int l = 0; l < NDIM; ++l)
        for (int k = 0; k < NDIM; ++k)
          value += invJ[i][l] * gkscon[l][k] * invJ[j][k];
      gcon[i][j] = value;
    }

  dual2 g[NDIM][NDIM];
  for (int i = 0; i < NDIM; ++i)
    for (int j = 0; j < NDIM; ++j)
      g[i][j] = d2_const(0.);
  g[0][0] = g00;
  g[0][1] = g[1][0] = g01;
  g[0][3] = g[3][0] = g03;
  g[1][1] = g11;
  g[1][2] = g[2][1] = g12;
  g[1][3] = g[3][1] = g13;
  g[2][2] = g22;
  g[3][3] = g33;

  for (int i = 0; i < NDIM; ++i)
    for (int j = 0; j < NDIM; ++j) {
      gcov[i][j] = g[i][j].v;
      dg[0][i][j] = 0.;
      dg[1][i][j] = g[i][j].d1;
      dg[2][i][j] = g[i][j].d2;
      dg[3][i][j] = 0.;
    }
}

/* Return Gamma^i_{j l} K^l directly.  Geodesic and polarization evolution
 * only use the connection in this contracted form, so constructing all 64
 * Christoffel symbols is unnecessary.  The implemented coordinate systems
 * are stationary and axisymmetric, leaving only r and theta derivatives. */
void get_connection_k(double X[NDIM], double Kcon[NDIM],
                      double out[NDIM][NDIM])
{
  double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
  double dg[NDIM][NDIM][NDIM] = {{{0.}}};
  double Xh[NDIM], Xl[NDIM], gh[NDIM][NDIM], gl[NDIM][NDIM];
  double contracted[NDIM][NDIM];

  int gcon_ready = 0;
  if (metric == METRIC_FMKS && !use_eKS_internal) {
    fmks_metric_and_derivs(X, gcov, gcon, dg);
    gcon_ready = 1;
  } else {
    gcov_func(X, gcov);
    for (int k = 1; k <= 2; ++k) {
      for (int l = 0; l < NDIM; ++l) {
        Xh[l] = X[l];
        Xl[l] = X[l];
      }
      Xh[k] += 1.e-7;
      Xl[k] -= 1.e-7;
      gcov_func(Xh, gh);
      gcov_func(Xl, gl);
      for (int i = 0; i < NDIM; ++i)
        for (int j = i; j < NDIM; ++j) {
          const double d = (gh[i][j] - gl[i][j]) / 2.e-7;
          dg[k][i][j] = d;
          dg[k][j][i] = d;
        }
    }
  }
  if (!gcon_ready && gcon_func(gcov, gcon)) {
    fprintf(stderr, "Encountered singluar gcov when getting connection!\n");
    print_vector("X", X);
    print_matrix("gcov", gcov);
  }

  /* B[m][j] = (d_j g_ml + d_l g_mj - d_m g_jl) K^l / 2. */
  for (int m = 0; m < NDIM; ++m)
    for (int j = 0; j < NDIM; ++j) {
      double term1 = 0., term2 = 0., term3 = 0.;
      if (j == 1 || j == 2)
        for (int l = 0; l < NDIM; ++l)
          term1 += dg[j][m][l] * Kcon[l];
      for (int l = 1; l <= 2; ++l)
        term2 += dg[l][m][j] * Kcon[l];
      if (m == 1 || m == 2)
        for (int l = 0; l < NDIM; ++l)
          term3 += dg[m][j][l] * Kcon[l];
      contracted[m][j] = 0.5 * (term1 + term2 - term3);
    }

  for (int i = 0; i < NDIM; ++i)
    for (int j = 0; j < NDIM; ++j) {
      double value = 0.;
      for (int m = 0; m < NDIM; ++m)
        value += gcon[i][m] * contracted[m][j];
      out[i][j] = value;
    }
}

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

/*
 * Normalize input vector so that |v . v| = target
 * Overwrites input
 */
void normalize_to(double vcon[NDIM], double Gcov[NDIM][NDIM], double target)
{
  int k, l;
  double norm;

  norm = 0.;
  for (k = 0; k < 4; k++)
    for (l = 0; l < 4; l++)
      norm += vcon[k] * vcon[l] * Gcov[k][l];

  norm = sqrt(fabs(norm));
  for (k = 0; k < 4; k++)
    vcon[k] *= target/norm;

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
int invert_matrix(double Am[NDIM][NDIM], double Aminv[NDIM][NDIM])
{

  int i, j;
  int n = NDIM;
  int permute[NDIM];
  double dxm[NDIM], Amtmp[NDIM][NDIM];

  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++) {
      Amtmp[i][j] = Am[i][j];
    }

  // Get the LU matrix:
  if (LU_decompose(Amtmp, permute) != 0) {
    //fprintf(stderr, "invert_matrix(): singular matrix encountered! \n");
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
int LU_decompose(double A[NDIM][NDIM], int permute[NDIM])
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
      //fprintf(stderr, "LU_decompose(): row-wise singular matrix!\n");
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
void LU_substitution(double A[NDIM][NDIM], double B[NDIM], int permute[NDIM])
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
