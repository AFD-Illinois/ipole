
/*
 * Debugging-specific utilities gathered from around ipole
 * Printing and sanity checks for tetrads
 */

#include "geometry.h"
#include "decs.h"

// The world needed these
// Maybe not in this form
void print_matrix(char *name, double g[NDIM][NDIM])
{
  // Print a name and a matrix
  fprintf(stderr, "%s:\n", name);
  fprintf(stderr,
          "%g\t%g\t%g\t%g\n%g\t%g\t%g\t%g\n%g\t%g\t%g\t%g\n%g\t%g\t%g\t%g\n",
          g[0][0], g[0][1], g[0][2], g[0][3], g[1][0], g[1][1], g[1][2],
          g[1][3], g[2][0], g[2][1], g[2][2], g[2][3], g[3][0], g[3][1],
          g[3][2], g[3][3]);
  // Additionally kill things if/when we hit NaNs
  MUNULOOP
    if (isnan(g[mu][nu]))
      exit(-1);
}
void print_matrix_c(char *name, double complex g[NDIM][NDIM])
{
  fprintf(stderr, "%s:\n", name);
  fprintf(stderr,
          "%g+i%g\t%g+i%g\t%g+i%g\t%g+i%g\n%g+i%g\t%g+i%g\t%g+i%g\t%g+i%g\n%g+i%g\t%g+i%g\t%g+i%g\t%g+i%g\n%g+i%g\t%g+i%g\t%g+i%g\t%g+i%g\n",
          creal(g[0][0]), cimag(g[0][0]), creal(g[0][1]), cimag(g[0][1]), creal(g[0][2]), cimag(g[0][2]), creal(g[0][3]), cimag(g[0][3]),
          creal(g[1][0]), cimag(g[1][0]), creal(g[1][1]), cimag(g[1][1]), creal(g[1][2]), cimag(g[1][2]), creal(g[1][3]), cimag(g[1][3]),
          creal(g[2][0]), cimag(g[2][0]), creal(g[2][1]), cimag(g[2][1]), creal(g[2][2]), cimag(g[2][2]), creal(g[2][3]), cimag(g[2][3]),
          creal(g[3][0]), cimag(g[3][0]), creal(g[3][1]), cimag(g[3][1]), creal(g[3][2]), cimag(g[3][2]), creal(g[3][3]), cimag(g[3][3])
          );
  MUNULOOP
    if ( isnan(creal(g[mu][nu])) || isnan(cimag(g[mu][nu])) )
      exit(-1);
}
void print_vector(char *name, double v[NDIM])
{
  fprintf(stderr, "%s: ", name);
  fprintf(stderr,
          "%g\t%g\t%g\t%g\n",
          v[0], v[1], v[2], v[3]);
  MUNULOOP
    if (isnan(v[nu]))
      exit(-1);
}

void check_ortho(double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM])
{
  MUNULOOP {
    double sum = 0. ;
    for(int m=0;m<NDIM;m++) sum += Econ[mu][m]*Ecov[nu][m];
    // TODO gauge normality as harshly as orhogonality?
    // TODO flag under DEBUG
    //if (sum > 1e-10 && mu != nu) fprintf(stderr,"non-orthonormal Econ,Ecov %d %d: %g\n", mu, nu, sum);
  }
}

void check_u(double Ucon[NDIM], double Ucov[NDIM]) {
  double U = 0;
  MULOOP U += Ucon[mu] * Ucov[mu];
  if (isnan(U) || fabs(fabs(U) - 1) > 1e-6) printf("U is wrong: u.u = %f\n", U);
}

void check_N(double complex N[NDIM][NDIM],
    double Kcon[NDIM], double gcov[NDIM][NDIM])
{
  double complex dot;
  double Kcov[NDIM];
  int i, j;

  fprintf(stderr, "enter check_N\n");

  /* compute k . N */
  flip_index(Kcon, gcov, Kcov);
  fprintf(stderr, "(k . N\n");
  /* first one way */
  for (i = 0; i < 4; i++) {
    dot = 0. + I * 0.;
    for (j = 0; j < 4; j++)
    dot += Kcov[j] * N[j][i];
    fprintf(stderr, "%d %g + i %g\n", i, creal(dot), cimag(dot));
  }
  /* then the other */
  for (i = 0; i < 4; i++) {
    dot = 0. + I * 0.;
    for (j = 0; j < 4; j++)
    dot += Kcov[j] * N[i][j];
    fprintf(stderr, "%d %g + i %g\n", i, creal(dot), cimag(dot));
  }
  fprintf(stderr, "k . N)\n");

  /* check for hermiticity */
  fprintf(stderr, "(herm:\n");
  for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++)
  fprintf(stderr, "%d %d %g + i %g\n", i, j,
      creal(N[i][j] - conj(N[j][i])),
      cimag(N[i][j] - conj(N[j][i]))
  );
  fprintf(stderr, "herm)\n");

  /* check invariants */
  double complex Nud[NDIM][NDIM];
  void complex_lower(double complex N[NDIM][NDIM],
      double gcov[NDIM][NDIM], int low1, int low2,
      double complex Nl[NDIM][NDIM]);
  complex_lower(N, gcov, 0, 1, Nud);
  for (i = 0; i < 4; i++)
  fprintf(stderr, "N: %d %g + i %g\n", i, creal(N[i][i]),
      cimag(N[i][i]));
  for (i = 0; i < 4; i++)
  fprintf(stderr, "Nud: %d %g + i %g\n", i, creal(Nud[i][i]),
      cimag(Nud[i][i]));
  dot = 0. + I * 0.;
  for (i = 0; i < 4; i++)
  dot += Nud[i][i];
  fprintf(stderr, "I: %g + i %g\n", creal(dot), cimag(dot));

  double complex Ndd[NDIM][NDIM];
  complex_lower(N, gcov, 1, 1, Ndd);
  dot = 0. + I * 0.;
  for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++)
  dot +=
  2. * 0.25 * (N[i][j] + N[j][i]) * (Ndd[i][j] + Ndd[j][i]);
  fprintf(stderr, "IQUsq: %g + i %g\n", creal(dot), cimag(dot));

  dot = 0. + I * 0.;
  for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++)
  dot +=
  -2. * 0.25 * (N[i][j] - N[j][i]) * (Ndd[i][j] - Ndd[j][i]);
  fprintf(stderr, "Vsqsq: %g + i %g\n", creal(dot), cimag(dot));

  fprintf(stderr, "leave check_N\n");
}
