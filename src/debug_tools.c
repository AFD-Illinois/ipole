

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
