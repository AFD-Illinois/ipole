#ifndef DECS_H
#define DECS_H

//
#include "constants.h"
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

#define NDIM	4

#define xstr(s) str(s)
#define str(s) #s

#define STRLEN (2048)

#define NIMG (5) // Stokes parameters plus faraday depth

#define sign(x) (((x) < 0) ? -1 : ((x) > 0))
#define delta(x, y) ( ((x) == (y)) ? 1 : 0 )

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

#define MULOOP for(int mu=0;mu<NDIM;mu++)
#define MUNULOOP for(int mu=0;mu<NDIM;mu++) for(int nu=0;nu<NDIM;nu++)

// TODO phase out for print_vector from debug_tools?
#define P4VEC(S,X) fprintf(stderr, "%s: %g %g %g %g\n", S, X[0], X[1], X[2], X[3]);

#define MAXNSTEP 50000

// used for slow light
struct of_image {
  int nstep;
  double intensity;
  double tau;
  double tauF;
  double complex N_coord[NDIM][NDIM];
};

#endif // DECS_H
