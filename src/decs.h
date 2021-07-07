#ifndef DECS_H
#define DECS_H

#include "constants.h"
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

// Strings and string tools
#define VERSION_STRING "ipole-beta-1.5"
#define xstr(s) str(s)
#define str(s) #s
#define STRLEN (2048)

// Dimensions
#define NDIM	4
#define NIMG (5) // Stokes parameters plus faraday depth

// "functions"
#define sign(x) (((x) < 0) ? -1 : ((x) > 0))
#define delta(x, y) ((x) == (y))

// MAX and MIN were here but are both unsafe & slow. Use fmax/fmin

#define MULOOP for(int mu=0;mu<NDIM;mu++)
#define MUNULOOP for(int mu=0;mu<NDIM;mu++) for(int nu=0;nu<NDIM;nu++)

struct of_traj {
  int nturns;
  double dl;
  double X[NDIM];
  double Kcon[NDIM];
  double Xhalf[NDIM];
  double Kconhalf[NDIM];
};

// used for slow light
struct of_image {
  int nstep;
  double intensity;
  double tau;
  double tauF;
  double complex N_coord[NDIM][NDIM];
};

#endif // DECS_H
