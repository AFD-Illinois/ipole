
#ifndef MODEL_PARAMS_H
#define MODEL_PARAMS_H

#include "decs.h"

// This could be fun with a moving object sometime
#define SLOW_LIGHT (0)
// There is really only one problem where we should record Stokes parameters per step
// it makes no sense in a relativistic context
#define INTEGRATOR_TEST (1)

// Model parameters (TODO eliminate if we can these are unused)
extern double rmax_geo;
extern double model_dl;

void record_stokes_parameters(double SI, double SQ, double SU, double SV, double lam);

#endif /* MODEL_PARAMS_H */
