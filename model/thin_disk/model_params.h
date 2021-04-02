
#ifndef MODEL_PARAMS_H
#define MODEL_PARAMS_H

#include "decs.h"

// No sense using slow light with a static model
#define SLOW_LIGHT (0)

#define THIN_DISK (1)

// Model parameters
extern double rmax_geo;

// New definitions needed for problem defined with boundary condition
// Uses are #ifdef'd off on the THIN_DISK flag above
int thindisk_region(double Xi[NDIM], double Xf[NDIM]);
void get_model_i(double X[NDIM], double K[NDIM], double *SI);
void get_model_stokes(double X[NDIM], double K[NDIM], double *SI, double *SQ, double *SU, double *SV);


#endif /* MODEL_PARAMS_H */
