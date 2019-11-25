/*
 * model_radiation.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef MODEL_RADIATION_H
#define MODEL_RADIATION_H

#include "decs.h"

/* transfer coefficients in tetrad frame */
void jar_calc(double X[NDIM], double Kcon[NDIM], double *jI, double *jQ,
              double *jU, double *jV, double *aI, double *aQ, double *aU,
              double *aV, double *rQ, double *rU, double *rV);

double jnu_synch(double nu, double Ne, double Thetae, double B, double theta);
void get_jkinv(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv);

#endif /* MODEL_RADIATION_H */
