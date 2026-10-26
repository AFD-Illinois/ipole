/*
 * model_radiation.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef MODEL_RADIATION_H
#define MODEL_RADIATION_H

#include "decs.h"
#include "par.h"

void try_set_radiation_parameter(const char *word, const char *value);

/* transfer coefficients in tetrad frame */
void jar_calc(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV, Params *params);

/* As above, also return the fluid state already sampled for the coefficients. */
void jar_calc_with_fluid(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV,
    double Ucon[NDIM], double Ucov[NDIM],
    double Bcon[NDIM], double Bcov[NDIM], int *has_fluid, Params *params);

void get_jkinv(double X[NDIM], double Kcon[NDIM], double *jI, double *aI, Params *params);

#endif /* MODEL_RADIATION_H */
