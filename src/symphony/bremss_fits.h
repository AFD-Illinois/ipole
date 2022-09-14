#ifndef SYMPHONY_BREMSS_FITS_H_
#define SYMPHONY_BREMSS_FITS_H_
#include "params.h"
#include <gsl/gsl_sf_bessel.h>

/* Fits */

// Initialize necessary cached values
void init_bremss_spline();

/* Emissivities */
double bremss_I(struct parameters * params, int bremss_type);

/* Absorptivities */
double bremss_I_abs(struct parameters * params, int bremss_type);

#endif /* SYMPHONY_MAXWELL_JUETTNER_H_ */
