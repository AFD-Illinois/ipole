#ifndef SYMPHONY_MAXWELL_JUETTNER_H_
#define SYMPHONY_MAXWELL_JUETTNER_H_
#include "params.h"
//#include "distribution_function_common_routines.h"
#include <gsl/gsl_sf_bessel.h>

double maxwell_juettner_f(double gamma, struct parameters * params);
double differential_of_maxwell_juettner(double gamma, struct parameters * params);
double maxwell_juettner_n_peak(struct parameters * params);

/* Fits */

/* Emissivities */
double maxwell_juettner_I(struct parameters * params);
double maxwell_juettner_Q(struct parameters * params);
double maxwell_juettner_V(struct parameters * params);

double planck_func(struct parameters * params);

/* Absorptivities */
double maxwell_juettner_I_abs(struct parameters * params);
double maxwell_juettner_Q_abs(struct parameters * params);
double maxwell_juettner_V_abs(struct parameters * params);

/* Faraday rotation/conversion coefficients */
double maxwell_juettner_rho_Q(struct parameters * params);
double maxwell_juettner_rho_V(struct parameters * params);

#endif /* SYMPHONY_MAXWELL_JUETTNER_H_ */
