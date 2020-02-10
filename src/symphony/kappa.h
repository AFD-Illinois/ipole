#ifndef SYMPHONY_KAPPA_H_
#define SYMPHONY_KAPPA_H_
#include "params.h"
//#include "distribution_function_common_routines.h"
#include "gsl/gsl_sf_hyperg.h"

double kappa_to_be_normalized(double gamma, void * paramsInput);
double kappa_f(double gamma, struct parameters * params);
double differential_of_kappa(double gamma, struct parameters * params);

double kappa_I(struct parameters * params);
double kappa_Q(struct parameters * params);
double kappa_V(struct parameters * params);

double kappa_I_abs(struct parameters * params);
double kappa_Q_abs(struct parameters * params);
double kappa_V_abs(struct parameters * params);

#endif /* SYMPHONY_KAPPA_H_ */
