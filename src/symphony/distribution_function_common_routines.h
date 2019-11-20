//#ifndef SYMPHONY_DISTRIBUTION_FUNCTION_COMMON_ROUTINES_H_
//#define SYMPHONY_DISTRIBUTION_FUNCTION_COMMON_ROUTINES_H_

#include <gsl/gsl_integration.h>
#include "params.h"
#include "maxwell_juettner.h"
#include "power_law.h"
#include "kappa.h"

double normalize_f(double (*distribution)(double, void *),
                   struct parameters * params
                  );

double numerical_differential_of_f(double gamma, struct parameters * params);
double analytic_differential_of_f(double gamma, struct parameters * params);

//#endif /* SYMPHONY_DISTRIBUTION_FUNCTION_COMMON_ROUTINES_H_ */
