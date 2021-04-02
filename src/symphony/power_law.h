#ifndef SYMPHONY_POWER_LAW_H_
#define SYMPHONY_POWER_LAW_H_
#include "params.h"
//#include "distribution_function_common_routines.h"

double power_law_to_be_normalized(double gamma, void * paramsInput);
double power_law_f(double gamma, struct parameters * params); 
double differential_of_power_law(double gamma, struct parameters * params);

double power_law_I(struct parameters * params);
double power_law_Q(struct parameters * params);
double power_law_V(struct parameters * params);

double power_law_I_abs(struct parameters * params);
double power_law_Q_abs(struct parameters * params);
double power_law_V_abs(struct parameters * params);

#endif /* SYMPHONY_POWER_LAW_H_ */
