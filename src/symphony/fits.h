#ifndef SYMPHONY_FITS_H_
#define SYMPHONY_FITS_H_

#include "params.h"

/* Maxwell-Juettner emissivity fits */
double maxwell_juettner_I(struct parameters * params);
double maxwell_juettner_Q(struct parameters * params);
double maxwell_juettner_V(struct parameters * params);

/* Maxwell-Juettner absorptivity fits */
double maxwell_juettner_I_abs(struct parameters * params);
double maxwell_juettner_Q_abs(struct parameters * params);
double maxwell_juettner_V_abs(struct parameters * params);

/* Maxwell-Juettner Faraday rotation fits */
double maxwell_juettner_rho_Q(struct parameters * params);
double maxwell_juettner_rho_V(struct parameters * params);

/* Power-law emissivity fits */
double power_law_I(struct parameters * params);
double power_law_Q(struct parameters * params);
double power_law_V(struct parameters * params);

/* Power-law absorptivity fits */
double power_law_I_abs(struct parameters * params);
double power_law_Q_abs(struct parameters * params);
double power_law_V_abs(struct parameters * params);

/* Kappa emissivity fits */
double kappa_I(struct parameters * params);
double kappa_Q(struct parameters * params);
double kappa_V(struct parameters * params);

/* Kappa absorptivity fits */
double kappa_I_abs(struct parameters * params);
double kappa_Q_abs(struct parameters * params);
double kappa_V_abs(struct parameters * params);

/* Kappa Faraday Rotation fits */
double kappa35_rho_Q(struct parameters * params);
double kappa4_rho_Q(struct parameters * params);
double kappa45_rho_Q(struct parameters * params);
double kappa5_rho_Q(struct parameters * params);
double kappa35_rho_V(struct parameters * params);
double kappa4_rho_V(struct parameters * params);
double kappa45_rho_V(struct parameters * params);
double kappa5_rho_V(struct parameters * params);

void check_for_errors(struct parameters * params);
double j_nu_fit(struct parameters * params, int polarization);
double alpha_nu_fit(struct parameters * params, int polarization);
double rho_nu_fit(struct parameters * params, int polarization);
#endif /* SYMPHONY_FITS_H_ */

