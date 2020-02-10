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

double check_for_errors(struct parameters * params);
double j_nu_fit(double nu,
                double magnetic_field,
                double electron_density,
                double observer_angle,
                int distribution,
                int polarization,
                double theta_e,
                double power_law_p,
                double gamma_min,
                double gamma_max,
                double gamma_cutoff,
                double kappa,
                double kappa_width
                );
double alpha_nu_fit(double nu,
                    double magnetic_field,
                    double electron_density,
                    double observer_angle,
                    int distribution,
                    int polarization,
                    double theta_e,
                    double power_law_p,
                    double gamma_min,
                    double gamma_max,
                    double gamma_cutoff,
                    double kappa,
                    double kappa_width
                    );

double rho_nu_fit(double nu,
                  double magnetic_field,
                  double electron_density,
                  double observer_angle,
                  int distribution,
                  int polarization,
                  double theta_e,
                  double power_law_p,
                  double gamma_min,
                  double gamma_max,
                  double gamma_cutoff,
                  double kappa,
                  double kappa_width
		  );
#endif /* SYMPHONY_FITS_H_ */

