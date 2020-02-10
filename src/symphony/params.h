#ifndef SYMPHONY_PARAMS_H_
#define SYMPHONY_PARAMS_H_

#include <math.h>

struct parameters
{
  /*parameters of calculation*/
  /*we use Gaussian CGS units*/
  double pi;        
  double mass_electron;
  double plancks_constant;
  double speed_light;
  double electron_charge;
  double n_max;
  int    C;       
  /*Keys for the distributions*/
  int    MAXWELL_JUETTNER;
  int    POWER_LAW;
  int    KAPPA_DIST;
  /*Keys for the polarization modes*/
  int    STOKES_I;
  int    STOKES_Q;
  int    STOKES_U;
  int    STOKES_V;
  /*Keys for the mode: absorptivity or emissivity*/
  int    ABSORPTIVITY;
  int    EMISSIVITY;
  /*Keys for alpha_nu computation method*/
  int SYMPHONY_METHOD;
  int SUSCEPT_METHOD;

  /*USER PARAMS:*/
  double nu;               /* GHz */
  double magnetic_field;   /* Gauss */
  double electron_density; /* 1/cc */
  double observer_angle;   /* rad */  
  int    distribution;     
  int    polarization; 
  int    mode;             /*Emissivity or Absorptivity*/
  double gamma_cutoff;
  /* Options for fits */
  int approximate;         /*Use approximate Bessel functions for speed*/
  int dexter_fit;

  /*Thermal distribution parameters*/
  double theta_e;

  /*power law parameters*/
  double power_law_p;
  double gamma_min;
  double gamma_max;

  /*kappa distribution parameters*/
  double kappa;
  double kappa_width;

  /*Choose if n-space peak is known, or if it must be found adaptively */
  int use_n_peak;
  double (*n_peak)(struct parameters *);

  /*Set distribution_function */
  double (*distribution_function)(double gamma, struct parameters *);

  /*analytic_differential_of_f, which can be used as a test of the 
    numerical differential_of_f */
  double (*analytic_differential)(double gamma, struct parameters *);

  int stokes_v_switch;

  char *error_message; /* if not NULL, records source of error in current calculation */

  /*susceptibility tensor paramsS */
  double gamma;
  double epsilon0;
  double epsilon;
  double omega;
  double omega_c;
  double omega_p;
  double real;
  double (*tau_integrand)(double, void * parameters);
  double (*gamma_integrand)(double, void * parameters);
};

struct parametersGSL
{
  struct parameters params;
  double n;
};

void setConstParams(struct parameters *params);
double get_nu_c(struct parameters params);
double get_omega_p(struct parameters params);

#endif /* SYMPHONY_PARAMS_H_ */
