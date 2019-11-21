#include "fits.h"

#include <stdlib.h>
#include <stdio.h>

/*Wrappers for the fitting formulae*/

/*j_nu_fit: wrapper for the emissivity fitting formulae.  Takes in the
 *          same parameters as j_nu(), populates the struct of 
 *          parameters, and passes them to the fitting formulae.  The
 *          fitting formulae for each distribution function are located
 *          in their corresponding folders, so for example the
 *          KAPPA_DIST fitting formulae are located in the kappa folder,
 *          in the file kappa_fits.c
 *
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa, 
 *         kappa_width 
 *@returns: the corresponding fitting formula (based on the distribution
 *          function) evaluated for the input parameters.
*/
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
                double kappa_width)
{
/*fill the struct with values*/
  struct parameters params;
  setConstParams(&params);
  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.polarization       = polarization;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;


//  check_for_errors(&params);
  
  if(params.distribution == params.MAXWELL_JUETTNER)
  {
    if     (params.polarization == params.STOKES_I) return maxwell_juettner_I(&params); 
    else if(params.polarization == params.STOKES_Q) return maxwell_juettner_Q(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return maxwell_juettner_V(&params);
  }

  else if(params.distribution == params.POWER_LAW)
  {
    if     (params.polarization == params.STOKES_I) return power_law_I(&params);
    else if(params.polarization == params.STOKES_Q) return power_law_Q(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return power_law_V(&params);
  }

  else if(params.distribution == params.KAPPA_DIST)
  {
    if     (params.polarization == params.STOKES_I) return kappa_I(&params);
    else if(params.polarization == params.STOKES_Q) return kappa_Q(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return kappa_V(&params);
  }

  return 0.;
}

/*alpha_nu_fit: wrapper for the absorptivity fitting formulae.  Takes in the
 *              same parameters as alpha_nu(), populates the struct of 
 *              parameters, and passes them to the fitting formulae.  The
 *              fitting formulae for each distribution function are located
 *              in their corresponding folders, so for example the
 *              KAPPA_DIST fitting formulae are located in the kappa folder,
 *              in the file kappa_fits.c
 *
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa, 
 *         kappa_width 
 *@returns: the corresponding fitting formula (based on the distribution
 *          function) evaluated for the input parameters.
*/

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
                    double kappa_width)
{
/*fill the struct with values*/
  struct parameters params;
  setConstParams(&params);
  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.polarization       = polarization;
  params.mode               = params.ABSORPTIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;


  if(params.distribution == params.MAXWELL_JUETTNER)
  {
    if     (params.polarization == params.STOKES_I) return maxwell_juettner_I_abs(&params);
    else if(params.polarization == params.STOKES_Q) return maxwell_juettner_Q_abs(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return maxwell_juettner_V_abs(&params);
  }

  else if(params.distribution == params.POWER_LAW)
  {
    if     (params.polarization == params.STOKES_I) return power_law_I_abs(&params);
    else if(params.polarization == params.STOKES_Q) return power_law_Q_abs(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return power_law_V_abs(&params);
  }

  else if(params.distribution == params.KAPPA_DIST)
  {
    if     (params.polarization == params.STOKES_I) return kappa_I_abs(&params);
    else if(params.polarization == params.STOKES_Q) return kappa_Q_abs(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return kappa_V_abs(&params);
  }

  return 0.;
}

/*rho_nu_fit: Fits to Faraday rotation/conversion coefficients; right now only has
 *            formulae for Maxwell-Juettner distribution, from Dexter (2016)
 * 
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa, 
 *         kappa_width 
 *@returns: the corresponding Faraday coefficient fitting formula 
 *          (based on the distribution function) evaluated for 
 *          the input parameters.
 */
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
                  double kappa_width)
{
/*fill the struct with values*/
  struct parameters params;
  setConstParams(&params);
  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.polarization       = polarization;
  params.mode               = params.ABSORPTIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

//  check_for_errors(&params);

  if(params.polarization == params.STOKES_I)
  {
    printf("No Faraday rotation of total intensity");
    return 0.;
  }

  if(params.distribution == params.MAXWELL_JUETTNER)
  {
    if     (params.polarization == params.STOKES_Q) return maxwell_juettner_rho_Q(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return maxwell_juettner_rho_V(&params);
  }

/*Faraday rotation coefficients for Power-law and Kappa distributions will be added later*/
//  else if(params.distribution == params.POWER_LAW)
//  {
//    if     (params.polarization == params.STOKES_Q) return power_law_rho_Q(&params);
//    else if(params.polarization == params.STOKES_U) return 0.;
//    else if(params.polarization == params.STOKES_V) return power_law_rho_V(&params);
//  }
//
//  else if(params.distribution == params.KAPPA_DIST)
//  {
//    if     (params.polarization == params.STOKES_Q) return kappa_rho_Q(&params);
//    else if(params.polarization == params.STOKES_U) return 0.;
//    else if(params.polarization == params.STOKES_V) return kappa_rho_V(&params);
//  }

  return 0.;
}

/*check_for_errors: takes in struct of parameters and checks the provided
 *                  user inputs for values that are not physical or are 
 *                  outside of the limits of validity for the fitting 
 *                  formulae.
 *
 *@params: struct of parameters params
 *@returns: prints error messages or quits if disallowed values
 *          (such as magnitude of magnetic field < 0) are entered.
 */
double check_for_errors(struct parameters * params)
{

  double nu_c = get_nu_c(*params);

  /* catch potential errors */
  if(params->nu/nu_c > 3e10 && 
     (params->polarization == params->STOKES_Q 
      || params->polarization == params->STOKES_V))
  {
    printf("\n WARNING: nu high; Stokes Q and V fits may be inaccurate \n");
  }
  if(params->magnetic_field < 0)
  {
    printf("\n ERROR: B magnitude cannot be negative \n");
    exit(0);
  }
  if(params->electron_density < 0)
  {
    printf("\n ERROR: cannot have negative electron number density \n");
    exit(0);
  }
  if(params->kappa < 2.5 || params->kappa > 7.5)
  {
    printf("\n WARNING: kappa out of range of fitting formula \n");
  }
  if(params->kappa_width < 3 || params->kappa_width > 200)
  {
    printf("\n WARNING: w out of range; fitting formula may be inaccurate\n");
  }
  if(params->gamma_min < 1)
  {
    printf("\n ERROR: gamma_min < 1\n");
    exit(0);
  }
  if(params->observer_angle < 5.*(params->pi)/180. 
     || params->observer_angle == 90.*(params->pi)/180.)
  {
    printf("\n WARNING: theta out of range; fitting formula may be inaccurate \n");
  }

  return 0.;
}
