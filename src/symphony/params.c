#include "params.h"
#include <stddef.h> /* for NULL */

/*setConstParams: sets values of constant parameters used throughout
 *                the calculation, such as pi or the charge of an 
 *                electron (in CGS units).  Also sets values of
 *                the keys for the distribution, emissivity/
 *                absorptivity, and Stokes parameter.
 *
 *@params: struct of parameters params
 *@returns: sets values of constant parameters.
 */
void setConstParams(struct parameters *params)
{
  params->pi               = 3.1415926535897932384;
  params->mass_electron    = 9.1093826e-28;
  params->plancks_constant = 6.6260693e-27;
  params->speed_light      = 2.99792458e10;
  params->electron_charge  = 4.80320680e-10;
  params->n_max            = 30;
  params->C                = 10;
  /* Keys for the distributions */
  params->MAXWELL_JUETTNER = 0;
  params->POWER_LAW        = 1;
  params->KAPPA_DIST       = 2;
  /* Keys for the polarization parameter */
  params->STOKES_I         = 15;
  params->STOKES_Q         = 16;
  params->STOKES_U         = 17;
  params->STOKES_V         = 18;
  /* Keys for the mode: absorptivity or emissivity */
  params->ABSORPTIVITY     = 10;
  params->EMISSIVITY       = 11;
  /* Default: find n-space peak adaptively */
  params->use_n_peak       = 0;
  params->error_message    = NULL;
  /* Key for approach to compute coefficients */
  params->SYMPHONY_METHOD  = 20;
  params->SUSCEPT_METHOD   = 21;
  /*Default options for fits*/
  params->approximate = 1;
  params->dexter_fit = 0;

  params->epsilon0  = 1./(4. * params->pi); //permittivity of free space, CGS units
  params->epsilon   = -1.;            //sign of electron charge
  
}

/*get_nu_c: takes in values of electron_charge, magnetic_field, mass_electron,
 *          and speed_light, and returns the cyclotron frequency nu_c.  
 * 
 *@params:  struct parameters params, contains parameters mentioned above
 *@returns: cyclotron frequency, nu_c, for the provided parameters
 */
double get_nu_c(struct parameters params)
{
  return  (params.electron_charge * params.magnetic_field)
        / (2. * params.pi * params.mass_electron * params.speed_light);
}

double get_omega_p(struct parameters params)
{
  return sqrt(params.electron_density * params.electron_charge
              * params.electron_charge 
              / (params.mass_electron * params.epsilon0));
}
