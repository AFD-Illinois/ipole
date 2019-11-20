#include "distribution_function_common_routines.h"

/*normalize_f: normalizes the distribution function using GSL's 
 *             QAGIU integrator.   
 *
 *@params: struct of parameters to pass to distribution function
 *@returns: 1 over the normalization constant for the chosen distribution
 */
double normalize_f(double (*distribution)(double, void *),
                   struct parameters * params
                  )
{
  
  /*set GSL QAGIU integrator parameters */
  double lower_bound = 1.;
  double absolute_error = 0.;
  double relative_error = 1e-8;
  int limit  = 1000;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
  double result, error;
  gsl_function F;

  F.function = distribution;

  F.params = params;
  
  gsl_integration_qagiu(&F, lower_bound, absolute_error, 
                        relative_error, limit, w, &result, &error
                       );


  gsl_integration_workspace_free(w);
  
  return result;
}


/*analytic_differential_of_f: The absorptivity integrand ([1] eq. 12) 
 *                   has a term dependent on a differential operator 
 *                   ([1] eq. 13) applied to the distribution function. 
 *                   This function calls each individual differential, 
 *                   evaluated analytically, for the 3 distributions 
 *                   studied. This function is provided for speed and as
 *                   a check of the approximate numerical differential
 *                   used below. numerical_differential_of_f(), below,
 *                   works with any gyrotropic distribution function.
 * 
 *@params: Lorentz factor gamma, struct of parameters params
 *@returns: differential of the distribution function term in the gamma
 *          integrand.  
 */
double analytic_differential_of_f(double gamma, struct parameters * params) 
{
  /*described in Section 2 of [1] */

  double Df = 0.;

  /*all of the distribution functions used are independent of gyrophase
    phi, so integrate out dphi to get 2*pi */
  double gyrophase_indep =  2. * params->pi;

  /*need to convert d^3p to dgamma dcos(xi) by multiplying by
    a factor of m^3 c^3 */
  double d3p_to_dgamma   =  pow(params->mass_electron, 3.)
                          * pow(params->speed_light, 3.);

  /*prefactor from [1] eq. 13, multiplied by 2*pi from integrating
    out phi and m^3 c^3 from changing from d3p to dgamma dcos(xi) */
  double prefactor = (2. * params->pi * params->nu
                      / (params->mass_electron
                         *params->speed_light*params->speed_light))
                     * gyrophase_indep
                     * d3p_to_dgamma;

  Df = params->analytic_differential(gamma, params);

  return prefactor * Df;
}

double numerical_differential_of_f(double gamma, struct parameters * params)
{
  /*this is "d^3p Df" from [1] eq. 12 and 13.*/ 

  double Df = 0.;
  double epsilon = 3e-4;

  /*all of the distribution functions used are independent of gyrophase
    phi, so integrate out dphi to get 2*pi */
  double gyrophase_indep =  2. * params->pi;

  /*need to convert d^3p to dgamma dcos(xi) by multiplying by
    a factor of m^3 c^3 */
  double d3p_to_dgamma   =  pow(params->mass_electron, 3.)
                          * pow(params->speed_light, 3.);

  /*prefactor from [1] eq. 13, multiplied by 2*pi from integrating
    out phi and m^3 c^3 from changing from d3p to dgamma dcos(xi) */ 
  double prefactor = (2. * params->pi * params->nu
                      / (params->mass_electron
                         *params->speed_light*params->speed_light))
                     * gyrophase_indep
                     * d3p_to_dgamma;

  /*The if statements below are necessary because for some values of
    gamma, the quantity distribution_function(gamma+epsilon) or
    distribution_function(gamma-epsilon) is complex, and returns
    NaN.  The if statements use a one-sided approximation to avoid
    these regions. */
  if(isnan(params->distribution_function(gamma+epsilon, params)) != 0)
  {
    Df =  (params->distribution_function(gamma, params)
           - params->distribution_function(gamma-epsilon, params))
          / (epsilon);
  }
  else if(isnan(params->distribution_function(gamma-epsilon, params)) != 0)
  {
    Df =  (params->distribution_function(gamma+epsilon, params)
           - params->distribution_function(gamma, params))
          / (epsilon);
  }
  else
  {
    Df =  (params->distribution_function(gamma+epsilon, params)
           - params->distribution_function(gamma-epsilon, params))
          / (2. * epsilon);
  }

  return prefactor * Df;
}

