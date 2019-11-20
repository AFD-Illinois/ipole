#include "power_law.h"

/*power_law_to_be_normalized: Unnormalized power-law distribution function
 *                            (eq. 14 and 17 of [1]).  Analytical normalization
 *                            for this distribution function is used here, but
 *                            it is often useful to append an exponential
 *                            cutoff to the distribution at some Lorentz factor
 *                            gamma.  Because of this exponential cutoff, the 
 *                            power-law must still be normalized (done
 *                            numerically by normalize_f()).  One can 
 *                            "turn off" the exponential cutoff by setting the
 *                            parameter gamma_cutoff to be very large, in which
 *                            case the cutoff term is approximately equal to 1.
 *
 *@params: Lorentz factor gamma, void pointer to struct of parameters
 *@returns: unnormalized power-law distribution function with exponential
 *          cutoff.
 */
double power_law_to_be_normalized(double gamma, void * paramsInput) 
{
  struct parameters * params = (struct parameters*) paramsInput;

  double norm_term = 4. * params->pi;

/*TODO: the discontinuity in PL is making NANs pop up in Df; fix this */
//  if (gamma < params->gamma_min || gamma > params->gamma_max)
//      return 0;

  double prefactor = (params->power_law_p - 1.) / 
                     (pow(params->gamma_min, 1. - params->power_law_p) 
                      - pow(params->gamma_max, 1. 
                            - params->power_law_p));

  double body = pow(gamma, -params->power_law_p) 
                * exp(- gamma / params->gamma_cutoff);

  double ans = norm_term * prefactor * body;

  return ans;
}

/*power_law_f: normalized power-law distribution function with exponential
 *             cutoff.  This is normalized by calling normalize_f(), in the
 *             file integrate.c, which uses GSL's QAGIU integrator.
 *
 *@params: Lorentz factor gamma, struct of parameters params
 *@returns: normalized power-law distribution function with
 *          exponential cutoff.
 */
double power_law_f(double gamma, struct parameters * params) 
{

  /*The distribution function only needs to be normalized once; the
    static variable declaration and if statement below ensure that
    it is not being normalized with each function call to power_law_f;
    it does need to be normalized again if any of the distribution
    function parameters change, though. */
  static double norm = 0.;
  static double previous_power_law_p  = 0.;
  static double previous_gamma_min    = 0.;
  static double previous_gamma_max    = 0.;
  static double previous_gamma_cutoff = 0.;
  if(norm == 0. || previous_power_law_p  != params->power_law_p
                || previous_gamma_min    != params->gamma_min
                || previous_gamma_max    != params->gamma_max
                || previous_gamma_cutoff != params->gamma_cutoff)
  {
    norm = 1./normalize_f(&power_law_to_be_normalized, params);
    previous_power_law_p  = params->power_law_p;
    previous_gamma_min    = params->gamma_min;
    previous_gamma_max    = params->gamma_max;
    previous_gamma_cutoff = params->gamma_cutoff;
  }

/*TODO: the discontinuity in PL is making NANs pop up in Df; fix this */
//  if (gamma < params->gamma_min || gamma > params->gamma_max)
//      return 0;

  double beta = sqrt(1. - 1./(gamma*gamma));

  double prefactor = params->electron_density * (params->power_law_p - 1.) 
                     / (pow(params->gamma_min, 1. - params->power_law_p) 
                        - pow(params->gamma_max, 1. - params->power_law_p));

  double body = pow(gamma, -params->power_law_p) 
                * exp(- gamma / params->gamma_cutoff);

  double ans = norm * prefactor * body 
               * 1./(pow(params->mass_electron, 3.) 
               * pow(params->speed_light, 3.) 
               * gamma*gamma * beta);

  return ans;

}

/*differential_of_power_law: The integrand for the absorptivity calculation 
 *                           ([1] eq. 12) depends on a differential of the 
 *                           distribution function ([1] eq. 13).  For the 
 *                           power-law distribution, this is evaluated
 *                           analytically for speed and accuracy. 
 *
 *@params: Lorentz factor gamma, struct of parameters params
 *@returns: the differential of the power-law distribution function
 *          ([1] eq. 12, 13, and 17) for the alpha_nu() calculation.
 *
 */

double differential_of_power_law(double gamma, struct parameters * params)
{
  /* There's a discontinuity in `f` at gamma_{min,max} that we're not going
   * to handle correctly, so it seems best to insist that any integration
   * not go beyond those bounds. */
  if (gamma <= params->gamma_min || gamma >= params->gamma_max)
      return NAN;

  double pl_norm = 1./(normalize_f(&power_law_to_be_normalized, params));

  double d3p_to_dgamma = 1./( pow(params->mass_electron, 3.)
                             * pow(params->speed_light, 3.));

  double prefactor = (params->electron_density*(params->power_law_p-1.))
                     / ( (  pow(params->gamma_min, 1.-params->power_law_p) 
                          - pow(params->gamma_max, 1.-params->power_law_p)
                         )
                       );

  double term1 = ((-params->power_law_p-1.)*exp(-gamma/params->gamma_cutoff)
             *pow(gamma,-params->power_law_p-2.)/(sqrt(gamma*gamma - 1.)));

  double term2 = (exp(-gamma/params->gamma_cutoff) 
                  * pow(gamma,(-params->power_law_p-1.))
                  /(params->gamma_cutoff * sqrt(gamma*gamma - 1.)));

  double term3 = (exp(-gamma/params->gamma_cutoff) 
                  * pow(gamma,-params->power_law_p))
                 /pow((gamma*gamma - 1.), (3./2.));

  double Df = pl_norm * d3p_to_dgamma * prefactor * (term1 - term2 - term3);

  return Df;
}
