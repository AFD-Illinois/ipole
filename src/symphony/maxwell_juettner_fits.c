#include "maxwell_juettner.h"

/*maxwell_juettner_I: fitting formula for the emissivity (polarized in Stokes I)
 *                    produced by a Maxwell-Juettner (relativistic thermal) 
 *                    distribution of electrons. (Eq. 29, 31 of [1])
 *
 *@params: struct of parameters params
 *@returns: fit to the emissivity, polarized in Stokes I, for the given 
 *          parameters for a Maxwell-Juettner distribution.
 */
double maxwell_juettner_I(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_s = (2./9.)*nu_c*sin(params->observer_angle)*params->theta_e
                *params->theta_e;

  double X = params->nu/nu_s;

  double prefactor = (params->electron_density 
                      * pow(params->electron_charge, 2.) 
                      * nu_c)/params->speed_light;

  double term1 = sqrt(2.)*params->pi/27. * sin(params->observer_angle);

  double term2 = pow(pow(X, 0.5)+pow(2., 11./12.)*pow(X, 1./6.), 2.);
  
  double term3 = exp(-pow(X, 1./3.));
  
  double ans = prefactor * term1 * term2 * term3;

  return ans;
}

/*maxwell_juettner_Q: fitting formula for the emissivity (polarized in Stokes Q)
 *                    produced by a Maxwell-Juettner (relativistic thermal) 
 *                    distribution of electrons. (Eq. 29, 31 of [1])
 *
 *@params: struct of parameters params
 *@returns: fit to the emissivity, polarized in Stokes Q, for the given 
 *          parameters for a Maxwell-Juettner distribution.
 */
double maxwell_juettner_Q(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_s = (2./9.)*nu_c*sin(params->observer_angle)
                *params->theta_e*params->theta_e;

  double X = params->nu/nu_s;

  double prefactor = (params->electron_density 
                      * pow(params->electron_charge, 2.) 
                      * nu_c)/params->speed_light;

  double term1 = sqrt(2.)*params->pi/27. * sin(params->observer_angle);

  double term2 = (7.*pow(params->theta_e, 24./25.)+35.)
		/(10.*pow(params->theta_e, 24./25.)+75.);

  double term3 = pow(pow(X, 0.5)+term2*pow(2., 11./12.)*pow(X, 1./6.), 2.);

  double ans = prefactor*term1*term3*exp(-pow(X, 1./3.));

  return -ans;
}

/*maxwell_juettner_V: fitting formula for the emissivity (polarized in Stokes V)
 *                    produced by a Maxwell-Juettner (relativistic thermal) 
 *                    distribution of electrons. (Eq. 29, 31 of [1])
 *
 *@params: struct of parameters params
 *@returns: fit to the emissivity, polarized in Stokes V, for the given 
 *          parameters for a Maxwell-Juettner distribution.
 */
double maxwell_juettner_V(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_s = (2./9.)*nu_c*sin(params->observer_angle)*params->theta_e
                *params->theta_e;

  double X = params->nu/nu_s;

  double prefactor = (params->electron_density 
                      * pow(params->electron_charge, 2.) 
                      * nu_c)/params->speed_light;

  double term1 = (37.-87.*sin(params->observer_angle-28./25.))
                /(100.*(params->theta_e+1.));

  double term2 = pow(1.+(pow(params->theta_e, 3./5.)/25.+7./10.)
		*pow(X, 9./25.), 5./3.);

  double ans = prefactor*term1*term2*exp(-pow(X, 1./3.));

  /*NOTE: Sign corrected; the sign in Leung et al. (2011)
    and Pandya et al. (2016) for Stokes V transfer coefficients
    does not follow the convention the papers describe (IEEE/IAU);
    the sign has been corrected here.*/
  return ans;
}

/*planck_func: The Planck function (used in eq. 25 of [1]) can be used to
 *             obtain alpha_nu() fitting formulae from the j_nu() fitting 
 *             formulae for the Maxwell-Juettner (relativistic thermal)
 *             distribution.
 *
 *@params: struct of parameters params
 *@returns: Planck function evaluated for the supplied parameters
 */
double planck_func(struct parameters * params)
{
  double term1 = (2.*params->plancks_constant*pow(params->nu, 3.))
                /pow(params->speed_light, 2.);

  double term2 = (exp(params->plancks_constant*params->nu
                  /(params->theta_e*params->mass_electron
                    *pow(params->speed_light, 2.)))-1.);

  double ans = term1 / term2;

  return ans;
}

/*maxwell_juettner_I_abs: Fitting formula for the absorptivity, polarized in
 *                        Stokes I, for a Maxwell-Juettner electron momentum 
 *                        distribution.  Uses eq. 30, 31, 32 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to the absorptivity, in Stokes I, for a Maxwell-
 *          Juettner distribution of electrons.
 */
double maxwell_juettner_I_abs(struct parameters * params)
{
  double ans = maxwell_juettner_I(params)/planck_func(params);
  return ans;
}

/*maxwell_juettner_Q_abs: Fitting formula for the absorptivity, polarized in
 *                        Stokes Q, for a Maxwell-Juettner electron momentum 
 *                        distribution.  Uses eq. 30, 31, 32 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to the absorptivity, in Stokes Q, for a Maxwell-
 *          Juettner distribution of electrons.
 */
double maxwell_juettner_Q_abs(struct parameters * params)
{
  double ans = maxwell_juettner_Q(params)/planck_func(params);
  return ans;
}

/*maxwell_juettner_V_abs: Fitting formula for the absorptivity, polarized in
 *                        Stokes V, for a Maxwell-Juettner electron momentum 
 *                        distribution.  Uses eq. 30, 31, 32 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to the absorptivity, in Stokes V, for a Maxwell-
 *          Juettner distribution of electrons.
 */
double maxwell_juettner_V_abs(struct parameters * params)
{
  double ans = maxwell_juettner_V(params)/planck_func(params);
  return ans;
}

/*maxwell_juettner_rho_Q: Fitting formula for Faraday conversion coefficient
 *                        rho_Q from Dexter (2016) for a Maxwell-Juettner
 *                        distribution.  This formula comes from his
 *                        equations B4, B6, B8, and B13.
 *
 *@params: struct of parameters params
 *@returns: fitting formula to the Faraday conversion coefficient
 *          for a Maxwell-Juettner distribution of electrons. 
 */
double maxwell_juettner_rho_Q(struct parameters * params)
{
  double omega0 = params->electron_charge*params->magnetic_field
                  / (params->mass_electron*params->speed_light);

  double wp2 = 4. * params->pi * params->electron_density 
	       * pow(params->electron_charge, 2.) / params->mass_electron;

  /* argument for function f(X) (called jffunc) below */
  double x = params->theta_e * sqrt(sqrt(2.) * sin(params->observer_angle)
                 * (1.e3*omega0 / (2. * params->pi * params->nu)));

  /* this is the definition of the modified factor f(X) from Dexter (2016) */
  double extraterm = (.011*exp(-x/47.2) - pow(2., (-1./3.)) / pow(3., (23./6.))
                      * params->pi * 1.e4 * pow((x + 1.e-16), (-8./3.)))
                     * (0.5 + 0.5 * tanh((log(x) - log(120.)) / 0.1));

  double jffunc = 2.011 * exp(-pow(x, 1.035)/4.7) - cos(x/2.) 
                  * exp(-pow(x, 1.2)/2.73) - .011 * exp(-x / 47.2) + extraterm;

  double k1 = gsl_sf_bessel_Kn(1, 1./params->theta_e);
  double k2 = gsl_sf_bessel_Kn(2, 1./params->theta_e);
  double k_ratio = (k2 > 0) ? k1/k2 : 1;

  double eps11m22 = jffunc * wp2 * pow(omega0, 2.) 
                    / pow(2.*params->pi * params->nu, 4.)
                    * (k_ratio + 6. * params->theta_e)
                    * pow(sin(params->observer_angle), 2.);

  double rhoq = 2. * params->pi * params->nu /(2. * params->speed_light) 
                * eps11m22;

  return rhoq;
}

/*maxwell_juettner_rho_V: Fitting formula for Faraday rotation coefficient
 *                        rho_V from Dexter (2016) for a Maxwell-Juettner
 *                        distribution.  This formula comes from his
 *                        equations B7, B8, B14, and B15.
 *
 *@params: struct of parameters params
 *@returns: fitting formula to the Faraday rotation coefficient
 *          for a Maxwell-Juettner distribution of electrons. 
 */
double maxwell_juettner_rho_V(struct parameters * params)
{
  double omega0 = params->electron_charge*params->magnetic_field
                  / (params->mass_electron*params->speed_light);

  double wp2 = 4. * params->pi * params->electron_density 
                  * pow(params->electron_charge, 2.) / params->mass_electron;

  /* argument for function g(X) (called shgmfunc) below */
  double x = params->theta_e * sqrt(sqrt(2.) * sin(params->observer_angle)
                 * (1.e3*omega0 / (2. * params->pi * params->nu)));

  /* Approximate the Bessel functions if allowed */
  double k2=0, k0=0;
  if (params->approximate && params->theta_e > 5) {
    k0 = -log(1 / (2. * params->theta_e)) - 0.5772;
    k2 = 2. * params->theta_e * params->theta_e;
  } else {
    k0 = gsl_sf_bessel_Kn(0, 1./params->theta_e);
    k2 = gsl_sf_bessel_Kn(2, 1./params->theta_e);
  }

  /* There are several fits of rho_V phrased as functions of x */
  // TODO add straight Bessel-approx -> constant extrapolation?
  double fit_factor = 0;
  if (params->dexter_fit && k2 > 0) { // TODO Further limit the usage here to match grtrans
    // Jason Dexter (2016) fits using the modified difference factor g(X)
    double shgmfunc = 0.43793091 * log(1. + 0.00185777 * pow(x, 1.50316886));
    fit_factor = (k0 - shgmfunc) / k2;
  } else {
    // Shcherbakov fits.  Good to the smallest Thetae at high freq but questionable for low frequencies
    double shgmfunc = 1 - 0.11*log(1 + 0.035*x);
    double k_ratio = (k2 > 0) ? k0/k2 : 1;
    fit_factor = k_ratio * shgmfunc;
  }

  double eps12 = wp2 * omega0 / pow((2. * params->pi * params->nu), 3.)
                     * fit_factor * cos(params->observer_angle);

  return 2. * params->pi * params->nu / params->speed_light * eps12;
}
