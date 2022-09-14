#include "maxwell_juettner.h"

#include "constants.h"
#include "radiation.h"

#include <stdio.h>

// Local functions for dispatched fits
double I_I(double x);
double I_Q(double x);
double I_V(double x);
double maxwell_juettner_dexter_I(struct parameters *params);
double maxwell_juettner_dexter_Q(struct parameters *params);
double maxwell_juettner_dexter_V(struct parameters *params);
double maxwell_juettner_leung_I(struct parameters *params);

/********* PANDYA FITS *********/

/*maxwell_juettner_I: fitting formula for the emissivity (polarized in Stokes I)
 *                    produced by a Maxwell-Juettner (relativistic thermal) 
 *                    distribution of electrons. (Eq. 29, 31 of [1])
 *
 *@params: struct of parameters params
 *@returns: fit to the emissivity, polarized in Stokes I, for the given 
 *          parameters for a Maxwell-Juettner distribution.
 */
double maxwell_juettner_I(struct parameters *params)
{
  if (params->dexter_fit == 2) {
    return maxwell_juettner_leung_I(params);
  } else if (params->dexter_fit) {
    return maxwell_juettner_dexter_I(params);
  }
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
double maxwell_juettner_Q(struct parameters *params)
{
  if (params->dexter_fit) {
    return maxwell_juettner_dexter_Q(params);
  }
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
double maxwell_juettner_V(struct parameters *params)
{
  if (params->dexter_fit) {
    return maxwell_juettner_dexter_V(params);
  }
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
double planck_func(struct parameters *params)
{
  return Bnu_inv(params->nu, params->theta_e) * pow(params->nu, 3);
}

/*maxwell_juettner_I_abs: Fitting formula for the absorptivity, polarized in
 *                        Stokes I, for a Maxwell-Juettner electron momentum 
 *                        distribution.  Uses eq. 30, 31, 32 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to the absorptivity, in Stokes I, for a Maxwell-
 *          Juettner distribution of electrons.
 */
double maxwell_juettner_I_abs(struct parameters *params)
{
  double Bnu = planck_func(params);
  double ans = 0.;
  if (Bnu > 0.) {
    ans = maxwell_juettner_I(params)/Bnu;
  }
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
double maxwell_juettner_Q_abs(struct parameters *params)
{
  double Bnu = planck_func(params);
  double ans = 0.;
  if (Bnu > 0.) {
    ans = maxwell_juettner_Q(params)/Bnu;
  }
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
double maxwell_juettner_V_abs(struct parameters *params)
{
  double Bnu = planck_func(params);
  double ans = 0.;
  if (Bnu > 0.) {
    ans = maxwell_juettner_V(params)/Bnu;
  }
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
double maxwell_juettner_rho_Q(struct parameters *params)
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

  double k0 = gsl_sf_bessel_Kn(0, 1./params->theta_e);
  double k2 = gsl_sf_bessel_Kn(2, 1./params->theta_e);

  /* There are several fits of rho_V phrased as functions of x */
  // TODO add straight Bessel-approx -> constant extrapolation?
  double fit_factor = 0;
  if (params->dexter_fit && k2 > 0) { // TODO Further limit the usage here to match grtrans
    // Jason Dexter (2016) fits using the modified difference factor g(X)
    double shgmfunc = 0.43793091 * log(1. + 0.00185777 * pow(x, 1.50316886));
    fit_factor = (k0 - shgmfunc) / k2;  // TODO might be unstable way to phrase
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

/********* DEXTER FITS *********/
double maxwell_juettner_dexter_I(struct parameters *params)
{
  double Ne = params->electron_density;
  double nu = params->nu;
  double Thetae = params->theta_e;
  double B = params->magnetic_field;
  double theta = params->observer_angle;
  // Synchrotron emissivity
  double nus = 3.0 * EE * B * sin(theta) / 4.0 / M_PI / ME / CL * Thetae * Thetae + 1.0;
  double x = nu / nus;

  return Ne * EE * EE * nu / 2. / sqrt(3) / CL / Thetae / Thetae * I_I(x); // [g/s^2/cm = ergs/s/cm^3]
}

double maxwell_juettner_dexter_Q(struct parameters *params)
{
  double Ne = params->electron_density;
  double nu = params->nu;
  double Thetae = params->theta_e;
  double B = params->magnetic_field;
  double theta = params->observer_angle;
  double nus = 3.0 * EE * B * sin(theta) / 4.0 / M_PI / ME / CL * Thetae * Thetae + 1.0;
  double x = nu / nus;
  return -Ne * EE * EE * nu / 2. / sqrt(3) / CL / Thetae / Thetae * I_Q(x); // Translate to Symphony convention
}

double maxwell_juettner_dexter_V(struct parameters *params)
{
  double Ne = params->electron_density;
  double nu = params->nu;
  double Thetae = params->theta_e;
  double B = params->magnetic_field;
  double theta = params->observer_angle;
  double nus = 3.0 * EE * B * sin(theta) / 4.0 / M_PI / ME / CL * Thetae * Thetae + 1.0;
  double x = nu / nus;
  double ans = 2. * Ne * EE * EE * nu / tan(theta) / 3. / sqrt(3) / CL / Thetae / Thetae / Thetae * I_V(x);
#if DEBUG
  if (isnan(ans)) {
    fprintf(stderr, "NaN in Dexter jV. Ne, nu, Thetae, B, theta, nus, x: %g %g %g %g %g %g %g\n",
                    Ne, nu, Thetae, B, theta, nus, x);
  }
#endif
  return ans;
}

// Supporting functions

double I_I(double x)
{
  return 2.5651 * (1 + 1.92 * pow(x, -1. / 3.) +
      0.9977 * pow(x, -2. / 3.)) * exp(-1.8899 * pow(x,
        1. /
        3.));
}

double I_Q(double x)
{
  return 2.5651 * (1 + 0.93193 * pow(x, -1. / 3.) +
      0.499873 * pow(x, -2. / 3.)) * exp(-1.8899 * pow(x,
        1. /
        3.));
}

double I_V(double x)
{
  return (1.81384 / x + 3.42319 * pow(x, -2. / 3.) +
      0.0292545 * pow(x, -0.5) + 2.03773 * pow(x,
        -1. / 3.)) *
    exp(-1.8899 * pow(x, 1. / 3.));
}

/********* LEUNG FIT *********/

/*
 * thermal synchrotron emissivity
 *
 * Interpolates between Petrosian limit and
 * classical thermal synchrotron limit
 * Good for Thetae >~ 1
 * See Leung+ 2011, restated Pandya+ 2016
 */
double maxwell_juettner_leung_I(struct parameters *params)
{
  double Ne = params->electron_density;
  double nu = params->nu;
  double Thetae = params->theta_e;
  double B = params->magnetic_field;
  double theta = params->observer_angle;

  double K2 = fmax(gsl_sf_bessel_Kn(2,1./Thetae), SMALL);

  double nuc = EE*B/(2.*M_PI*ME*CL);
  double nus = (2./9.)*nuc*Thetae*Thetae*sin(theta);

  if(nu > 1.e12*nus)
    return 0.;

  double x = nu/nus ;
  double f = pow( pow(x,1./2.) + pow(2.,11./12.)*pow(x,1./6.), 2 );
  double j = (sqrt(2.)*M_PI*EE*EE*Ne*nus/(3.*CL*K2)) * f * exp(-pow(x,1./3.));

#if DEBUG
  if (isnan(j) || isinf(j)) {
    fprintf(stderr, "j nan in Leung fit: j %g f %g x %g nu %g nus %g nuc %g K2 %g Thetae %g\n", j, f, x, nu, nus, nuc, K2, Thetae);
  }
#endif

  return j;
}
