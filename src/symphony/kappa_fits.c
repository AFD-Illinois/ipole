#include "kappa.h"

#include "stdio.h"

#define SMALL 1e-40

/**
 * Stabilized version of the hypergeometric fn:
 * Use GSL's version on small domain of validity,
 * use an extension to maintain stability outside that
 *
 * Final form courtesy of Angelo Ricarte
 */
static inline double stable_hyp2f1(double a, double b, double c, double z)
{
  if (z > -1 && z < 1) {
    return gsl_sf_hyperg_2F1(a, b, c, z);
  } else {
    /*GSL 2F1 only works for |z| < 1; had to apply a hypergeometric function
    identity because in our case z = -kappa*w, so |z| > 1 */
    return pow(1.-z, -a) * tgamma(c) * tgamma(b-a)
                 / (tgamma(b)*tgamma(c-a))
                 * gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))
                 + pow(1.-z, -b) * tgamma(c) * tgamma(a-b)
                 / (tgamma(a) * tgamma(c-b))
                 * gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));
  }
}

/*kappa_I: fitting formula to the emissivity, in Stokes I, produced by a
 *         kappa distribution of electrons (without any exponential
 *         cutoff).  Uses eq. 29, 35, 36, 37, 38 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to kappa emissivity polarized in Stokes I
 */
double kappa_I(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_w = pow(params->kappa_width * params->kappa, 2.) *nu_c
                *sin(params->observer_angle);

  double X_k = params->nu/nu_w;

  double prefactor = (params->electron_density 
                      * pow(params->electron_charge, 2.) * nu_c 
                      * sin(params->observer_angle))
	             /params->speed_light;

  double Nlow = 4. * params->pi * tgamma(params->kappa-4./3.) 
                / (pow(3., 7./3.) * tgamma(params->kappa-2.));

  double Nhigh = (1./4.) * pow(3., (params->kappa-1.)/2.) 
                 * (params->kappa-2.) * (params->kappa-1.)
		 * tgamma(params->kappa/4.-1./3.) 
                 * tgamma(params->kappa/4.+4./3.);

  double x = 3. * pow(params->kappa, -3./2.);

  double ans = prefactor * Nlow * pow(X_k, 1./3.) 
               * pow(1.+pow(X_k, x * (3. * params->kappa-4.)/6.)
                     * pow(Nlow/Nhigh, x), -1./x);

  return ans;

}

/*kappa_Q: fitting formula to the emissivity, in Stokes Q, produced by a
 *         kappa distribution of electrons (without any exponential
 *         cutoff).  Uses eq. 29, 35, 36, 37, 38 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to kappa emissivity polarized in Stokes Q
 */
double kappa_Q(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_w = pow(params->kappa_width*params->kappa, 2.) * nu_c 
                * sin(params->observer_angle);

  double X_k = params->nu/nu_w;

  double prefactor = (params->electron_density 
                      * pow(params->electron_charge, 2.) 
                      * nu_c * sin(params->observer_angle))
	             /params->speed_light;

  double Nlow = -(1./2.) * 4. * params->pi * tgamma(params->kappa-4./3.)
                 /(pow(3., 7./3.) * tgamma(params->kappa-2.));

  double Nhigh = -(pow(4./5., 2)+params->kappa/50.) * (1./4.) 
                 * pow(3., (params->kappa-1.)/2.) * (params->kappa-2.) 
                 * (params->kappa-1.) * tgamma(params->kappa/4.-1./3.)
                 * tgamma(params->kappa/4.+4./3.);
 
  double x = (37./10.)*pow(params->kappa, -8./5.);

  double ans = prefactor * Nlow * pow(X_k, 1./3.)
              * pow(1. + pow(X_k, x * (3. * params->kappa-4.)/6.)
              * pow(Nlow/Nhigh, x), -1./x);

  return ans;

}

/*kappa_V: fitting formula to the emissivity, in Stokes V, produced by a
 *         kappa distribution of electrons (without any exponential
 *         cutoff).  Uses eq. 29, 35, 36, 37, 38 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to kappa emissivity polarized in Stokes V
 */
double kappa_V(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_w = pow(params->kappa_width * params->kappa, 2.) 
                * nu_c * sin(params->observer_angle);

  double X_k = params->nu/nu_w;

  double prefactor = (params->electron_density 
                     * pow(params->electron_charge, 2.)
                     * nu_c * sin(params->observer_angle))
                    /params->speed_light;

  double Nlow = -pow(3./4., 2.) 
                * pow(pow(sin(params->observer_angle), -12./5.)-1., 12./25.)
                * (pow(params->kappa, -66./125.) / params->kappa_width) 
                * pow(X_k, -7./20.) * 4. * params->pi 
                * tgamma(params->kappa-4./3.)
                / (pow(3., 7./3.)*tgamma(params->kappa-2.));

  double Nhigh = -pow(7./8., 2.) 
                 * pow(pow(sin(params->observer_angle), -5./2.)-1., 11./25.)
                 * (pow(params->kappa, -11./25.) / params->kappa_width) 
                 * pow(X_k, -1./2.) * (1./4.) * pow(3., (params->kappa-1.)/2.) 
                 * (params->kappa-2.) * (params->kappa-1.) 
                 * tgamma(params->kappa/4.-1./3.) 
                 * tgamma(params->kappa/4.+4./3.);

  double x = 3.*pow(params->kappa, -3./2.);

  double ans = (Nhigh < SMALL*SMALL) ? 0:
               prefactor * Nlow * pow(X_k, 1./3.)
               * pow(1.+pow(X_k, x * (3.*params->kappa-4.)/6.)
               * pow(Nlow/Nhigh, x), -1./x);

  /*The Stokes V absorption coefficient changes sign at observer_angle
    equals 90deg, but this formula does not.  This discrepancy is a 
    bug in this formula, and is patched by the term below.*/
  double sign_bug_patch = cos(params->observer_angle) /
                          fabs(cos(params->observer_angle));

  /*NOTE: Sign corrected; the sign in Leung et al. (2011)
    and Pandya et al. (2016) for Stokes V transfer coefficients
    does not follow the convention the papers describe (IEEE/IAU);
    the sign has been corrected here.*/

//  if (isnan(ans) || isnan(sign_bug_patch)) {
//    fprintf(stderr, "NaN in kappa. X_k, prefactor, Nlow, Nhigh, x, ans, sign: %g %g %g %g %g %g %g",
//                     X_k, prefactor, Nlow, Nhigh, x, ans, sign_bug_patch);
//  }

  return -ans * sign_bug_patch;

}

/*kappa_I_abs: fitting formula to the absorptivity, in Stokes I, from
 *             by a kappa distribution of electrons (without any 
 *             exponential cutoff).  Uses eq. 30, 39, 40, 41, 42 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to kappa absorptivity polarized in Stokes I
 */
double kappa_I_abs(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_w = pow(params->kappa_width * params->kappa, 2.) 
                * nu_c * sin(params->observer_angle);

  double X_k = params->nu/nu_w;

  double prefactor = params->electron_density * params->electron_charge 
                     / (params->magnetic_field * sin(params->observer_angle));

  double a = params->kappa - 1./3.;
 
  double b = params->kappa + 1.;

  double c = params->kappa + 2./3.;

  double z = -params->kappa*params->kappa_width;

  double hyp2f1 = stable_hyp2f1(a, b, c, z);

  double Nlow = pow(3., 1./6.) * (10./41.) * pow(2. * params->pi, 2.)
               / pow(params->kappa_width * params->kappa, 16./3.-params->kappa)
               * (params->kappa-2.) * (params->kappa-1.) * params->kappa
               / (3.*params->kappa-1.) * tgamma(5./3.) * hyp2f1;

  double Nhigh = 2. * pow(params->pi, 5./2.)/3. * (params->kappa-2.) 
                 * (params->kappa-1.) * params->kappa
                 / pow(params->kappa_width * params->kappa, 5.) 
                 * (2 * tgamma(2. + params->kappa/2.)
                 / (2.+params->kappa)-1.) 
                 * (pow(3./params->kappa, 19./4.) + 3./5.);

  double x = pow(-7./4. + 8. * params->kappa/5., -43./50.);

  double ans = prefactor * Nlow * pow(X_k, -5./3.)
               * pow(1. + pow(X_k, x * (3. * params->kappa-1.)/6.)
               * pow(Nlow/Nhigh, x), -1./x);

  return ans;
}

/*kappa_Q_abs: fitting formula to the absorptivity, in Stokes Q, from
 *             by a kappa distribution of electrons (without any 
 *             exponential cutoff).  Uses eq. 30, 39, 40, 41, 42 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to kappa absorptivity polarized in Stokes Q
 */
double kappa_Q_abs(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_w = pow(params->kappa_width*params->kappa, 2.) 
                * nu_c * sin(params->observer_angle);

  double X_k = params->nu/nu_w;

  double prefactor = params->electron_density * params->electron_charge
                     / (params->magnetic_field * sin(params->observer_angle));

  double a = params->kappa - 1./3.;
 
  double b = params->kappa + 1.;

  double c = params->kappa + 2./3.;

  double z = -params->kappa * params->kappa_width;

  double hyp2f1 = stable_hyp2f1(a, b, c, z);

  double Nlow = -(25./48.) * pow(3., 1./6.) * (10./41.) 
                * pow(2. * params->pi, 2.)
                / pow(params->kappa_width*params->kappa, 16./3.-params->kappa)
                * (params->kappa-2.) * (params->kappa-1.)
                * params->kappa / (3. * params->kappa-1.) 
                * tgamma(5./3.) * hyp2f1;
  double Nhigh = -(pow(21., 2.) * pow(params->kappa, -144./25.) + 11./20.) 
                 * 2. * pow(params->pi, 5./2.)/3. * (params->kappa-2.) 
                 * (params->kappa-1.) * params->kappa 
                 / pow(params->kappa_width * params->kappa, 5.) 
                 * (2 * tgamma(2. + params->kappa/2.)
                 / (2. + params->kappa)-1.);

  double x = (7./5.) * pow(params->kappa, -23./20.);

  double ans = prefactor * Nlow * pow(X_k, -5./3.)
              * pow(1. + pow(X_k, x * (3. * params->kappa-1.) / 6.)
              * pow(Nlow/Nhigh, x), -1./x);

  return ans;
}

/*kappa_V_abs: fitting formula to the absorptivity, in Stokes V, from
 *             by a kappa distribution of electrons (without any 
 *             exponential cutoff).  Uses eq. 30, 39, 40, 41, 42 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to kappa absorptivity polarized in Stokes V
 */
double kappa_V_abs(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_w = pow(params->kappa_width * params->kappa, 2.) 
               * nu_c * sin(params->observer_angle);

  double X_k = params->nu/nu_w;

  double prefactor = params->electron_density * params->electron_charge 
                     / (params->magnetic_field * sin(params->observer_angle));

  double a = params->kappa - 1./3.;
 
  double b = params->kappa + 1.;

  double c = params->kappa + 2./3.;

  double z = -params->kappa * params->kappa_width;

  double hyp2f1 = stable_hyp2f1(a, b, c, z);

  double Nlow = -(77./(100. * params->kappa_width)) 
                * pow(pow(sin(params->observer_angle), -114./50.)
                -1., 223./500.) * pow(X_k, -7./20.) 
                * pow(params->kappa, -7./10) * pow(3., 1./6.) 
                * (10./41.) * pow(2. * params->pi, 2.)
                / pow(params->kappa_width*params->kappa, 16./3.-params->kappa)
                * (params->kappa-2.) * (params->kappa-1.)
                * params->kappa / (3. * params->kappa-1.)
                * tgamma(5./3.) * hyp2f1;

  double Nhigh = -(143./10. * pow(params->kappa_width, -116./125.))
                 * pow(pow(sin(params->observer_angle), -41./20.)-1., 1./2.)
                 * (13.*13. * pow(params->kappa, -8.) + 13./(2500.) 
                 * params->kappa - 263./5000. + 47.
                 / (200.*params->kappa)) * pow(X_k, -1./2.) 
                 * 2. * pow(params->pi, 5./2.) / 3.
                 * (params->kappa-2.) * (params->kappa-1.) * params->kappa 
                 / pow(params->kappa_width*params->kappa, 5.)
                 * (2 * tgamma(2. + params->kappa/2.)
                 / (2. + params->kappa) - 1.);

  double x = (61./50.)*pow(params->kappa, -142./125.)+7./1000.;

  double ans = prefactor * Nlow * pow(X_k, -5./3.)
               * pow(1.+pow(X_k, x * (3. * params->kappa-1.) / 6.)
               * pow(Nlow/Nhigh, x), -1./x);

  /*The Stokes V absorption coefficient changes sign at observer_angle
    equals 90deg, but this formula does not.  This discrepancy is a 
    bug in this formula, and is patched by the term below.*/
  double sign_bug_patch = cos(params->observer_angle) /
                          fabs(cos(params->observer_angle));

  /*NOTE: Sign corrected; the sign in Leung et al. (2011)
    and Pandya et al. (2016) for Stokes V transfer coefficients
    does not follow the convention the papers describe (IEEE/IAU);
    the sign has been corrected here.*/
  return -ans * sign_bug_patch;
}
