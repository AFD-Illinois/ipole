#include "power_law.h"

/*power_law_I: fitting formula to the emissivity, in Stokes I, produced by a
 *             power-law distribution of electrons (without any exponential
 *             cutoff).  Uses eq. 29, 33 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law emissivity polarized in Stokes I
 */
double power_law_I(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double prefactor = (params->electron_density*pow(params->electron_charge,2.)
                      *nu_c)/params->speed_light;

  double term1 = pow(3., params->power_law_p/2.)*(params->power_law_p-1.)
                 *sin(params->observer_angle);

  double term2 = 2.*(params->power_law_p+1.)
                 *(pow(params->gamma_min, 1.-params->power_law_p)
                   -pow(params->gamma_max, 1.-params->power_law_p));

  double term3 = tgamma((3.*params->power_law_p-1.)/12.)
                *tgamma((3.*params->power_law_p+19.)/12.);

  double term4 = pow(params->nu/(nu_c*sin(params->observer_angle)), 
                     -(params->power_law_p-1.)/2.);

  double ans = prefactor*term1/term2*term3*term4;

  return ans;
}

/*power_law_Q: fitting formula to the emissivity, in Stokes Q, produced by a
 *             power-law distribution of electrons (without any exponential
 *             cutoff).  Uses eq. 29, 33 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law emissivity polarized in Stokes Q
 */
double power_law_Q(struct parameters * params)
{
  double p_term = -(params->power_law_p + 1.)/(params->power_law_p + 7./3.);

  double ans = p_term * power_law_I(params);

  return ans;
}

/*power_law_V: fitting formula to the emissivity, in Stokes V, produced by a
 *             power-law distribution of electrons (without any exponential
 *             cutoff).  Uses eq. 29, 33 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law emissivity polarized in Stokes V
 */
double power_law_V(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double term1 = -(171./250.)*pow(params->power_law_p, 49./100.);

  double term2 = 1./tan(params->observer_angle) 
                 * pow(params->nu
                       /(3.*nu_c*sin(params->observer_angle)), -1./2.);
  
  double ans = term1*term2*power_law_I(params);
 
  /*NOTE: Sign corrected; the sign in Leung et al. (2011)
    and Pandya et al. (2016) for Stokes V transfer coefficients
    does not follow the convention the papers describe (IEEE/IAU);
    the sign has been corrected here.*/ 
  return -ans;
}

/*power_law_I_abs: fitting formula to the absorptivity, in Stokes I, from
 *                 by a power-law distribution of electrons (without any 
 *                 exponential cutoff).  Uses eq. 30, 34 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law absorptivity polarized in Stokes I
 */
double power_law_I_abs(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double prefactor = (params->electron_density*pow(params->electron_charge,2.))
                    /(params->nu*params->mass_electron*params->speed_light);

  double term1 = pow(3., (params->power_law_p+1.)/2.)*(params->power_law_p-1.);

  double term2 = 4.*(pow(params->gamma_min, 1.-params->power_law_p)
                -pow(params->gamma_max, 1.-params->power_law_p));

  double term3 = tgamma((3.*params->power_law_p+2.)/12.)
                *tgamma((3.*params->power_law_p+22.)/12.);

  double term4 = pow(params->nu/(nu_c*sin(params->observer_angle)), 
                     -(params->power_law_p+2.)/2.);

  double ans = prefactor*term1/term2*term3*term4;

  return ans;
}

/*power_law_Q_abs: fitting formula to the absorptivity, in Stokes Q, from
 *                 by a power-law distribution of electrons (without any 
 *                 exponential cutoff).  Uses eq. 30, 34 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law absorptivity polarized in Stokes Q
 */
double power_law_Q_abs(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double prefactor = (params->electron_density*pow(params->electron_charge,2.))
                    /(params->nu*params->mass_electron*params->speed_light);

  double term1 = pow(3., (params->power_law_p+1.)/2.)*(params->power_law_p-1.);

  double term2 = 4.*(pow(params->gamma_min, 1.-params->power_law_p)
                -pow(params->gamma_max, 1.-params->power_law_p));

  double term3 = tgamma((3.*params->power_law_p+2.)/12.)
                *tgamma((3.*params->power_law_p+22.)/12.);

  double term4 = pow(params->nu/(nu_c*sin(params->observer_angle)), 
                     -(params->power_law_p+2.)/2.);

  double term5 = -pow((17./500.)*params->power_law_p - 43./1250., 43./500.);

  double ans = prefactor*term1/term2*term3*term4*term5;

  return ans;
}

/*power_law_V_abs: fitting formula to the absorptivity, in Stokes V, from
 *                 by a power-law distribution of electrons (without any 
 *                 exponential cutoff).  Uses eq. 30, 34 of [1], but
 *                 is modified slightly for increased accuracy at the
 *                 cost of a more complicated function.
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law absorptivity polarized in Stokes V
 */
double power_law_V_abs(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double prefactor = (params->electron_density*pow(params->electron_charge,2.))
                    /(params->nu*params->mass_electron*params->speed_light);

  double term1 = pow(3., (params->power_law_p+1.)/2.)*(params->power_law_p-1.);

  double term2 = 4.*(pow(params->gamma_min, 1.-params->power_law_p)
                -pow(params->gamma_max, 1.-params->power_law_p));

  double term3 = tgamma((3.*params->power_law_p+2.)/12.)
                *tgamma((3.*params->power_law_p+22.)/12.);

  double term4 = pow(params->nu/(nu_c*sin(params->observer_angle)), 
                     -(params->power_law_p+2.)/2.);

  double term5 = -pow((71./100.)*params->power_law_p+22./625.,197./500.);

  double term6 = pow((31./10.)*pow(sin(params->observer_angle),-48./25)-31./10.,
                     64./125.);

  double term7 = pow(params->nu/(nu_c*sin(params->observer_angle)), -1./2.);
 
  double ans = prefactor*term1/term2*term3*term4*term5*term6*term7;

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


/*power_law_rho_Q: Fitting formula for Faraday conversion coefficient
 *                        rho_Q from Dexter (2016) for a power-law
 *                        distribution.  This formula comes from his
 *                        equations B1, B2, and B3.
 *
 *@params: struct of parameters params
 *@returns: fitting formula to the Faraday conversion coefficient
 *          for a power-law distribution of electrons. 
 */
double power_law_rho_Q(struct parameters *params)
{
  double echargehere = params->electron_charge;
  double mehere = params->mass_electron;
  double chere = params->speed_light;
  double pihere = params->speed_light;
  double Bhere = params->magnetic_field;
  double thetaB = params->observer_angle;
  double nehere = params->electron_density;
  double nuhere = params->nu;
  double phere = params->power_law_p;
  double gammamin = params->gamma_min;
  double gammamax = params->gamma_max;

  double nuB = echargehere*Bhere / (mehere*chere) / (2.0 * pihere);
  double numin = 1.5*sin(thetaB)*nuB*gammamin*gammamin;
  double rhoperp = nehere * pow(echargehere, 2.0) / (mehere * chere * nuB * sin(thetaB)) * (phere - 1.0) / (pow(gammamin, 1.0-phere)-pow(gammamax, 1.0-phere));
  double rhoQ = -rhoperp*pow(nuB/nuhere,3.0)*pow(gammamin, 2.0-phere)*(1.0-pow(2.0/3.0*numin/nuhere, phere/2.0-1.0))/(phere/2.0-1.0); //note: there's an ambiguity about whether to use the factor of 2/3 (see arXiv vs. published version of Dexter and the Jones+Odell 1977 paper)

  return rhoQ;
}

/*power_law_rho_V: Fitting formula for Faraday rotation coefficient
 *                        rho_V from Dexter (2016) for a power-law
 *                        distribution.  This formula comes from his
 *                        equations B2 and B3.
 *
 *@params: struct of parameters params
 *@returns: fitting formula to the Faraday rotation coefficient
 *          for a power-law distribution of electrons. 
 */
double power_law_rho_V(struct parameters *params)
{
  double echargehere = params->electron_charge;
  double mehere = params->mass_electron;
  double chere = params->speed_light;
  double pihere = params->speed_light;
  double Bhere = params->magnetic_field;
  double thetaB = params->observer_angle;
  double nehere = params->electron_density;
  double nuhere = params->nu;
  double phere = params->power_law_p;
  double gammamin = params->gamma_min;
  double gammamax = params->gamma_max;

  double nuB = echargehere*Bhere / (mehere*chere) / (2.0 * pihere);
  double rhoperp = nehere * pow(echargehere, 2.0) / (mehere * chere * nuB * sin(thetaB)) * (phere - 1.0) / (pow(gammamin, 1.0-phere)-pow(gammamax, 1.0-phere));
  double rhoV = 2.0*(phere+2.0)/(phere+1.0)*rhoperp*pow(nuB*sin(thetaB)/nuhere,2.0)*pow(gammamin, -(phere+1.0))*log(gammamin)/tan(thetaB);

  return rhoV;
}


