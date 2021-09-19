#include "fits.h"

#include <stdio.h>
#include <stdlib.h>

/*Wrappers for the fitting formulae, ipole-specific calling convention*/

/*j_nu_fit: wrapper for the emissivity fitting formulae.  Takes in the
 *          same parameters as j_nu(), and passes them to the fitting formulae.
 *
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa, 
 *         kappa_width 
 *@returns: the corresponding fitting formula (based on the distribution
 *          function) evaluated for the input parameters.
*/
double j_nu_fit(struct parameters *params, int polarization)
{
  params->polarization = polarization;

#if DEBUG
  check_for_errors(params);
#endif
  
  if(params->distribution == params->MAXWELL_JUETTNER)
  {
    if     (params->polarization == params->STOKES_I) return maxwell_juettner_I(params); 
    else if(params->polarization == params->STOKES_Q) return maxwell_juettner_Q(params);
    else if(params->polarization == params->STOKES_U) return 0.;
    else if(params->polarization == params->STOKES_V) return maxwell_juettner_V(params);
  }

  else if(params->distribution == params->POWER_LAW)
  {
    if     (params->polarization == params->STOKES_I) return power_law_I(params);
    else if(params->polarization == params->STOKES_Q) return power_law_Q(params);
    else if(params->polarization == params->STOKES_U) return 0.;
    else if(params->polarization == params->STOKES_V) return power_law_V(params);
  }

  else if(params->distribution == params->KAPPA_DIST)
  {
    if (params->polarization == params->STOKES_I) {
      if (params->kappa < params->kappa_interp_begin) {
        return kappa_I(params);
      } else if (params->kappa >= params->kappa_interp_begin && params->kappa < params->kappa_interp_end) {
        return ((params->kappa_interp_end - params->kappa) * maxwell_juettner_I(params) +
                (params->kappa - params->kappa_interp_begin) * kappa_I(params)) /
                (params->kappa_interp_end - params->kappa_interp_begin);
      } else {
        return maxwell_juettner_I(params);
      }
    } else if (params->polarization == params->STOKES_Q) {
      if (params->kappa < params->kappa_interp_begin) {
        return kappa_Q(params);
      } else if (params->kappa >= params->kappa_interp_begin && params->kappa < params->kappa_interp_end) {
        return ((params->kappa_interp_end - params->kappa) * maxwell_juettner_Q(params) +
                (params->kappa - params->kappa_interp_begin) * kappa_Q(params)) /
                (params->kappa_interp_end - params->kappa_interp_begin);
      } else {
        return maxwell_juettner_Q(params);
      }
    } else if (params->polarization == params->STOKES_U) {
      return 0.;
    } else if (params->polarization == params->STOKES_V) {
      if (params->kappa < params->kappa_interp_begin) {
        return kappa_V(params);
      } else if (params->kappa >= params->kappa_interp_begin && params->kappa < params->kappa_interp_end) {
        return ((params->kappa_interp_end - params->kappa) * maxwell_juettner_V(params) +
                (params->kappa - params->kappa_interp_begin) * kappa_V(params)) /
                (params->kappa_interp_end - params->kappa_interp_begin);
      } else {
        return maxwell_juettner_V(params);
      }
    }
  }

  return 0.;
}

/*alpha_nu_fit: wrapper for the absorptivity fitting formulae.  Takes in the
 *              same parameters as alpha_nu(), and passes them to the fitting formulae.
 *
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa, 
 *         kappa_width 
 *@returns: the corresponding fitting formula (based on the distribution
 *          function) evaluated for the input parameters.
*/

double alpha_nu_fit(struct parameters *params, int polarization)
{
  params->polarization = polarization;

#if DEBUG
  check_for_errors(params);
#endif

  if(params->distribution == params->MAXWELL_JUETTNER)
  {
    if     (params->polarization == params->STOKES_I) return maxwell_juettner_I_abs(params);
    else if(params->polarization == params->STOKES_Q) return maxwell_juettner_Q_abs(params);
    else if(params->polarization == params->STOKES_U) return 0.;
    else if(params->polarization == params->STOKES_V) return maxwell_juettner_V_abs(params);
  }

  else if(params->distribution == params->POWER_LAW)
  {
    if     (params->polarization == params->STOKES_I) return power_law_I_abs(params);
    else if(params->polarization == params->STOKES_Q) return power_law_Q_abs(params);
    else if(params->polarization == params->STOKES_U) return 0.;
    else if(params->polarization == params->STOKES_V) return power_law_V_abs(params);
  }

  else if(params->distribution == params->KAPPA_DIST)
  {
    if     (params->polarization == params->STOKES_I) {
      if (params->kappa < params->kappa_interp_begin) {
        return kappa_I_abs(params);
      } else if (params->kappa >= params->kappa_interp_begin && params->kappa < params->kappa_interp_end) {
        return ((15 - params->kappa) * maxwell_juettner_I_abs(params) +
                (params->kappa - params->kappa_interp_begin) * kappa_I_abs(params)) /
                (params->kappa_interp_end - params->kappa_interp_begin);
      } else {
        return maxwell_juettner_I_abs(params);
      }
    } else if(params->polarization == params->STOKES_Q) {
      if (params->kappa < params->kappa_interp_begin) {
        return kappa_Q_abs(params);
      } else if (params->kappa >= params->kappa_interp_begin && params->kappa < params->kappa_interp_end) {
        return ((params->kappa_interp_end - params->kappa) * maxwell_juettner_Q_abs(params) +
                (params->kappa - params->kappa_interp_begin) * kappa_Q_abs(params)) /
                (params->kappa_interp_end - params->kappa_interp_begin);
      } else {
        return maxwell_juettner_Q_abs(params);
      }
    } else if(params->polarization == params->STOKES_U) {
      return 0.;
    } else if(params->polarization == params->STOKES_V) {
      if (params->kappa < params->kappa_interp_begin) {
        return kappa_V_abs(params);
      } else if (params->kappa >= params->kappa_interp_begin && params->kappa < params->kappa_interp_end) {
        return ((params->kappa_interp_end - params->kappa) * maxwell_juettner_V_abs(params) +
                (params->kappa - params->kappa_interp_begin) * kappa_V_abs(params)) /
                (params->kappa_interp_end - params->kappa_interp_begin);
      } else {
        return maxwell_juettner_V_abs(params);
      }
    }
  }

  return 0.;
}

/*rho_nu_fit: Fits to Faraday rotation/conversion coefficients; right now only has
 *            formulae for Maxwell-Juettner distribution, from Dexter (2016),
 *            and the Kappa distribution, from Marszewski+ (2021)
 * 
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa, 
 *         kappa_width 
 *@returns: the corresponding Faraday coefficient fitting formula 
 *          (based on the distribution function) evaluated for 
 *          the input parameters.
 */
double rho_nu_fit(struct parameters *params, int polarization)
{
  params->polarization = polarization;

#if DEBUG
  check_for_errors(params);
  if(params->polarization == params->STOKES_I)
  {
    printf("No Faraday rotation of total intensity");
    return 0.;
  }
#endif

  if(params->distribution == params->MAXWELL_JUETTNER)
  {
    if     (params->polarization == params->STOKES_Q) return maxwell_juettner_rho_Q(params);
    else if(params->polarization == params->STOKES_U) return 0.;
    else if(params->polarization == params->STOKES_V) return maxwell_juettner_rho_V(params);
  
  } else if(params->distribution == params->KAPPA_DIST) {
#if DEBUG
    if(params->nu/(2.8e6 * params->magnetic_field) < 100 || 
       params->nu/(pow(params->kappa_width * params->kappa, 2.) * (2.8e6 * params->magnetic_field) * sin(params->observer_angle)) < 0.1)
    {
      printf("\n WARNING: nu and/or X_kappa low; rho kappa fits may be inaccurate \n");
    }
#endif
  if (params->polarization == params->STOKES_Q)
  {
    if (params->kappa < 3.5)
      return kappa35_rho_Q(params);
    else if (params->kappa >= 3.5 && params->kappa < 4.0)
      return ((4.0 - params->kappa) * kappa35_rho_Q(params) + (params->kappa - 3.5) * kappa4_rho_Q(params)) / 0.5;
    else if (params->kappa >= 4.0 && params->kappa < 4.5)
      return ((4.5 - params->kappa) * kappa4_rho_Q(params) + (params->kappa - 4.0) * kappa45_rho_Q(params)) / 0.5;
    else if (params->kappa >= 4.5 && params->kappa < 5.0)
      return ((5.0 - params->kappa) * kappa45_rho_Q(params) + (params->kappa - 4.5) * kappa5_rho_Q(params)) / 0.5;
    else if (params->kappa >= 5.0 && params->kappa <= 8.0)
      return ((8.0 - params->kappa) * kappa5_rho_Q(params) + (params->kappa - 5.0) * maxwell_juettner_rho_Q(params)) / 5.0;
    else if (params->kappa > 8.0)
      return maxwell_juettner_rho_Q(params);
  }
  else if(params->polarization == params->STOKES_U) return 0.;
  else if(params->polarization == params->STOKES_V)
  {
    if(params->kappa < 3.5)
      return kappa35_rho_V(params);
    else if(params->kappa >= 3.5 && params->kappa < 4.0)
      return ((4.0 - params->kappa) * kappa35_rho_V(params) + (params->kappa - 3.5) * kappa4_rho_V(params)) / 0.5;
    else if(params->kappa >= 4.0 && params->kappa < 4.5)
      return ((4.5 - params->kappa) * kappa4_rho_V(params) + (params->kappa - 4.0) * kappa45_rho_V(params)) / 0.5;
    else if(params->kappa >= 4.5 && params->kappa <= 5.0)
      return ((5.0 - params->kappa) * kappa45_rho_V(params) + (params->kappa - 4.5) * kappa5_rho_V(params)) / 0.5;
    else if (params->kappa >= 5.0 && params->kappa <= 8.0)
      return ((8.0 - params->kappa) * kappa5_rho_V(params) + (params->kappa - 5.0) * maxwell_juettner_rho_V(params)) / 5.0;
    else if (params->kappa > 8.0)
      return maxwell_juettner_rho_V(params);
  }
} else if(params->distribution == params->POWER_LAW) {
  fprintf(stderr, "\nNo rotativity fits are implemented for power-law diestributions!!\n");
  exit(-1);
}

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
void check_for_errors(struct parameters *params)
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
  if(params->distribution == params->KAPPA_DIST) {
    // This is taken care of by switching fits, see above
    // if(params->kappa < 2.5 || params->kappa > 7.5)
    // {
    //   printf("\n WARNING: kappa out of range of fitting formula \n");
    // }
    //if(params->kappa_width < 3 || 
    if (params->kappa_width > 200)
    {
      printf("\n WARNING: w out of range; fitting formula may be inaccurate\n");
    }
  }
  if(params->distribution == params->POWER_LAW) {
    if(params->gamma_min < 1)
    {
      printf("\n ERROR: gamma_min < 1\n");
      exit(0);
    }
  }
  // if(params->observer_angle < 5.*(params->pi)/180. 
  //   || params->observer_angle == 90.*(params->pi)/180.)
  // {
  //   printf("\n WARNING: theta out of range; fitting formula may be inaccurate \n");
  // }
}