/**********************************************************/
/*** all you need to make a polarized radiative transfer***/
/***** used in ipole to evolve complex tensor N ***********/
/***** along with standard evolution for I scalar *********/
/**********************************************************/
/**** written by Monika Moscibrodzka on 09 July 2014 ******/
/************ @ Eindhoven Airport *************************/
/****************and then rewritten by C.Gammie ***********/
/**********************************************************/

#include "ipolarray.h"

#include "decs.h"
#include "coordinates.h"
#include "geometry.h"
#include "model.h"
#include "model_radiation.h"
#include "model_tetrads.h"
#include "debug_tools.h"
#include <complex.h>

// These are related mostly to where exp() and 1/x overflow
// desired accuracy.
#define CUT_HIGH_ABS 500
#define CUT_LOW_ABS SMALL

// Sub-functions
void push_polar(double Xi[NDIM], double Xm[NDIM], double Xf[NDIM],
    double Ki[NDIM], double Km[NDIM], double Kf[NDIM],
    complex double Ni[NDIM][NDIM],
    complex double Nm[NDIM][NDIM],
    complex double Nf[NDIM][NDIM], double dlam);

/* tensor tools */
void complex_lower(double complex N[NDIM][NDIM], double gcov[NDIM][NDIM],
    int low1, int low2, double complex Nl[NDIM][NDIM]);
void stokes_to_tensor(double fI, double fQ, double fU, double fV,
    double complex f_tetrad[NDIM][NDIM]);
void tensor_to_stokes(double complex f_tetrad[NDIM][NDIM], double *fI,
    double *fQ, double *fU, double *fV);
void complex_coord_to_tetrad_rank2(double complex T_coord[NDIM][NDIM],
    double Ecov[NDIM][NDIM],
    double complex T_tetrad[NDIM][NDIM]);
void complex_tetrad_to_coord_rank2(double complex T_tetrad[NDIM][NDIM],
    double Econ[NDIM][NDIM],
    double complex T_coord[NDIM][NDIM]);



/************************PRIMARY FUNCTION*******************************/
/**
 * Calculate the emission produced/absorbed/rotated along the given trajectory traj
 * of total length nsteps.
 * Return arguments of intensity, total optical depth, total Faraday depth, and complex polarized emission tensor N^alpha^beta
 * 
 * Returns flag indicating at least one step either used a questionable tetrad, or produced a NaN value
 */
int integrate_emission(struct of_traj *traj, int nsteps,
                    double *Intensity, double *Tau, double *tauF,
                    double complex N_coord[NDIM][NDIM], Params *params)
{
  //fprintf(stderr, "Begin integrate emission\n");
  // Initialize
  MUNULOOP N_coord[mu][nu] = 0.0 + I * 0.0;
  *tauF = 0.;
  // Unpolarized
  *Intensity = 0.;
  *Tau = 0.;
  // Error flag
  int oddflag = 0;

  // Integrate the transfer equation (& parallel transport) forwards along trajectory
  for (int nstep=nsteps; nstep > 0; --nstep) {
    int sflag = 0;
    struct of_traj ti = traj[nstep];
    struct of_traj tf = traj[nstep-1];

    // Parallel transport polarization vector if necessary
    //fprintf(stderr, "Push Polar\n");
    if (!params->only_unpolarized) {
      double complex Nh[NDIM][NDIM];
      push_polar(ti.X, ti.X, ti.Xhalf, ti.Kcon, ti.Kcon, ti.Kconhalf, N_coord, N_coord, Nh, 0.5 * ti.dl);
      push_polar(ti.X, ti.Xhalf, tf.X, ti.Kcon, ti.Kconhalf, tf.Kcon, N_coord, Nh, N_coord, ti.dl);
    }

#if THIN_DISK
    if (thindisk_region(ti.X, tf.X)) {
      // The thin disk problem emits nowhere but uses a boundary condition region defined by thindisk_region
      // There we just get a starting value for intensity with get_model_i
      get_model_i(ti.X, ti.Kcon, Intensity);

      if (!params->only_unpolarized) {
        // For polarized emission, it sets all Stokes parameters according to the function get_model_stokes
        // We also have to transform into & out of fluid frame

        // Make a tetrad
        double gcov[NDIM][NDIM];
        gcov_func(tf.X, gcov);
        double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
        get_model_fourv(tf.X, tf.Kcon, Ucon, Ucov, Bcon, Bcov);
        double Ecov[NDIM][NDIM], Econ[NDIM][NDIM];
        sflag |= make_plasma_tetrad(Ucon, tf.Kcon, Bcon, gcov, Econ, Ecov);

        // Get the Stokes parameters
        double SI, SQ, SU, SV;
        get_model_stokes(tf.X, tf.Kcon, &SI, &SQ, &SU, &SV);

        // make N_tetrad, and transform
        double complex N_tetrad[NDIM][NDIM];
        complex_coord_to_tetrad_rank2(N_coord, Ecov, N_tetrad);
        stokes_to_tensor(SI, SQ, SU, SV, N_tetrad);
        complex_tetrad_to_coord_rank2(N_tetrad, Econ, N_coord);
      }
    }
#endif

    if (radiating_region(tf.X)) {

      int ZERO_EMISSION = 0;
      if (params->target_nturns >= 0 && ti.nturns != params->target_nturns) {
        ZERO_EMISSION = 1; 
      }

      // Solve unpolarized transport
      double ji, ki, jf, kf;
      get_jkinv(ti.X, ti.Kcon, &ji, &ki, params);
      get_jkinv(tf.X, tf.Kcon, &jf, &kf, params);

      if (ZERO_EMISSION) {
        *Intensity = approximate_solve(*Intensity, 0, ki, 0, kf, ti.dl, Tau);
      } else {
        *Intensity = approximate_solve(*Intensity, ji, ki, jf, kf, ti.dl, Tau);
      }
      //fprintf(stderr, "Unpolarized transport\n");

      // Solve polarized transport
      if (!params->only_unpolarized) {
        sflag |= evolve_N(ti.X, ti.Kcon,
                        ti.Xhalf, ti.Kconhalf,
                        tf.X, tf.Kcon,
                        ti.dl, N_coord, tauF, 
                        ZERO_EMISSION, params);
        //fprintf(stderr, "Polarized transport\n");
      }
    }

    // smoosh together all the flags we hit along a geodesic
    oddflag |= sflag;



    // Cry immediately on bad tetrads, even if we're not debugging
    if (sflag & 1) {
      fprintf(stderr, "that's odd: no orthonormal tetrad found at\n");
      fprintf(stderr, "nstep: %d\n", nstep);
      double r, th;
      bl_coord(tf.X, &r, &th);
      fprintf(stderr, "X: %g %g %g %g\n", tf.X[0], tf.X[1], tf.X[2], tf.X[3]);
      fprintf(stderr, "r,th: %g %g\n", r, th);
      double gcov[NDIM][NDIM];
      gcov_func(tf.X, gcov);
      double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
      get_model_fourv(tf.X, tf.Kcon, Ucon, Ucov, Bcon, Bcov);
      print_vector("Ucon", Ucon);
      print_vector("Kcon", tf.Kcon);
      print_vector("Bcon", Bcon);

      double ucov[4];
      flip_index(Ucon, gcov, ucov);
      double bsq = 0., udotu = 0., udotb = 0., kdotu = 0., kdotb = 0.;
      MULOOP {
        bsq += Bcon[mu] * Bcov[mu];
        udotu += Ucon[mu] * ucov[mu];
        udotb += Bcon[mu] * ucov[mu];
        kdotu += tf.Kcon[mu] * ucov[mu];
        kdotb += tf.Kcon[mu] * Bcov[mu];
      }
      double bsq_reported = get_model_b(tf.X);
      fprintf(stderr, "If all of the following are nonzero, file a bug:\n");
      fprintf(stderr, "bsq = %g bsq_reported = %g u.u = %g  u.b = %g k.u = %g k.b = %g\n",
                      bsq, bsq_reported, udotu, udotb, kdotu, kdotb);
      // exit(-1);
    }
    // Same if there was something in gcov
    if (sflag & 16) {
      fprintf(stderr, "Matrix inversion failed in tetrad check, step %d:\n", nstep);
      // TODO
    }

    // TODO pull more relevant stuff back out here
#if DEBUG
    // Cry on bad tauF
    if (sflag & 2) {
      printf("tauF = %e dlam = %e\n", *tauF, ti.dl);
      fprintf(stderr, "nstep: %d\n", nstep);
      exit(-1);
    }

    // Cry on bad N
    if (sflag & 4) {
      fprintf(stderr, "\nNaN in N00!\n");
      fprintf(stderr, "nstep: %d\n", nstep);
      print_vector("Xi", ti.X);
      print_vector("Xf", tf.X);
      print_vector("Kconi", ti.Kcon);
      print_vector("Kconf", tf.Kcon);
      //printf("Stokes initial: [%e %e %e %e]\n", SI0, SQ0, SU0, SV0);
      //printf("Stokes final: [%e %e %e %e] dlam: %e\n", SI, SQ, SU, SV, dlam);
      //printf("Coefficients: j: [%e %e %e %e] a: [%e %e %e %e] rho: [%e %e %e]\n", jI, jQ, jU, jV, aI, aQ, aU, aV, rQ, rU, rV);
      //MUNULOOP printf("Econ[%i][%i] = %e Ncoord = %e Ntet = %e\n", mu, nu, Econ[mu][nu], creal(N_coord[mu][nu]), creal(N_tetrad[mu][nu]));
      exit(-1);
    }
#endif
  //fprintf(stderr, "End Loop\n");
  }

  //fprintf(stderr, "End integrate emission\n");
  // Otherwise propagate the full flag so caller can yell or record
  return oddflag;
}

/***************************MAIN FUNCTIONS******************************/
/*
 * parallel transport N over dl
 */
void push_polar(double Xi[NDIM], double Xm[NDIM], double Xf[NDIM],
    double Ki[NDIM], double Km[NDIM], double Kf[NDIM],
    complex double Ni[NDIM][NDIM],
    complex double Nm[NDIM][NDIM],
    complex double Nf[NDIM][NDIM], double dlam)
{
#if INTEGRATOR_TEST
  double dl = dlam;
#else
  double dl = dlam / (L_unit * HPL / (ME * CL * CL));
#endif

  /* find the connection */
  double lconn[NDIM][NDIM][NDIM];
  get_connection(Xm, lconn);
  int i, j, k, l;

  /* push N */
  for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++)
  Nf[i][j] = Ni[i][j];

  for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++)
  for (k = 0; k < 4; k++)
  for (l = 0; l < 4; l++)
  Nf[i][j] += -(lconn[i][k][l] * Nm[k][j] * Km[l] +
      lconn[j][k][l] * Nm[i][k] * Km[l]
  ) * dl;

}

/*
 * Updates N for one step on geodesics, using the previous step N
 * here we compute new right-hand side of the equation
 * and somehow rotate this along the geodesics knowing
 * first point and last point X and K
 * 
 * Return an error flag indicating any singular matrix, bad tetrad, etc.
 */
int evolve_N(double Xi[NDIM], double Kconi[NDIM],
    double Xhalf[NDIM], double Kconhalf[NDIM],
    double Xf[NDIM], double Kconf[NDIM],
    double dlam, double complex N_coord[NDIM][NDIM], double *tauF, 
    int ZERO_EMISSION, Params *params)
{
  // TODO might be useful to split this into flat-space S->S portion and transformations to/from N
  double gcov[NDIM][NDIM];
  double Ucon[NDIM],Bcon[NDIM];
  double Ucov[NDIM],Bcov[NDIM];
  double Ecov[NDIM][NDIM], Econ[NDIM][NDIM];
  double complex N_tetrad[NDIM][NDIM];
  double B;
  double jI, jQ, jU, jV;
  double aI, aQ, aU, aV;
  double rV, rU, rQ;
  double rho2, rho, rdS;
  double SI = 0, SQ = 0, SU = 0, SV = 0;
  double SI0, SQ0, SU0, SV0;
  double SI1, SQ1, SU1, SV1;
  double SI2, SQ2, SU2, SV2;
  int oddflag = 0;

  // get fluid parameters at Xf
  get_model_fourv(Xf, Kconf, Ucon, Ucov, Bcon, Bcov);

  // evaluate transport coefficients
  jar_calc(Xf, Kconf, &jI, &jQ, &jU, &jV,
      &aI, &aQ, &aU, &aV, &rQ, &rU, &rV, params);
  if (ZERO_EMISSION) {
    jI = 0.;
    jQ = 0.;
    jU = 0.;
    jV = 0.;
  }

  // Guess B if we *absolutely must*
  // Note get_model_b (rightly) returns 0 outside the domain,
  // but we can cling to the 4-vectors a bit longer
  B = 0.;
  MULOOP B += Bcon[mu] * Bcov[mu];
  if (B <= 0.) {
    Bcon[0] = 0.;
    Bcon[1] = 1.;
    Bcon[2] = 1.;
    Bcon[3] = 1.;
  }

  // make plasma tetrad
  gcov_func(Xf, gcov);
  oddflag |= make_plasma_tetrad(Ucon, Kconf, Bcon, gcov, Econ, Ecov);

  // TODO If B is 0, just keep guessing
  //int exhausted = 0;
  //while (oddflag && B <= 0. && !exhausted) { etc }

  /* convert N to Stokes */
  complex_coord_to_tetrad_rank2(N_coord, Ecov, N_tetrad);
  tensor_to_stokes(N_tetrad, &SI0, &SQ0, &SU0, &SV0);

  /* apply the Faraday rotation solution for a half step */
  double x = dlam * 0.5;

  rdS = rQ * SQ0 + rU * SU0 + rV * SV0;
  rho2 = rQ * rQ + rU * rU + rV * rV;
  rho = sqrt(rho2);
  double c, s, sh;
  c = cos(rho * x);
  s = sin(rho * x);
  sh = sin(0.5 * rho * x);
  if (rho2 > 0) {
    SI1 = SI0;
    SQ1 = SQ0 * c + 2 * rQ * rdS / rho2 * sh * sh + (rU * SV0 - rV * SU0) / rho * s;
    SU1 = SU0 * c + 2 * rU * rdS / rho2 * sh * sh + (rV * SQ0 - rQ * SV0) / rho * s;
    SV1 = SV0 * c + 2 * rV * rdS / rho2 * sh * sh + (rQ * SU0 - rU * SQ0) / rho * s;
  } else {
    SI1 = SI0;
    SQ1 = SQ0;
    SU1 = SU0;
    SV1 = SV0;
  }
  /* done rotation solution half step */

  /* apply full absorption/emission step */
  x = dlam;
  double aI2 = aI * aI;
  double aP2 = aQ * aQ + aU * aU + aV * aV;
  double aP = sqrt(aP2);
  double ads0 = aQ * SQ1 + aU * SU1 + aV * SV1;
  double adj = aQ * jQ + aU * jU + aV * jV;

  if (aP*x > CUT_HIGH_ABS || aI*x > CUT_HIGH_ABS) {
    // Solution assuming aI ~ aP >> 1, the worst case.
    // This covers the case aI >> aP by sending exp(aP-aI)->0, which is safe for all terms
    double expDiffx = exp((aP-aI) * x)/2;

    SI2 = (SI1 * expDiffx
        - (ads0 / aP) * expDiffx
        + adj / (aI2 - aP2) * (-1 + (aI * expDiffx + aP * expDiffx) / aP)
        + aI * jI / (aI2 - aP2) * (1 - (aI * expDiffx + aP * expDiffx) / aI));

    SQ2 = (ads0 * aQ / aP2 * expDiffx
        - aQ / aP * SI1 * expDiffx
        + jQ / aI
        + adj * aQ / (aI * (aI2 - aP2)) *
        (1 - aI / aP2 * (aI * expDiffx + aP * expDiffx))
        + jI * aQ / (aP * (aI2 - aP2)) * (-aP + (aP * expDiffx + aI * expDiffx)));

    SU2 = (ads0 * aU / aP2 * expDiffx
        - aU / aP * SI1 * expDiffx
        + jU / aI
        + adj * aU / (aI * (aI2 - aP2)) *
        (1 - aI / aP2 * (aI * expDiffx + aP * expDiffx))
        + jI * aU / (aP * (aI2 - aP2)) * (-aP + (aP * expDiffx + aI * expDiffx)));

    SV2 = (ads0 * aV / aP2 * expDiffx
        - aV / aP * SI1 * expDiffx
        + jV / aI
        + adj * aV / (aI * (aI2 - aP2)) * (1 -
            aI / aP2 * (aI * expDiffx + aP * expDiffx))
        + jI * aV / (aP * (aI2 - aP2)) *
        (-aP + (aP * expDiffx + aI * expDiffx)));
  } else if (aP*x > CUT_LOW_ABS) { /* full analytic solution has trouble if polarized absorptivity is small */
    double expaIx = exp(-aI * x);
    double sinhaPx = sinh(aP * x);
    double coshaPx = cosh(aP * x);

    SI2 = (SI1 * coshaPx * expaIx
        - (ads0 / aP) * sinhaPx * expaIx
        + adj / (aI2 - aP2) * (-1 + (aI * sinhaPx + aP * coshaPx) / aP * expaIx)
        + aI * jI / (aI2 - aP2) * (1 - (aI * coshaPx + aP * sinhaPx) / aI * expaIx));

    SQ2 = (SQ1 * expaIx
        + ads0 * aQ / aP2 * (-1 + coshaPx) * expaIx
        - aQ / aP * SI1 * sinhaPx * expaIx
        + jQ * (1 - expaIx) / aI
        + adj * aQ / (aI * (aI2 - aP2)) * (1 - (1 - aI2 / aP2) * expaIx
            - aI / aP2 * (aI * coshaPx + aP * sinhaPx) * expaIx)
        + jI * aQ / (aP * (aI2 - aP2)) * (-aP + (aP * coshaPx + aI * sinhaPx) * expaIx));

    SU2 = (SU1 * expaIx
        + ads0 * aU / aP2 * (-1 + coshaPx) * expaIx
        - aU / aP * SI1 * sinhaPx * expaIx
        + jU * (1 - expaIx) / aI
        + adj * aU / (aI * (aI2 - aP2)) *
        (1 - (1 - aI2 / aP2) * expaIx -
            aI / aP2 * (aI * coshaPx +
                aP * sinhaPx) * expaIx)
        + jI * aU / (aP * (aI2 - aP2)) *
        (-aP + (aP * coshaPx + aI * sinhaPx) * expaIx));

    SV2 = (SV1 * expaIx
        + ads0 * aV / aP2 * (-1 + coshaPx) * expaIx
        - aV / aP * SI1 * sinhaPx * expaIx
        + jV * (1 - expaIx) / aI
        + adj * aV / (aI * (aI2 - aP2)) * (1 -
            (1 - aI2 / aP2) * expaIx -
            aI / aP2 * (aI * coshaPx +
                aP * sinhaPx) * expaIx)
        + jI * aV / (aP * (aI2 - aP2)) *
        (-aP + (aP * coshaPx + aI * sinhaPx) * expaIx));

  } else {
    // Still account for aI which may be >> aP, e.g. simulating unpolarized transport
    // Should still make this an expansion in aP as well
    double tau_fake = 0;
    SI2 = approximate_solve(SI1, jI, aI, jI, aI, x, &tau_fake);
    SQ2 = approximate_solve(SQ1, jQ, aI, jQ, aI, x, &tau_fake);
    SU2 = approximate_solve(SU1, jU, aI, jU, aI, x, &tau_fake);
    SV2 = approximate_solve(SV1, jV, aI, jV, aI, x, &tau_fake);
  }
  /* done absorption/emission full step */

  /* apply second rotation half-step */
  x = dlam * 0.5;
  rdS = rQ * SQ2 + rU * SU2 + rV * SV2;
  rho2 = rQ * rQ + rU * rU + rV * rV;
  rho = sqrt(rho2);
  c = cos(rho * x);
  s = sin(rho * x);
  sh = sin(0.5 * rho * x);
  if (rho2 > 0) {
    SI = SI2;
    SQ = SQ2 * c + 2 * rQ * rdS / rho2 * sh * sh + (rU * SV2 - rV * SU2) / rho * s;
    SU = SU2 * c + 2 * rU * rdS / rho2 * sh * sh + (rV * SQ2 - rQ * SV2) / rho * s;
    SV = SV2 * c + 2 * rV * rdS / rho2 * sh * sh + (rQ * SU2 - rU * SQ2) / rho * s;
  } else {
    SI = SI2;
    SQ = SQ2;
    SU = SU2;
    SV = SV2;
  }
  /* done second rotation half-step */

  *tauF += dlam*fabs(rV); //*sqrt(SQ*SQ + SU*SU);

  if (params->stokes_floors) {
    // Correct the resulting Stokes parameters to guarantee:
    // 1. I > 0
    // 2. sqrt(Q^2 + U^2 + V^2) < I
    if (SI < 0) {
      SI = 0;
      SQ = 0;
      SU = 0;
      SV = 0;
    } else {
      double pol_frac = sqrt(SQ*SQ + SU*SU + SV*SV) / SI;
      if ( pol_frac > 1. ) {
        SQ /= pol_frac;
        SU /= pol_frac;
        SV /= pol_frac;
      }
    }
  }

  /* re-pack the Stokes parameters into N */
  stokes_to_tensor(SI, SQ, SU, SV, N_tetrad);
  complex_tetrad_to_coord_rank2(N_tetrad, Econ, N_coord);

  // Record Stokes parameters iff we're doing an integrator test
#if INTEGRATOR_TEST
  static double lam = 0;
  lam += dlam;
  record_stokes_parameters(SI, SQ, SU, SV, lam);
#endif


  // Flag if something is wrong
  if (*tauF > 1.e100 || *tauF < -1.e100 || isnan(*tauF)) oddflag |= 2;
  if (isnan(creal(N_tetrad[0][0])) || isnan(creal(N_coord[0][0]))) oddflag |= 4;

  /* SOURCE STEP DONE */
  return oddflag;
}

/* converts tensor N to Stokes parameters detected at the camera*/
void project_N(double X[NDIM], double Kcon[NDIM],
    double complex N_coord[NDIM][NDIM], double *Stokes_I,
    double *Stokes_Q, double *Stokes_U, double *Stokes_V,
    double rotcam)
{
  double complex N_tetrad[NDIM][NDIM];
  double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];
  double Q,U;

  make_camera_tetrad(X, Econ, Ecov);

  complex_coord_to_tetrad_rank2(N_coord, Ecov, N_tetrad);

  tensor_to_stokes(N_tetrad, Stokes_I, &Q, &U, Stokes_V);

  // rotate Stokes Q, U to be oriented correctly according to
  // a rotated camera
  rotcam *= -2.;
  *Stokes_Q = Q*cos(rotcam) - U*sin(rotcam);
  *Stokes_U = Q*sin(rotcam) + U*cos(rotcam);
}

/*
 * must be a stable, approximate solution to radiative transfer
 * that runs between points w/ initial intensity I, emissivity
 * ji, opacity ki, and ends with emissivity jf, opacity kf.
 *
 * return final intensity
 */
double approximate_solve(double Ii, double ji, double ki, double jf,
    double kf, double dl, double *tau)
{
  double efac, If, javg, kavg, dtau;

  javg = (ji + jf) / 2.;
  kavg = (ki + kf) / 2.;

  dtau = dl * kavg;
  *tau += dtau;

  if (dtau < 1.e-3) {
    If = Ii + (javg - Ii * kavg) * dl * (1. -
        (dtau / 2.) * (1. -
          dtau / 3.));
  } else {
    efac = exp(-dtau);
    If = Ii * efac + (javg / kavg) * (1. - efac);
  }

  return If;
}

/*************************SUPPORTING FUNCTIONS******************************/

void complex_lower(double complex N[NDIM][NDIM],
    double gcov[NDIM][NDIM],
    int low1, int low2, double complex Nl[NDIM][NDIM])
{
  int i, j, k, l;

  if (!low1 && !low2)
  return;

  if (low1 && low2) {
    for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) {
      Nl[i][j] = 0. + I * 0.;
      for (k = 0; k < 4; k++)
      for (l = 0; l < 4; l++) {
        Nl[i][j] += N[k][l] * gcov[k][i] * gcov[l][j];
      }
    }
    return;
  }

  if (low1) {
    for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) {
      Nl[i][j] = 0. + I * 0.;
      for (k = 0; k < 4; k++)
      Nl[i][j] += N[k][j] * gcov[k][i];
    }
    return;
  }

  for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++) {
    Nl[i][j] = 0. + I * 0.;
    for (l = 0; l < 4; l++)
    Nl[i][j] += N[i][l] * gcov[l][j];
  }
  return;

}

void stokes_to_tensor(double fI, double fQ, double fU, double fV,
    double complex f_tetrad[NDIM][NDIM])
{
  /*notice that I swapped sign of the imaginary part [2][3] in [3][2] - which one is correct? */
  f_tetrad[1][1] = (fI + fQ + 0. * I);
  f_tetrad[1][2] = (fU - I * fV);
  f_tetrad[2][1] = (fU + I * fV);
  f_tetrad[2][2] = (fI - fQ + 0. * I);

}

void tensor_to_stokes(double complex f_tetrad[NDIM][NDIM],
    double *fI, double *fQ, double *fU, double *fV)
{

  /*here I divide by two to agree with above */
  *fI = creal(f_tetrad[1][1] + f_tetrad[2][2]) / 2;
  *fQ = creal(f_tetrad[1][1] - f_tetrad[2][2]) / 2;
  *fU = creal(f_tetrad[1][2] + f_tetrad[2][1]) / 2;
  *fV = cimag(f_tetrad[2][1] - f_tetrad[1][2]) / 2;

}

void complex_coord_to_tetrad_rank2(double complex T_coord[NDIM][NDIM],
    double Ecov[NDIM][NDIM],
    double complex T_tetrad[NDIM][NDIM])
{
  int i, j, k, l;

  for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++)
  T_tetrad[i][j] = 0. + I * 0.;

  for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++)
  for (k = 0; k < 4; k++)
  for (l = 0; l < 4; l++)
  T_tetrad[i][j] +=
  T_coord[k][l] * Ecov[i][k] * Ecov[j][l];

  return;
}

void complex_tetrad_to_coord_rank2(double complex T_tetrad[NDIM][NDIM],
    double Econ[NDIM][NDIM],
    double complex T_coord[NDIM][NDIM])
{
  int i, j, k, l;

  for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++)
  T_coord[i][j] = 0. + I * 0.;

  for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++)
  for (k = 0; k < 4; k++)
  for (l = 0; l < 4; l++)
  T_coord[i][j] +=
  T_tetrad[k][l] * Econ[k][i] * Econ[l][j];

  return;
}
