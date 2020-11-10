/**********************************************************/
/*** all you need to make a polarized radiative transfer***/
/***** used in ipole to evolve complex tensor N ***********/
/***** along with standard evolution for I scalar *********/
/**********************************************************/
/**** written by Monika Moscibrodzka on 09 July 2014 ******/
/************ @ Eindhoven Airport *************************/
/****************and then rewritten by C.Gammie ***********/
/**********************************************************/

#include "decs.h"
#include "coordinates.h"
#include "geometry.h"
#include "model.h"
#include "model_radiation.h"
#include "model_tetrads.h"
#include "debug_tools.h"
#include <complex.h>

/* tensor tools*/
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

/***************************MAIN FUNCTIONS******************************/
/* initialize tensor N in the coordinate frame at the beginning of the *
 geodesics integration = it is zero */
void init_N(double X[NDIM], double Kcon[NDIM],
    double complex N_coord[NDIM][NDIM])
{
  MUNULOOP N_coord[mu][nu] = 0.0 + I * 0.0;
  return;
}

/*

 parallel transport N over dl

 */
void push_polar(double Xi[NDIM], double Xm[NDIM], double Xf[NDIM],
    double Ki[NDIM], double Km[NDIM], double Kf[NDIM],
    complex double Ni[NDIM][NDIM],
    complex double Nm[NDIM][NDIM],
    complex double Nf[NDIM][NDIM], double dl)
{

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
  ) * dl / (L_unit * HPL / (ME * CL * CL));

  return;
}

/* updates N for one step on geodesics, using the previous step N*/
/* here we compute new right-hand side of the equation */
/* and somehow rotate this along the geodesics knowing */
/* first point and last point X and K*/

void evolve_N(double Xi[NDIM], double Kconi[NDIM],
    double Xhalf[NDIM], double Kconhalf[NDIM],
    double Xf[NDIM], double Kconf[NDIM],
    double dlam, double complex N_coord[NDIM][NDIM], double *tauF)
{
  double gcov[NDIM][NDIM];
  double Ucon[NDIM],Bcon[NDIM];
  double Ucov[NDIM],Bcov[NDIM];
  double Ecov[NDIM][NDIM], Econ[NDIM][NDIM];
  double complex Nh[NDIM][NDIM];
  double complex N_tetrad[NDIM][NDIM];
  double B;
  double jI, jQ, jU, jV;
  double aI, aQ, aU, aV;
  double rV, rU, rQ;
  double rho2, rho, rdS;
  double SI, SQ, SU, SV;
  double SI0, SQ0, SU0, SV0;
  double SI1, SQ1, SU1, SV1;
  double SI2, SQ2, SU2, SV2;

  /* parallel transport N by a half, and then full, step */
  push_polar(Xi, Xi, Xhalf, Kconi, Kconi, Kconhalf, N_coord, N_coord, Nh, 0.5 * dlam);
  push_polar(Xi, Xhalf, Xf, Kconi, Kconhalf, Kconf, N_coord, Nh, N_coord, dlam);

  /* absorption/emission/rotation step.  only complete if radiating_region condition is satisfied */
  if ( radiating_region(Xf) ) {

    // get fluid parameters at Xf
    get_model_fourv(Xf, Ucon, Ucov, Bcon, Bcov);

    /* evaluate transport coefficients */
    gcov_func(Xf, gcov);
    jar_calc(Xf, Kconf, &jI, &jQ, &jU, &jV,
        &aI, &aQ, &aU, &aV, &rQ, &rU, &rV);

    if (counterjet == 1) { // Emission from X[2] > 0.5 only
      if (Xf[2] < 0.5) {
        jI = jQ = jU = jV = 0.;
      }
    } else if (counterjet == 2) { // Emission from X[2] < 0.5 only
      if (Xf[2] > 0.5) {
        jI = jQ = jU = jV = 0.;
      }
    }

    //ARR:  Modify according to positron ratio
    jI *= (1 + positronRatio);
    jQ *= (1 + positronRatio);
    jU *= (1 + positronRatio);
    jV *= (1 - positronRatio);
    aI *= (1 + positronRatio);
    aQ *= (1 + positronRatio);
    aU *= (1 + positronRatio);
    aV *= (1 - positronRatio);
    rQ *= (1 + positronRatio);
    rU *= (1 + positronRatio);
    rV *= (1 - positronRatio);

    /* make plasma tetrad */
    B = get_model_b(Xf); /* field in G */
    if (B < 0.) {
      Bcon[0] = 0.;
      Bcon[1] = 1.;
      Bcon[2] = 1.;
      Bcon[3] = 1.;
    }

    make_plasma_tetrad(Ucon, Kconf, Bcon, gcov, Econ, Ecov);

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

    *tauF += dlam*fabs(rV);
    if (*tauF > 1.e100 || *tauF < -1.e100 || isnan(*tauF)) {
      printf("tauF = %e dlam = %e rV = %e\n", *tauF, dlam, rV);
      exit(-1);
    }

    if (aP > SMALL) { /* full analytic solution has trouble if polarized absorptivity is small */
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

    } else { /* this should really be a series expansion in aP */
      SI2 = SI1 + x * jI;
      SQ2 = SQ1 + x * jQ;
      SU2 = SU1 + x * jU;
      SV2 = SV1 + x * jV;
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

    /* re-pack the Stokes parameters into N */
    stokes_to_tensor(SI, SQ, SU, SV, N_tetrad);
    complex_tetrad_to_coord_rank2(N_tetrad, Econ, N_coord);
    /*if (isnan(creal(N_tetrad[0][0])) || isnan(creal(N_coord[0][0]))) {
     printf("N_tet = %e N_coord = %e Econ = %e\n", creal(N_tetrad[0][0]),
     creal(N_coord[0][0]), Econ[0][0]);
     MUNULOOP printf("Econ[%i][%i] = %e Ntet = %e\n", mu,nu,Econ[mu][nu],
     creal(N_tetrad[mu][nu]));
     }*/

  }

#if THIN_DISK
  // The thin disk problem emits nowhere but uses a boundary condition region defined by thindisk_region
  // In this region, it sets all Stokes parameters according to the function get_model_stokes

  //check_N(N[NDIM][NDIM], Kcon[NDIM], gcov[NDIM][NDIM]);

  // If we're in exactly the thin disk...
  if (thindisk_region(Xi, Xf)) {
    // get fluid parameters at Xf -- B is set to grtrans' "polarization direction"
    double gcov[NDIM][NDIM];
    gcov_func(Xf, gcov);
    get_model_fourv_K(Xf, Kconf, Ucon, Ucov, Bcon, Bcov);

    make_plasma_tetrad(Ucon, Kconf, Bcon, gcov, Econ, Ecov);
    // Otherwise N_tetrad is uninitialized memory.  Probably fine but meh
    complex_coord_to_tetrad_rank2(N_coord, Ecov, N_tetrad);
    tensor_to_stokes(N_tetrad, &SI0, &SQ0, &SU0, &SV0);

    get_model_stokes(Xf, Kconf, &SI, &SQ, &SU, &SV);

    stokes_to_tensor(SI, SQ, SU, SV, N_tetrad);
    complex_tetrad_to_coord_rank2(N_tetrad, Econ, N_coord);
  }

  //check_N(N[NDIM][NDIM], Kcon[NDIM], gcov[NDIM][NDIM]);
#endif

  /* SOURCE STEP DONE */

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

  return;

}

/***************************END MAIN FUNCTIONS******************************/

/*************************SUPPORTING FUNCTIONS******************************/

/*
 *
 */
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
  int i, j;

  for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++)
  f_tetrad[i][j] = 0. + I * 0.;
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
