
/**********************************************************/

/*** all you need to make a polarized radiative transfer***/
/***** used in ipole to evolve complex tensor N *******/
/***** along with standard evolution for I scalar *********/
/**********************************************************/
/**** written by Monika Moscibrodzka on 09 July 2014 ******/
/************ @ Eindhoven Airport *************************/
/************  last update: 9 May 2017   ******************/

/****************and then rewritten by C.Gammie ***********/

/**********************************************************/

#include "decs.h"
#include "defs.h"

/* the following definitions are used only locally */
#define MNLOOP   for(m=0;m<NDIM;m++)for(n=0;n<NDIM;n++)

/* transfer coefficients in tetrad frame */
void jar_calc(double X[NDIM], double Kcon[NDIM],
	      double *jI, double *jQ, double *jU, double *jV,
	      double *aI, double *aQ, double *aU, double *aV,
	      double *rQ, double *rU, double *rV);

/* tensor tools*/
void check_N(double complex N[NDIM][NDIM], double Kcon[NDIM],
	     double gcov[NDIM][NDIM]);
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
/* initialize tensor N in the coordinate frame at the bening of the *
geodesics integration = it is zero */
void init_N(double X[NDIM], double Kcon[NDIM],
	    double complex N_coord[NDIM][NDIM])
{
    int m, n;

    MNLOOP N_coord[m][n] = 0.0 + I * 0.0;

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
    int k;
    double gcov[NDIM][NDIM];
    double Ucon[NDIM],Bcon[NDIM];
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
    int radiating_region(double X[4]);

    /* parallel transport N by a half, and then full, step */
    push_polar(Xi, Xi, Xhalf, Kconi, Kconi, Kconhalf, N_coord, N_coord, Nh, 0.5 * dlam);
    push_polar(Xi, Xhalf, Xf, Kconi, Kconhalf, Kconf, N_coord, Nh, N_coord, dlam);

    /* absorption/emission/rotation step.  only complete if radiating_region condition is satisfied */
    if ( radiating_region(Xf) ) {

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

	/* make plasma tetrad */
	get_model_ucon(Xf, Ucon);
	B = get_model_b(Xf);	/* field in G */
	if (B > 0.) {
	    get_model_bcon(Xf, Bcon);
	}
	else {
	    Bcon[0] = 0.;
	    for (k = 1; k < NDIM; k++)
		Bcon[k] = 1.;
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

	if (aP > SMALL) {  /* full analytic solution has trouble if polarized absorptivity is small */
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

	} else {  /* this should really be a series expansion in aP */
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

    }

    /* SOURCE STEP DONE */

}


/* converts tensor N to Stokes parameters detected at the camera*/
void project_N(double X[NDIM], double Kcon[NDIM],
	       double complex N_coord[NDIM][NDIM], double *Stokes_I,
	       double *Stokes_Q, double *Stokes_U, double *Stokes_V)
{
    double complex N_tetrad[NDIM][NDIM];
    double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];

    make_camera_tetrad(X, Econ, Ecov);

    complex_coord_to_tetrad_rank2(N_coord, Ecov, N_tetrad);

    tensor_to_stokes(N_tetrad, Stokes_I, Stokes_Q, Stokes_U, Stokes_V);

    return;

}

/***************************END MAIN FUNCTIONS******************************/


/*************************SUPPORTING FUNCTIONS******************************/

/*

    call this function if you want to check
    that the coherency tensor N satisfies certain
    basic properties:
    k . N = N . k = 0
    hermitian
    evaluate the invariants: I, Q^2 + U^2, V^2

*/

void check_N(double complex N[NDIM][NDIM],
	     double Kcon[NDIM], double gcov[NDIM][NDIM])
{
    double complex dot;
    double Kcov[NDIM];
    int i, j;

    fprintf(stderr, "enter check_N\n");

    /* compute k . N */
    lower(Kcon, gcov, Kcov);
    fprintf(stderr, "(k . N\n");
    /* first one way */
    for (i = 0; i < 4; i++) {
	dot = 0. + I * 0.;
	for (j = 0; j < 4; j++)
	    dot += Kcov[j] * N[j][i];
	fprintf(stderr, "%d %g + i %g\n", i, creal(dot), cimag(dot));
    }
    /* then the other */
    for (i = 0; i < 4; i++) {
	dot = 0. + I * 0.;
	for (j = 0; j < 4; j++)
	    dot += Kcov[j] * N[i][j];
	fprintf(stderr, "%d %g + i %g\n", i, creal(dot), cimag(dot));
    }
    fprintf(stderr, "k . N)\n");

    /* check for hermiticity */
    fprintf(stderr, "(herm:\n");
    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    fprintf(stderr, "%d %d %g + i %g\n", i, j,
		    creal(N[i][j] - conj(N[j][i])),
		    cimag(N[i][j] - conj(N[j][i]))
		);
    fprintf(stderr, "herm)\n");

    /* check invariants */
    double complex Nud[NDIM][NDIM];
    void complex_lower(double complex N[NDIM][NDIM],
		       double gcov[NDIM][NDIM], int low1, int low2,
		       double complex Nl[NDIM][NDIM]);
    complex_lower(N, gcov, 0, 1, Nud);
    for (i = 0; i < 4; i++)
	fprintf(stderr, "N: %d %g + i %g\n", i, creal(N[i][i]),
		cimag(N[i][i]));
    for (i = 0; i < 4; i++)
	fprintf(stderr, "Nud: %d %g + i %g\n", i, creal(Nud[i][i]),
		cimag(Nud[i][i]));
    dot = 0. + I * 0.;
    for (i = 0; i < 4; i++)
	dot += Nud[i][i];
    fprintf(stderr, "I: %g + i %g\n", creal(dot), cimag(dot));

    double complex Ndd[NDIM][NDIM];
    complex_lower(N, gcov, 1, 1, Ndd);
    dot = 0. + I * 0.;
    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    dot +=
		2. * 0.25 * (N[i][j] + N[j][i]) * (Ndd[i][j] + Ndd[j][i]);
    fprintf(stderr, "IQUsq: %g + i %g\n", creal(dot), cimag(dot));

    dot = 0. + I * 0.;
    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    dot +=
		-2. * 0.25 * (N[i][j] - N[j][i]) * (Ndd[i][j] - Ndd[j][i]);
    fprintf(stderr, "Vsqsq: %g + i %g\n", creal(dot), cimag(dot));

    fprintf(stderr, "leave check_N\n");
}

void complex_lower(double complex N[NDIM][NDIM],
		   double gcov[NDIM][NDIM],
		   int low1, int low2, double complex Nl[NDIM][NDIM]
    )
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

#undef MNLOOP
