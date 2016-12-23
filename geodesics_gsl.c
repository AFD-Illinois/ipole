
#include "decs.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

int geodesic_deriv(double t, const double y[], double dy[], void *params)
{
	double X[NDIM],Kcon[NDIM],lconn[NDIM][NDIM][NDIM] ;
	int i,j,k ;

	for(i=0;i<NDIM;i++) X[i] = y[i] ;
	for(i=0;i<NDIM;i++) Kcon[i] = y[i+NDIM] ;

	get_connection(X, lconn) ;
	
	for(i=0;i<NDIM;i++) dy[i] = Kcon[i] ;

	for (k = 0; k < 4; k++) dy[k+NDIM] = 0. ;

	for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	for (k = 0; k < 4; k++)
		dy[k+NDIM] += -lconn[k][i][j] * Kcon[i] * Kcon[j];

	return(GSL_SUCCESS) ;
}

/*

this is the main photon orbit integrator 

*/

#define EPSACC 1.e-5

void push_photon_gsl(double X[NDIM], double Kcon[NDIM], double dl)
{

	int i ;
	static int firstc = 1 ;
	double XK[2*NDIM] ;
	double l ;
	double lf ;
	double dl0 ;
	static gsl_odeiv_step    *s ;
	static gsl_odeiv_control *c ;
	static gsl_odeiv_evolve  *e ;
	gsl_odeiv_system sys = {geodesic_deriv, NULL, 2*NDIM, NULL} ;

	const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd ;

	if(firstc) {
		firstc = 0 ;
		s = gsl_odeiv_step_alloc(T, 2*NDIM) ;
		c = gsl_odeiv_control_y_new (EPSACC, EPSACC) ;
		e = gsl_odeiv_evolve_alloc (2*NDIM) ;
	}

	/* set initial state */
	for(i=0;i<NDIM;i++) XK[i] = X[i] ;
	for(i=0;i<NDIM;i++) XK[i+NDIM] = Kcon[i] ;

	l = 0. ;
	lf = dl ;
	dl0 = dl ;
	while(fabs(l - lf) > 0.) {
		int status = gsl_odeiv_evolve_apply(e,c,s,&sys,&l,lf,&dl0, XK) ;


		if(status != GSL_SUCCESS) {
			fprintf(stderr,"integrator failure\n") ;
			exit(100) ;
		}

		//fprintf(stderr,"%g %g %g\n",l, lf, dl0) ;
		//fprintf(stderr,"%g %g\n",XK[1],XK[2]) ;
	}

	for(i=0;i<NDIM;i++) X[i] = XK[i] ;
	for(i=0;i<NDIM;i++) Kcon[i] = XK[i+NDIM] ;

	//fprintf(stderr,"%g %g %g\n",X[1],X[2],X[3]) ;

	/* done! */
}
