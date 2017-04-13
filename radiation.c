
/* 

model-independent
radiation-related utilities.

*/

#include "decs.h"



double Bnu_inv(double nu, double Thetae)
{

	double x;

	x = HPL*nu/(ME*CL*CL*Thetae);

	//original
	
	//	if(x < 1.e-3)	/* Taylor expand */
	//  return ((2.*HPL/(CL*CL))/( x/24. * (24. + x*(12. + x*(4. + x)))));
	//else
	//  return((2.*HPL/(CL*CL))/(exp(x) - 1.)) ;
	
	//same as jasons
	if(x < 1.e-6)	
	  return ((2.*HPL/(CL*CL))/x);
	else
	  return((2.*HPL/(CL*CL))/(exp(x) - 1.)) ;//bnu/nu^3
	
	/*grtrans
	where(h*nu/k/T.lt.1d-6)
         bnu=2d0*nu*nu*k*T/c2
      elsewhere
	  bnu = 2d0*h*nu/c2*nu*nu/(exp(h*nu/k/T)-1d0)
      endwhere
	*/



}


/*------------------------------------------------------------------------------------------*/
double jnu_inv(double nu, double Thetae, double Ne, double B, double theta, double Be)
{
  double j ;
#if(SYN)
  j = jnu_synch(nu,Ne,Thetae,B,theta,Be) ;
#endif

#if(FF)
  j = jnu_ff(nu,Ne,Thetae,B,theta) ;
#endif

  return(j/(nu*nu)) ;
}
double anu_inv(double nu, double Thetae, double Ne, double B, double theta,double Be)
{
  double a ;
  a = anu_synch(nu,Ne,Thetae,B,theta,Be) ;
  return(a*nu) ;

}
/*------------------------------------------------------------------------------------------*/

/*******************************************************************************************/
double jnu_inv_nth(double nu, double Thetae, double Ne, double B, double theta, double Be)
{
  double j ;
  j = jnu_synch_nth(nu,Ne,Thetae,B,theta,Be) ;
  return(j/(nu*nu)) ;
}

double anu_inv_nth(double nu, double Thetae, double Ne, double B, double theta, double Be)
{
  double a ;
  a = anu_synch_nth(nu,Ne,Thetae,B,theta,Be) ;
  return(a*nu) ;
}
/*******************************************************************************************/


/* get frequency in fluid frame, in Hz */
double get_fluid_nu(double Kcon[NDIM], double Ucov[NDIM])
{
  //	int i;
	double ener, mec2_h;
	double nu ;

	mec2_h = ME*CL*CL/HPL;

	/* this is the energy in electron rest-mass units */
	ener = -( Kcon[0] * Ucov[0] + 
		  Kcon[1] * Ucov[1] + 
		  Kcon[2] * Ucov[2] + 
		  Kcon[3] * Ucov[3] );


	nu = ener;//*mec2_h;

	if(nu < 0.) {
	  nu = 1. ;
	  //		fprintf(stdout,"nu=1 get_fluid_nu, K: %g %g %g %g\n", Kcon[0],Kcon[1],Kcon[2],Kcon[3]) ;
	  //	fprintf(stdout,"nu=1 get_fluid_nu, U: %g %g %g %g\n", Ucov[0],Ucov[1],Ucov[2],Ucov[3]) ;
	}
	//fprintf(stderr,"nu = %g\n", nu) ;

	if(isnan(ener)) {
		fprintf(stderr,"isnan get_fluid_nu, K: %g %g %g %g\n",
		Kcon[0],Kcon[1],Kcon[2],Kcon[3]) ;
		//fprintf(stderr,"isnan get_fluid_nu, X: %g %g %g %g\n", X[0],X[1],X[2],X[3]) ;
		fprintf(stderr,"isnan get_fluid_nu, U: %g %g %g %g\n", Ucov[0],Ucov[1],Ucov[2],Ucov[3]) ;
		exit(1);
	}

	return(nu) ;
}

/* return angle between magnetic field and wavevector */
double get_bk_angle(double X[NDIM], double Kcon[NDIM], double Ucov[NDIM])
{
	double Bcon[NDIM],Bcov[NDIM] ;
	double B,k,mu ;

	get_model_bcov(X, Bcov) ;
	get_model_bcon(X, Bcon) ;

	B = sqrt(fabs(Bcon[0]*Bcov[0] + Bcon[1]*Bcov[1] + Bcon[2]*Bcov[2] + Bcon[3]*Bcov[3])) ;

	if(B == 0.) return(M_PI/2.) ;

	k = fabs(Kcon[0]*Ucov[0] + Kcon[1]*Ucov[1] + Kcon[2]*Ucov[2] + Kcon[3]*Ucov[3]) ;

	mu = ( Kcon[0]*Bcov[0] + Kcon[1]*Bcov[1] + Kcon[2]*Bcov[2] + Kcon[3]*Bcov[3])/(k*B) ;

	if(fabs(mu) > 1. ) mu /= fabs(mu) ;

	if(isnan(mu)) fprintf(stderr,"isnan get_bk_angle\n") ;

 	return( acos(mu) );

}

