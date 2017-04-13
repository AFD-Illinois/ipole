
#include "decs.h"

/*  needed for bremsstrahlung  */
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf.h>
/******************************/


/* 
synchrotron "mixed" emissivity formula 
interpolates between Petrosian limit and
classical thermal synchrotron limit
good for Thetae > 1
*/
double jnu_synch(double nu, double Ne, double Thetae, double B, double theta,double Ber)
{
        double K2,nuc,nus,x,f,j,sth ;

	if(nu==1.) return 0;
	if(Thetae<0.1) return 0;
	if(Thetae<5.0)	 K2 = gsl_sf_bessel_Kn(2,1./Thetae) ;
	else K2 = 2.*Thetae*Thetae ;
	
        nuc = EE*B/(2.*M_PI*ME*CL) ;

	
	sth = sin(theta) ; 
        nus = (2./9.)*nuc*Thetae*Thetae*sth ;
	if(nu > 1.e12*nus) return(0.) ;
        x = nu/nus ;
        f = pow( pow(x,1./2.) + pow(2.,11./12.)*pow(x,1./6.), 2 ) ;
	j = (sqrt(2.)*M_PI*EE*EE*Ne*nus/(3.*CL*K2)) * f * exp(-pow(x,1./3.)) ;
	
	//Jasons emissivity Stokes I
	x=nu/Thetae/Thetae/(nuc*1.5*sth);                                                                                                                                     
	j=Ne*EE*EE*nu/2./1.7/CL/Thetae/Thetae*2.5651*(1+1.92*pow(x,-1./3.)+0.9977*pow(x,-2./3.))*exp(-1.8899*pow(x,1./3.)); 
	
	
	//insert here the kappa function to simulate particle accelration k=3.5 w=theta
	/*
	double w=Thetae;
	double kappa=Ber; 
	double nu_w = pow(w * kappa, 2.) *nuc *sin(theta);
	double X_k = nu/nu_w;
	double prefactor = (Ne * pow(EE, 2.) * nuc * sin(theta)) /CL;
	
	double Nlow = 4. * M_PI * tgamma(kappa-4./3.) / (pow(3., 7./3.) * tgamma(kappa-2.));
	double Nhigh = (1./4.) * pow(3., (kappa-1.)/2.)
	  * (kappa-2.) * (kappa-1.)
	  * tgamma(kappa/4.-1./3.)
	  * tgamma(kappa/4.+4./3.);
	x = 3. * pow(kappa, -3./2.);
	double ans = prefactor * Nlow * pow(X_k, 1./3.)
         	   * pow(1.+pow(X_k, x * (3. * kappa-4.)/6.)
			 * pow(Nlow/Nhigh, x), -1./x);
	
	j=ans;
	*/	
	
	//thermal again
	/*
	Thetae=10.;
	double nu_s = (2./9.)*nuc*sin(theta)*Thetae*Thetae;
	double X = nu/nu_s;
	double prefactor = (Ne
			    * pow(EE, 2.)
			    * nuc)/CL;

	double term1 = sqrt(2.)*M_PI/27. * sin(theta);
	double term2 = pow(pow(X, 0.5)+pow(2., 11./12.)*pow(X, 1./6.), 2.);
	double term3 = exp(-pow(X, 1./3.));
	double ans = prefactor * term1 * term2 * term3;
	j=ans;
	*/

        return(j) ;
}

//synchrotron absorptivity
double anu_synch(double nu, double Ne, double Thetae, double B, double theta, double Ber)
{
  double  nuc = EE*B/(2.*M_PI*ME*CL) ;

  double w=Thetae;
  double kappa=Ber;

  //insert here the kappa function to simulate particle accelration
  double nu_w = pow(w * kappa, 2.) * nuc * sin(theta);

  double X_k = nu/nu_w;

  double prefactor = Ne * EE / (B * sin(theta));

  double a = kappa - 1./3.;
  double b = kappa + 1.;
  double c = kappa + 2./3.;
  double z = -kappa*w;
  /*GSL 2F1 only works for |z| < 1; had to apply a hypergeometric function                                                                                                                      
    identity because in our case z = -kappa*w, so |z| > 1 */
  double hyp2f1 = pow(1.-z, -a) * tgamma(c) * tgamma(b-a)
    / (tgamma(b)*tgamma(c-a))
    * gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))
    + pow(1.-z, -b) * tgamma(c) * tgamma(a-b)
    / (tgamma(a) * tgamma(c-b))
    * gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));

  double Nlow = pow(3., 1./6.) * (10./41.) * pow(2. * M_PI, 2.)
    / pow(w * kappa, 16./3.-kappa)
    * (kappa-2.) * (kappa-1.) * kappa
    / (3.*kappa-1.) * tgamma(5./3.) * hyp2f1;

  double Nhigh = 2. * pow(M_PI, 5./2.)/3. * (kappa-2.)
    * (kappa-1.) * kappa
    / pow(w * kappa, 5.)
    * (2 * tgamma(2. + kappa/2.)
       / (2.+kappa)-1.)
    * (pow(3./kappa, 19./4.) + 3./5.);

  double x = pow(-7./4. + 8. * kappa/5., -43./50.);

  double ans = prefactor * Nlow * pow(X_k, -5./3.)
    * pow(1. + pow(X_k, x * (3. * kappa-1.)/6.)
	  * pow(Nlow/Nhigh, x), -1./x);

	
  return(ans) ;
}




struct my_f_params {double Thetae; double w; };
/*bremsstrahlung integrant for e-p interactions for Maxwellian distribution*/
/* see below */
double f1 (double gamma, void * p){

  struct my_f_params * params = (struct my_f_params *)p;
  double thetae1 = (params->Thetae);
  double w1 = (params->w);
  double beta,beta2,gamma2;
  double K2;
  double w_dsigma_dw(double w1,double gamma);

  K2 = gsl_sf_bessel_Kn(2,1./thetae1) ;
  gamma2=gamma*gamma;
  beta=sqrt(gamma2-1.)/gamma;
  beta2=beta*beta;
  double f1=w_dsigma_dw(w1,gamma)*beta2*gamma2*exp(-gamma/thetae1)/thetae1/K2;
  return f1;

}


/* isotropic free-free transition emissivity function */ 
/* uses GSL library to integrate (Gauss quadratures) emisivity for e-i collisions*/
/* the emissivity is isotropic and does not depend on B field */
double jnu_ff(double nu, double Ne, double Thetae, double B, double theta)
{
	
        double j;
        double gamma_min,gamma_max;
        double integral_result, error;
        double j_nu_ee,j_nu_ei,omega;
        gsl_integration_workspace * w;


	/* emissivity from e-p for all temperatures */
	omega=nu*HPL_MECL2; //photon energy in electron rest-mass units
	
	//parameters that go into function that we need to integrate over gammas
	struct my_f_params alpha = {Thetae, omega};
	gsl_function F1;
	F1.function = &f1;
	F1.params = &alpha;

	gamma_min=1.+omega ;
	gamma_max=1e7; // thoretically integration to infinity 
	
	/* Gaussian quadratures */
	//w  = gsl_integration_workspace_alloc(1000);
	//gsl_integration_qag(&F1, gamma_min, gamma_max, 0, 1e-7, 1000, GSL_INTEG_GAUSS61, w, &integral_result, &error);
	//gsl_integration_workspace_free(w);
	//j_nu_ei=integral_result*Ne*Ne*CL*HPL/4./M_PI ;
	//we assume gaunt=1;
	//this is more or less ok for rel and ultrarel plasma
	j_nu_ei=6.8e-38*Ne*Ne/sqrt(Thetae*ME*CL*CL/KBOL)*exp(-omega/Thetae);


	/* emissivity from e-e for all temperatures */
	//this process domiantes for thetae>1 anyway
	j_nu_ee=jnu_ff_ee(nu, Ne, Thetae);

	/* total emissivity ee and ei*/
	j=j_nu_ei+j_nu_ee ;

        return(j) ;
}


/* synchortron emissivity and absorptivity for power-law electrons */
/* function needs precomputed tables, that are calculates by a separate program */
/* the tables should be precomputed at the begining of ray tracer*/
double jnu_synch_nth(double nu, double Ne, double Thetae, double B, double theta, double Be)
{

  int ind;
  double j,nuc,sth,nee,nub ;
  double gmax,gmin,p;
  double lx_min,lx_max,dlnx;
  double xmin,xmax,Gfunmin,Gfunmax;    
  double b2,b4,b6,b8;
  b2=Be;
  b4=b2*b2;
  b6=b2*b4;
  b8=b4*b4;

  //Ne=1.e6;
  //B=10;
  //Thetae=5.;

  p=3.5;
  if(Thetae>1./4.){gmin=4.*Thetae;}
  else{gmin=1.0;}
  gmax=1.e8;
  delta_acc=0.2;
  nee=delta_acc*Ne;


  if(theta <= 0. || theta >= M_PI) return 0.0;
  sth=sin(theta);
  j=0.0;

  if(nee>0.0){


    nub=EE*B/2.0/M_PI/ME/CL;
    nuc=1.5*nub*sth;

    lx_min=log(1.e-10);
    lx_max=log(702.);
    dlnx=(lx_max-lx_min)/(201.);
    xmin=nu/(nuc*gmin*gmin);
    xmax=nu/(nuc*gmax*gmax);
      
    if(xmin <= exp(xtab[0])) Gfunmin=exp(Gtab[0]);
    if(xmin >= exp(xtab[200]) ) Gfunmin=0.0;
    if( xmin > exp(xtab[0]) && xmin < exp(xtab[200])){
      ind=(int) ((log(xmin)-lx_min)/dlnx);
      Gfunmin=exp(Gtab[ind]+(Gtab[ind+1]-Gtab[ind])*(log(xmin)-xtab[ind])/(xtab[ind+1]-xtab[ind]));
    }
      
    if(xmax <= exp(xtab[0])) Gfunmax=exp(Gtab[0]);
    if(xmax >= exp(xtab[200])) Gfunmax=0.0;
    if( xmax > exp(xtab[0]) && xmax < exp(xtab[200])){
      ind=(int) ((log(xmax)-lx_min)/dlnx);
      Gfunmax=exp(Gtab[ind]+(Gtab[ind+1]-Gtab[ind])*(log(xmax)-xtab[ind])/(xtab[ind+1]-xtab[ind]));
    }
    
    j=nee*EE*EE*(p-1.)*nuc/2./sqrt(3.)/CL/(pow(gmin,1.-p)-pow(gmax,1.-p))*pow(nu/nuc,0.5-0.5*p)*(Gfunmax-Gfunmin);
  }

  return(j) ;
}

/*absorption for non-thermal electrons*/
double anu_synch_nth(double nu, double Ne, double Thetae, double B, double theta, double Be)
{
  int ind;
  double aaa,nuc,sth,nee,nub ;
  double gmax,gmin,p;
  double lx_min,lx_max,dlnx;
  double xmin,xmax,Gafunmin,Gafunmax;    
  double b2,b4,b6,b8;
  b2=Be*Be;
  b4=b2*b2;
  b6=b2*b4;
  b8=b4*b4;

  p=3.5;
  if(Thetae>1./4.){gmin=4.*Thetae;}
  else{gmin=1.0;}
  gmax=1.e8;
  delta_acc=0.2;
  nee=delta_acc*Ne;


  if(theta <= 0. || theta >= M_PI) return 0.0;
  sth=sin(theta);
  aaa=0.0;

  
  if(nee>0.0){
    nub=EE*B/2.0/M_PI/ME/CL;
    nuc=1.5*nub*sth;
    
    lx_min=log(1.e-10);
    lx_max=log(702.);
    dlnx=(lx_max-lx_min)/(201.);
    xmin=nu/(nuc*gmin*gmin);
    xmax=nu/(nuc*gmax*gmax);
    
    if(xmin <= exp(xatab[0])) Gafunmin=exp(Gatab[0]);
    if(xmin >= exp(xatab[200])) Gafunmin=0.0;
    if(xmin > exp(xatab[0]) && xmin < exp(xatab[200])){
      ind=(int) ((log(xmin)-lx_min)/dlnx);
      Gafunmin=exp(Gatab[ind]+(Gatab[ind+1]-Gatab[ind])*(log(xmin)-xatab[ind])/(xatab[ind+1]-xatab[ind]));
    }
    
    
    if(xmax <= exp(xatab[0])) Gafunmax=exp(Gatab[0]);
    if(xmax >= exp(xatab[200])) Gafunmax=0.0;
    if(xmax > exp(xatab[0]) && xmax < exp(xatab[200])){
      ind=(int) ((log(xmax)-lx_min)/dlnx);
      Gafunmax=exp(Gatab[ind]+(Gatab[ind+1]-Gatab[ind])*(log(xmax)-xatab[ind])/(xatab[ind+1]-xatab[ind]));
    }
    
    aaa=nee*EE*EE*(p-1.)*(p+2)/(4.*sqrt(3.)*ME*CL*nuc)/(pow(gmin,1.-p)-pow(gmax,1.-p)) 
      *pow(nu/nuc,-p/2.-2.)*(Gafunmax-Gafunmin);
    
  }

  return(aaa) ;
}

