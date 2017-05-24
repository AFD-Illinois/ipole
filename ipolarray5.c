/*** all you need to make a polarized radiative transfer***/
/***** used in mibothros to evolve complex tensor N *******/
/***** along with standard evolution for I scalar *********/
/**********************************************************/
/**** written by Monika Moscibrodzka on 09 July 2014 ******/
/************ @ Eindhoven Airport *************************/
/************  last update: 9 May 2017   ******************/

/****************and then rewritten by C.Gammie ***********/

/**********************************************************/
/******** compile with: gcc -lm -lgsl -lgslcblas **********/
/**********************************************************/
/**********************************************************/
#include "decs.h"
#include "defs.h"

/* the following definitions are used only locally */
#define S2     (1.4142135637) //sqrt(2)
#define S3     (1.73205080757) //sqrt(3)
#define ILOOP  for(i=0;i<NDIM;i++)
#define JLOOP  for(j=0;j<NDIM;j++)
#define KLOOP  for(k=0;k<NDIM;k++)
#define LLOOP  for(l=0;l<NDIM;l++)
#define IJLOOP for(i=0;i<NDIM;i++)for(j=0;j<NDIM;j++)
#define IJKLLOOP for(i=0;i<NDIM;i++)for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)for(l=0;l<NDIM;l++)
#define MNLOOP   for(m=0;m<NDIM;m++)for(n=0;n<NDIM;n++)
#define MNKLLOOP for(m=0;m<NDIM;m++)for(n=0;n<NDIM;n++)for(k=0;k<NDIM;k++)for(l=0;l<NDIM;l++)
#define KLLOOP   for(k=0;k<NDIM;k++)for(l=0;l<NDIM;l++)
#define sign(x) (((x) < 0) ? -1 : ((x) > 0))

/* declatations of functions that are used locally defined functions , not visible from outside*/
/*thermal plasma emissivity, absorptivity and Faraday conversion and rotation */
double g(double Xe);
double h(double Xe);
double Je(double Xe);
double jffunc(double Xe);
double I_I(double x);
double I_Q(double x);
double I_V(double x);
double Bnu(double nu, double Thetae);
double besselk_asym(int n, double x);

/*complex tensor in tetrad frame */
void jar_calc(const double X[NDIM],const double Kcon[NDIM], 
		    double *jI,double *jQ,double *jU,double *jV,
		    double *aI,double *aQ,double *aU,double *aV,
	      double *rQ,double *rU,double *rV,
	      double Op[NDIM][NDIM],double xx);

/* tensor tools*/
void check_N(double complex N[NDIM][NDIM], double Kcon[NDIM], double gcov[NDIM][NDIM]);
void complex_lower(double complex N[NDIM][NDIM], 
	double gcov[NDIM][NDIM], 
	int low1, int low2,
	double complex Nl[NDIM][NDIM]
		   );
void stokes_to_tensor( double fI, double fQ, double fU, double fV,
		       double complex f_tetrad[NDIM][NDIM]);
void tensor_to_stokes(double complex f_tetrad[NDIM][NDIM],
		      double *fI, double *fQ, double *fU, double *fV);
void complex_coord_to_tetrad_rank2(double complex T_coord[NDIM][NDIM],
				   double Ecov[NDIM][NDIM], double complex T_tetrad[NDIM][NDIM]);
void complex_tetrad_to_coord_rank2(double complex T_tetrad[NDIM][NDIM],
				   double Econ[NDIM][NDIM], double complex T_coord[NDIM][NDIM]);



/***************************MAIN FUNCTIONS******************************/
/* initialize tensor N in the coordinate frame at the bening of the *
geodesics integration = it is zero */
void init_N(const double X[NDIM],const double Kcon[NDIM],double complex N_coord[NDIM][NDIM])
{
  int m,n,k,i,j,l;
  double SI,SQ,SU,SV,SI0,SQ0,SU0,SV0;
  double complex N_tetrad[NDIM][NDIM];
  double gcov[NDIM][NDIM],gcon[NDIM][NDIM];
  double Ucon[NDIM],Ucov[NDIM],Kcov[NDIM];
  double Ecov[NDIM][NDIM],Econ[NDIM][NDIM];
  double trial1[NDIM],trial2[NDIM];
  
  MNLOOP N_coord[m][n] = 0.0+I*0.0;

  return 0;

  
  //if( exp(X[1]) < 1.06*Rh ) return 0;

  // or this, in case if you want to test transport of Stokes parameters, and test transport part
  SI=4.0;
  SQ=3.0;
  SU=2.0;
  SV=1.0;

  stokes_to_tensor(SI, SQ, SU, SV, N_tetrad);
  
  gcov_func(X,gcov);   
  gcon_func(gcov, gcon);                                                                                             
  get_model_ucon(X, Ucon);                                                                                
  get_model_ucov(X, Ucov);   
  
  lower(Kcon, gcov, Kcov);
  double dot = 0.;
  for(i=0;i<4;i++) dot -= Ucon[i]*Kcov[i];
  trial1[0] = Kcon[0]/dot - Ucon[0];
  trial1[1] = Kcon[1]/dot - Ucon[1];
  trial1[2] = Kcon[2]/dot - Ucon[2];
  trial1[3] = Kcon[3]/dot - Ucon[3];
  trial2[0] = 0. ; trial2[1] = 0. ; trial2[2] = 1. ; trial2[3] = 0. ;

  make_tetrad(Ucon, trial1, trial2, gcov, Econ, Ecov);
  complex_tetrad_to_coord_rank2(N_tetrad,Econ,N_coord);

}


void push_polar(double Xi[NDIM], double Xm[NDIM], double Xf[NDIM],
                double Ki[NDIM], double Km[NDIM], double Kf[NDIM],
                complex double Ni[NDIM][NDIM], complex double Nm[NDIM][NDIM],
                complex double Nf[NDIM][NDIM], double dl)
{
  /* find the connection */
  double lconn[NDIM][NDIM][NDIM];
  get_connection(Xm, lconn);
  int i,j,k,l;

  //no need below, this is done outside and recoreded
  /* push X */
  /*
  for(i=0;i<4;i++)
    Xf[i] = Xi[i] + Km[i]*dl;
  */
  /* push k */
  /*
  for(i=0;i<4;i++)
    Kf[i] = Ki[i];
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      for(k=0;k<4;k++)
	Kf[i] += -lconn[i][j][k]*Km[j]*Km[k]*dl;
  */
  /* push N */
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Nf[i][j] = Ni[i][j];
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      for(k=0;k<4;k++)
        for(l=0;l<4;l++)
	  Nf[i][j] += -(
                        lconn[i][k][l]*Nm[k][j]*Km[l] +
                        lconn[j][k][l]*Nm[i][k]*Km[l]
                        )*dl/(L_unit*HPL/(ME*CL*CL));

  return;
}
/* updates N for one step on geodesics, using the previous step N*/
/* here we compute new rigth-hand side of the equation */
/* and somehow rotate this along the geodesics knowing*/
/* first point and last point X and K*/

void evolve_N(const double Xi[NDIM],const double Kconi[NDIM],
	      const double Xhalf[NDIM], const double Kconhalf[NDIM],
              const double Xf[NDIM],const double Kconf[NDIM],
	      const double dlam,
	      double complex N_coord[NDIM][NDIM])
{
    int i,j,k,l,m,n;
    double gcov[NDIM][NDIM],gcon[NDIM][NDIM],gdet;
    double Kcovi[NDIM],Kcov[NDIM],Kcon[NDIM],Kcovf[NDIM];
    double Ucon[NDIM],Ucov[NDIM];
    double Bcon[NDIM],Bcov[NDIM];
    double Ecov[NDIM][NDIM],Econ[NDIM][NDIM];
    double lconn[NDIM][NDIM][NDIM];
    double complex H[NDIM][NDIM][NDIM][NDIM];
    double complex Nh[NDIM][NDIM];
    double complex Ndd[NDIM][NDIM];
    double complex Jtetrad[NDIM][NDIM],Jcoord[NDIM][NDIM];
    double complex Atetrad[NDIM][NDIM],Acoord[NDIM][NDIM];
    double complex Rtetrad[NDIM][NDIM],Rcoord[NDIM][NDIM];
    double trial1[NDIM];
    double trial2[NDIM];
    double Bhatcon[NDIM],B;
    double complex N_tetrad[NDIM][NDIM];
    double jI,jQ,jU,jV;
    double aI,aQ,aU,aV;
    double rV,rU,rQ,rho2,rho,rdS,dlamrot;
    double SI,SQ,SU,SV;
    double SI0,SQ0,SU0,SV0;
    double Op[NDIM][NDIM];

    /*void push_polar(double Xi[NDIM], double Xm[NDIM], double Xf[NDIM],
                double Ki[NDIM], double Km[NDIM], double Kf[NDIM],
                complex double Ni[NDIM][NDIM], complex double Nm[NDIM][NDIM],
                complex double Nf[NDIM][NDIM], double dl)*/

    
    push_polar(Xi, Xi, Xhalf,
	       Kconi, Kconi, Kconhalf,
	       N_coord, N_coord, Nh,
	       0.5*dlam) ;
    push_polar(Xi, Xhalf, Xf,
	       Kconi, Kconhalf, Kconf,
	       N_coord, Nh, N_coord,
	       dlam) ;
    
    
    //if no source steps, then uncomment below, transports Stokes perfectly even for larger EPS
    //return 0;

    /**************************** SOURCE STEP at Xf *******************************************/
    double Ne = get_model_ne(Xf) ; 

    if(Ne>0.0){

      gcov_func(Xf,gcov);                                                                                                

      jar_calc(Xf,Kconf,&jI,&jQ,&jU,&jV,
	       &aI,&aQ,&aU,&aV,
	       &rQ,&rU,&rV,
	       Op,dlam);

      B = get_model_b(Xf) ;            /* field in G */
      get_model_bcon(Xf, Bcon);
      if (B > 0.) { 
	for (k = 0; k < NDIM; k++) Bhatcon[k] = Bcon[k] / (B / B_unit); 
      }else{
	for (k = 0; k < NDIM; k++) Bhatcon[k] = 0.0; 
	Bhatcon[2] = 1.;
      }
      
      get_model_ucon(Xf, Ucon);                                                                                 
      make_tetrad(Ucon, Kconf, Bhatcon, gcov, Econ, Ecov);
      
      complex_coord_to_tetrad_rank2(N_coord, Ecov, N_tetrad);
      tensor_to_stokes(N_tetrad, &SI0, &SQ0, &SU0, &SV0);
      
      double sum1,sum2;

      /*** All steps (rot/abs/emis) analytically at once, assuming that conditions dont change over a step ***/
      // operator Op is returned by jar_calc, together with jIQUV
    
      sum1=Op[0][0]*jI  + Op[0][1]*jQ  + Op[0][2]*jU  + Op[0][3]*jV;
      sum2=Op[0][0]*SI0 + Op[0][1]*SQ0 + Op[0][2]*SU0 + Op[0][3]*SV0;
      SI = dlam*sum1 + sum2;

      sum1=Op[1][0]*jI  + Op[1][1]*jQ  + Op[1][2]*jU  + Op[1][3]*jV;
      sum2=Op[1][0]*SI0 + Op[1][1]*SQ0 + Op[1][2]*SU0 + Op[1][3]*SV0;
      SQ = dlam*sum1 + sum2;

      sum1=Op[2][0]*jI  + Op[2][1]*jQ  + Op[2][2]*jU  + Op[2][3]*jV;
      sum2=Op[2][0]*SI0 + Op[2][1]*SQ0 + Op[2][2]*SU0 + Op[2][3]*SV0;
      SU = dlam*sum1 + sum2;

      sum1=Op[3][0]*jI  + Op[3][1]*jQ  + Op[3][2]*jU  + Op[3][3]*jV;
      sum2=Op[3][0]*SI0 + Op[3][1]*SQ0 + Op[3][2]*SU0 + Op[3][3]*SV0;
      SV = dlam*sum1 + sum2;

      //pack IQUV back to N_tetrad and transform tensor N from fluid/tetrad to coordinate frame for transport
      stokes_to_tensor(SI, SQ, SU, SV, N_tetrad);
      complex_tetrad_to_coord_rank2(N_tetrad, Econ, N_coord);

    }

    /*************************** SOURCE STEP DONE **************************/	
 	

}


/*converts tensor N to tokes parameters detected at the camera*/
void project_N(const double X[NDIM],const double Kcon[NDIM],const double Ucam[NDIM],const double complex N_coord[NDIM][NDIM],double *Stokes_I,double *Stokes_Q,double *Stokes_U,double *Stokes_V)
{
  int i,j,l,k,m,n;
  double gcov[NDIM][NDIM];
  double Kcov[NDIM];
  double complex N_tetrad[NDIM][NDIM];
  double Econ[NDIM][NDIM],Ecov[NDIM][NDIM] ;
  double trial1[NDIM],trial2[NDIM] ;
  double SI,SQ,SU,SV;
  
  gcov_func(X,gcov);    
  lower(Kcon,gcov,Kcov);
  
  double dot = 0.;
  for(i=0;i<4;i++) dot -= Ucam[i]*Kcov[i];
  trial1[0] = Kcon[0]/dot - Ucam[0];
  trial1[1] = Kcon[1]/dot - Ucam[1];
  trial1[2] = Kcon[2]/dot - Ucam[2];
  trial1[3] = Kcon[3]/dot - Ucam[3];
  trial2[0] = 0. ;
  trial2[1] = 0. ;
  trial2[2] = 1. ;
  trial2[3] = 0. ;
  make_tetrad(Ucam, trial1, trial2, gcov, Econ, Ecov) ;
  
  complex_coord_to_tetrad_rank2(N_coord, Ecov, N_tetrad);
  tensor_to_stokes(N_tetrad, &SI, &SQ, &SU, &SV);
  
  *Stokes_I=SI;
  *Stokes_Q=SQ;
  *Stokes_U=SU;
  *Stokes_V=SV;

  return;

}

/***************************END MAIN FUNCTIONS******************************/


/*************************SUPPORTING FUNCTIONS******************************/

/*invariant plasma emissivities/abs/rot in tetrad frame */
void jar_calc(const double X[NDIM], const double Kcon[NDIM], 
		    double *jI, double *jQ, double *jU, double *jV, 
		    double *aI, double *aQ, double *aU, double *aV, 
	            double *rQ, double *rU, double *rV,
                    double Op[NDIM][NDIM],double xx)
{
  int m,n,i,j,k,l,ii,jj,iii,jjj;
  double nu,Thetae,Ne,B,theta,en,nu2;
  double x,Xe,omega0,omega,nuc,nub;
  double Bnu1,Binv;
  double Ucov[NDIM];
  double Thetaer, wp2;
  double M1[NDIM][NDIM];
  double M2[NDIM][NDIM];
  double M3[NDIM][NDIM];
  double M4[NDIM][NDIM];

  Ne = get_model_ne(X) ; 

  get_model_ucov(X, Ucov);
  theta = get_bk_angle(X,Kcon,Ucov) ;     /* angle between k & b  */

  if(theta <= 0. || theta >= M_PI) {  /* no emission along field  */

    //if no absorption then operator O=unity;
    KLLOOP Op[k][l]=0.0;
    Op[0][0]=1.0;
    Op[1][1]=1.0;
    Op[2][2]=1.0;  
    Op[3][3]=1.0;

    *jI=0.0;
    *jQ=0.0;
    *jU=0.0;
    *jV=0.0;
  
    *aI=0.0;
    *aQ=0.0;
    *aU=0.0;
    *aV=0.0;
    
    *rQ=0.0;
    *rV=0.0;
    *rU=0.0;
    
    return;
    
  } else {
    
    nu = get_fluid_nu(Kcon,Ucov);// * freqcgs1 ;  /* freq in Hz          */
    nu2=nu*nu;
    en = HPL* nu;                        /* erg s Hz                  */ 
    B = get_model_b(X) ;            /* field in G                */
    Thetae = get_model_thetae(X) ;  /* temp in e rest-mass units */
    Thetaer=1./Thetae;
    
    omega0=EE*B/ME/CL;
    wp2=4.*M_PI*Ne*EE*EE/ME;

    /*Faraday rotativities for thermal plasma */
    Xe=Thetae*sqrt(S2*sin(theta)*(1.e3*omega0/2./M_PI/nu));
    
      /*below expressions are correct but use gsl functions*/
      /*
	    *rV=2.0*M_PI*nu/CL *
	      wp2*omega0/pow(2.*M_PI*nu,3)*
	      (gsl_sf_bessel_Kn(0,Thetaer)-Je(Xe))/gsl_sf_bessel_Kn(2,Thetaer)*cos(theta);

	    *rQ=2.*M_PI*nu/2./CL*jffunc(Xe)*wp2*omega0*omega0/pow(2*M_PI*nu,4)*
	      (gsl_sf_bessel_Kn(1,Thetaer)/gsl_sf_bessel_Kn(2,Thetaer)+6.*Thetae)*sin(theta)*sin(theta);

	    *rU=0.0;
	    
	    *rQ *= nu;
	    *rV *= nu;
    */

      /* Here I use approximate bessel functions to match rhoqv with grtrans */
    *rQ=2.*M_PI*nu/2./CL*
      wp2*omega0*omega0/pow(2*M_PI*nu,4)* 
      jffunc(Xe)*
      (besselk_asym(1,Thetaer)/besselk_asym(2,Thetaer)+6.*Thetae)*sin(theta)*sin(theta);
      
    *rU=0.0; 
      
    *rV=2.0*M_PI*nu/CL *                                                                                               
      wp2*omega0/pow(2.*M_PI*nu,3)*
      (besselk_asym(0,Thetaer)-Je(Xe))/besselk_asym(2,Thetaer)*cos(theta);                                                                                                           
    /* invariant rotativities, because I update invatiant stokes parameters using analytical solution */
    *rQ *= nu;                                                                                          
    *rV *= nu; 
      
      //*rQ=0.0; //uncomment to turn off term responsible for Faraday conversion
      //*rV=0.0; //uncomment to turn off term responsible for Faraday rotation
      
      /*synchrotron emissivity*/
    nuc=3.0*EE*B*sin(theta)/4.0/M_PI/ME/CL*Thetae*Thetae+1.0;
    x=nu/nuc;

    *jI=Ne*EE*EE*nu/2./S3/CL/Thetae/Thetae*I_I(x) ; // [g/s^2/cm]
    *jQ=Ne*EE*EE*nu/2./S3/CL/Thetae/Thetae*I_Q(x) ;
    *jU=0.0;    // just a convention 
    *jV=2.*Ne*EE*EE*nu/tan(theta)/3./S3/CL/Thetae/Thetae/Thetae*I_V(x) ;

    /* Planck function */
    Bnu1=Bnu(nu,Thetae);
    
    /* invariant synchrotron absopritivty */
    *aI=*jI/(Bnu1)*nu;
    *aQ=*jQ/(Bnu1)*nu;
    *aU=*jU/(Bnu1)*nu;
    *aV=*jV/(Bnu1)*nu;
      
    /* invaraint emissivity */
    *jI/=nu2;
    *jQ/=nu2;
    *jU/=nu2;
    *jV/=nu2;
    
    /*****************************************************************/

    //all additional coefficients used to compute operator O, page 432
      double alpha2=(*aQ)*(*aQ)+(*aU)*(*aU)+(*aV)*(*aV);
      double rho2=(*rQ)*(*rQ)+(*rU)*(*rU)+(*rV)*(*rV);
      double alphadrho=(*aQ)*(*rQ)+(*aU)*(*rU)+(*aV)*(*rV);
      double sigm=sign(alphadrho);
      double th=2.*sqrt(  pow(alpha2 - rho2,2)/4.+ alphadrho*alphadrho );
      double ith=1./th;
      double Lam1=sqrt( sqrt(  pow(alpha2 - rho2,2)/4. + alphadrho*alphadrho ) + (alpha2-rho2)/2.);
      double Lam2=sqrt( sqrt(  pow(alpha2 - rho2,2)/4. + alphadrho*alphadrho ) - (alpha2-rho2)/2.);


      //M1
      KLLOOP M1[k][l]=0.0;
      M1[0][0]=1.0;
      M1[1][1]=1.0;
      M1[2][2]=1.0;
      M1[3][3]=1.0;

      //M2
      //row 0
      M2[0][0]=0.0;
      M2[0][1]=ith*(Lam2*(*aQ)-sigm*Lam1*(*rQ));
      M2[0][2]=ith*(Lam2*(*aU)-sigm*Lam1*(*rU));
      M2[0][3]=ith*(Lam2*(*aV)-sigm*Lam1*(*rV));
    
      //row 1
      M2[1][0]=ith*(Lam2*(*aQ)-sigm*Lam1*(*rQ));
      M2[1][1]=0.0;
      M2[1][2]=ith*(sigm*Lam1*(*aV)+Lam2*(*rV));
      M2[1][3]=ith*(-sigm*Lam1*(*aU)-Lam2*(*rU));
      
      M2[2][0]=ith*(Lam2*(*aU)-sigm*Lam1*(*rU));
      M2[2][1]=ith*(-sigm*Lam1*(*aV)-Lam2*(*rV));
      M2[2][2]=0.0;
      M2[2][3]=ith*(sigm*Lam1*(*aQ)+Lam2*(*rQ));

      M2[3][0]=ith*(Lam2*(*aV)-sigm*Lam1*(*rV));
      M2[3][1]=ith*(sigm*Lam1*(*aU)+Lam2*(*rU));
      M2[3][2]=ith*(-sigm*Lam1*(*aQ)-Lam2*(*rQ));
      M2[3][3]=0.0;


      //M3
      //row 0
      M3[0][0]=0.0;
      M3[0][1]=ith*(Lam1*(*aQ)+sigm*Lam2*(*rQ));
      M3[0][2]=ith*(Lam1*(*aU)+sigm*Lam2*(*rU));//jason paper has mistake here, I follow book, grtrans is OK. same as here
      M3[0][3]=ith*(Lam1*(*aV)+sigm*Lam2*(*rV));

      //row 1
      M3[1][0]=ith*(Lam1*(*aQ)+sigm*Lam2*(*rQ));
      M3[1][1]=0.0;
      M3[1][2]=ith*(-sigm*Lam2*(*aV)+Lam1*(*rV));
      M3[1][3]=ith*(sigm*Lam2*(*aU)-Lam1*(*rU));

      M3[2][0]=ith*(Lam1*(*aU)+sigm*Lam2*(*rU));
      M3[2][1]=ith*(sigm*Lam2*(*aV)-Lam1*(*rV));
      M3[2][2]=0.0;
      M3[2][3]=ith*(-sigm*Lam2*(*aQ)+Lam1*(*rQ));

      M3[3][0]=ith*(Lam1*(*aV)+sigm*Lam2*(*rV));
      M3[3][1]=ith*(-sigm*Lam2*(*aU)+Lam1*(*rU));
      M3[3][2]=ith*(sigm*Lam2*(*aQ)-Lam1*(*rQ));
      M3[3][3]=0.0;

      //M4
      //row 0 in M4 matrix
      M4[0][0]=2.*ith*(alpha2+rho2)/2.;
      M4[0][1]=2.*ith*((*aV)*(*rU)-(*aU)*(*rV));
      M4[0][2]=2.*ith*((*aQ)*(*rV)-(*aV)*(*rQ));
      M4[0][3]=2.*ith*((*aU)*(*rQ)-(*aQ)*(*rU));

      //row 1 etc
      M4[1][0]=2.*ith*((*aU)*(*rV)-(*aV)*(*rU));
      M4[1][1]=2.*ith*((*aQ)*(*aQ)+(*rQ)*(*rQ) - (alpha2+rho2)/2.);
      M4[1][2]=2.*ith*((*aQ)*(*aU)-(*rQ)*(*rU));
      M4[1][3]=2.*ith*((*aV)*(*aQ)+(*rV)*(*rQ));

      M4[2][0]=2.*ith*((*aV)*(*rQ)-(*aQ)*(*rV));
      M4[2][1]=2.*ith*((*aQ)*(*aU)+(*rQ)*(*rU));
      M4[2][2]=2.*ith*((*aU)*(*aU)+(*rU)*(*rU) - (alpha2+rho2)/2.);
      M4[2][3]=2.*ith*((*aU)*(*aV)+(*rU)*(*rV));

      M4[3][0]=2.*ith*((*aQ)*(*rU)-(*aU)*(*rQ));
      M4[3][1]=2.*ith*((*aV)*(*aQ)+(*rV)*(*rQ));
      M4[3][2]=2.*ith*((*aU)*(*aV)+(*rU)*(*rV));
      M4[3][3]=2.*ith*((*aV)*(*aV)+(*rV)*(*rV) - (alpha2+rho2)/2.);

      //finally operator O
      
      KLLOOP Op[k][l]=exp(- (*aI)*xx)*( 0.5*(cosh(Lam1*xx)+cos(Lam2*xx))*M1[k][l]  //here typo in book cos instead of cosh
					- sin(Lam2*xx)*M2[k][l]
					- sinh(Lam1*xx)*M3[k][l]
					+ 0.5*(cosh(Lam1*xx)-cos(Lam2*xx))*M4[k][l]
      			    );
      

      
      /******************************************************************/      
      
    }
  
  

 
}

/*emissivity functions and functions used for Faraday conversion and rotation*/
/*from J. Dexter PhD thesis (checked with Leung harmony program, and Huang & Shcherbakov 2011*/
double g(double Xe){ 
return 1.-0.11*log(1+0.035*Xe);
}
double h(double Xe){
  return 2.011*exp(-pow(Xe,1.035)/4.7)-cos(Xe*0.5)*exp(-pow(Xe,1.2)/2.73)-0.011*exp(-Xe/47.2);
}
double Je(double Xe){ 
  return 0.43793091*log(1.+0.00185777*pow(Xe,1.50316886));
}
double jffunc(double Xe){
  double extraterm;
  extraterm = (0.011*exp(-Xe/47.2)-pow(2.,-1./3.)/pow(3.,23./6.)*M_PI*1e4* pow(Xe+1e-16,-8./3.))*(0.5+0.5*tanh((log(Xe)-log(120.))/0.1));
  return 2.011*exp(-pow(Xe,1.035)/4.7) - cos(Xe*0.5)*exp(-pow(Xe,1.2)/2.73) - 0.011*exp(-Xe/47.2) + extraterm;
}
double I_I(double x){
    return 2.5651*(1+1.92*pow(x,-1./3.)+0.9977*pow(x,-2./3.))*exp(-1.8899*pow(x,1./3.));
}
double I_Q(double x){
    return 2.5651*(1+0.93193*pow(x,-1./3.)+0.499873*pow(x,-2./3.))*exp(-1.8899*pow(x,1./3.));
}
double I_V(double x){
    return (1.81348/x+3.42319*pow(x,-2./3.)+0.0292545*pow(x,-0.5)+2.03773*pow(x,-1./3.))*exp(-1.8899*pow(x,1./3.));
}

/* Planck funtion: Bnu Bnu_inv*nu**3*/
double Bnu(double nu, double Thetae)
{
    double x;
    x = HPL*nu/(ME*CL*CL*Thetae);

    //    if(x < 1.e-3)   /* Taylor expand */
	    // return ((2.*HPL/(CL*CL))/( x/24. * (24. + x*(12. + x*(4. + x))))*nu*nu*nu);
    if(x < 1.e-6)   /* Taylor expand */
      return 2*nu*nu*ME*Thetae; //this way to match grtrans
    else
      return((2.*HPL/(CL*CL))/(exp(x) - 1.)*nu*nu*nu) ;
    
}

double besselk_asym(int n, double x){

  if(n==0) return -log(x/2.)-0.5772;

  if(n==1) return 1./x;

  if(n==2) return 2./x/x;


}

/*end of emissivity functions*/


void check_N(double complex N[NDIM][NDIM], 
		double Kcon[NDIM], 
		double gcov[NDIM][NDIM]) {
	double complex dot;
	double Kcov[NDIM];
	int i,j;

	fprintf(stderr, "enter check_N\n");

	/* compute k . N */
	lower(Kcon, gcov, Kcov);
	fprintf(stderr,"(k . N\n");
	/* first one way */
	for(i=0;i<4;i++) {
		dot = 0. + I*0.;
		for(j=0;j<4;j++) dot += Kcov[j]*N[j][i];
		fprintf(stderr,"%d %g + i %g\n",i,creal(dot),cimag(dot));
	}
	/* then the other */
	for(i=0;i<4;i++) {
		dot = 0. + I*0.;
		for(j=0;j<4;j++) dot += Kcov[j]*N[i][j];
		fprintf(stderr,"%d %g + i %g\n",i,creal(dot),cimag(dot));
	}
	fprintf(stderr,"k . N)\n");

#if 1
	/* check for hermiticity */
	fprintf(stderr,"(herm:\n");
	for(i=0;i<4;i++) 
	for(j=0;j<4;j++) fprintf(stderr,"%d %d %g + i %g\n",i,j,
		creal(N[i][j] - conj(N[j][i]) ),
		cimag(N[i][j] - conj(N[j][i]) )
		) ;
	fprintf(stderr,"herm)\n");

	/* check invariants */
	double complex Nud[NDIM][NDIM];
	void complex_lower(double complex N[NDIM][NDIM], double gcov[NDIM][NDIM], int low1, int low2,
		double complex Nl[NDIM][NDIM]);
	complex_lower(N, gcov, 0, 1, Nud) ;
	for(i=0;i<4;i++) fprintf(stderr,"N: %d %g + i %g\n",i, creal(N[i][i]),cimag(N[i][i]));
	for(i=0;i<4;i++) fprintf(stderr,"Nud: %d %g + i %g\n",i, creal(Nud[i][i]),cimag(Nud[i][i]));
	dot = 0. + I*0.;
	for(i=0;i<4;i++) dot += Nud[i][i];
	fprintf(stderr,"I: %g + i %g\n",creal(dot),cimag(dot));

	double complex Ndd[NDIM][NDIM];
	complex_lower(N, gcov, 1, 1, Ndd);
	dot = 0. + I*0.;
	for(i=0;i<4;i++) 
	for(j=0;j<4;j++) 
		dot += 2. * 0.25* (N[i][j] + N[j][i])* (Ndd[i][j] + Ndd[j][i]);
	fprintf(stderr,"IQUsq: %g + i %g\n",creal(dot),cimag(dot));

	dot = 0. + I*0.;
	for(i=0;i<4;i++) 
	for(j=0;j<4;j++) 
		dot += -2. * 0.25* (N[i][j] - N[j][i])* (Ndd[i][j] - Ndd[j][i]);
	fprintf(stderr,"Vsqsq: %g + i %g\n",creal(dot),cimag(dot));
#endif

	fprintf(stderr, "leave check_N\n");
}



void complex_lower(double complex N[NDIM][NDIM], 
	double gcov[NDIM][NDIM], 
	int low1, int low2,
	double complex Nl[NDIM][NDIM]
	)
{
	int i,j,k,l;

	if(!low1 && !low2) return;

	if(low1 && low2) {
		for(i=0;i<4;i++)
		for(j=0;j<4;j++) {
			Nl[i][j] = 0. + I*0.;
			for(k=0;k<4;k++)
			for(l=0;l<4;l++) {
				Nl[i][j] += N[k][l]*gcov[k][i]*gcov[l][j] ;
			}
		}
		return;
	}

	if(low1) {
		for(i=0;i<4;i++)
		for(j=0;j<4;j++) {
			Nl[i][j] = 0. + I*0.;
			for(k=0;k<4;k++)
				Nl[i][j] += N[k][j]*gcov[k][i];
		}
		return;
	}

	for(i=0;i<4;i++)
	for(j=0;j<4;j++) {
		Nl[i][j] = 0. + I*0.;
		for(l=0;l<4;l++) 
			Nl[i][j] += N[i][l]*gcov[l][j];
	}
	return;
	
}

void stokes_to_tensor( double fI, double fQ, double fU, double fV,
		       double complex f_tetrad[NDIM][NDIM])
{
  int i,j;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++) f_tetrad[i][j] = 0. + I*0.;
  /*notice that I swapped sign of the imaginary part [2][3] in [3][2] - which one is correct? */
  f_tetrad[2][2] = (fI + fQ + 0.*I);
  f_tetrad[2][3] = (fU - I*fV);
  f_tetrad[3][2] = (fU + I*fV);
  f_tetrad[3][3] = (fI - fQ + 0.*I);

}
void tensor_to_stokes(double complex f_tetrad[NDIM][NDIM],
		      double *fI, double *fQ, double *fU, double *fV)
{

  /*here I divide by two to agree with above*/
  *fI = creal(f_tetrad[2][2] + f_tetrad[3][3])/2;
  *fQ = creal(f_tetrad[2][2] - f_tetrad[3][3])/2;
  *fU = creal(f_tetrad[2][3] + f_tetrad[3][2])/2;
  *fV = cimag(f_tetrad[3][2] - f_tetrad[2][3])/2;

}

void complex_coord_to_tetrad_rank2(double complex T_coord[NDIM][NDIM],
				   double Ecov[NDIM][NDIM], double complex T_tetrad[NDIM][NDIM])
{
  int i,j,k,l;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      T_tetrad[i][j] = 0. + I*0.;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      for(k=0;k<4;k++)
        for(l=0;l<4;l++)
	  T_tetrad[i][j] += T_coord[k][l]*Ecov[i][k]*Ecov[j][l] ;

  return;
}

void complex_tetrad_to_coord_rank2(double complex T_tetrad[NDIM][NDIM],
				   double Econ[NDIM][NDIM], double complex T_coord[NDIM][NDIM])
{
  int i,j,k,l;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++) T_coord[i][j] = 0. + I*0.;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      for(k=0;k<4;k++)
        for(l=0;l<4;l++) T_coord[i][j] += T_tetrad[k][l]*Econ[k][i]*Econ[l][j];

  return;
}


//these definitions are only used in ipolarray.c
#undef S2     
#undef S3     
#undef KLOOP
#undef ILOOP
#undef JLOOP
#undef LLOOP
#undef IJLOOP
#undef MNLOOP
#undef IJKLLOOP
#undef MNKLLOOP






