/*Author: Monika Moscibrodzka */
/*Last update 1.05.2015*/
/*Contains all functions needed to calculate free-free emission*/

#include "decs.h"
#include <gsl/gsl_sf_bessel.h>



/*cross-section for p-e bremsstrahlung in the Born approximation by Heitler (1954)*/
/*this must be in units of SIGMA_THOMSON*/
double w_dsigma_dw(double w, double gamma)
{
  double w_dsigma_dwp ;
  double gammap,gamma2,gammap2;

  gammap=gamma-w;
  gamma2=gamma*gamma;  
  gammap2=gammap*gammap;

  /*for ultrarelativistic limit Jauch&Rohlich - from Stepney thesis*/
  w_dsigma_dwp=3.*ALPHAF/2./M_PI * (1.-2.*gammap/3./gamma+gammap2/gamma2)*(log(2.*gamma*gammap/w)-0.5);
  /*need one in the Born approximation from Heitler 1954*/

  return(w_dsigma_dwp*SIGMA_THOMSON) ;

}


/* bremsstrahlung emissivity for e-e interactions, ergs/s/cm3/Hz/str */                                                            
/* analytical functions for thermal distribution of electrons*/
double jnu_ff_ee(double nu, double Ne, double Thetae){

  double x,j;

  x=HPL_MECL2*nu/Thetae;

  /*for low temperatures theta<0.02 (50 keV) - Gould 1980,81, nonrelativistic*/
  if(Thetae < 0.1){
    j=2./M_PI*Bl(x)*sqrt(2.*Thetae/M_PI)*exp(-x/2.)*gsl_sf_bessel_Kn(0,x/2.);
  }

  /*this is for temperaures 0.01-100 (50keV-1MeV) - Stepney and Guilbert 1983*/
  if(Thetae >= 0.1 && Thetae <= 2.){
    j=Thetae*exp(-x)*G_ff(x,Thetae);
  }

  /*for high temepratures theta>2 (1 MeV ) - Alexian 1968, ultrarelativistic*/
  if(Thetae > 2.){
    j=3./4./M_PI*exp(-x)*
      (28./3. + 2.*x + x*x/2. + 2.*(8./3.+4./3.*x+x)*(log(2.*Thetae)-0.5772)-
       exp(-x)*Ei(-x)*(8./3. - 43.*x + x*x)) ;
  }

  return j/4./M_PI*Ne*Ne*SIGMA_THOMSON*CL*HPL*ALPHAF;
 
}


/*Approximation from Maxon 1972*/
double Ei(double xx){
    double x,x2,Eii;
    x=fabs(xx);	
    x2=x*x;

    if (x<=1.){
	Eii=log(x)+0.577-x+0.25*x2-0.055*x2*x;
    }else if(x>1.){
	Eii=-exp(-x)/x*(x2+2.33*x+0.25)/(x2+3.33*x+1.68);
    }
    return Eii;
}


double Bl(double x){
   
    return 0.85+1.35*sqrt(x)+0.38*x;

}


/*calculates e-e spectrum using tables and interpolations*/
double G_ff(double x,double Thetae){

    int i,ii;
    double en;
    double f1A,f1B,f1C,f1D;
    double f2A,f2B,f2C,f2D;
    double x1,x2;
    double AA,BB,CC,DD;
    double alphaa,betab,gammag,deltad;
    double Gfun1;
    double energy[]={50,75,100,150,200,300,400,500,600,700,800,900,1000}; //keV
  
  /* low energy data for different temperatures*/

  double A[] = { 1.584, 1.357, 1.197, 1.023, 0.883, 
		 0.700, 0.572, 0.484, 0.417, 0.361, 
		 0.322, 0.286, 0.259 };
  
  double B[] = { 0.578, 0.437, 0.291, 0.204, 0.0835, 
		 -0.0494, -0.139, -0.181, -0.209, -0.204, 
		 -0.244, -0.257, -0.258 };
  
  double C[] = {4.565, 3.842, 3.506, 3.036, 2.831, 
		2.545, 2.352, 2.175, 2.028, 1.914, 
		1.795, 1.705, 1.617 };
  
  double D[] = {2.091, 1.855, 1.672, 1.593, 1.487, 
		1.364, 1.254, 1.179, 1.108, 1.030, 
		0.982, 0.923, 0.879 };
  

    /* high energy data for different temperatures*/
  double alpha[] = {0.0387,0.0633, 0.0862, 0.128, 0.159, 
		    0.208, 0.234, 0.245, 0.248, 0.247, 
		    0.243, 0.239, 0.235};
  
  double beta[] = {0.523, 0.540, 0.569, 0.596, 0.658, 
		   0.633, 0.643, 0.695, 0.729, 0.756, 
		   0.763, 0.755, 0.735};
    
  double gamma[] = {5.319, 4.412, 3.897, 3.383, 2.974, 
		    2.738, 2.424, 2.025, 1.716, 1.457, 
		    1.271, 1.140, 1.060};

  double delta[] = {0.782, 0.689, 0.633, 0.523, 0.532, 
		    0.326, 0.302, 0.394, 0.453, 0.500, 
		    0.515, 0.508, 0.478};

  en=Thetae*ME*CL*CL/KEV;

  /* find ABCD or alpha,beta,gamma,detla for a given temperature*/

  
  if(x > 0.05 && x < 1.1){
      
      if(en <  50){
        x1 = energy[0];
        x2 = energy[1];
	f1A = A[0];
	f2A = A[1];
	f1B = B[0];
	f2B = B[1];
	f1C = C[0];
	f2C = C[1];
	f1D = D[0];
	f2D = D[1];
      }else if(en > 1000){
	x1 = energy[11];
	x2 = energy[12];
	f1A = A[11];
	f2A = A[12];
	f1B = B[11];
	f2B = B[12];
	f1C = C[11];
	f2C = C[12];
	f1D = D[11];
	f2D = D[12];
      }else{

	  for(ii=0;ii<12;ii++){
	      if( en >= energy[ii] && en <= energy[ii+1]){
		  i=ii;
		  x1 = energy[i] ;
		  x2 = energy[i+1];
		  //fprintf(stdout,"1:i=%d ii=%d en=%g x1=%g x2=%g\n",i,ii,en,x1,x2);
		  break;
	      }
	  }
	  
	  
	  f1A = A[i];
	  f2A = A[i + 1];
	  f1B = B[i];
	  f2B = B[i + 1];
	  f1C = C[i];
	  f2C = C[i + 1];
	  f1D = D[i];
	  f2D = D[i + 1];
      }
      AA = (f2A - f1A) / (x2 - x1) * en + (f1A * x2 - x1 * f2A) / (x2 - x1);
      BB = (f2B - f1B) / (x2 - x1) * en + (f1B * x2 - x1 * f2B) / (x2 - x1);
      CC = (f2C - f1C) / (x2 - x1) * en + (f1C * x2 - x1 * f2C) / (x2 - x1);
      DD = (f2D - f1D) / (x2 - x1) * en + (f1D * x2 - x1 * f2D) / (x2 - x1);

      /*interpolation or extrapolation*/
      Gfun1=(AA+BB*x)*log(1./x)+CC+DD*x;
      
    }else if(x > 1.1 && x < 10.){
	    
      if(en < 50){
	x1 = energy[0];
	x2 = energy[1];
	f1A = alpha[0];
	f2A = alpha[1];
	f1B = beta[0];
	f2B = beta[1];
	f1C = gamma[0];
	f2C = gamma[1];
	f1D = delta[0];
	f2D = delta[1];
      }else if(en > 1000){
	x1 = energy[11];
	x2 = energy[12];
	f1A = alpha[11];
	f2A = alpha[12];
	f1B = beta[11];
	f2B = beta[12];
	f1C = gamma[11];
	f2C = gamma[12];
	f1D = delta[11];
	f2D = delta[12];
      }else{
	  for(ii=0;ii<12;ii++){
	      if(en >= energy[ii] && en <= energy[ii+1]){
		  i=ii;
		  x1 = energy[i] ;
		  x2 = energy[i+1];
		  //fprintf(stdout,"2:i=%d ii=%d en=%g x1=%g x2=%g\n",i,ii,en,x1,x2);
		  break;
	      }
	  }
	f1A = alpha[i];
	f2A = alpha[i + 1];
	f1B = beta[i];
	f2B = beta[i + 1];
	f1C = gamma[i];
	f2C = gamma[i + 1];
	f1D = delta[i];
	f2D = delta[i + 1];
      }
      alphaa =(f2A - f1A) / (x2 - x1) * en + (f1A * x2 - x1 * f2A) / (x2 - x1);
      betab  =(f2B - f1B) / (x2 - x1) * en + (f1B * x2 - x1 * f2B) / (x2 - x1);
      gammag =(f2C - f1C) / (x2 - x1) * en + (f1C * x2 - x1 * f2C) / (x2 - x1);
      deltad =(f2D - f1D) / (x2 - x1) * en + (f1D * x2 - x1 * f2D) / (x2 - x1);

      Gfun1=alphaa*x*x+betab*x+gammag+deltad/x;
    }else{

      Gfun1=0.0;       /*for other x'es it is zero???*/

    }
    
    return Gfun1;
    
}

