
#include "decs.h"

/* 

"mixed" emissivity formula 

interpolates between Petrosian limit and
classical thermal synchrotron limit

good for Thetae >~ 1

*/

double jnu_synch(double nu, double Ne, double Thetae, double B, double theta)
{
        double K2,nuc,nus,x,f,j,sth ;

        //K2 = gsl_sf_bessel_Kn(2,1./Thetae) ;
        K2 = 2.*Thetae*Thetae ;

        nuc = EE*B/(2.*M_PI*ME*CL) ;
	sth = sin(theta) ;
        nus = (2./9.)*nuc*Thetae*Thetae*sth ;
	if(nu > 1.e12*nus) return(0.) ;
        x = nu/nus ;
        f = pow( pow(x,1./2.) + pow(2.,11./12.)*pow(x,1./6.), 2 ) ;
        j = (sqrt(2.)*M_PI*EE*EE*Ne*nus/(3.*CL*K2)) * f * exp(-pow(x,1./3.)) ;

        return(j) ;
}

