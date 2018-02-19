/*

	model-dependent geometry routines:
		gcov
		get_connection

*/

#include "decs.h"

/* 
   Current metric: modified Kerr-Schild, squashed in theta
   to give higher resolution at the equator 
*/

// !!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!! //
//                                       //
// MUST CHANGE th(X[2]) IN THREE PLACES! //
//                                       //
// !!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!! //

void gcov_func(double *X, double gcov[][NDIM])
{
	int k, l;
	double sth, cth, s2, rho2;
	double r, th;
	double tfac, rfac, hfac, pfac;

	DLOOP gcov[k][l] = 0.;

        r = exp(X[1]) + R0 ;

        //th =  M_PI*X[2] ;
        th = M_PI * X[2] + (1.-hslope)/2.*sin(2. * M_PI * X[2]);

	cth = cos(th);
	sth = fabs(sin(th));
	if (sth < SMALL)
	  sth = SMALL;

	s2 = sth * sth;
	rho2 = r * r + a * a * cth * cth;

	/* transformation for Kerr-Schild -> modified Kerr-Schild */
	tfac = 1.;
	rfac = r - R0;
	hfac = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]);//M_PI;
	pfac = 1.;

	gcov[0][0] = (-1. + 2. * r / rho2) * tfac * tfac;
	gcov[0][1] = (2. * r / rho2) * tfac * rfac;
	gcov[0][3] = (-2. * a * r * s2 / rho2) * tfac * pfac;

	gcov[1][0] = gcov[0][1];
	gcov[1][1] = (1. + 2. * r / rho2) * rfac * rfac;
	gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2)) * rfac * pfac;

	gcov[2][2] = rho2 * hfac * hfac;

	gcov[3][0] = gcov[0][3];
	gcov[3][1] = gcov[1][3];
	gcov[3][3] =
	    s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2)) * pfac * pfac;

}

/* 

   connection calculated analytically for modified Kerr-Schild
   	coordinates 


   this gives the connection coefficient
	\Gamma^{i}_{j,k} = conn[..][i][j][k]
   where i = {1,2,3,4} corresponds to, e.g., {t,ln(r),theta,phi}
*/

void get_connection(double X[4], double lconn[4][4][4])
{
  //get_connection_num(X, lconn);
  //return;
  
	double r1,r2,r3,r4,sx,cx;
	double th,dthdx2,dthdx22,d2thdx22,sth,cth,sth2,cth2,sth4,cth4,s2th,c2th;
	double a2,a3,a4,rho2,irho2,rho22,irho22,rho23,irho23,irho23_dthdx2;
	double fac1,fac1_rho23,fac2,fac3,a2cth2,a2sth2,r1sth2,a4cth4;

	r1 = exp(X[1]);
	r2 = r1*r1;
	r3 = r2*r1;
	r4 = r3*r1;

	sx = sin(2.*M_PI*X[2]);
	cx = cos(2.*M_PI*X[2]);

	/* HARM-2D MKS */
	th = M_PI*X[2] + 0.5*(1-hslope)*sx;
	dthdx2 = M_PI*(1.+(1-hslope)*cx);
	d2thdx22 = -2.*M_PI*M_PI*(1-hslope)*sx;

	/* HARM-3D MKS */
	//th = th_beg + th_len*X[2] + hslope*sx;
	//th = th_len*X[2] + hslope*sx;
	//dthdx2 = th_len + 2.*M_PI*hslope*cx;
	//d2thdx22 = -4.*M_PI*M_PI*hslope*sx;	/* d^2(th)/dx2^2 */

	//this is needed for radiative transfer which is independend of the grid used in the harm3d code
	//so unless the grid is squezzed in theta direction (our grid is not, it is just trimmed)
	//it should be same as in the old one
	/*new from hotakarun*/
	//is this for geodesics? should be independent of mhd model
  
	//th =  M_PI*X[2];
	//dthdx2 = M_PI;
	//d2thdx22 = 0.0; 

	dthdx22 = dthdx2*dthdx2;

	sth = sin(th);
	cth = cos(th);
	sth2 = sth*sth;
	r1sth2 = r1*sth2;
	sth4 = sth2*sth2;
	cth2 = cth*cth;
	cth4 = cth2*cth2;
	s2th = 2.*sth*cth;
	c2th = 2*cth2 - 1.;

	a2 = a*a;
	a2sth2 = a2*sth2;
	a2cth2 = a2*cth2;
	a3 = a2*a;
	a4 = a3*a;
	a4cth4 = a4*cth4;

	rho2 = r2 + a2cth2;
	rho22 = rho2*rho2;
	rho23 = rho22*rho2;
	irho2 = 1./rho2;
	irho22 = irho2*irho2;
	irho23 = irho22*irho2;
	irho23_dthdx2 = irho23/dthdx2;

	fac1 = r2 - a2cth2;
	fac1_rho23 = fac1*irho23;
	fac2 = a2 + 2*r2 + a2*c2th;
	fac3 = a2 + r1*(-2. + r1);

	lconn[0][0][0] = 2.*r1*fac1_rho23;
	lconn[0][0][1] = r1*(2.*r1+rho2)*fac1_rho23;
	lconn[0][0][2] = -a2*r1*s2th*dthdx2*irho22;
	lconn[0][0][3] = -2.*a*r1sth2*fac1_rho23;

	lconn[0][1][0] = lconn[0][0][1];
	lconn[0][1][1] = 2.*r2*(r4 + r1*fac1 - a4cth4)*irho23;
	lconn[0][1][2] = -a2*r2*s2th*dthdx2*irho22;
	lconn[0][1][3] = a*r1*(-r1*(r3 + 2*fac1) + a4cth4)*sth2*irho23;

	lconn[0][2][0] = lconn[0][0][2];
	lconn[0][2][1] = lconn[0][1][2];
	lconn[0][2][2] = -2.*r2*dthdx22*irho2;
	lconn[0][2][3] = a3*r1sth2*s2th*dthdx2*irho22;

	lconn[0][3][0] = lconn[0][0][3];
	lconn[0][3][1] = lconn[0][1][3];
	lconn[0][3][2] = lconn[0][2][3];
	lconn[0][3][3] = 2.*r1sth2*(-r1*rho22 + a2sth2*fac1)*irho23;

	lconn[1][0][0] = fac3*fac1/(r1*rho23);
	lconn[1][0][1] = fac1*(-2.*r1 + a2sth2)*irho23;
	lconn[1][0][2] = 0.;
	lconn[1][0][3] = -a*sth2*fac3*fac1/(r1*rho23);

	lconn[1][1][0] = lconn[1][0][1];
	lconn[1][1][1] = (r4*(-2. + r1)*(1. + r1) + a2*(a2*r1*(1. + 3.*r1)*cth4 + a4cth4*cth2 + r3*sth2 + r1*cth2*(2.*r1 + 3.*r3 - a2sth2)))*irho23;
	lconn[1][1][2] = -a2*dthdx2*s2th/fac2;
	lconn[1][1][3] = a*sth2*(a4*r1*cth4 + r2*(2*r1 + r3 - a2sth2) + a2cth2*(2.*r1*(-1. + r2) + a2sth2))*irho23;

	lconn[1][2][0] = lconn[1][0][2];
	lconn[1][2][1] = lconn[1][1][2];
	lconn[1][2][2] = -fac3*dthdx22*irho2;
	lconn[1][2][3] = 0.;

	lconn[1][3][0] = lconn[1][0][3];
	lconn[1][3][1] = lconn[1][1][3];
	lconn[1][3][2] = lconn[1][2][3];
	lconn[1][3][3] = -fac3*sth2*(r1*rho22 - a2*fac1*sth2)/(r1*rho23);

	lconn[2][0][0] = -a2*r1*s2th*irho23_dthdx2;
	lconn[2][0][1] = r1*lconn[2][0][0];
	lconn[2][0][2] = 0.;
	lconn[2][0][3] = a*r1*(a2+r2)*s2th*irho23_dthdx2;

	lconn[2][1][0] = lconn[2][0][1];
	lconn[2][1][1] = r2*lconn[2][0][0];
	lconn[2][1][2] = r2*irho2;
	lconn[2][1][3] = (a*r1*cth*sth*(r3*(2. + r1) + a2*(2.*r1*(1. + r1)*cth2 + a2*cth4 + 2*r1sth2)))*irho23_dthdx2;

	lconn[2][2][0] = lconn[2][0][2];
	lconn[2][2][1] = lconn[2][1][2];
	lconn[2][2][2] = -a2*cth*sth*dthdx2*irho2 + d2thdx22/dthdx2;
	lconn[2][2][3] = 0.;

	lconn[2][3][0] = lconn[2][0][3];
	lconn[2][3][1] = lconn[2][1][3];
	lconn[2][3][2] = lconn[2][2][3];
	lconn[2][3][3] = -cth*sth*(rho23 + a2sth2*rho2*(r1*(4. + r1) + a2cth2) + 2.*r1*a4*sth4)*irho23_dthdx2;

	lconn[3][0][0] = a*fac1_rho23;
	lconn[3][0][1] = r1*lconn[3][0][0];
	lconn[3][0][2] = -2.*a*r1*cth*dthdx2/(sth*rho22);
	lconn[3][0][3] = -a2sth2*fac1_rho23;

	lconn[3][1][0] = lconn[3][0][1];
	lconn[3][1][1] = a*r2*fac1_rho23;
	lconn[3][1][2] = -2*a*r1*(a2 + 2*r1*(2. + r1) + a2*c2th)*cth*dthdx2/(sth*fac2*fac2);
	lconn[3][1][3] = r1*(r1*rho22 - a2sth2*fac1)*irho23;

	lconn[3][2][0] = lconn[3][0][2];
	lconn[3][2][1] = lconn[3][1][2];
	lconn[3][2][2] = -a*r1*dthdx22*irho2;
	lconn[3][2][3] = dthdx2*(0.25*fac2*fac2*cth/sth + a2*r1*s2th)*irho22;

	lconn[3][3][0] = lconn[3][0][3];
	lconn[3][3][1] = lconn[3][1][3];
	lconn[3][3][2] = lconn[3][2][3];
	lconn[3][3][3] = (-a*r1sth2*rho22 + a3*sth4*fac1)*irho23;

}

