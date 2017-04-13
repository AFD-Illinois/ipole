/*
	HARM model specification routines 
*/

#include "decs.h"

/*
	HARM 3d grid functions
*/
double ****bcon;
double ****bcov;
double ****ucon;
double ****ucov;
double ****p;
double ***ne;
double ***thetae;
double ***b;
double ***Ber;

void interp_fourv(double X[NDIM], double ****fourv, double Fourv[NDIM]) ;
double interp_scalar(double X[NDIM], double ***var) ;

/* 
   Current metric: modified Kerr-Schild, squashed in theta
   to give higher resolution at the equator 
*/

/* mnemonics for dimensional indices */
#define TT      0
#define RR      1
#define TH      2
#define PH      3

void gcov_func(double *X, double gcov[][NDIM])
{
	int k, l;
	double sth, cth, s2, rho2;
	double r, th;
	double tfac, rfac, hfac, pfac;

	DLOOP gcov[k][l] = 0.;
	//	bl_coord(X, &r, &th);

        r = exp(X[1]) + R0 ;
        th =  M_PI*X[2] ;

	cth = cos(th);
	sth = fabs(sin(th));
	if (sth < SMALL)
	  sth = SMALL;

	s2 = sth * sth;
	rho2 = r * r + a * a * cth * cth;

	/* transformation for Kerr-Schild -> modified Kerr-Schild */
	tfac = 1.;
	rfac = r - R0;
	//hfac = th_len;
	hfac = M_PI;
	//hfac = th_len + hslope*2.*M_PI*cos(2. * M_PI * X[2]);
	pfac = 1.;

	gcov[TT][TT] = (-1. + 2. * r / rho2) * tfac * tfac;
	gcov[TT][1] = (2. * r / rho2) * tfac * rfac;
	gcov[TT][3] = (-2. * a * r * s2 / rho2) * tfac * pfac;

	gcov[1][TT] = gcov[TT][1];
	gcov[1][1] = (1. + 2. * r / rho2) * rfac * rfac;
	gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2)) * rfac * pfac;

	gcov[2][2] = rho2 * hfac * hfac;

	gcov[3][TT] = gcov[TT][3];
	gcov[3][1] = gcov[1][3];
	gcov[3][3] =
	    s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2)) * pfac * pfac;
}



void gcov_func_rec(double *X, double gcov[][NDIM])
{
	int k, l;
	double sth, cth, s2, rho2;
	double r, th;
	double tfac, rfac, hfac, pfac;

	DLOOP gcov[k][l] = 0.;

	// for scn data files                                                                                                                                                      
        r = exp(X[1]) + R0 ;
        th = th_beg + th_len *X[2] ;// + hslope*sin(2. * M_PI * X[2]);                                                                                                         
               
	cth = cos(th);
	sth = fabs(sin(th));
	if (sth < SMALL)
	  sth = SMALL;

	s2 = sth * sth;
	rho2 = r * r + a * a * cth * cth;

	/* transformation for Kerr-Schild -> modified Kerr-Schild */
	tfac = 1.;
	rfac = r - R0;
	hfac = th_len;
	//hfac = th_len + hslope*2.*M_PI*cos(2. * M_PI * X[2]);
 	pfac = 1.;

	gcov[TT][TT] = (-1. + 2. * r / rho2) * tfac * tfac;
	gcov[TT][1] = (2. * r / rho2) * tfac * rfac;
	gcov[TT][3] = (-2. * a * r * s2 / rho2) * tfac * pfac;

	gcov[1][TT] = gcov[TT][1];
	gcov[1][1] = (1. + 2. * r / rho2) * rfac * rfac;
	gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2)) * rfac * pfac;

	gcov[2][2] = rho2 * hfac * hfac;

	gcov[3][TT] = gcov[TT][3];
	gcov[3][1] = gcov[1][3];
	gcov[3][3] =
	    s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2)) * pfac * pfac;
}




#undef TT
#undef RR
#undef TH
#undef PH

/* 

   connection calculated analytically for modified Kerr-Schild
   	coordinates 


   this gives the connection coefficient
	\Gamma^{i}_{j,k} = conn[..][i][j][k]
   where i = {1,2,3,4} corresponds to, e.g., {t,ln(r),theta,phi}
*/

void get_connection(double X[4], double lconn[4][4][4])
{

  
	double r1,r2,r3,r4,sx,cx;
	double th,dthdx2,dthdx22,d2thdx22,sth,cth,sth2,cth2,sth4,cth4,s2th,c2th,c4th;
	double a2,a3,a4,rho2,irho2,rho22,irho22,rho23,irho23,irho23_dthdx2;
	double fac1,fac1_rho23,fac2,fac3,a2cth2,a2sth2,r1sth2,a4cth4;

	r1 = exp(X[1]);
	r2 = r1*r1;
	r3 = r2*r1;
	r4 = r3*r1;

	sincos(2.*M_PI*X[2], &sx, &cx);

	/* HARM-2D MKS */
	//th = M_PI*X[2] + 0.5*(1-hslope)*sx;
	//dthdx2 = M_PI*(1.+(1-hslope)*cx);
	//d2thdx22 = -2.*M_PI*M_PI*(1-hslope)*sx;

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
  
	th =  M_PI*X[2];
	dthdx2 = M_PI;
	d2thdx22 = 0.0; 

	dthdx22 = dthdx2*dthdx2;

	sincos(th, &sth, &cth);
	sth2 = sth*sth;
	r1sth2 = r1*sth2;
	sth4 = sth2*sth2;
	cth2 = cth*cth;
	cth4 = cth2*cth2;
	s2th = 2.*sth*cth;
	c2th = 2*cth2 - 1.;
	c4th = cth4 - 6.*cth2*sth2 + sth4;

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

/* Sets the spatial discretization in numerical derivatives : */
#define DEL 1.e-7

void conn_func(double X[NDIM], double conn[NDIM][NDIM][NDIM])
{
  int i, j, k, l;
  double tmp[NDIM][NDIM][NDIM];
  double Xh[NDIM], Xl[NDIM];
  double gcon[NDIM][NDIM];
  double gcov[NDIM][NDIM];
  double gh[NDIM][NDIM];
  double gl[NDIM][NDIM];
  double delta,deltah,deltal;

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  for (k = 0; k < NDIM; k++) {
    for (l = 0; l < NDIM; l++)   Xh[l] = X[l];
    for (l = 0; l < NDIM; l++)   Xl[l] = X[l];
    Xh[k] += DEL;
    Xl[k] -= DEL;
    gcov_func(Xh, gh);
    gcov_func(Xl, gl);

    for (i = 0; i < NDIM; i++){
      for (j = 0; j < NDIM; j++){
        conn[i][j][k] =  (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
      }
    }
  }

  /* now rearrange to find \Gamma_{ijk} */
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      for (k = 0; k < NDIM; k++)
        tmp[i][j][k] =  0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);


  /* finally, raise index */
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      for (k = 0; k < NDIM; k++) {
        conn[i][j][k] = 0.;
        for (l = 0; l < NDIM; l++)   conn[i][j][k] += gcon[i][l] * tmp[l][j][k];
      }


  /* done! */

}
#undef DEL


#define EPS	0.03

double stepsize(double X[NDIM], double Kcon[NDIM])
{
	int i;
	double dl, dlx1, dlx2, dlx3;
	double idlx1,idlx2,idlx3 ;

	dlx1 = EPS / (fabs(Kcon[1]) + SMALL*SMALL) ;
	dlx2 = EPS * GSL_MIN(X[2], 1. - X[2]) / (fabs(Kcon[2]) + SMALL*SMALL) ;
	dlx3 = EPS / (fabs(Kcon[3]) + SMALL*SMALL) ;

	idlx1 = 1./(fabs(dlx1) + SMALL*SMALL) ;
	idlx2 = 1./(fabs(dlx2) + SMALL*SMALL) ;
	idlx3 = 1./(fabs(dlx3) + SMALL*SMALL) ;

	dl = 1. / (idlx1 + idlx2 + idlx3) ;

	return (dl);
}


void init_model(char *args[])
{
	int i,j,k,l;
	double X[NDIM],del[NDIM];
	void init_harm3d_grid(char *);
	void init_harm3d_data(char *);

	fprintf(stderr, "\nreading data header...");
	/* Read in header and allocate space for grid data */
	init_harm3d_grid(args[3]);
	fprintf(stderr, "success\n");


	/* find dimensional quantities from black hole
		mass and its accretion rate */
	set_units(args[4]);

	fprintf(stderr, "reading data...");
	/* Read in the grid data */
	init_harm3d_data(args[3]);
	fprintf(stderr, "success\n");


	/* pre-compute densities, field strengths, etc. */
	init_physical_quantities() ;

	Rh = 1 + sqrt(1. - a * a) ;

}

/* 

	these supply basic model data to grmonty

*/

void get_model_ucov(double X[NDIM], double Ucov[NDIM])
{
	double gcov[NDIM][NDIM], Ucon[NDIM] ;

	gcov_func(X, gcov);

	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   
	   	/* sensible default value */
		Ucov[0] = -1./sqrt(-gcov[0][0]) ;
		Ucov[1] = 0. ;
		Ucov[2] = 0. ;
		Ucov[3] = 0. ;

		return ;
	}

	//get_model_ucon(X, Ucon);
	//lower(Ucon, gcov, Ucov);

	interp_fourv(X, ucov, Ucov) ;

}

void get_model_ucon(double X[NDIM], double Ucon[NDIM])
{

	double gcov[NDIM][NDIM] ;
	double gcon[NDIM][NDIM] ;
	double tmp[NDIM] ;

	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   	/* sensible default value */
	   	gcov_func(X, gcov) ;

		tmp[0] = -1./sqrt(-gcov[0][0]) ;
		tmp[1] = 0. ;
		tmp[2] = 0. ;
		tmp[3] = 0. ;

	   	gcon_func(gcov, gcon) ;
		Ucon[0] = 
			tmp[0]*gcon[0][0] +
			tmp[1]*gcon[0][1] +
			tmp[2]*gcon[0][2] +
			tmp[3]*gcon[0][3] ;
		Ucon[1] = 
			tmp[0]*gcon[1][0] +
			tmp[1]*gcon[1][1] +
			tmp[2]*gcon[1][2] +
			tmp[3]*gcon[1][3] ;
		Ucon[2] = 
			tmp[0]*gcon[2][0] +
			tmp[1]*gcon[2][1] +
			tmp[2]*gcon[2][2] +
			tmp[3]*gcon[2][3] ;
		Ucon[3] = 
			tmp[0]*gcon[3][0] +
			tmp[1]*gcon[3][1] +
			tmp[2]*gcon[3][2] +
			tmp[3]*gcon[3][3] ;
	
		return ;
	}
	   
	interp_fourv(X, ucon, Ucon) ;
}

void get_model_bcov(double X[NDIM], double Bcov[NDIM])
{
	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {

	   	Bcov[0] = 0. ;
	   	Bcov[1] = 0. ;
	   	Bcov[2] = 0. ;
	   	Bcov[3] = 0. ;

		return ;
	}
	interp_fourv(X, bcov, Bcov) ;
}

void get_model_bcon(double X[NDIM], double Bcon[NDIM])
{
	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {

	   	Bcon[0] = 0. ;
	   	Bcon[1] = 0. ;
	   	Bcon[2] = 0. ;
	   	Bcon[3] = 0. ;

		return ;
	}
	interp_fourv(X, bcon, Bcon) ;
}

double get_model_thetae(double X[NDIM])
{
	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   	return(0.) ;
	}
	
	return(interp_scalar(X, thetae)) ;
}

//b field strength in Gauss
double get_model_b(double X[NDIM])
{
  int i,j,k,del[3];

	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   	return(0.) ;
	}

	return(interp_scalar(X, b)) ;
	

}

double get_model_ne(double X[NDIM])
{
	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   	return(0.) ;
	}
	
	return(interp_scalar(X, ne)) ;
}



double get_model_Ber(double X[NDIM])
{
	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   	return(0.) ;
	}
	
	return(interp_scalar(X, Ber)) ;
}

