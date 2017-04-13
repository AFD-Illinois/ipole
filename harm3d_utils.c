
#include "decs.h"

/** HARM utilities **/

extern double ****bcon;
extern double ****bcov;
extern double ****ucon;
extern double ****ucov;
extern double ****p;
extern double ***ne;
extern double ***thetae;
extern double ***b;
extern double ***Ber;

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]) ;

/********************************************************************

				Interpolation routines

 ********************************************************************/

/* return fluid four-vector in simulation units */
void interp_fourv(double X[NDIM], double ****fourv, double Fourv[NDIM]){
	double del[NDIM],b1,b2,b3,d1,d2,d3,d4;
	int i, j, k, ip1, jp1, kp1;

	/* find the current zone location and offsets del[0], del[1] */
	Xtoijk(X, &i, &j, &k, del);

	ip1 = i + 1;
	jp1 = j + 1;
	kp1 = k + 1;
	
	//conditions at the x3 periodic boundary (at the last active zone)
	if(k==(N3-1)) kp1=0;   

	b1 = 1.-del[1];
	b2 = 1.-del[2];
	b3 = 1.-del[3];

	d1 = b1*b2;
	d3 = del[1] * b2;
	d2 = del[2] * b1;
	d4 = del[1] * del[2];


	/* Interpolate along x1,x2 first */
	Fourv[0] = d1*fourv[i][j][k][0] + d2*fourv[i][jp1][k][0] + d3*fourv[ip1][j][k][0] + d4*fourv[ip1][jp1][k][0];
	Fourv[1] = d1*fourv[i][j][k][1] + d2*fourv[i][jp1][k][1] + d3*fourv[ip1][j][k][1] + d4*fourv[ip1][jp1][k][1];
	Fourv[2] = d1*fourv[i][j][k][2] + d2*fourv[i][jp1][k][2] + d3*fourv[ip1][j][k][2] + d4*fourv[ip1][jp1][k][2];
	Fourv[3] = d1*fourv[i][j][k][3] + d2*fourv[i][jp1][k][3] + d3*fourv[ip1][j][k][3] + d4*fourv[ip1][jp1][k][3];

	/* Now interpolate above in x3 */
	Fourv[0] = b3*Fourv[0] + del[3]*(d1*fourv[i][j][kp1][0] + d2*fourv[i][jp1][kp1][0] + d3*fourv[ip1][j][kp1][0] + d4*fourv[ip1][jp1][kp1][0]);
	Fourv[1] = b3*Fourv[1] + del[3]*(d1*fourv[i][j][kp1][1] + d2*fourv[i][jp1][kp1][1] + d3*fourv[ip1][j][kp1][1] + d4*fourv[ip1][jp1][kp1][1]);
	Fourv[2] = b3*Fourv[2] + del[3]*(d1*fourv[i][j][kp1][2] + d2*fourv[i][jp1][kp1][2] + d3*fourv[ip1][j][kp1][2] + d4*fourv[ip1][jp1][kp1][2]);
	Fourv[3] = b3*Fourv[3] + del[3]*(d1*fourv[i][j][kp1][3] + d2*fourv[i][jp1][kp1][3] + d3*fourv[ip1][j][kp1][3] + d4*fourv[ip1][jp1][kp1][3]);
	//new

	//no interpolation of vectors at all
	/*
	Fourv[0]=fourv[i][j][k][0];
	Fourv[1]=fourv[i][j][k][1];
	Fourv[2]=fourv[i][j][k][2];
	Fourv[3]=fourv[i][j][k][3];
	*/
}

/* return	 scalar in cgs units */
double interp_scalar(double X[NDIM], double ***var)
{
	double del[NDIM],b1,b2,i1,i2,i3,i4,i5,i6,i7,i8,interp;
	int i, j, k, ip1, jp1, kp1;

	/* find the current zone location and offsets del[0], del[1] */
	Xtoijk(X, &i, &j, &k, del);

	ip1 = i+1;
	jp1 = j+1;
	kp1 = k+1;
	if(k==(N3-1)) kp1=0;   	

	b1 = 1.-del[1];
	b2 = 1.-del[2];

	/* Interpolate in x1,x2 first */

	interp = var[i][j][k]*b1*b2 + 
	  var[i][jp1][k]*b1*del[2] + 
	  var[ip1][j][k]*del[1]*b2 + 
	  var[ip1][jp1][k]*del[1]*del[2];


	/* Now interpolate above in x3 */

	interp = (1.-del[3])*interp + 
     	  del[3]*(var[i  ][j  ][kp1]*b1*b2 +
		  var[i  ][jp1][kp1]*del[2]*b1 +
		  var[ip1][j  ][kp1]*del[1]*b2 +
		  var[ip1][jp1][kp1]*del[1]*del[2]);
	
	//new, no interpolations what so ever
	//	interp=var[i][j][k];
	/* use bilinear interpolation to find rho; piecewise constant
	   near the boundaries */
	
	return(interp);

}

/***********************************************************************************

					End interpolation routines

 ***********************************************************************************/


void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM])
{
  double phi,b3;

	/* Map X[3] into sim range, assume startx[3] = 0 */
	phi = fmod(X[3], stopx[3]);
	//fold it to be positive and find index
	if(phi < 0.0) phi = stopx[3]+phi;
	
	//give index of a zone - zone index is moved to the grid zone center/
	//to account for the fact that all variables are reconstrucuted at zone centers?
	*i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
	*j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
	*k = (int) ((phi  - startx[3]) / dx[3] - 0.5 + 1000) - 1000;	

	//this makes sense, interpolate with outflow condition
	if(*i < 0) {
	  *i = 0 ;
	  del[1] = 0. ;
	}
	else if(*i > N1-2) { //OK because del[1]=1 and only terms with ip1=N1-1 will be important in interpolation
	  *i = N1-2 ;
	  del[1] = 1. ;
	}
	else {
	  del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
	}
	
	if(*j < 0) {
          *j = 0 ;
          del[2] = 0. ;
        }
        else if(*j > N2-2) {
          *j = N2-2 ;
          del[2] = 1. ;
        }
        else {
          del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
        }

        if(*k < 0) {
          *k = 0 ;
          del[3] = 0. ;
        }
        else if(*k > N3-1) {
          *k = N3-1;
          del[3] = 1. ;
        }
        else {
          del[3] = (phi - ((*k + 0.5) * dx[3] + startx[3])) / dx[3];
        }

	return;
}

//#define SINGSMALL (1.E-20)
/* return boyer-lindquist coordinate of point */
void bl_coord(double *X, double *r, double *th)
{

	// for cfg data files
	//*r = Rin * exp(X[1]);
	// for scn data files
	*r = exp(X[1]) + R0 ;
	//*th = th_beg + th_len *X[2] + hslope*sin(2. * M_PI * X[2]);  //2D
	//hotaka run, hslope=0?
	//*th = M_PI * X[2] + hslope*sin(2. * M_PI * X[2]);
	*th = th_beg + M_PI*X[2]  ;

	//fix coord
	/*
	if (fabs(*th) < SINGSMALL) {
	  if ((*th) >= 0)
	    *th = SINGSMALL;
	  if ((*th) < 0)
	    *th = -SINGSMALL;
        }
        if (fabs(M_PI - (*th)) < SINGSMALL) {
	  if ((*th) >= M_PI)
	    *th = M_PI + SINGSMALL;
	  if ((*th) < M_PI)
	    *th = M_PI - SINGSMALL;
        }
		
	*/
	return;
}

void coord(int i, int j, int k, double *X)
{

	/* returns zone-centered values for coordinates */
	X[0] = startx[0];
	X[1] = startx[1] + (i + 0.5) * dx[1];
	X[2] = startx[2] + (j + 0.5) * dx[2];
	X[3] = startx[3] + (k + 0.5) * dx[3];

	return;
}


void set_units(char *munitstr)
{
	double MBH ;
	FILE *fp ;

	fp = fopen("model_param.dat","r") ;
	if(fp == NULL) {
		fprintf(stderr,"Can't find model_param.dat\n") ;
		exit(1) ;
	}
	fscanf(fp,"%lf",&MBH) ;
	fclose(fp) ;

	sscanf(munitstr,"%lf",&M_unit) ;

	/** input parameters appropriate to Sgr A* **/
	MBH *= MSUN ;

	/** from this, calculate units of length, time, mass,
	    and derivative units **/
	L_unit = GNEWT * MBH / (CL * CL);
	T_unit = L_unit / CL;
	fprintf(stderr,"UNITS\n") ;
	fprintf(stderr,"L,T,M: %g %g %g\n",L_unit,T_unit,M_unit) ;

	RHO_unit = M_unit / pow(L_unit, 3);
	U_unit = RHO_unit * CL * CL;
	B_unit = CL * sqrt(4.*M_PI*RHO_unit);
	fprintf(stderr,"rho,u,B: %g %g %g\n",RHO_unit,U_unit,B_unit) ;

}



void init_physical_quantities(void)
{
	int i, j, k;
	double gcov[NDIM][NDIM] ;
        double gcon[NDIM][NDIM] ;
        double X[NDIM],Be,lor,betaf,bsq;
	double THETAE_MAX=500.;
	double r,th,thmax,thmin,two_temp_gam,Thetae_unit,sigma_m,beta,b2,beta_trans;
	double trat;
	double Kin,Res,Epoynting;

	for (i = 0; i < N1; i++) {
	  X[1] = startx[1] + ( i + 0.5)*dx[1];
		for (j = 0; j < N2; j++) {
		  X[2] = startx[2] + (j+0.5)*dx[2];
			for (k = 0; k < N3; k++) {
			  ne[i][j][k] = p[KRHO][i][j][k] * RHO_unit/(MP+ME) ;

			  bsq= bcon[i][j][k][0] * bcov[i][j][k][0] +
			       bcon[i][j][k][1] * bcov[i][j][k][1] +
			       bcon[i][j][k][2] * bcov[i][j][k][2] +
			       bcon[i][j][k][3] * bcov[i][j][k][3] ;

			  b[i][j][k] = sqrt(bsq)*B_unit ;
			  sigma_m=bsq/p[KRHO][i][j][k] ;

			  //Hydro-bernoulli parameter
			  Be = -(1.+ p[UU][i][j][k]/p[KRHO][i][j][k]*gam)*ucov[i][j][k][0]; //>1 
			  Ber[i][j][k]=Be;                           
			  
			  //for large scale jet old prescription
			  //two_temp_gam = 0.5 * ((1. + 2. / 3. * (trat_d + 1.) / (trat_d + 2.)) + gam);             
			  //Thetae_unit = (two_temp_gam - 1.) * (MP / ME) / (1. + trat_d);  
			  //simplified

			  /*
			  Thetae_unit = (gam - 1.) * (MP / ME) / trat_d;
			  thetae[i][j][k] = (p[UU][i][j][k]/p[KRHO][i][j][k])* Thetae_unit;
			  if(Be>=1.02) thetae[i][j][k] = theta_j;
			  */

			  // beta presciption
			  
			  beta=p[UU][i][j][k]*(gam-1.)/0.5/bsq;
			  b2=pow(beta,2);
			  trat = trat_d * b2/(1. + b2) + trat_j /(1. + b2);
			  Thetae_unit = (gam - 1.) * (MP / ME) / trat;
			  thetae[i][j][k] = (p[UU][i][j][k]/p[KRHO][i][j][k])* Thetae_unit;
			  

			  /*p=kappa-1*/
			  Ber[i][j][k] = 7.5 * b2/(1. + b2) + 2.5 /(1. + b2);
			  
		  
			  //strongly magnetized = empty, no shiny spine
			  if(sigma_m > 2.0) ne[i][j][k]=0.0;

			}
		}
	}

	return ;
}


void *malloc_rank1(int n1, int size)
{
	void *A;

	if((A = malloc(n1*size)) == NULL){
		fprintf(stderr,"malloc failure in malloc_rank1\n");
		exit(123);
	}

	return A;
}

double **malloc_rank2(int n1, int n2)
{

	double **A;
	double *space;
	int i;

	space = malloc_rank1(n1*n2, sizeof(double));

	A = malloc_rank1(n1, sizeof(double *));

	for(i = 0; i < n1; i++) A[i] = &(space[i*n2]);

	return A;
}


double ***malloc_rank3(int n1, int n2, int n3)
{

	double ***A;
	double *space;
	int i,j;

	space = malloc_rank1(n1*n2*n3, sizeof(double));

	A = malloc_rank1(n1, sizeof(double *));

	for(i = 0; i < n1; i++){
		A[i] = malloc_rank1(n2,sizeof(double *));
		for(j = 0; j < n2; j++){
			A[i][j] = &(space[n3*(j + n2*i)]);
		}
	}

	return A;
}


double ****malloc_rank4(int n1, int n2, int n3, int n4)
{

	double ****A;
	double *space;
	int i,j,k;

	space = malloc_rank1(n1*n2*n3*n4, sizeof(double));

	A = malloc_rank1(n1, sizeof(double *));
	
	for(i=0;i<n1;i++){
		A[i] = malloc_rank1(n2,sizeof(double *));
		for(j=0;j<n2;j++){
			A[i][j] = malloc_rank1(n3,sizeof(double *));
			for(k=0;k<n3;k++){
				A[i][j][k] = &(space[n4*(k + n3*(j + n2*i))]);
			}
		}
	}

	return A;
}

double *****malloc_rank5(int n1, int n2, int n3, int n4, int n5)
{

	double *****A;
	double *space;
	int i,j,k,l;

	space = malloc_rank1(n1*n2*n3*n4*n5, sizeof(double));

	A = malloc_rank1(n1, sizeof(double *));
	
	for(i=0;i<n1;i++){
		A[i] = malloc_rank1(n2, sizeof(double *));
		for(j=0;j<n2;j++){
			A[i][j] = malloc_rank1(n3, sizeof(double *));
			for(k=0;k<n3;k++){
				A[i][j][k] = malloc_rank1(n4, sizeof(double *));
				for(l=0;l<n4;l++){
					A[i][j][k][l] = &(space[n5*(l + n4*(k + n3*(j + n2*i)))]);
				}
			}
		}
	}

	return A;
}

void init_storage(void)
{
	int i;

	bcon = malloc_rank4(N1,N2,N3,NDIM);
	bcov = malloc_rank4(N1,N2,N3,NDIM);
	ucon = malloc_rank4(N1,N2,N3,NDIM);
	ucov = malloc_rank4(N1,N2,N3,NDIM);
	p = (double ****)malloc_rank1(NPRIM,sizeof(double *));
	for(i = 0; i < NPRIM; i++) p[i] = malloc_rank3(N1,N2,N3);
	ne = malloc_rank3(N1,N2,N3);
	thetae = malloc_rank3(N1,N2,N3);
	b = malloc_rank3(N1,N2,N3);
	Ber = malloc_rank3(N1,N2,N3);

	return;
}


