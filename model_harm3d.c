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

void interp_fourv(double X[NDIM], double ****fourv, double Fourv[NDIM]) ;
double interp_scalar(double X[NDIM], double ***var) ;

void init_model(char *args[])
{
	void init_harm3d_grid(char *);
	void init_harm3d_data(char *);

	fprintf(stderr, "reading data header...\n");
	/* Read in header and allocate space for grid data */
	init_harm3d_grid(args[3]);
	fprintf(stderr, "success\n");

	/* find dimensional quantities from black hole
		mass and its accretion rate */
	set_units(args[4]);

	fprintf(stderr, "reading data...\n");
	/* Read in the grid data */
	init_harm3d_data(args[3]);
	fprintf(stderr, "success\n");

	/* pre-compute densities, field strengths, etc. */
	init_physical_quantities() ;

	/* horizon radius */
	Rh = 1 + sqrt(1. - a * a) ;

}

/* 

	these supply basic model data to grmonty

*/

void get_model_ucov(double X[NDIM], double Ucov[NDIM])
{
	double gcov[NDIM][NDIM];

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


/** HARM utilities **/

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
	double del[NDIM],b1,b2,interp;
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
  double phi;

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
	RHO_unit = M_unit / pow(L_unit, 3);
	U_unit = RHO_unit * CL * CL;
	B_unit = CL * sqrt(4.*M_PI*RHO_unit);

	fprintf(stderr,"L,T,M units: %g [cm] %g [g] %g [sec]\n",L_unit,T_unit,M_unit) ;
	fprintf(stderr,"rho,u,B units: %g [g cm^-3] %g [g cm^-1 sec^-2] %g [G] \n",RHO_unit,U_unit,B_unit) ;
}



void init_physical_quantities(void)
{
	int i, j, k;
        double bsq,Thetae_unit,sigma_m,beta,b2,trat;

	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			for (k = 0; k < N3; k++) {
			  ne[i][j][k] = p[KRHO][i][j][k] * RHO_unit/(MP+ME) ;

			  bsq= bcon[i][j][k][0] * bcov[i][j][k][0] +
			       bcon[i][j][k][1] * bcov[i][j][k][1] +
			       bcon[i][j][k][2] * bcov[i][j][k][2] +
			       bcon[i][j][k][3] * bcov[i][j][k][3] ;

			  b[i][j][k] = sqrt(bsq)*B_unit ;
			  sigma_m=bsq/p[KRHO][i][j][k] ;

			  // beta presciption
			  beta=p[UU][i][j][k]*(gam-1.)/0.5/bsq;
			  b2=pow(beta,2);
			  trat = trat_d * b2/(1. + b2) + trat_j /(1. + b2);
			  Thetae_unit = (gam - 1.) * (MP / ME) / trat;
			  thetae[i][j][k] = (p[UU][i][j][k]/p[KRHO][i][j][k])* Thetae_unit;
			  
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

	return;
}

/* HDF5 v1.6 API */
//#include <H5LT.h>

/* HDF5 v1.8 API */
#include <hdf5.h>
#include <hdf5_hl.h>

/* Harm3d globals */

extern double ****bcon;
extern double ****bcov;
extern double ****ucon;
extern double ****ucov;
extern double ****p;
extern double ***ne;
extern double ***thetae;
extern double ***b;

void init_harm3d_grid(char *fname)
{
	hid_t file_id;
	double th_end,th_cutout;

	file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(file_id < 0){
		fprintf(stderr,"file %s does not exist, aborting...\n",fname);
		exit(1234);
	}

	H5LTread_dataset_int(file_id,		"/Header/Grid/totalsize1", 	&N1);
	H5LTread_dataset_int(file_id,		"/Header/Grid/totalsize2", 	&N2);
	H5LTread_dataset_int(file_id,		"/Header/Grid/totalsize3", 	&N3);
	H5LTread_dataset_double(file_id, 	"/Header/Grid/dx1", 		&(dx[1]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/dx2", 		&(dx[2]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/dx3", 		&(dx[3]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/startx0", 	&(startx[0]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/startx1", 	&(startx[1]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/startx2", 	&(startx[2]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/startx3", 	&(startx[3]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/a",		&a);
	H5LTread_dataset_double(file_id,	"/Header/Grid/gam",		&gam);
	H5LTread_dataset_double(file_id,	"/Header/Grid/R0",		&R0);
	H5LTread_dataset_double(file_id,	"/Header/Grid/Rin",		&Rin);
	H5LTread_dataset_double(file_id, 	"/Header/Grid/Rout",		&Rout);
	H5LTread_dataset_double(file_id, 	"/Header/Grid/h_slope",		&hslope);
	//H5LTread_dataset_double(file_id, 	"/Header/Grid/X1_slope",	&X1slope);
	//H5LTread_dataset_double(file_id,	"/Header/Grid/X1_0",		&X10);
	H5LTread_dataset_double(file_id,	"/Header/Grid/th_cutout",  	&th_cutout);
	//H5LTread_dataset_double(file_id,	"/Header/Grid/th_beg",		&th_beg);
	//H5LTread_dataset_double(file_id,	"/Header/Grid/th_end",		&th_end);

	/* Account for 3 ghost zones */

	
	startx[1] += 3*dx[1];
	startx[2] += 3*dx[2];
        startx[3] += 3*dx[3];
	
	fprintf(stdout,"start: %g %g %g \n",startx[1],startx[2],startx[3]);
	//below is equivalent to the above
	/*
	startx[1] = log(Rin-R0);       
        startx[2] = th_cutout/M_PI ; 
	*/

	fprintf(stdout,"th_cutout: %g  %d x %d x %d\n",th_cutout,N1,N2,N3);
	

	th_beg=th_cutout;
	th_end=M_PI-th_cutout;
	th_len = th_end-th_beg;

	stopx[0] = 1.;
	stopx[1] = startx[1]+N1*dx[1];
	stopx[2] = startx[2]+N2*dx[2];
	stopx[3] = startx[3]+N3*dx[3];

	fprintf(stdout,"stop: %g %g %g \n",stopx[1],stopx[2],stopx[3]);


	init_storage();

	H5Fclose(file_id);
}

void init_harm3d_data(char *fname)
{
	hid_t file_id;
	int i,j,k,l,m;
	double X[NDIM],UdotU,ufac,udotB;
	double gcov[NDIM][NDIM],gcon[NDIM][NDIM], g;
	double dMact, Ladv, MBH;
	double r,th;
	FILE *fp;

	fp = fopen("model_param.dat","r") ;
        if(fp == NULL) {
	  fprintf(stderr,"Can't find model_param.dat\n") ;
	  exit(1) ;
        }
        fscanf(fp,"%lf",&MBH) ;
        fclose(fp) ;
	MBH=MBH*MSUN;
        L_unit = GNEWT * MBH / (CL * CL);
        T_unit = L_unit / CL;

	
	file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(file_id < 0){
		fprintf(stderr,"file %s does not exist, aborting...\n",fname);
		exit(12345);
	}


	//	fprintf(stderr,"data incoming...");
	H5LTread_dataset_double(file_id, "rho",	&(p[KRHO][0][0][0]));
	H5LTread_dataset_double(file_id, "uu",	&(p[UU][0][0][0]));
	H5LTread_dataset_double(file_id, "v1",	&(p[U1][0][0][0]));
	H5LTread_dataset_double(file_id, "v2",	&(p[U2][0][0][0]));
	H5LTread_dataset_double(file_id, "v3",	&(p[U3][0][0][0]));
	H5LTread_dataset_double(file_id, "B1",	&(p[B1][0][0][0]));
	H5LTread_dataset_double(file_id, "B2",	&(p[B2][0][0][0]));
	H5LTread_dataset_double(file_id, "B3",	&(p[B3][0][0][0]));

	H5Fclose(file_id);
	X[0] = 0.;
	X[3] = 0.;

	//fprintf(stderr,"reconstructing 4-vectors...\n");
	dMact = Ladv = 0.;

	//reconstruction of variables at the zone center!
	for(i = 0; i < N1; i++){
	  X[1] = startx[1] + ( i + 0.5)*dx[1];
	  for(j = 0; j < N2; j++){
	    X[2] = startx[2] + (j+0.5)*dx[2];
	    gcov_func(X, gcov); // in system with cut off
	    gcon_func(gcov, gcon);
	    g = gdet_func(gcov);

	    bl_coord(X, &r, &th);

	    for(k = 0; k < N3; k++){
	      UdotU = 0.;
	      X[3] = startx[3] + (k+0.5)*dx[3];
	      
	      //the four-vector reconstruction should have gcov and gcon and gdet using the modified coordinates
	      //interpolating the four vectors to the zone center !!!!
	      for(l = 1; l < NDIM; l++) for(m = 1; m < NDIM; m++) UdotU += gcov[l][m]*p[U1+l-1][i][j][k]*p[U1+m-1][i][j][k];
	      ufac = sqrt(-1./gcon[0][0]*(1 + fabs(UdotU)));
	      ucon[i][j][k][0] = -ufac*gcon[0][0];
	      for(l = 1; l < NDIM; l++) ucon[i][j][k][l] = p[U1+l-1][i][j][k] - ufac*gcon[0][l];
	      lower(ucon[i][j][k], gcov, ucov[i][j][k]);

	      //reconstruct the magnetic field three vectors
	      udotB = 0.;
	      for(l = 1; l < NDIM; l++) udotB += ucov[i][j][k][l]*p[B1+l-1][i][j][k];
	      bcon[i][j][k][0] = udotB;
	      for(l = 1; l < NDIM; l++) bcon[i][j][k][l] = (p[B1+l-1][i][j][k] + ucon[i][j][k][l]*udotB)/ucon[i][j][k][0];
	      lower(bcon[i][j][k], gcov, bcov[i][j][k]);

	      if(i <= 20) dMact += g * p[KRHO][i][j][k] * ucon[i][j][k][1] ;
	      if(i >= 20 && i < 40) Ladv += g * p[UU][i][j][k] * ucon[i][j][k][1] * ucov[i][j][k][0] ;

	    }
	  }
	}
	
	dMact *= dx[3]*dx[2] ;
	dMact /= 21. ;
	Ladv *= dx[3]*dx[2] ;
	Ladv /= 21. ;

	fprintf(stderr,"dMact: %g [code]\n",dMact) ;
	fprintf(stderr,"Ladv: %g [code]\n",Ladv) ;
	fprintf(stderr,"Mdot: %g [g/s] \n",-dMact*M_unit/T_unit) ;
	fprintf(stderr,"Mdot: %g [MSUN/YR] \n",-dMact*M_unit/T_unit/(MSUN / YEAR)) ;
        double Mdotedd=4.*M_PI*GNEWT*MBH*MP/CL/0.1/SIGMA_THOMSON;
        fprintf(stderr,"Mdot: %g [Mdotedd]\n",-dMact*M_unit/T_unit/Mdotedd) ;
        fprintf(stderr,"Mdotedd: %g [g/s]\n",Mdotedd) ;
        fprintf(stderr,"Mdotedd: %g [MSUN/YR]\n",Mdotedd/(MSUN/YEAR)) ;
	
}

double root_find(double th)
{
    int i;
    double X2a, X2b, X2c, tha, thb, thc;
    double dthdX2, dtheta_func(double y), theta_func(double y);

    if (th < M_PI / 2.) {
	X2a = 0. - SMALL;
	X2b = 0.5 + SMALL;
    } else {
	X2a = 0.5 - SMALL;
	X2b = 1. + SMALL;
    }

    tha = theta_func(X2a);
    thb = theta_func(X2b);

    /* bisect for a bit */
    for (i = 0; i < 10; i++) {
	X2c = 0.5 * (X2a + X2b);
	thc = theta_func(X2c);

	if ((thc - th) * (thb - th) < 0.)
	    X2a = X2c;
	else
	    X2b = X2c;
    }

    /* now do a couple of newton-raphson strokes */
    tha = theta_func(X2a);
    for (i = 0; i < 2; i++) {
	dthdX2 = dtheta_func(X2a);
	X2a -= (tha - th) / dthdX2;
	tha = theta_func(X2a);
    }

    return (X2a);
}


/*this does not depend on theta cut-outs there is no squizzing*/
double theta_func(double x)
{
    //2D 
    //return (M_PI * x + 0.5 * (1. - hslope) * sin(2. * M_PI * x));
    //3D new
    // return th_len * x + th_beg +  hslope * sin(2 * M_PI * x);
    //3D run 
    //return (M_PI * x + th_beg);
    return (M_PI * x);
}

double dtheta_func(double x)
{
    //2D 
    //return (M_PI * (1. + (1. - hslope) * cos(2. * M_PI * x)));
    //3D new
    //return th_len + 2. * M_PI * hslope * cos(2 * M_PI * x);
    //3D run 
    return M_PI;
    //return th_len;
}

