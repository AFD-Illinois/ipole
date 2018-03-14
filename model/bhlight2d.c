// bhlight2d model specification routines

#include "decs.h"
#include <unistd.h>
#include <stdarg.h>

// Macros
#define ZSLOOP(istart, istop, jstart, jstop) \
  for (int i = istart; i <= istop; i++) \
  for (int j = jstart; j <= jstop; j++)
#define PLOOP for(int k = 0; k < nprim; k++)

// bhlight 2d grid functions
double ***bcon;
double ***bcov;
double ***ucon;
double ***ucov;
double ***p;
double **ne;
double **thetae;
double **b;

void interp_fourv(double X[NDIM], double ***fourv, double Fourv[NDIM]) ;
double interp_scalar(double X[NDIM], double **var) ;
void safe_fscanf(FILE *stream, const char *format, ...) ;

void parse_input(int argc, char *argv[])
{
  if (argc != 5) {
    fprintf(stderr, "ERROR format is\n");
    fprintf(stderr, "  ipole theta[deg] freq[cgs] filename counterjet\n");
    exit(-1);
  }

  sscanf(argv[1], "%lf", &thetacam);
  sscanf(argv[2], "%lf", &freqcgs);
  strcpy(fnam, argv[3]);
  sscanf(argv[4], "%d",  &counterjet);
}

void update_data() {
}

double game, gamp, Thetae_unit;
int WITH_ELECTRONS, nprim;

void init_model(char *args[])
{
	void init_bhlight2d_data(char *);

	fprintf(stderr, "reading data...\n");

  /* Read in header and allocate space for grid data */
	init_bhlight2d_data(args[3]);
	fprintf(stderr, "success\n");

	//fprintf(stderr, "reading data...\n");
	/* Read in the grid data */
	//init_bhlight2d_data(args[3]);
	//fprintf(stderr, "success\n");

	/* pre-compute densities, field strengths, etc. */
	init_physical_quantities() ;

	/* horizon radius */
	Rh = 1 + sqrt(1. - a*a) ;
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

	  // sensible default value
		Ucov[0] = -1./sqrt(-gcov[0][0]);
		Ucov[1] = 0.;
		Ucov[2] = 0.;
		Ucov[3] = 0.;

    return;
	}

	interp_fourv(X, ucov, Ucov);
}

void get_model_ucon(double X[NDIM], double Ucon[NDIM])
{

	double gcov[NDIM][NDIM];
	double gcon[NDIM][NDIM];
	double tmp[NDIM];

	if(X[1] < startx[1] ||
	   X[1] > stopx[1]  ||
	   X[2] < startx[2] ||
	   X[2] > stopx[2]) {
	   	// Sensible default value
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

  double te = interp_scalar(X, thetae);
  if (te < 1.e-10) {
    printf("startx[] = %e %e %e\n", startx[1],startx[2],startx[3]);
    printf("stopx[] = %e %e %e\n", stopx[1],stopx[2],stopx[3]);
    printf("\n\nBAD!! te = %e X = %e %e %e %e\n", te, X[0],X[1],X[2],X[3]);
    exit(-1);
  }

	return interp_scalar(X, thetae);
}

// Magnetic field strength in Gauss
double get_model_b(double X[NDIM])
{

	if(X[1] < startx[1] ||
	   X[1] > stopx[1]  ||
	   X[2] < startx[2] ||
	   X[2] > stopx[2]) {

     return 0.;
	}

	return interp_scalar(X, b);
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

void Xtoij(double X[NDIM], int *i, int *j, double del[NDIM]);

/********************************************************************

				Interpolation routines

 ********************************************************************/

/* return fluid four-vector in simulation units */
void interp_fourv(double X[NDIM], double ***fourv, double Fourv[NDIM]){
	double del[NDIM],b1,b2,/*b3,*/d1,d2,d3,d4;
	int i, j, ip1, jp1;

	/* find the current zone location and offsets del[0], del[1] */
	Xtoij(X, &i, &j, del);

	ip1 = i + 1;
	jp1 = j + 1;

	b1 = 1.-del[1];
	b2 = 1.-del[2];
	//b3 = 1.-del[3];

	d1 = b1*b2;
	d3 = del[1] * b2;
	d2 = del[2] * b1;
	d4 = del[1] * del[2];


	// Interpolate along X1, X2
	Fourv[0] = d1*fourv[i][j][0] + d2*fourv[i][jp1][0] + d3*fourv[ip1][j][0] +
             d4*fourv[ip1][jp1][0];
	Fourv[1] = d1*fourv[i][j][1] + d2*fourv[i][jp1][1] + d3*fourv[ip1][j][1] +
             d4*fourv[ip1][jp1][1];
	Fourv[2] = d1*fourv[i][j][2] + d2*fourv[i][jp1][2] + d3*fourv[ip1][j][2] +
             d4*fourv[ip1][jp1][2];
	Fourv[3] = d1*fourv[i][j][3] + d2*fourv[i][jp1][3] + d3*fourv[ip1][j][3] +
             d4*fourv[ip1][jp1][3];
}

/* return	 scalar in cgs units */
double interp_scalar(double X[NDIM], double **var)
{
	double del[NDIM],b1,b2,interp;
	int i, j, ip1, jp1;

	/* find the current zone location and offsets del[0], del[1] */
	Xtoij(X, &i, &j, del);

	ip1 = i+1;
	jp1 = j+1;

	b1 = 1.-del[1];
	b2 = 1.-del[2];

	// Interpolate in X1, X2
	interp = var[i][j]*b1*b2 +
	  var[i][jp1]*b1*del[2] +
	  var[ip1][j]*del[1]*b2 +
	  var[ip1][jp1]*del[1]*del[2];

  return interp;
}

/***********************************************************************************

					End interpolation routines

 ***********************************************************************************/


void Xtoij(double X[NDIM], int *i, int *j, double del[NDIM])
{
	//give index of a zone - zone index is moved to the grid zone center/
	//to account for the fact that all variables are reconstrucuted at zone centers?
	// BRR: why?
  *i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
	*j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;

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

	return;
}


void coord(int i, int j, double *X)
{

	/* returns zone-centered values for coordinates */
	X[0] = startx[0];
	X[1] = startx[1] + (i + 0.5) * dx[1];
	X[2] = startx[2] + (j + 0.5) * dx[2];
	X[3] = 0.;

	return;
}

void init_physical_quantities(void)
{
	int i, j;
  double bsq, sigma_m;///, beta, b2, trat;

	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
		 ne[i][j] = p[KRHO][i][j]*RHO_unit/(MP + ME);

		 bsq= bcon[i][j][0] * bcov[i][j][0] +
		      bcon[i][j][1] * bcov[i][j][1] +
		      bcon[i][j][2] * bcov[i][j][2] +
		      bcon[i][j][3] * bcov[i][j][3] ;

		 b[i][j] = sqrt(bsq)*B_unit;
		 sigma_m = bsq/p[KRHO][i][j];

		 // beta presciption
		 //beta=p[UU][i][j]*(gam-1.)/0.5/bsq;
		 //b2=pow(beta,2);
		 //trat = trat_d * b2/(1. + b2) + trat_j /(1. + b2);
		 //Thetae_unit = (gam - 1.) * (MP / ME) / trat;
		 //thetae[i][j] = (p[UU][i][j]/p[KRHO][i][j])* Thetae_unit;

		 //strongly magnetized = empty, no shiny spine
		 if(sigma_m > 2.0) ne[i][j]=0.0;
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

void init_storage(void)
{
	bcon = malloc_rank3(N1,N2,NDIM);
	bcov = malloc_rank3(N1,N2,NDIM);
	ucon = malloc_rank3(N1,N2,NDIM);
	ucov = malloc_rank3(N1,N2,NDIM);
	p = (double ***)malloc_rank1(nprim,sizeof(double *));
	for(int i = 0; i < nprim; i++) p[i] = malloc_rank2(N1,N2);
	ne = malloc_rank2(N1,N2);
	thetae = malloc_rank2(N1,N2);
	b = malloc_rank2(N1,N2);
}

/* bhlight2d globals */
void init_bhlight2d_data(char *fname)
{
	int idum;
  double fdum;

  //FILE *fp = fopen(fname, "r");
  FILE *fp = fopen(fnam, "r");
  printf("fname = %s\n", fnam);
  if (fp == NULL) {
    fprintf(stderr, "file %s does not exist, aborting...\n", fnam);
    exit(1234);
  }

  printf("about to read file\n");

  double MBH;
  safe_fscanf(fp, "%lf", &fdum); // t
  safe_fscanf(fp, "%d", &N1);
  safe_fscanf(fp, "%d", &N2);
  safe_fscanf(fp, "%d", &idum); // "N3" = 1
  safe_fscanf(fp, "%lf", &(startx[1]));
  safe_fscanf(fp, "%lf", &(startx[2]));
  safe_fscanf(fp, "%lf", &(startx[3]));
  safe_fscanf(fp, "%lf", &(dx[1]));
  safe_fscanf(fp, "%lf", &(dx[2]));
  safe_fscanf(fp, "%lf", &(dx[3]));
  safe_fscanf(fp, "%lf", &fdum); // tf
  safe_fscanf(fp, "%d", &idum); // nstep
  safe_fscanf(fp, "%lf", &MBH); // in cgs
  safe_fscanf(fp, "%lf", &a);
  safe_fscanf(fp, "%lf", &L_unit);
  safe_fscanf(fp, "%lf", &T_unit);
  safe_fscanf(fp, "%lf", &M_unit);
  safe_fscanf(fp, "%lf", &Thetae_unit);
  safe_fscanf(fp, "%lf", &gam);
  safe_fscanf(fp, "%lf", &game);
  safe_fscanf(fp, "%lf", &gamp);
  safe_fscanf(fp, "%d", &idum); // RADMODEL
  safe_fscanf(fp, "%lf", &fdum); // tp_over_te
  safe_fscanf(fp, "%lf", &fdum); // cour
  safe_fscanf(fp, "%lf", &DTd); // DTd
  safe_fscanf(fp, "%lf", &fdum); // DTl
  safe_fscanf(fp, "%lf", &fdum); // DTi
  safe_fscanf(fp, "%d", &idum); // DTr
  safe_fscanf(fp, "%d", &idum); // dump_cnt
  safe_fscanf(fp, "%d", &idum); // image_cnt
  safe_fscanf(fp, "%d", &idum); // rdump_cnt
  safe_fscanf(fp, "%lf", &fdum); // dt
  safe_fscanf(fp, "%d", &idum); // lim
  safe_fscanf(fp, "%d", &idum); // failed
  safe_fscanf(fp, "%lf", &fdum); // Rin
  safe_fscanf(fp, "%lf", &Rout); // Rout
  safe_fscanf(fp, "%lf", &hslope); // hslope
  safe_fscanf(fp, "%lf", &fdum); // R0
  while ( (fgetc(fp)) != '\n' ) ;
  //safe_fscanf(fp, "%d", &WITH_ELECTRONS); // WITH_ELECTRONS
  //safe_fscanf(fp, "%d", &idum); // SPEC_THETABINS
  //safe_fscanf(fp, "%d", &idum); // SPEC_FREQBINS
  //safe_fscanf(fp, "%lf", &fdum); // SPEC_NUMIN
  //safe_fscanf(fp, "%lf", &fdum); // SPEC_NUMAX
  //safe_fscanf(fp, "%d", &idum); // MONIKA_TPTE

  // Finish setting units
  RHO_unit = M_unit / pow(L_unit, 3);
	U_unit = RHO_unit * CL * CL;
	B_unit = CL * sqrt(4.*M_PI*RHO_unit);
  double MdotEdd = 4.*M_PI*GNEWT*MBH*MP/CL/0.1/SIGMA_THOMSON;

  printf("L T M U B RHO = %e %e %e %e %e %e\n", L_unit, T_unit, M_unit, U_unit, B_unit, RHO_unit);

  nprim = NPRIM;

  /*if (WITH_ELECTRONS == 0) {
    nprim = 8;
  } else {
    nprim = 13;
  }*/

	fprintf(stdout,"start: %g %g %g \n",startx[1],startx[2],startx[3]);

	//stopx[0] = 1.;
	stopx[1] = startx[1] + N1*dx[1];
	stopx[2] = startx[2] + N2*dx[2];
	stopx[3] = startx[3] + N3*dx[3];
  rmax = MIN(50., Rout);
  th_beg = 0.0174;

	fprintf(stdout,"stop: %g %g %g \n",stopx[1],stopx[2],stopx[3]);

  // Allocate memory for fluid data
	init_storage();

  // Read grid data
  ZSLOOP(0, N1 - 1, 0, N2 - 1) {
    safe_fscanf(fp, "%lf", &fdum); // X[1]
    safe_fscanf(fp, "%lf", &fdum); // X[2]
    safe_fscanf(fp, "%lf", &fdum); // r
    safe_fscanf(fp, "%lf", &fdum); // th
    PLOOP safe_fscanf(fp, "%lf", &(p[k][i][j]));
    safe_fscanf(fp, "%lf", &fdum); // divb
    safe_fscanf(fp, "%lf", &ucon[i][j][0]); // ucon[0]
    safe_fscanf(fp, "%lf", &ucon[i][j][1]); // ucon[1]
    safe_fscanf(fp, "%lf", &ucon[i][j][2]); // ucon[2]
    safe_fscanf(fp, "%lf", &ucon[i][j][3]); // ucon[3]
    safe_fscanf(fp, "%lf", &ucov[i][j][0]); // ucov[0]
    safe_fscanf(fp, "%lf", &ucov[i][j][1]); // ucov[1]
    safe_fscanf(fp, "%lf", &ucov[i][j][2]); // ucov[2]
    safe_fscanf(fp, "%lf", &ucov[i][j][3]); // ucov[3]
    safe_fscanf(fp, "%lf", &bcon[i][j][0]); // bcon[0]
    safe_fscanf(fp, "%lf", &bcon[i][j][1]); // bcon[1]
    safe_fscanf(fp, "%lf", &bcon[i][j][2]); // bcon[2]
    safe_fscanf(fp, "%lf", &bcon[i][j][3]); // bcon[3]
    safe_fscanf(fp, "%lf", &bcov[i][j][0]); // bcov[0]
    safe_fscanf(fp, "%lf", &bcov[i][j][1]); // bcov[1]
    safe_fscanf(fp, "%lf", &bcov[i][j][2]); // bcov[2]
    safe_fscanf(fp, "%lf", &bcov[i][j][3]); // bcov[3]
    while ( (fgetc(fp)) != '\n' ) ;
    /*safe_fscanf(fp, "%lf", &fdum); // vmin
    safe_fscanf(fp, "%lf", &fdum); // vmax
    safe_fscanf(fp, "%lf", &fdum); // vmin
    safe_fscanf(fp, "%lf", &fdum); // vmax
    safe_fscanf(fp, "%lf", &fdum); // g
    safe_fscanf(fp, "%lf", &fdum); // <G[0]>
    safe_fscanf(fp, "%lf", &fdum); // <G[1]>
    safe_fscanf(fp, "%lf", &fdum); // <G[2]>
    safe_fscanf(fp, "%lf", &fdum); // <G[3]>
    safe_fscanf(fp, "%lf", &fdum); // qud
    safe_fscanf(fp, "%lf", &fdum); // qvisc
    safe_fscanf(fp, "%lf", &fdum); // qcoul
    safe_fscanf(fp, "%lf", &fdum); // N_esuper
    safe_fscanf(fp, "%lf", &fdum); // N_esuper_electron
    if (WITH_ELECTRONS == 0) {
      safe_fscanf(fp, "%lf", &thetae[i][j]); // Thetae
    } else {
      thetae[i][j] = p[KELCOND][i][j]*pow(p[KRHO][i][j],game-1.)*Thetae_unit;
      safe_fscanf(fp, "%lf", &thetae[i][j]); // Thetae
    }*/
    thetae[i][j] = p[KELCOND][i][j]*pow(p[KRHO][i][j],game-1.)*Thetae_unit;

    // Need a floor on thetae to avoid NANs
    thetae[i][j] = MAX(thetae[i][j], 1.e-4);
  }

  printf("DONE!\n");

  double dMact = 0., Ladv = 0.;
  double r, th, X[NDIM], gcov[NDIM][NDIM], gcon[NDIM][NDIM], g;
  ZSLOOP(0, N1 - 1, 0, N2 - 1) {
    X[1] = startx[1] + (i + 0.5)*dx[1];
    X[2] = startx[2] + (j + 0.5)*dx[2];
    gcov_func(X, gcov);
    gcon_func(gcov, gcon);
    g = gdet_func(gcov);
    bl_coord(X, &r, &th);

    if (i <= 20) dMact += g*p[KRHO][i][j]*ucon[i][j][1];
    if (i >= 20 && i < 40) Ladv += g*p[UU][i][j]*ucon[i][j][1]*ucov[i][j][0];
  }

  dMact *= dx[2]*dx[3]/21.;
  Ladv *= dx[2]*dx[3]/21.;

  fprintf(stdout, "MBH:   %e Msolar\n", MBH/MSUN);
  fprintf(stdout, "a:     %g\n", a);
  fprintf(stdout, "Mdot:  %e g/s\n", -dMact*M_unit/T_unit);
  fprintf(stdout, "Mdot:  %e Msolar/yr\n", -dMact*M_unit/T_unit/(MSUN/YEAR));
  fprintf(stdout, "Mdot:  %e MdotEdd\n", -dMact*M_unit/T_unit/MdotEdd);

  //exit(-1);
}

/*
void init_bhlight2d_data(char *fname)
{
	hid_t file_id;
	int i,j,k,l,m;
	double X[NDIM],UdotU,ufac,udotB;
	double gcov[NDIM][NDIM],gcon[NDIM][NDIM], g;
	double dMact, Ladv, MBH;
	double r,th;
	FILE *fp;

  // Set MBH, L_unit, T_unit

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

}*/

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
    return (M_PI * x + 0.5 * (1. - hslope) * sin(2. * M_PI * x));
    //3D new
    // return th_len * x + th_beg +  hslope * sin(2 * M_PI * x);
    //3D run
    //return (M_PI * x + th_beg);
    //return (M_PI * x);
}

double dtheta_func(double x)
{
    //2D
    return (M_PI * (1. + (1. - hslope) * cos(2. * M_PI * x)));
    //3D new
    //return th_len + 2. * M_PI * hslope * cos(2 * M_PI * x);
    //3D run
    //return M_PI;
    //return th_len;
}

void safe_fscanf(FILE *stream, const char *format, ...)
{
  va_list args;
  va_start(args, format);
  int vfscanfReturn = vfscanf(stream, format, args);
  va_end(args);
  if (vfscanfReturn == -1) {
    fprintf(stderr, "fscanf() call failed! Exiting!\n");
    exit(-1);
  }
}

// COORDINATES
void gcov_func(double X[NDIM], double gcov[NDIM][NDIM])
{
  MUNULOOP gcov[mu][nu] = 0.;

  double sth, cth, s2, rho2;
  double r, th;

  bl_coord(X, &r, &th);

  cth = cos(th);
  sth = sin(th);

  s2 = sth*sth;
  rho2 = r*r + a*a*cth*cth;

  // KS -> MKS transformation
  double tfac, rfac, hfac, pfac;
  tfac = 1.;
  rfac = r - R0;
  hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
  pfac = 1.;

  gcov[0][0] = (-1. + 2.*r/rho2)*tfac*tfac;
  gcov[0][1] = (2.*r/rho2)*tfac*rfac;
  gcov[0][3] = (-2.*a*r*s2/rho2)*tfac*pfac;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = (1. + 2.*r/rho2)*rfac*rfac;
  gcov[1][3] = (-a*s2*(1. + 2.*r/rho2))*rfac*pfac;

  gcov[2][2] = (rho2)*hfac*hfac;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = (s2*(rho2 + a*a*s2*(1. + 2.*r/rho2)))*pfac*pfac;
}

void bl_coord(double *X, double *r, double *th)
{
  *r = exp(X[1]) + R0;
  *th = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
}

void get_connection(double X[4], double lconn[4][4][4])
{
  get_connection_num(X, lconn);
}

int radiating_region(double X[NDIM])
{
  if (X[1] < log(rmax) && X[2]>th_beg/M_PI && X[2]<(1.-th_beg/M_PI) ) {
    return 1;
  } else {
    return 0;
  }
}

