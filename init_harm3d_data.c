#include "decs.h"

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
extern double ***Ber;

/*emissivity functions and functions used for Faraday conversion and rotation*/
/*from Dexter PhD thesis (checked with Leung harmony program, and Huang & Shcherbakov 2011*/
double gfun(double Xe){
  return 1.-0.11*log(1+0.035*Xe);
}
double hfun(double Xe){
  return 2.011*exp(-pow(Xe,1.035)/4.7)-cos(Xe*0.5)*exp(-pow(Xe,1.2)/2.73)-0.011*exp(-Xe/47.2);
}


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

	fprintf(stdout,"th_cutout: %g  \n",th_cutout);
	

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
	double BSQ;
	double r,th;
	double Ne,beta,b2,trat,Thetae_unit,Thetae,theta,B,delr,frel,Omega0,omega,S2,Xe,rhov,RM,RMacc,kcon[NDIM],nu,ucovrm[NDIM];

	FILE *fp;
	FILE *fp1;
	FILE *fp2;

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
	RMacc=0.0;
	//null normalized radial wavevector
	kcon[0]=1;
	kcon[1]=1;
	kcon[2]=0;
	kcon[3]=0;

	//reconstruction of variables at the zone center!
	for(i = 0; i < N1; i++){
	  X[1] = startx[1] + ( i + 0.5)*dx[1];
	  for(j = 0; j < N2; j++){
	    X[2] = startx[2] + (j+0.5)*dx[2];
	    gcov_func_rec(X, gcov); // in system with cut off
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

	      Ber[i][j][k] = -(1.+ p[UU][i][j][k]/p[KRHO][i][j][k]*gam)*ucov[i][j][k][0];

	      BSQ=bcon[i][j][k][0]*bcov[i][j][k][0] +
		bcon[i][j][k][1]*bcov[i][j][k][1] +
		bcon[i][j][k][2]*bcov[i][j][k][2] +
		bcon[i][j][k][3]*bcov[i][j][k][3] ;
	      
	      if(i <= 20) dMact += g * p[KRHO][i][j][k] * ucon[i][j][k][1] ;
	      if(i >= 20 && i < 40) Ladv += g * p[UU][i][j][k] * ucon[i][j][k][1] * ucov[i][j][k][0] ;

	      
	      /* chose th and phi angles, 20 degrees (192 points in theta dx[2]) and 0 degrees */
	      /* only within 50 Rg */
	      /*
	      if(j==17 && k==0 ){
		ucovrm[0]=ucov[i][j][k][0];
		ucovrm[1]=ucov[i][j][k][1];
		ucovrm[2]=ucov[i][j][k][2];
		ucovrm[3]=ucov[i][j][k][3];
		Ne= p[KRHO][i][j][k]* RHO_unit/(MP+ME);
		beta=p[UU][i][j][k]*(gam-1.)/0.5/BSQ;
		b2=beta*beta;
		trat = trat_d * b2/(1. + b2) + trat_j /(1. + b2);
		Thetae_unit = (gam - 1.) * (MP / ME) / trat;
		Thetae= (p[UU][i][j][k]/p[KRHO][i][j][k])* Thetae_unit;
		theta=get_bk_angle(X,kcon,ucovrm);
		B=sqrt(BSQ)*B_unit;
		nu = get_fluid_nu(kcon,ucovrm);
		delr=r*log(Rout/Rin)/N1;
		frel=log(Thetae)/pow(Thetae,2)/2. * pow(Thetae,4)/(1+pow(Thetae,4)) + 1./(1.+pow(Thetae,4));
		Omega0=EE*B/ME/CL;
		omega=2.*M_PI*nu;
		S2=1.4142; //sqrt(2)                                                                                                                                                       
  		Xe=Thetae*sqrt(S2*sin(theta)*(1e3*Omega0/omega));
		rhov=Ne*4.*M_PI*EE*EE*Omega0/ME/CL/pow(omega,2)*
		  gsl_sf_bessel_Kn(0,1./Thetae)/gsl_sf_bessel_Kn(2,1./Thetae)*cos(theta)*gfun(Xe);
		RM=frel*Ne*B*cos(theta)*(delr*L_unit*2.6312e-13);
		RMacc=RM+RMacc;		  
		fprintf(stdout," %g %g %g %g %g %g %g %g\n",r,delr,RMacc,Ne,B*cos(theta),Thetae,frel,rhov);
		}*/
	    

	    }
	  }
	}
	
	//	exit(1);
	
	dMact *= dx[3]*dx[2] ;
	dMact /= 21. ;
	Ladv *= dx[3]*dx[2] ;
	Ladv /= 21. ;

	fprintf(stderr,"dMact: %g\n",dMact) ;
	fprintf(stderr,"Ladv: %g\n",Ladv) ;
        fprintf(stderr,"Lunit: %g %g %g\n",L_unit,T_unit,M_unit) ;
        fprintf(stderr,"Mdot: %g x M_unit [MSUN/YR] \n",-dMact / T_unit / (MSUN / YEAR));
	fprintf(stderr,"Mdot: %g [MSUN/YR] \n",-dMact*M_unit/T_unit/(MSUN / YEAR)) ;
        double Medd=4.*M_PI*GNEWT*MBH*MP/CL/0.1/SIGMA_THOMSON;
	fprintf(stderr,"Mdot in Eddington units") ;
        fprintf(stderr,"Mdotedd: %g [g/s]\n",Medd) ;
        fprintf(stderr,"Mdotedd: %g [MSUN/YR]\n",Medd/MSUN*YEAR) ;
        fprintf(stderr,"Mdot: %g [Medd]\n",-dMact*M_unit/T_unit/Medd) ;
	
	/*read in the non-theraml emissivity tables*/
	/*	
        fp1 = fopen("Gtab.3.5.dat", "r");
        if (fp1 == NULL) {
	  fprintf(stderr,"can't open Gtab.dat data file\n");
	  exit(1);
        } else {
	  fprintf(stderr,"successfully opened %s\n","Gtab.dat");
        }
	for(i=0;i<202;i++) {
	  fscanf(fp1,"%lf %lf\n", &xtab[i], &Gtab[i]) ;
	  }
	fclose(fp1);


        fp2 = fopen("Gatab.3.5.dat", "r");
        if (fp2 == NULL) {
	  fprintf(stderr,"can't open Gatab.dat data file\n");
	  exit(1);
        } else {
	  fprintf(stderr,"successfully opened %s\n","Gatab.dat");
        }
	for(i=0;i<202;i++) {
	  fscanf(fp2,"%lf %lf\n", &xatab[i], &Gatab[i]) ;
	  }
	fclose(fp2);
	*/

}

