#include "decs.h"
#include "defs.h"
#include <omp.h>

#define MAXNSTEP 	50000

#define MAX(a,b) (((a)>(b))?(a):(b))
#define sign(x) (((x) < 0) ? -1 : ((x) > 0))

struct of_traj {
	double dl ;
	double X[NDIM] ;
	double Kcon[NDIM] ;
	double Xhalf[NDIM] ;
	double Kconhalf[NDIM] ;
} traj[MAXNSTEP] ;


extern double image[NX][NY];
double image[NX][NY];
extern double imageS[NX][NY][NDIM];
double imageS[NX][NY][NDIM];


int main(int argc, char *argv[]) 
{
	double X[NDIM],Kcon[NDIM];
	double Xhalf[NDIM],Kconhalf[NDIM];
	double dl,Intensity,tau,dlsubstep;
	double DX,DY,fovx,fovy ;
	double thetacam,phicam,rcam,Xcam[NDIM],Kcam[NDIM],Ucam[NDIM] ;
	double Gcov[NDIM][NDIM] ;
	double Gcon[NDIM][NDIM] ;

	double freq,freqcgs ;
	double Ftot, Dsource ;
	int i,j,k,l,nstep ;
	double Xi[NDIM],Xf[NDIM],Kconi[NDIM],Kconf[NDIM],ji,ki,jf,kf;
	double scale;
	double root_find(double th);
      	double jtau,ktau,RMtau;
	int imax,jmax;
	double Imax,Iavg;//,dOmega;
	double dltot;
	double Stokes_I,Stokes_Q,Stokes_U,Stokes_V;
	double complex N_coord[NDIM][NDIM];
	double gdet;

	if(argc < 6) {
	  // fprintf(stderr,"usage: mibothros theta freq filename Munit theta_j trat_d\n") ;
	  fprintf(stderr,"usage: mibothros theta freq filename Munit trat_j trat_d\n") ;
	  exit(0) ;
	}

	sscanf(argv[1],"%lf",&thetacam) ;
	sscanf(argv[2],"%lf",&freqcgs) ;
	sscanf(argv[4],"%lf",&M_unit) ;
	sscanf(argv[5],"%lf",&trat_j) ;
 	sscanf(argv[6],"%lf",&trat_d) ;

	init_model(argv) ;
	/* normalize frequency to electron rest-mass energy */
	freq = freqcgs*HPL/(ME*CL*CL) ;
	freqcgs1=freqcgs;

	/* initialize local parameters */
	/* fix camera worldline */
	rcam = 100.;
	phicam = 0.0;// -1.5708;//corresponds to +1.57 in grtrans
	Xcam[0] = 0.0 ;
	Xcam[1] = log(rcam);
	Xcam[2] = root_find(thetacam/180.*M_PI);
	Xcam[3] = phicam;

	fprintf(stdout,"cam_th_cal=%g [deg] th_beg=%g th_len=%g a=%g R0=%g hslope=%g\n",(th_beg+th_len*Xcam[2])*180./M_PI,th_beg,th_len,a,R0,hslope);
	
	R0=0.0;
	Ucam[0] = 1. ;
	Ucam[1] = 0. ;
	Ucam[2] = 0. ;
	Ucam[3] = 0. ;
	

	gcov_func(Xcam,Gcov) ; /* uses harm3D coordinates - grid with cut-offs, so
 				the k four-vector is set in the coorindates that cut off the polar regions, hence connections should be cut like that as well...
			       , but what happens when the four vector gets into the forbidden zone? */
	normalize(Ucam,Gcov) ;
	
	/* fix camera field of view */
	/* size of field of view in the plane of the black hole, in units of GM/c^2 */
	DX = 26.0 ;	
	DY = 26.0 ;

	fovx = DX/rcam ;
	fovy = DY/rcam ;

	Dsource = DM87 ;
	scale=(DX*L_unit/NX)*(DY*L_unit/NY)/(Dsource*Dsource)/JY ;
	fprintf(stderr,"scale=%g D=%g cm FOVx=%g FOVy=%g Lunit=%g %g\n",scale,Dsource,DX,DY,L_unit,4*M_PI*Dsource*Dsource/1e23) ;

	for(i=0;i<NX;i++) {
	  //for(i=84;i<85;i++) {
	  fprintf(stderr," %d",i);
#pragma omp parallel for default(none) private(j,k,l,ki,kf,ktau,ji,jf,jtau,nstep,dlsubstep,dl,X,Xhalf,Kcon,Kconhalf,Xi,Xf,Kconi,Kconf,traj,Intensity,tau,dltot,N_coord,Stokes_I,Stokes_Q,Stokes_U,Stokes_V,Gcov,Gcon,gdet,Kcam) shared(i,Xcam,Ucam,fovx,fovy,freq,freqcgs,image,imageS,L_unit,stderr,stdout,th_beg,th_len,hslope) schedule(static,1) 
	  //for(j=130;j<131;j++) {
	          	    for(j=0;j<NY;j++) {

                  init(i,j,Xcam,Ucam,fovx,fovy,X,Kcon) ;

		  /* integrate backwards along trajectory */
		  nstep = 0;
		  
		  while(!stop_backward_integration(X,Kcon,Xcam) ) {

		    /* This stepsize function can be troublesome inside of R = 2M,and should be used cautiously in this region. */
		    dl = stepsize(X,Kcon) ;
		    
		    traj[nstep].dl = dl;
		    for(k=0;k<NDIM;k++) traj[nstep].X[k] = X[k] ;
		    for(k=0;k<NDIM;k++) traj[nstep].Kcon[k] = Kcon[k] ;

		    /* move photon */
		    push_photon(X,Kcon,-dl,Xhalf,Kconhalf) ;

		    nstep++ ;
		    for(k=0;k<NDIM;k++) traj[nstep].Xhalf[k] = Xhalf[k] ;
		    for(k=0;k<NDIM;k++) traj[nstep].Kconhalf[k] = Kconhalf[k] ;

		    if(nstep > MAXNSTEP-2) {
		      fprintf(stderr,"MAXNSTEP exceeded on j=%d i=%d\n",j,i) ;
		      exit(1);
		    }
		    }
		    nstep-- ;
		    /* final step violated the "stop" condition,
				   so don't record it */

		    /* integrate forwards along trajectory, including
		       radiative transfer equation */
		    //init N, Intensity
		    for(l=0;l<NDIM;l++) {
		      Xi[l] = traj[nstep].X[l] ;
		      Kconi[l] = traj[nstep].Kcon[l] ;
		    }
		    init_N(Xi,Kconi,N_coord);
		    Intensity = 0.0 ;


		    while(nstep > 0) {
		      /* find start point emissivity */
		  
		      for(l=0;l<NDIM;l++) {
			Xi[l] = traj[nstep].X[l] ;
			Kconi[l] = traj[nstep].Kcon[l] ;
			Xhalf[l] = traj[nstep].Xhalf[l] ;
			Kconhalf[l] = traj[nstep].Kconhalf[l] ;
		      }
		      get_jkinv(Xi,Kconi,&ji,&ki);
		
			/* loop over substeps.  stepsize really should
			   be set adaptively */
		
		      dlsubstep = traj[nstep-1].dl;

		      	for(l=0;l<NDIM;l++) {
			  Xf[l]=traj[nstep-1].X[l];
			  Kconf[l]=traj[nstep-1].Kcon[l];  
			}

			get_jkinv(Xf,Kconf,&jf,&kf);

			Intensity = approximate_solve(Intensity,ji,ki,jf,kf,dlsubstep*L_unit/freqcgs) ;
		   
			evolve_N(Xi,Kconi,Xhalf,Kconhalf,Xf,Kconf,dlsubstep,N_coord);
		    
			/* swap start and finish */
			ji = jf ;
			ki = kf ;
	    
		      nstep-- ;
		    }
		    
		/* deposit intensity, and Stokes parameters in pixels */
		image[i][j] = Intensity*pow(freqcgs,3);
		project_N(Xf,Kconf,Ucam,N_coord,&Stokes_I,&Stokes_Q,&Stokes_U,&Stokes_V);
		imageS[i][j][0]=Stokes_I*pow(freqcgs,3) ;
		imageS[i][j][1]=Stokes_Q*pow(freqcgs,3) ;
		imageS[i][j][2]=Stokes_U*pow(freqcgs,3) ;
		imageS[i][j][3]=Stokes_V*pow(freqcgs,3) ;

	      }
#pragma omp barrier
	}


	/* printing out to files and on stderr*/
 	scale=(DX*L_unit/NX)*(DY*L_unit/NY)/(Dsource*Dsource)/JY ;
	fprintf(stderr,"\nscale=%g\n",scale) ;
	
	Ftot = 0. ;
	Imax=0.0;
	Iavg=0.0;
	for(i=0;i<NX;i++)
	  for(j=0;j<NY;j++) {
	    //image[i][j] *= scale;
	    Ftot += image[i][j]*scale ;//if this is Jansky per pix^2 then just adding up should give a total flux
	    Iavg+=image[i][j];
	    if(image[i][j]>Imax)
	      {
		imax=i;
		jmax=j;
		Imax=image[i][j];
	      }
	  }

	fprintf(stderr,"imax=%d jmax=%d Imax=%g Iavg=%g\n",imax,jmax,Imax,Iavg/(NX*NY));
	fprintf(stderr,"Ftot: %g %g scale=%g\n",freqcgs,Ftot,scale) ;
        fprintf(stderr,"nuLnu = %g\n",Ftot*Dsource*Dsource*JY*freqcgs) ;

	/* image, dump result */
	make_ppm(image, freq, "mibothros_fnu.ppm") ;
	dump(image,imageS, "mibothros.dat", scale ) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) image[i][j] = log(image[i][j] + 1.e-50) ;
	make_ppm(image, freq, "mibothros_lfnu.ppm") ;

	/* done! */


}


void dump(double image[NX][NY], double imageS[NX][NY][NDIM], char *fname, double scale) 
{
	FILE *fp ;
       	int i,j ;
	//double x,y,dx,dy,xstart,ystart;
	double xx,yy,dum_r;
       	double sum_rm,sum_i;


	fp = fopen(fname,"w") ;
	if(fp == NULL) {
		fprintf(stderr,"unable to open %s\n",fname) ;
		exit(1) ;
	}


	sum_i=0.0;
	for(i=0;i<NX;i++){
	  for(j=0;j<NY;j++){ 
	    sum_i += image[i][j];
	    fprintf(fp,"%d %d %15.10g %15.10g %15.10g %15.10g %15.10g \n",i,j,image[i][j],
                    imageS[i][j][0],imageS[i][j][1],imageS[i][j][2],imageS[i][j][3]) ;
	  }
	    fprintf(fp,"\n") ;
	}
	fclose(fp) ;


	/*dump vtk file*/
        fp = fopen("mibothros.vtk","w") ;
        if(fp == NULL) {
	  fprintf(stderr,"unable to open %s\n","mibothros.vtk") ;
	  exit(1) ;
        }

        fprintf(fp,"# vtk DataFile Version 2.0\n") ;
        fprintf(fp,"Image Simulation Result\n") ;
        fprintf(fp,"ASCII\n") ;
        fprintf(fp,"DATASET STRUCTURED_GRID\n");
        fprintf(fp,"DIMENSIONS %d %d %d\n", NX+1, NY+1,1);
        fprintf(fp,"POINTS %d float\n", (NX+1)*(NY+1));

        for(j=0;j<=NY;j++){
	  for(i=0;i<=NX;i++){
	    xx=i*1.0;//-x_size/2.+i*stepx;                                                                                                                 
	    yy=j*1.0;//-y_size/2.+j*stepy;                                                                                                                 
	    fprintf(fp,"%g %g %g\n",xx,yy,0.0);
	  }
        }
        fprintf(fp,"\nPOINT_DATA %d\n",(NX+1)*(NY+1));
        fprintf(fp,"SCALARS Intensity float\n");
        fprintf(fp,"LOOKUP_TABLE default\n");

        for(j=0;j<=NY;j++){
	  for(i=0;i<=NX;i++){
	    dum_r=(image[i][j]);
	    fprintf(fp,"%g\n",dum_r);
	  }
        }

        fclose(fp);



}


void init(int i, int j, 
	double Xcam[4], double Ucam[4], 
	double fovx, double fovy,	/* field of view, in radians */
	double X[4], double Kcon[4]	/* position, wavevector */
) 
{
	double Gcov[NDIM][NDIM] ;
	double Econ[NDIM][NDIM] ;
	double Ecov[NDIM][NDIM] ;
	double Kcon_tetrad[NDIM] ;
	double trial1[NDIM] ;
	double trial2[NDIM] ;
	int k ;

	/* construct orthonormal tetrad.
		e^0 along Ucam
		e^1 inward along radius vector
		e^2 toward north pole of coordinate system
			("y" for the image plane)
		e^3 in the remaining direction
			("x" for the image plane)
		note: this is *modified* from the make_tetrad 
			routine used in grmonty

		this could easily be modified for a fly-through.
	*/
	/* set up trial vectors */
	trial1[0] = 0. ; trial1[1] = -1. ; trial1[2] = 0. ; trial1[3] = 0. ;
	trial2[0] = 0. ; trial2[1] = 0. ; trial2[2] = -1. ; trial2[3] = 0. ;
	gcov_func(Xcam, Gcov) ;

        make_tetrad(Ucam, trial1, trial2, Gcov, Econ, Ecov) ;

	/* construct *outgoing* wavevectors */
	Kcon_tetrad[0] = 0. ;
	Kcon_tetrad[1] = -1. ;
	Kcon_tetrad[2] = -(j/((double)NY) - 0.5)*fovy ;
	Kcon_tetrad[3] = -(i/((double)NX) - 0.5)*fovx ;
	
	/* normalize */
	null_normalize(Kcon_tetrad, 1.) ;

	/* translate into coordinate frame */
	tetrad_to_coordinate(Econ, Kcon_tetrad, Kcon) ;
       
	/* set position */
	for(k=0;k<NDIM;k++) X[k] = Xcam[k] ;

	/* done! */
}

/* normalize null vector in a tetrad frame */
void null_normalize(double Kcon[NDIM], double fnorm) 
{
	double inorm ;

	inorm = sqrt( Kcon[1]*Kcon[1] + Kcon[2]*Kcon[2] + Kcon[3]*Kcon[3]) ;

	Kcon[0] = fnorm ;
	Kcon[1] *= fnorm/inorm ;
	Kcon[2] *= fnorm/inorm ;
	Kcon[3] *= fnorm/inorm ;

}

/* 
   must be a stable, approximate solution to radiative transfer
   that runs between points w/ initial intensity I, emissivity
   ji, opacity ki, and ends with emissivity jf, opacity kf.

   Return final intensity
*/


double approximate_solve(double Ii,
	double ji,
	double ki,
	double jf,
	double kf,
	double dl)
{
	double efac,If,javg,kavg,dtau ;

	javg = (ji + jf)/2. ;
	kavg = (ki + kf)/2. ;
	
	dtau = dl*kavg ;

	if(dtau < 1.e-3) {
		If = Ii + (javg - Ii*kavg)*dl*(1. - (dtau/2.)*(1. - dtau/3.)) ;
	}
	else {
		efac = exp(-dtau) ;
		If = Ii*efac + (javg/kavg) * (1. - efac) ;
	}

	return(If) ;
}

/* condition for stopping the backwards-in-lambda
   integration of the photon geodesic */

#define LRMAX (log(1.1*Rout))
#define LRMIN (log(1.05*Rh))

int stop_backward_integration(
	double X[NDIM],
	double Kcon[NDIM],
	double Xcam[NDIM])
{

	if(
               (X[1] > LRMAX && Kcon[1] < 0.) || 	/* out far */
		X[1] < LRMIN				/* in deep */
	) return(1) ;
	else return(0) ;                                /* neither out far nor in deep */

}

/* get the invariant emissivity and opacity at a given position
   for a given wavevector */

void get_jkinv(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv)
{
  int ii,jj,kk;
  double del[NDIM];
  double nu,theta,B,Thetae,Ne,Bnuinv;
  double Ucov[NDIM],Bcov[NDIM];
  double Ucon[NDIM],Bcon[NDIM];
  double Kcov[NDIM],gcov[NDIM][NDIM];
  double Be=1;

	/* get fluid parameters */
	Ne = get_model_ne(X) ;	/* check to see if we're outside fluid model */

	if(Ne == 0.) {
		*jnuinv = 0. ;
		*knuinv = 0. ;
		return ;
	}

	/* get covariant four-velocity of fluid for use in get_bk_angle and get_fluid_nu */
	get_model_ucov(X, Ucov);
	get_model_bcov(X, Bcov);

	/*extra: print out stuff to test tetrads*/
	get_model_ucon(X, Ucon);
	get_model_bcon(X, Bcon);


	gcov_func(X,gcov);
	lower(Kcon, gcov, Kcov);

	//theta = M_PI/2.;//get_bk_angle(X,Kcon,Ucov) ;	/* angle between k & b */
	theta = get_bk_angle(X,Kcon,Ucov) ;	/* angle between k & b */
	if(theta <= 0. || theta >= M_PI) {	/* no emission along field */
		*jnuinv = 0. ;
		*knuinv = 0. ;
		return ;
	}

	nu = get_fluid_nu(Kcon,Ucov)*freqcgs1 ;	/* freq in Hz */

	B = get_model_b(X) ;		/* field in G */
	Thetae = get_model_thetae(X) ;	/* temp in e rest-mass units */

	/* assume emission is thermal */
	Bnuinv = Bnu_inv(nu,Thetae) ;
	*jnuinv = jnu_inv(nu,Thetae,Ne,B,theta,Be);

	if(Bnuinv < SMALL)
	  *knuinv = SMALL;
	else
	  *knuinv = *jnuinv/Bnuinv;

	if(isnan(*jnuinv) || isnan(*knuinv)) {
	  fprintf(stderr,"\nisnan get_jkinv\n") ;
	  fprintf(stderr,">> %g %g %g %g %g %g %g %g\n",*jnuinv,*knuinv,Ne,theta,nu,B,Thetae,Bnuinv) ;
	  fprintf(stderr,">>> %g %g",jnu_inv(nu,Thetae,Ne,B,theta,Be),anu_inv(nu,Thetae,Ne,B,theta,Be));
	}


	return ;
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


