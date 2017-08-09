#include "decs.h"
#include "defs.h"
#include <omp.h>

#define MAXNSTEP 	50000

struct of_traj {
    double dl;
    double X[NDIM];
    double Kcon[NDIM];
    double Xhalf[NDIM];
    double Kconhalf[NDIM];
} traj[MAXNSTEP];


extern double image[NX][NY];
double image[NX][NY];
extern double imageS[NX][NY][NDIM];
double imageS[NX][NY][NDIM];

int threadid,nthreads;

int main(int argc, char *argv[])
{
    double X[NDIM], Kcon[NDIM];
    double Xhalf[NDIM], Kconhalf[NDIM];
    double dl, Intensity;
    double DX, DY, fovx, fovy;
    double thetacam, phicam, rcam, Xcam[NDIM];

    double freq, freqcgs;
    double Ftot, Dsource;
    int i, j, k, l, nstep;
    double Xi[NDIM], Xf[NDIM], Kconi[NDIM], Kconf[NDIM], ji, ki, jf, kf;
    double scale;
    double root_find(double th);
    int imax, jmax;
    double Imax, Iavg;		//,dOmega;
    double Stokes_I, Stokes_Q, Stokes_U, Stokes_V;
    double complex N_coord[NDIM][NDIM];

#ifdef _OPENMP
#pragma omp parallel default(none) shared(nthreads) private(threadid)
    {
        threadid = omp_get_thread_num();
        printf("tid = %d\n", threadid);
        if(threadid==0) {
            nthreads = omp_get_num_threads();
            printf("nthreads = %d\n",nthreads);
        }
    }
#else
    printf("_OPENMP not defined: running in serial\n");
#endif

    if (argc < 6) {
	// fprintf(stderr,"usage: ipole theta freq filename Munit theta_j trat_d\n") ;
	fprintf(stderr,
		"usage: ipole theta freq filename Munit trat_j trat_d\n");
	exit(0);
    }

    sscanf(argv[1], "%lf", &thetacam);
    sscanf(argv[2], "%lf", &freqcgs);
    sscanf(argv[4], "%lf", &M_unit);
    sscanf(argv[5], "%lf", &trat_j);
    sscanf(argv[6], "%lf", &trat_d);

    init_model(argv);
    R0 = 0.0;

    /* normalize frequency to electron rest-mass energy */
    freq = freqcgs * HPL / (ME * CL * CL);

    /* fix camera location */
    rcam = 240.;
    phicam = 0.0;	
    Xcam[0] = 0.0;
    Xcam[1] = log(rcam);
    Xcam[2] = root_find(thetacam / 180. * M_PI);
    Xcam[3] = phicam;

    fprintf(stdout,
	    "cam_th_cal=%g [deg] th_beg=%g th_len=%g a=%g R0=%g hslope=%g\n",
	    (th_beg + th_len * Xcam[2]) * 180. / M_PI, th_beg, th_len, a,
	    R0, hslope);


    /* fix camera field of view */
    /* units = GM/c^2 in plane of the hole */
    DX = 40.0;
    DY = 40.0;
    fovx = DX / rcam;
    fovy = DY / rcam;

    //Dsource = DM87 ;
    Dsource = DSGRA;
    scale = (DX * L_unit / NX) * (DY * L_unit / NY) / (Dsource * Dsource) / JY;
    fprintf(stderr,"intensity [cgs] to flux per pixel [Jy] conversion: %g\n",scale);
    fprintf(stderr,"Dsource: %g [cm]\n",Dsource);
    fprintf(stderr,"Dsource: %g [kpc]\n",Dsource/(1.e3*PC));
    fprintf(stderr,"FOVx, FOVy: %g %g [GM/c^2]\n",DX,DY);
    fprintf(stderr,"FOVx, FOVy: %g %g [rad]\n",DX*L_unit/Dsource,DY*L_unit/Dsource);
    fprintf(stderr,"FOVx, FOVy: %g %g [muas]\n",DX*L_unit/Dsource * 2.06265e11 ,DY*L_unit/Dsource * 2.06265e11);

    int nprogress = 0;

#pragma omp parallel for \
default(none) \
schedule(static,NX*NY/nthreads) \
collapse(2) \
private(i,j,k,l,ki,kf,ji,jf,nstep,dl,X,Xhalf,Kcon,Kconhalf,\
   Xi,Xf,Kconi,Kconf,traj,Intensity,N_coord,Stokes_I,Stokes_Q,Stokes_U,\
   Stokes_V) \
shared(Xcam,fovx,fovy,freq,freqcgs,image,imageS,L_unit,stderr,stdout,\
   th_beg,th_len,hslope,nthreads,nprogress) 
    for (i = 0; i < NX; i++) {
	for (j = 0; j < NY; j++) {

	    init_XK(i, j, Xcam, fovx, fovy, X, Kcon);

	    for (k = 0; k < NDIM; k++)
		Kcon[k] *= freq;

	    /* integrate geodesic backwards along trajectory */
	    nstep = 0;
	    while (!stop_backward_integration(X, Kcon, Xcam)) {

		/* This stepsize function can be troublesome inside of R = 2M,
		   and should be used cautiously in this region. */
		dl = stepsize(X, Kcon);

		/* move photon one step backwards, the procecure updates X 
		   and Kcon full step and returns also values in the middle */
		push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
		nstep++;

		traj[nstep].dl = dl * L_unit * HPL / (ME * CL * CL);
		for (k = 0; k < NDIM; k++)
		    traj[nstep].X[k] = X[k];
		for (k = 0; k < NDIM; k++)
		    traj[nstep].Kcon[k] = Kcon[k];
		for (k = 0; k < NDIM; k++)
		    traj[nstep].Xhalf[k] = Xhalf[k];
		for (k = 0; k < NDIM; k++)
		    traj[nstep].Kconhalf[k] = Kconhalf[k];

		if (nstep > MAXNSTEP - 2) {
		    fprintf(stderr, "MAXNSTEP exceeded on j=%d i=%d\n", j,
			    i);
		    exit(1);
		}

	    }
	    nstep--; /* final step violated the "stop" condition,so don't record it */
	    /* DONE geodesic integration */

	    /* integrate forwards along trajectory, including radiative transfer equation */
	    // initialize N,Intensity; need X, K for this.
	    for (l = 0; l < NDIM; l++) {
		Xi[l] = traj[nstep].X[l];
		Kconi[l] = traj[nstep].Kcon[l];
	    }
	    init_N(Xi, Kconi, N_coord);
	    Intensity = 0.0;

	    while (nstep > 1) {	

		/* initialize X,K */
		for (l = 0; l < NDIM; l++) {
		    Xi[l]       = traj[nstep].X[l];
		    Kconi[l]    = traj[nstep].Kcon[l];
		    Xhalf[l]    = traj[nstep].Xhalf[l];
		    Kconhalf[l] = traj[nstep].Kconhalf[l];
		    Xf[l]       = traj[nstep - 1].X[l];
		    Kconf[l]    = traj[nstep - 1].Kcon[l];
		}

		/* solve total intensity equation alone */
		get_jkinv(Xi, Kconi, &ji, &ki);
		get_jkinv(Xf, Kconf, &jf, &kf);
		Intensity = approximate_solve(Intensity, ji, ki, jf, kf, traj[nstep].dl);

		/* solve polarized transport */
		evolve_N(Xi, Kconi, Xhalf, Kconhalf, Xf, Kconf, traj[nstep].dl, N_coord);

                /* swap start and finish */
		ji = jf;
		ki = kf;

		nstep--;
	    }

	    /* deposit intensity, and Stokes parameters in pixels */
	    image[i][j] = Intensity * pow(freqcgs, 3);
	    project_N(Xf, Kconf, N_coord, &Stokes_I, &Stokes_Q, &Stokes_U, &Stokes_V);
	    imageS[i][j][0] = Stokes_I * pow(freqcgs, 3);
	    imageS[i][j][1] = Stokes_Q * pow(freqcgs, 3);
	    imageS[i][j][2] = Stokes_U * pow(freqcgs, 3);
	    imageS[i][j][3] = Stokes_V * pow(freqcgs, 3);

	    if( (nprogress % NY) == 0 ) {
	        fprintf(stderr,"%d ",(nprogress/NY));
	    }

#pragma omp atomic	
            nprogress += 1;
	}
    }


    /* printing out to files and on stderr */
    Ftot = 0.;
    Imax = 0.0;
    Iavg = 0.0;
    imax = jmax = 0;
    for (i = 0; i < NX; i++)
	for (j = 0; j < NY; j++) {
	    Ftot += image[i][j] * scale;	
	    Iavg += image[i][j];
	    if (image[i][j] > Imax) {
		imax = i;
		jmax = j;
		Imax = image[i][j];
	    }
	}

    fprintf(stderr, "imax=%d jmax=%d Imax=%g Iavg=%g\n", imax, jmax, Imax,
	    Iavg / (NX * NY));
    fprintf(stderr, "Ftot: %g %g scale=%g\n", freqcgs, Ftot, scale);
    fprintf(stderr, "nuLnu = %g\n",
	    Ftot * Dsource * Dsource * JY * freqcgs);

    /* image, dump result */
    make_ppm(image, freq, "ipole_fnu.ppm");
    dump(image, imageS, "ipole.dat", scale);
    for (i = 0; i < NX; i++)
	for (j = 0; j < NY; j++)
	    image[i][j] = log(image[i][j] + 1.e-50);
    make_ppm(image, freq, "ipole_lfnu.ppm");

    /* done! */


}

void dump(double image[NX][NY], double imageS[NX][NY][NDIM], char *fname,
	  double scale)
{
    FILE *fp;
    int i, j;
    double xx, yy, dum_r;
    double sum_i;


    fp = fopen(fname, "w");
    if (fp == NULL) {
	fprintf(stderr, "unable to open %s\n", fname);
	exit(1);
    }


    sum_i = 0.0;
    for (i = 0; i < NX; i++) {
	for (j = 0; j < NY; j++) {
	    sum_i += image[i][j];
	    fprintf(fp, "%d %d %15.10g %15.10g %15.10g %15.10g %15.10g \n",
		    i, j, image[i][j], imageS[i][j][0], imageS[i][j][1],
		    imageS[i][j][2], imageS[i][j][3]);
	}
	fprintf(fp, "\n");
    }
    fclose(fp);


    /*dump vtk file */
    fp = fopen("ipole.vtk", "w");
    if (fp == NULL) {
	fprintf(stderr, "unable to open %s\n", "ipole.vtk");
	exit(1);
    }

    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Image Simulation Result\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", NX + 1, NY + 1, 1);
    fprintf(fp, "POINTS %d float\n", (NX + 1) * (NY + 1));

    for (j = 0; j <= NY; j++) {
	for (i = 0; i <= NX; i++) {
	    xx = i * 1.0;	//-x_size/2.+i*stepx;                                                                                                                 
	    yy = j * 1.0;	//-y_size/2.+j*stepy;                                                                                                                 
	    fprintf(fp, "%g %g %g\n", xx, yy, 0.0);
	}
    }
    fprintf(fp, "\nPOINT_DATA %d\n", (NX + 1) * (NY + 1));
    fprintf(fp, "SCALARS Intensity float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    for (j = 0; j <= NY; j++) {
	for (i = 0; i <= NX; i++) {
	    dum_r = (image[i][j]);
	    fprintf(fp, "%g\n", dum_r);
	}
    }

    fclose(fp);



}

void init_XK(int i, int j, double Xcam[4], double fovx, double fovy,	/* field of view, in radians */
	  double X[4], double Kcon[4]	/* position, wavevector */
    )
{
    double Econ[NDIM][NDIM];
    double Ecov[NDIM][NDIM];
    double Kcon_tetrad[NDIM];
    int k;

    /* construct orthonormal tetrad.
       e^0 along Ucam
       e^3 outward (!) along radius vector
       e^2 toward north pole of coordinate system
       ("y" for the image plane)
       e^1 in the remaining direction
       ("x" for the image plane)

    */

    make_camera_tetrad(Xcam, Econ, Ecov);

    /* construct *outgoing* wavevectors */
    Kcon_tetrad[0] = 0.;
    Kcon_tetrad[1] = (i / ((double) NX) - 0.5) * fovx;
    Kcon_tetrad[2] = (j / ((double) NY) - 0.5) * fovy;
    Kcon_tetrad[3] = 1.;

    /* normalize */
    null_normalize(Kcon_tetrad, 1.);

    /* translate into coordinate frame */
    tetrad_to_coordinate(Econ, Kcon_tetrad, Kcon);

    /* set position */
    for (k = 0; k < NDIM; k++)
	X[k] = Xcam[k];

    /* done! */
}

/* normalize null vector in a tetrad frame */
void null_normalize(double Kcon[NDIM], double fnorm)
{
    double inorm;

    inorm =
	sqrt(Kcon[1] * Kcon[1] + Kcon[2] * Kcon[2] + Kcon[3] * Kcon[3]);

    Kcon[0] = fnorm;
    Kcon[1] *= fnorm / inorm;
    Kcon[2] *= fnorm / inorm;
    Kcon[3] *= fnorm / inorm;

}

/* 

   must be a stable, approximate solution to radiative transfer
   that runs between points w/ initial intensity I, emissivity
   ji, opacity ki, and ends with emissivity jf, opacity kf.

   Return final intensity

*/

double approximate_solve(double Ii, double ji,
			 double ki, double jf, double kf, double dl)
{
    double efac, If, javg, kavg, dtau;

    javg = (ji + jf) / 2.;
    kavg = (ki + kf) / 2.;

    dtau = dl * kavg;

    if (dtau < 1.e-3) {
	If = Ii + (javg - Ii * kavg) * dl * (1. -
					     (dtau / 2.) * (1. -
							    dtau / 3.));
    } else {
	efac = exp(-dtau);
	If = Ii * efac + (javg / kavg) * (1. - efac);
    }

    return (If);
}

