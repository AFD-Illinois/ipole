#include "decs.h"
#include "defs.h"
#include <omp.h>

#define MAXNSTEP 	50000

//extern double image[NX][NY];
static double image[NX][NY];
//extern double imageS[NX][NY][NIMG];
static double imageS[NX][NY][NIMG];
static double complex Ngeo[NX][NY][NDIM][NDIM];
static double tauFgeo[NX][NY];

static double Xgeo[NX][NY][NDIM], Kcongeo[NX][NY][NDIM], t[NX][NY], tmin[NX][NY];

static double DX, DY, fovx, fovy;

int main(int argc, char *argv[])
{
    //omp_set_num_threads(1);
    //double X[NDIM], Kcon[NDIM];
    //double Xhalf[NDIM], Kconhalf[NDIM];
    //double dl, Intensity;
    //double DX, DY, fovx, fovy;
    double thetacam, phicam, rcam, Xcam[NDIM];

    set_levi_civita();

    double freq, freqcgs;
    double Ftot, Dsource;
    //int i, j, k, l, nstep;
    //double Xi[NDIM], Xf[NDIM], Kconi[NDIM], Kconf[NDIM], ji, ki, jf, kf;
    double scale;
    //double root_find(double th);
    double root_find(double x[NDIM]);
    int imax, jmax;
    double Imax, Iavg;		//,dOmega;
    //double Stokes_I, Stokes_Q, Stokes_U, Stokes_V;
    //double tauF;
    //double complex N_coord[NDIM][NDIM];

    #pragma omp parallel
    {
      if (omp_get_thread_num() == 0) {
        int nthreads = omp_get_num_threads();
        printf("nthreads = %d\n", nthreads);
      }
    }

    if (argc < 7) {
	// fprintf(stderr,"usage: ipole theta freq filename Munit theta_j trat_d\n") ;
	fprintf(stderr,
		"usage: ipole theta freq filename Munit trat_j trat_d counterjet\n");
	exit(0);
    }

    sscanf(argv[1], "%lf", &thetacam);
    sscanf(argv[2], "%lf", &freqcgs);
    sscanf(argv[4], "%lf", &M_unit);
    sscanf(argv[5], "%lf", &trat_j);
    sscanf(argv[6], "%lf", &trat_d);
    sscanf(argv[7], "%d", &counterjet);

    init_model(argv);
    R0 = 0.0;

    /* normalize frequency to electron rest-mass energy */
    freq = freqcgs * HPL / (ME * CL * CL);

    /* fix camera location */
    rcam = 240.;
    phicam = 0.0;
    Xcam[0] = 0.0;
    Xcam[1] = log(rcam);
    double x[NDIM] = {0., rcam, thetacam/180.*M_PI, phicam};
    Xcam[2] = root_find(x);
    Xcam[3] = phicam;

    printf("X[] = %e %e %e %e\n", Xcam[0], Xcam[1], Xcam[2], Xcam[3]);

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

    // Maximum radius of radiation interactions (GM/c^2)
    rmax = 50.;

    //Dsource = DM87 ;
    Dsource = DSGRA;
    scale = (DX * L_unit / NX) * (DY * L_unit / NY) / (Dsource * Dsource) / JY;
    printf("L_unit = %e DX = %e NX = %i Dsource = %e JY = %e\n", L_unit, DX, NX, Dsource,JY);
    fprintf(stderr,"intensity [cgs] to flux per pixel [Jy] conversion: %g\n",scale);
    fprintf(stderr,"Dsource: %g [cm]\n",Dsource);
    fprintf(stderr,"Dsource: %g [kpc]\n",Dsource/(1.e3*PC));
    fprintf(stderr,"FOVx, FOVy: %g %g [GM/c^2]\n",DX,DY);
    fprintf(stderr,"FOVx, FOVy: %g %g [rad]\n",DX*L_unit/Dsource,DY*L_unit/Dsource);
    fprintf(stderr,"FOVx, FOVy: %g %g [muas]\n",DX*L_unit/Dsource * 2.06265e11 ,DY*L_unit/Dsource * 2.06265e11);

  //  int nprogress = 0;

  printf("\nBACKWARD\n");
  #pragma omp parallel
  {
    int ngeomax = 0;
    double Xhalf[NDIM], Kconhalf[NDIM];
    #pragma omp for collapse(2)
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NY; j++) {
        if (j == 0) {printf("%i\n", i);}
        init_XK(i, j, Xcam, fovx, fovy, Xgeo[i][j], Kcongeo[i][j]);

        for (int mu = 0; mu < NDIM; mu++) Kcongeo[i][j][mu] *= freq;

        // Integrate geodesic backwards in time
        int ngeo = 0;
        int tmin_set = 0;
        double dl;
        while (!stop_backward_integration(Xgeo[i][j], Kcongeo[i][j], Xcam)) {
		      dl = stepsize(Xgeo[i][j], Kcongeo[i][j]);
          t[i][j] = Xgeo[i][j][0];
          if (Xgeo[i][j][1] < log(rmax) && tmin_set == 0) {
            tmin_set = 1;
            tmin[i][j] = t[i][j];
          }
          
          push_photon(Xgeo[i][j], Kcongeo[i][j], -dl, Xhalf, Kconhalf);
		      ngeo++;
        }

		    if (ngeo > MAXNSTEP - 2) {
		      fprintf(stderr, "MAXNSTEP exceeded on j=%d i=%d\n", j, i);
		      exit(1);
		    }

        if (ngeo > ngeomax) ngeomax = ngeo;
      }
    }
    printf("[%i] ngeomax = %i\n", omp_get_thread_num(), ngeomax);
  } // pragma omp parallel

  double tmax = 0.;
  //int imax, jmax;
  imax = jmax = 0;
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      if (t[i][j] < tmax) {
        tmax = t[i][j];
        imax = i;
        jmax = j;
      }
      //if (tmin
    }
  }
  printf("tmax = %e [%i %i]\n", tmax, imax, jmax);

  int done[NX][NY];//, done_thisslice[NX][NY];
  #pragma omp parallel for collapse(2)
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      //t[i][j] = tmax - t[i][j];
      done[i][j] = 0;
      image[i][j] = 0.;
	    init_N(Xgeo[i][j], Kcongeo[i][j], Ngeo[i][j]);
      tauFgeo[i][j] = 0.;
    }
  }

  DTd = 5.;
  //double DTd = 5.; // READ FROM FILES! ASSUME CONSTANT DTd!
  double tcurr = tmax;

  int nloop = 0;

  // Loop while all rays not yet finished
  while (tcurr < 0.) {

  ///// NOTE! dl MUST correspond to < DTd! Probably < DTd/2 for safety

    // Loop over rays, push to tcurr if necessary.
/*#pragma omp parallel for \
default(none) \
schedule(static,NX*NY/nthreads) \
collapse(2) \
private(i,j,k,l,ki,kf,ji,jf,nstep,dl,X,Xhalf,Kcon,Kconhalf,\
   Xi,Xf,Kconi,Kconf,traj,Intensity,N_coord,Stokes_I,Stokes_Q,Stokes_U,\
   Stokes_V,tauF) \
shared(Xcam,fovx,fovy,freq,freqcgs,image,imageS,L_unit,stderr,stdout,\
   th_beg,th_len,hslope,nthreads,nprogress,Kcongeo,Xgeo,t,tmin,tauFgeo,Ngeo,\
   rmax,tcurr,done)*/
   
    //#pragma omp parallel
    {

    //#pragma for collapse(2)
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NY; j++) {
        double Xhalfgeo[NDIM], Kconhalfgeo[NDIM];
        double Xprevgeo[NDIM], Kconprevgeo[NDIM];
        while (Xgeo[i][j][0] < tcurr && Xgeo[i][j][0] < 0.) {
          //printf("i,j = %i %i\n", i,j);
          // push
          //t[i][j] += DTd;

          // Integrate geodesic, retaining n, n+1/2, and n+1 X and K
          for (int mu = 0; mu < NDIM; mu++) {
            Xprevgeo[mu] = Xgeo[i][j][mu];
            Kconprevgeo[mu] = Kcongeo[i][j][mu];
          }
          double dl;
          dl = stepsize(Xgeo[i][j], Kcongeo[i][j]);
          push_photon(Xgeo[i][j], Kcongeo[i][j], dl, Xhalfgeo, Kconhalfgeo);

          double ji, jf, ki, kf;
          get_jkinv(Xprevgeo, Kconprevgeo, &ji, &ki);
		      get_jkinv(Xgeo[i][j], Kcongeo[i][j], &jf, &kf);
          dl *= L_unit*HPL/(ME*CL*CL);
		      image[i][j] = approximate_solve(image[i][j], ji, ki, jf, kf, dl);

          evolve_N(Xprevgeo, Kconprevgeo, Xhalfgeo, Kconhalfgeo, Xgeo[i][j],
            Kcongeo[i][j], dl, Ngeo[i][j], &tauFgeo[i][j]);
          /*if (i == 63 && j == 78) {
            printf("X[0] = %e Ia = %e N[0][0] = %e tauF = %e\n", Xgeo[i][j][0],
              image[i][j], Ngeo[i][j][0][0], tauFgeo[i][j]);
          }*/

          /*if (i == 65 && j == 71) {
            printf("X[] = %e %e %e %e\n", Xgeo[i][j][0], Xgeo[i][j][1], Xgeo[i][j][2], Xgeo[i][j][3]);
            printf("K[] = %e %e %e %e\n\n", Kcongeo[i][j][0], Kcongeo[i][j][1], Kcongeo[i][j][2], Kcongeo[i][j][3]);
          }*/

          if (tauFgeo[i][j] > 1.e100 || tauFgeo[i][j] < -1.e100) {
            printf("i,j = %i %i t = %e, X0 = %e\n", i,j, t[i][j], Xgeo[i][j][0]);
            for (int mu = 0; mu < NDIM; mu++) {
              printf("[%i] X = %e Kcon = %e\n", mu, Xgeo[i][j][mu], Kcongeo[i][j][mu]);
              printf("     Xprev = %e Kconprev = %e\n", Xprevgeo[mu], Kconprevgeo[mu]);
              printf("     Xhalf = %e Kconhalf = %e\n", Xhalfgeo[mu], Kconhalfgeo[mu]);
            }
            exit(-1);
          }


        }
        //printf("Xgeo[i][j][1] = %e log(rmax) = %e\n", Xgeo[i][j][1], log(rmax));
        if (Kcongeo[i][j][1] > 0. && Xgeo[i][j][1] > log(rmax)) {
          done[i][j] = 1;
        }
      }
    }

    printf("tcurr = %e\n", tcurr);

    int alldone = 1;
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NY; j++) {
        if (done[i][j] == 0)
          alldone = 0;
      }
    }

    if (!alldone) {
      printf("LOADING NEW DATA\n");
    }
    // load newest timeslice, rename two most recent. Don't do this if all superphotons outside of radiative region.

    tcurr += DTd;
    nloop++;
  }
  } // pragma omp for

  // Conversion factors
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      image[i][j] *= pow(freqcgs, 3);
      project_N(Xgeo[i][j], Kcongeo[i][j], Ngeo[i][j], &imageS[i][j][0],
        &imageS[i][j][1], &imageS[i][j][2], &imageS[i][j][3]);
      imageS[i][j][0] *= pow(freqcgs, 3);
      imageS[i][j][1] *= pow(freqcgs, 3);
      imageS[i][j][2] *= pow(freqcgs, 3);
      imageS[i][j][3] *= pow(freqcgs, 3);
      imageS[i][j][4] = tauFgeo[i][j];
    }
  }

  /* printing out to files and on stderr */
  Ftot = 0.;
  Imax = 0.0;
  Iavg = 0.0;
  imax = jmax = 0;
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      Ftot += image[i][j] * scale;
      Iavg += image[i][j];
      if (image[i][j] > Imax) {
        imax = i;
        jmax = j;
        Imax = image[i][j];
      }
    }
  } 

    fprintf(stderr, "imax=%d jmax=%d Imax=%g Iavg=%g\n", imax, jmax, Imax,
	    Iavg / (NX * NY));
    fprintf(stderr, "Ftot: %g %g scale=%g\n", freqcgs, Ftot, scale);
    fprintf(stderr, "nuLnu = %g\n",
	    4.*M_PI*Ftot * Dsource * Dsource * JY * freqcgs);

    /* image, dump result */
    //make_ppm(image, freq, "ipole_fnu.ppm");
    dump(image, imageS, "ipole.dat", scale);
    //for (int i = 0; i < NX; i++)
	//for (int j = 0; j < NY; j++)
	  IMLOOP image[i][j] = log(image[i][j] + 1.e-50);
    make_ppm(image, freq, "ipole_lfnu.ppm");
}

void dump(double image[NX][NY], double imageS[NX][NY][NIMG], char *fname,
	  double scale)
{
    FILE *fp;
    int i, j;
    double xx, yy;
    double sum_i;


    fp = fopen(fname, "w");
    if (fp == NULL) {
	fprintf(stderr, "unable to open %s\n", fname);
	exit(1);
    }
    fprintf(fp, "%d %d %e %e %e %e %e\n", NX, NY, DX, DY, scale, L_unit, M_unit);


    sum_i = 0.0;
    for (i = 0; i < NX; i++) {
	for (j = 0; j < NY; j++) {
    if (i == 50 && j == 50) {
      printf("I = %e\n", image[i][j]);
      //exit(-1);
    }
	    sum_i += image[i][j];
      xx = (i+0.5)*DX/NX;
      yy = (j+0.5)*DY/NY;
      fprintf(fp, "%d %d %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g \n",
		    i, j, xx, yy, image[i][j], imageS[i][j][0], imageS[i][j][1],
		    imageS[i][j][2], imageS[i][j][3], imageS[i][j][4]);
	}
	fprintf(fp, "\n");
    }
    fclose(fp);


    /*dump vtk file */
    /*fp = fopen("ipole.vtk", "w");
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
	    double dum_r = (image[i][j]);
	    fprintf(fp, "%g\n", dum_r);
	}
    }

    fclose(fp);*/



}

void init_XK(int i, int j, double Xcam[4], double fovx, double fovy,	/* field of view, in radians */
	  double X[4], double Kcon[4]	/* position, wavevector */
    )
{
    double Econ[NDIM][NDIM];
    double Ecov[NDIM][NDIM];
    double Kcon_tetrad[NDIM];
    //int k;

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
    for (int mu = 0; mu < NDIM; mu++)
	    X[mu] = Xcam[mu];
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

