#include "decs.h"
#include "defs.h"
#include <omp.h>

#define MAXNSTEP 50000

struct of_traj {
  double dl;
  double X[NDIM];
  double Kcon[NDIM];
  double Xhalf[NDIM];
  double Kconhalf[NDIM];
} traj[MAXNSTEP];
#pragma omp threadprivate(traj)

// global variables
double thetacam, freqcgs;
char fnam[STRLEN];

Params params = { 0 };

int main(int argc, char *argv[]) 
{
  // initialization
  double time = omp_get_wtime();
  set_levi_civita();

  double phicam, rcam, Xcam[NDIM];
  double freq, Dsource, scale;
  double DX, DY, fovx, fovy;

  double image[NX][NY];
  double imageS[NX][NY][NIMG];

  // "forward declarations"
  double root_find(double X[NDIM]);

#pragma omp parallel
  if (omp_get_thread_num() == 0) {
    fprintf(stderr, "nthreads = %d\n", omp_get_num_threads());
  }

  // load values from parameter file if passed. handle all
  // actual model parameter comprehension in model/* files
  for (int i=0; i<argc-1; ++i) {
    if ( strcmp(argv[i], "-par") == 0 ) {
      load_par(argv[i+1], &params);
    }
  }

  parse_input(argc, argv, &params);

  init_model(argv);
  R0 = 0.0;

  /* normalize frequency to electron rest-mass energy */
  freq = freqcgs * HPL / (ME * CL * CL);

  /* fix camera location */
  rcam = 240.;
  phicam = 0.0;

  if (params.loaded) {
    phicam = params.phicam;
  }

  Xcam[0] = 0.0;
  Xcam[1] = log(rcam);
  double x[NDIM] = {0., rcam, thetacam/180.*M_PI, phicam/180.*M_PI};
  Xcam[2] = root_find(x);
  Xcam[3] = phicam/180.*M_PI;

  printf("Xcam[] = %e %e %e %e\n", Xcam[0], Xcam[1], Xcam[2], Xcam[3]);

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
  //rmax = 50.;

  Dsource = DM87 ;
  //Dsource = DSGRA;
  scale = (DX * L_unit / NX) * (DY * L_unit / NY) / (Dsource * Dsource) / JY;
  printf("L_unit = %e DX = %e NX = %i Dsource = %e JY = %e\n", L_unit, DX, NX, Dsource,JY);
  fprintf(stderr,"intensity [cgs] to flux per pixel [Jy] conversion: %g\n",scale);
  fprintf(stderr,"Dsource: %g [cm]\n",Dsource);
  fprintf(stderr,"Dsource: %g [kpc]\n",Dsource/(1.e3*PC));
  fprintf(stderr,"FOVx, FOVy: %g %g [GM/c^2]\n",DX,DY);
  fprintf(stderr,"FOVx, FOVy: %g %g [rad]\n",DX*L_unit/Dsource,DY*L_unit/Dsource);
  fprintf(stderr,"FOVx, FOVy: %g %g [muas]\n",DX*L_unit/Dsource * 2.06265e11 ,DY*L_unit/Dsource * 2.06265e11);

  // slow light is implemented here. skipping
  if (0) {

  } 
  
  // fast light
  else {

    int nprogress = 0;

#pragma omp parallel for schedule(dynamic,1) collapse(2) shared(nprogress,image,imageS)
    for (int i=0; i<NX; ++i) {
      for (int j=0; j<NY; ++j) {
       
        double X[NDIM],Xhalf[NDIM],Xi[NDIM],Xf[NDIM],Kcon[NDIM],Kconhalf[NDIM],Kconi[NDIM],Kconf[NDIM];
        double dl, ji,ki, jf,kf;
        double complex N_coord[NDIM][NDIM];
        double Intensity, Stokes_I, Stokes_Q, Stokes_U, Stokes_V, tauF;

        init_XK(i,j, Xcam, fovx,fovy, X, Kcon); 
        for (int k=0; k<NDIM; ++ k) Kcon[k] *= freq;

        /* integrate geodesic backwards */
        int nstep = 0;

        while (!stop_backward_integration(X, Kcon, Xcam)) {
          /* This stepsize function can be troublesome inside of R = 2M,
             and should be used cautiously in this region. */
          dl = stepsize(X, Kcon);
          
          /* move photon one step backwards, the procecure updates X
             and Kcon full step and returns also values in the middle */
          push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
          nstep++;

          traj[nstep].dl = dl * L_unit * HPL / (ME * CL * CL);

          for (int k = 0; k < NDIM; k++)
            traj[nstep].X[k] = X[k];
          for (int k = 0; k < NDIM; k++)
            traj[nstep].Kcon[k] = Kcon[k];
          for (int k = 0; k < NDIM; k++)
            traj[nstep].Xhalf[k] = Xhalf[k];
          for (int k = 0; k < NDIM; k++)
            traj[nstep].Kconhalf[k] = Kconhalf[k];

          if (nstep > MAXNSTEP - 2) {
            fprintf(stderr, "MAXNSTEP exceeded on j=%d i=%d\n", j,i);
            exit(1);
          }
        }

        nstep--; // don't record final step because it violated "stop" condition
        /* DONE geodesic integration */
        
        // integrate forwards along trajectory, including radiative transfer equation
        // initialize N, Intensity -- need X, K for this.
        for (int l=0; l<NDIM; ++l) {
            Xi[l] = traj[nstep].X[l];
            Kconi[l] = traj[nstep].Kcon[l];
        }

        init_N(Xi, Kconi, N_coord);
        Intensity = 0.0;
        tauF = 0.;

        get_jkinv(traj[nstep].X, traj[nstep].Kcon, &ji, &ki);

        int flag = 0;

        while (nstep > 1) {

          // initialize X,K
          for (int l=0; l<NDIM; ++l) {
            Xi[l]       = traj[nstep].X[l];
            Kconi[l]    = traj[nstep].Kcon[l];
            Xhalf[l]    = traj[nstep].Xhalf[l]; 
            Kconhalf[l] = traj[nstep].Kconhalf[l];
            Xf[l]       = traj[nstep - 1].X[l];
            Kconf[l]    = traj[nstep - 1].Kcon[l];
          }

          // solve total intensity equation alone
          get_jkinv(Xf, Kconi, &jf, &kf);
          Intensity = approximate_solve(Intensity, ji, ki, jf, kf, traj[nstep].dl);

          // solve polarized transport
          evolve_N(Xi, Kconi, Xhalf, Kconhalf, Xf, Kconf, traj[nstep].dl, N_coord, &tauF);
          if (isnan(creal(N_coord[0][0]))) exit(-1);

          if (flag == 0) {
            project_N(Xf, Kconf, N_coord, &Stokes_I, &Stokes_Q, &Stokes_U, &Stokes_V);
            if (Stokes_I < 0) {
              void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]);
              double tdel[NDIM];
              int ti,tj,tk;
              Xtoijk(Xf, &ti, &tj, &tk, tdel);
              fprintf(stderr, "Stokes_I %g @ %d %d %d (%g %g %g %g)\n", Stokes_I, ti, tj, tk, Xf[0], Xf[1], Xf[2], Xf[3]);
              flag = 1;
            }
          }

          // swap start and finish
          ji = jf;
          ki = kf;

          nstep--;
        }

        // deposit intensity and Stokes parameter in pixels
        image[i][j] = Intensity * pow(freqcgs, 3);
        project_N(Xf, Kconf, N_coord, &Stokes_I, &Stokes_Q, &Stokes_U, &Stokes_V);
        imageS[i][j][0] = Stokes_I * pow(freqcgs, 3);
        imageS[i][j][1] = Stokes_Q * pow(freqcgs, 3);
        imageS[i][j][2] = Stokes_U * pow(freqcgs, 3);
        imageS[i][j][3] = Stokes_V * pow(freqcgs, 3);
        imageS[i][j][4] = tauF;
        if (isnan(imageS[i][j][0])) exit(-1);

        // print progress
        if ( (nprogress % NY) == 0 ) {
          fprintf(stderr, "%d ", nprogress/NY);
        }

#pragma omp atomic
        nprogress += 1;

      }
    }
  }

  // finish up
  printf("\n");
  printf("scale = %e\n", scale);

  double Ftot = 0.;
  double Ftot_unpol = 0.;
  double Imax = 0.0;
  double Iavg = 0.0;
  int imax = 0;
  int jmax = 0;
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      Ftot_unpol += image[i][j]*scale;
      Ftot += imageS[i][j][0] * scale;
      Iavg += imageS[i][j][0];
      if (imageS[i][j][0] > Imax) {
        imax = i;
        jmax = j;
        Imax = imageS[i][j][0];
      }
    }
  }
  
  fprintf(stderr, "imax=%d jmax=%d Imax=%g Iavg=%g\n", imax, jmax, Imax, Iavg/(NX*NY));
  fprintf(stderr, "freq: %g Ftot: %g (%g unpol) scale=%g\n", freqcgs, Ftot, Ftot_unpol, scale);
  fprintf(stderr, "nuLnu = %g\n", 4.*M_PI*Ftot * Dsource * Dsource * JY * freqcgs);

  // dump result. if parameters have been loaded, don't also
  // output image
  if (params.loaded) {
    dump(image, imageS, params.outf, scale, Dsource, Xcam, DX, DY, fovx, fovy);
  } else {
    dump(image, imageS, "ipole.dat", scale, Dsource, Xcam, DX, DY, fovx, fovy);
    IMLOOP image[i][j] = log(image[i][j] + 1.e-50);
    make_ppm(image, freq, "ipole_lfnu.ppm");
  }

  time = omp_get_wtime() - time;
  printf("Total wallclock time: %g s\n", time);

  return 0;
}


void dump(double image[NX][NY], double imageS[NX][NY][NIMG], const char *fname,
          double scale, double Dsource, double cam[NDIM], double DX, double DY,
          double fovx, double fovy)
{
  hid_t fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  if (fid < 0) {
    fprintf(stderr, "! unable to open/create hdf5 output file.\n");
    exit(-3);
  }

  h5io_add_attribute_str(fid, "/", "githash", xstr(VERSION));

  h5io_add_group(fid, "/header");
  h5io_add_data_dbl(fid, "/header/freqcgs", freqcgs);
  h5io_add_data_dbl(fid, "/header/scale", scale);
  h5io_add_data_dbl(fid, "/header/dsource", Dsource);

  h5io_add_group(fid, "/header/camera");
  h5io_add_data_int(fid, "/header/camera/nx", NX);
  h5io_add_data_int(fid, "/header/camera/ny", NY);
  h5io_add_data_dbl(fid, "/header/camera/dx", DX);
  h5io_add_data_dbl(fid, "/header/camera/dy", DY);
  h5io_add_data_dbl(fid, "/header/camera/fovx", fovx);
  h5io_add_data_dbl(fid, "/header/camera/fovy", fovy);
  h5io_add_data_dbl_1d(fid, "/header/camera/x", NDIM, cam);

  h5io_add_group(fid, "/header/units");
  h5io_add_data_dbl(fid, "/header/units/L_unit", L_unit);
  h5io_add_data_dbl(fid, "/header/units/M_unit", M_unit);
  h5io_add_data_dbl(fid, "/header/units/T_unit", T_unit);
  h5io_add_data_dbl(fid, "/header/units/Thetae_unit", Te_unit);

  // processing
  double Ftot_unpol=0., Ftot=0.;
  for (int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      Ftot_unpol += image[i][j] * scale;
      Ftot += imageS[i][j][0] * scale;
    }
  }
      
  // output stuff
  h5io_add_data_dbl(fid, "/Ftot", Ftot);
  h5io_add_data_dbl(fid, "/Ftot_unpol", Ftot_unpol);
  h5io_add_data_dbl(fid, "/nuLnu", 4. * M_PI * Ftot * Dsource * Dsource * JY * freqcgs);
  h5io_add_data_dbl(fid, "/nuLnu_unpol", 4. * M_PI * Ftot_unpol * Dsource * Dsource * JY * freqcgs);

  h5io_add_data_dbl_2d(fid, "/unpol", NX, NY, image);
  h5io_add_data_dbl_3d(fid, "/pol", NX, NY, 5, imageS);

  // allow model to output
  output_hdf5(fid);

  // housekeeping
  H5Fclose(fid);
}

/* construct orthonormal tetrad.
 *   e^0 along Ucam
 *   e^3 outward (!) along radius vector
 *   e^2 toward north pole of coordinate system
 *   ("y" for the image plane)
 *   e^1 in the remaining direction
 *   ("x" for the image plane)
 */
void init_XK(int i, int j, double Xcam[NDIM], double fovx, double fovy, 
             double X[NDIM], double Kcon[NDIM]) 
{
  double Econ[NDIM][NDIM];
  double Ecov[NDIM][NDIM];
  double Kcon_tetrad[NDIM];
  
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
  for (int mu = 0; mu < NDIM; mu++) X[mu] = Xcam[mu];

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

/* must be a stable, approximate solution to radiative transfer
 * that runs between points w/ initial intensity I, emissivity
 * ji, opacity ki, and ends with emissivity jf, opacity kf.
 *
 * return final intensity
 */
double approximate_solve(double Ii, double ji, double ki, double jf,
                         double kf, double dl)
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

  return If;
}

