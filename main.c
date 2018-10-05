#include "decs.h"
#include "defs.h"
#include <omp.h>

#define MAXNSTEP 50000

/* slow light globals here */
static double Xgeo[NX][NY][NDIM], Kcongeo[NX][NY][NDIM];
static double t[NX][NY], tmin[NX][NY], tauFgeo[NX][NY];
static double complex Ngeo[NX][NY][NDIM][NDIM];
static double X[NX][NY], Y[NX][NY];
#pragma omp threadprivate(Xgeo,Kcongeo,t,tmin,tauFgeo,Ngeo)


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
  counterjet = 0;

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
  rcam = 1000.;
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

  // slow light (untested 2018.09.25 gnw)
  if (0) {

    fprintf(stderr, "\nBACKWARDS\n");

    // get photon trajectories
#pragma omp parallel
    {
      int ngeomax = 0;
#pragma omp for collapse(2)
      for (int i=0; i<NX; ++i) {
        for (int j=0; j<NY; ++j) {

          double Xhalf[NDIM], Kconhalf[NDIM];

          if (!j) fprintf(stderr, "%i\n", i);

          init_XK(i, j, Xcam, fovx, fovy, Xgeo[i][j], Kcongeo[i][j]);
          for (int k=0; k<NDIM; ++k) Kcongeo[i][j][k] *= freq;

          // integrate geodesic backwards
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
    }  // omp parallel

    // get information about "age" of geodesics
    double tmax = 0;
    int imax = 0;
    int jmax = 0;
    for (int i=0; i<NX; ++i) {
      for (int j=0; j<NY; ++j) {
        if (t[i][j] < tmax) {
          tmax = t[i][j];
          imax = i;
          jmax = j;
        }
      }
    }
    for (int i=0; i<NX; ++i) {
      for (int j=0; j<NY; ++j) {
        Xgeo[i][j][0] = t0 + t[i][j] - tmax;
      }
    }
    fprintf(stderr, "tmax = %g @ pixel [%d %d]\n", tmax, imax, jmax);

    // radiative transfer below

    int done[NX][NY];

#pragma omp parallel for collapse(2)
    for (int i=0; i<NX; ++i) {
      for (int j=0; j<NY; ++j) {
        done[i][j] = 0;
        image[i][j] = 0;
        init_N(Xgeo[i][j], Kcongeo[i][j], Ngeo[i][j]);
        tauFgeo[i][j] = 0;
      }
    }

    // TODO: note there is the assumption that DTd is constant.

    double tcurr = t0;
    int nloop = 0;

    while (tcurr < t0 - tmax) {

      fprintf(stderr, "\ntcurr = %g\n\n", tcurr);

#pragma omp parallel
      {

        // first loop over all pixels as much as can be done given the time
        // range that the loaded fluid data corresponds to
#pragma omp for collapse(2)
        for (int i=0; i<NX; ++i) {
          for (int j=0; j<NY; ++j) {

            double Xhalfgeo[NDIM], Kconhalfgeo[NDIM];
            double Xprevgeo[NDIM], Kconprevgeo[NDIM];

            while (Xgeo[i][j][0] < tcurr + DTd) {

              for (int mu = 0; mu < NDIM; mu++) {
                Xprevgeo[mu] = Xgeo[i][j][mu];
                Kconprevgeo[mu] = Kcongeo[i][j][mu];
              }

              double dl;

              dl = stepsize(Xgeo[i][j], Kcongeo[i][j]);
              push_photon(Xgeo[i][j], Kcongeo[i][j], dl, Xhalfgeo, Kconhalfgeo);

              if (done[i][j] == 0) {

                double ji, jf, ki, kf;

                get_jkinv(Xprevgeo, Kconprevgeo, &ji, &ki);
                get_jkinv(Xgeo[i][j], Kcongeo[i][j], &jf, &kf);
                dl *= L_unit*HPL/(ME*CL*CL);
                image[i][j] = approximate_solve(image[i][j], ji, ki, jf, kf, dl);

                evolve_N(Xprevgeo, Kconprevgeo, Xhalfgeo, Kconhalfgeo, Xgeo[i][j],
                    Kcongeo[i][j], dl, Ngeo[i][j], &tauFgeo[i][j]);

                if (isnan(creal(Ngeo[i][j][0][0]))) {
                  printf("NGEO NAN! %e\n", creal(Ngeo[i][j][0][0]));
                  printf("[%i %i]\n", i,j);
                  printf("X[] = %e %e %e %e\n", Xgeo[i][j][0], Xgeo[i][j][1], Xgeo[i][j][2], Xgeo[i][j][3]);
                  printf("K[] = %e %e %e %e\n", Kcongeo[i][j][0], Kcongeo[i][j][1], Kcongeo[i][j][2], Kcongeo[i][j][3]);
                  exit(-1);
                }

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
            } // should not continuing this geodesic because we're after the valid set
            // of loaded fluid data. this closes the while.

            if (Kcongeo[i][j][1] > 0. && Xgeo[i][j][1] > log(rmax)) {
              done[i][j] = 1;
            }
          }
        }

#pragma omp single

        // check if we've finished every pixel
        fprintf(stderr, "tcurr = %g\n", tcurr);
        int alldone = 1;
        for (int i=0; i<NX; ++i) {
          for (int j=0; j<NY; ++j) {
            if (done[i][j] == 0) alldone = 0;
          }
        }

        // update fluid
        if (!alldone) {
#pragma omp single
          fprintf(stderr, "LOADING NEW DATA\n");
          update_data();
        }

#pragma omp single
        tcurr += DTd;
#pragma omp single
        nloop++;

      } 
    } // omp parallel

    // conversion factors
    for (int i=0; i<NX; ++i) {
      for (int j=0; j<NY; ++j) {
        fprintf(stderr, "Xgeo[%i][%i] = %e %e %e %e\n", i,j, Xgeo[i][j][0], Xgeo[i][j][1], Xgeo[i][j][2], Xgeo[i][j][3]);
        image[i][j] *= pow(freqcgs, 3);
        project_N(Xgeo[i][j], Kcongeo[i][j], Ngeo[i][j], &imageS[i][j][0], &imageS[i][j][1], &imageS[i][j][2], &imageS[i][j][3]);
        imageS[i][j][0] *= pow(freqcgs, 3);
        imageS[i][j][1] *= pow(freqcgs, 3);
        imageS[i][j][2] *= pow(freqcgs, 3);
        imageS[i][j][3] *= pow(freqcgs, 3);
        imageS[i][j][4] = tauFgeo[i][j];
      }
    }

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

        //if ( X[1] < (log(1.05*Rh)) )
        //fprintf(stderr, "Stopped at %g %g %g\n", X[1], X[2], X[3]);

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

          if (flag == 0 && 0) {
            project_N(Xf, Kconf, N_coord, &Stokes_I, &Stokes_Q, &Stokes_U, &Stokes_V);
            //            if (Stokes_I < 0 || fabs(Intensity - Stokes_I) > 1e-50 || 1) {
            if ( Stokes_I < 0 || 1 ) { // || 1 ) {
              //if ( i == 10 ) {
              void Xtoijk(double Xf[NDIM], int *i, int *j, int *k, double del[NDIM]);
              double tdel[NDIM];
              int ti,tj,tk;
              Xtoijk(Xf, &ti, &tj, &tk, tdel);
              fprintf(stderr, "Stokes_I %g (%g) @ %d %d %d (%g %g %g %g) [%d/%d]\n", Stokes_I, Intensity, ti, tj, tk, Xf[0], Xf[1], Xf[2], Xf[3], nstep, nstep);
              flag = 0;
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
      X[i][j] = (i - NX/2 + 0.5)*DX/NX;
      Y[i][j] = (j - NY/2 + 0.5)*DY/NY;
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

  // Grid
  h5io_add_data_dbl_2d(fid, "/X", NX, NY, X);
  h5io_add_data_dbl_2d(fid, "/Y", NX, NY, Y);

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

