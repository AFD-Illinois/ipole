#include "decs.h"

#include "model.h"
#include "model_geodesics.h"
#include "model_radiation.h"
#include "model_tetrads.h"

#include "radiation.h"
#include "coordinates.h"
#include "tetrads.h"
#include "geometry.h"
#include "geodesics.h"
#include "image.h"
#include "io.h"
#include "ipolarray.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <omp.h>

// global variables. TODO scope into main
static double tf = 0.;

Params params = { 0 };

// Functions defined & used only here.  TODO move. Introduce a camera.c or similar?
void init_XK (int i, int j, double Xcam[4], double fovx, double fovy,
              double X[4], double Kcon[4], double rotcam, double xoff,
              double yoff);

double approximate_solve (double Ii, double ji, double ki, double jf, double kf,
                   double dl, double *tau);

// Try to keep dynamic image sizes fast
static int nx, ny;
static inline int imgindex(int n, int i, int j) {return (n*nx + i)*ny + j;}

int main(int argc, char *argv[]) 
{
  // motd
  fprintf(stderr, "ipole. githash: %s\n", xstr(VERSION));
  fprintf(stderr, "notes: %s\n\n", xstr(NOTES));

  // initialization
  double time = omp_get_wtime();

  double tA, tB; // for slow light
  double rcam, thetacam, phicam, rotcam, Xcam[NDIM];
  double xoff, yoff;
  double freq, scale;
  double DX, DY, fovx, fovy;

  // Options not needed outside this monster function
  int quench_output = 0;
  int only_unpolarized = 0;

#pragma omp parallel
  if (omp_get_thread_num() == 0) {
    fprintf(stderr, "nthreads = %d\n", omp_get_num_threads());
  }

  // load values from parameter file. handle all actual
  // model parameter comprehension in the model/* files
  load_par_from_argv(argc, argv, &params);

  // figure out if we should run in a custom mode
  // TODO handle with other parameters
  for (int i=0; i<argc; ++i) {
    if ( strcmp(argv[i], "-quench") == 0 ) quench_output = 1;
    else if ( strcmp(argv[i], "-unpol") == 0 ) only_unpolarized = 1;
  }

  // now that we've loaded all parameters, tell our model about
  // them and use init_model to load the first dump
  init_model(&tA, &tB);

  // Then allocate now we know the image size
  nx = params.nx;
  ny = params.ny;
  double *taus = malloc(sizeof(*taus) * nx*ny);
  double *imageS = malloc(sizeof(*imageS) * nx*ny*NIMG);
  double *image = malloc(sizeof(*image) * nx*ny);

  // normalize frequency to electron rest-mass energy
  double freqcgs = params.freqcgs;
  freq = params.freqcgs * HPL / (ME * CL * CL);

  // TODO if spherical
  // Initialize the camera
  rcam = params.rcam;
  thetacam = params.thetacam;
  phicam = params.phicam;
  rotcam = params.rotcam*M_PI/180.;
  xoff = params.xoff;
  yoff = params.yoff;
  // translate to geodesic coordinates
  double x[NDIM] = {0., rcam, thetacam/180.*M_PI, phicam/180.*M_PI};
  Xcam[0] = 0.0;
  Xcam[1] = log(rcam);
  Xcam[2] = root_find(x);
  Xcam[3] = phicam/180.*M_PI;
  fprintf(stderr, "Xcam[] = %e %e %e %e\n", Xcam[0], Xcam[1], Xcam[2], Xcam[3]);
  fprintf(stderr, "a=%g R0=%g hslope=%g\n", a, R0, hslope);

  params.dsource *= PC;
  double Dsource = params.dsource; // Shorthand

  // set DX/DY using fov_dsource if possible, otherwise DX, otherwise old default
  double fov_to_d = Dsource / L_unit / 2.06265e11;

  if (params.fovx_dsource != 0.0) { // FOV was specified
    // Uncomment to be even more option lenient
    //if (params.fovy_dsource == 0.0) params.fovy_dsource = params.fovx_dsource;
  } else if (params.dx != 0.0) {
    //if (params.dy == 0.0) params.dy = params.dx;
    params.fovx_dsource = params.dx / fov_to_d;
    params.fovy_dsource = params.dy / fov_to_d;
  } else {
    fprintf(stderr, "No FOV was specified. Using default 160muas!\n");
    params.fovx_dsource = 160.;
    params.fovy_dsource = 160.;
  }
  DX = params.fovx_dsource * fov_to_d;
  DY = params.fovy_dsource * fov_to_d;
  params.dx = DX;
  params.dy = DY;

  // Set the *camera* fov values
  // We don't set these like other parameters, but output them for historical reasons
  fovx = DX / rcam;
  fovy = DY / rcam;

  scale = (DX * L_unit / nx) * (DY * L_unit / ny) / (Dsource * Dsource) / JY;
  fprintf(stderr,"L_unit = %e DX = %e NX = %i Dsource = %e JY = %e\n", L_unit, DX, nx, Dsource,JY);
  fprintf(stderr,"intensity [cgs] to flux per pixel [Jy] conversion: %g\n",scale);
  fprintf(stderr,"Dsource: %g [cm]\n",Dsource);
  fprintf(stderr,"Dsource: %g [kpc]\n",Dsource/(1.e3*PC));
  fprintf(stderr,"FOVx, FOVy: %g %g [GM/c^2]\n",DX,DY);
  fprintf(stderr,"FOVx, FOVy: %g %g [rad]\n",DX*L_unit/Dsource,DY*L_unit/Dsource);
  fprintf(stderr,"FOVx, FOVy: %g %g [muas]\n",DX*L_unit/Dsource * 2.06265e11 ,DY*L_unit/Dsource * 2.06265e11);

  // slow light
  if (SLOW_LIGHT) {

    // used to track when to write restart files
    double next_restart_after = -1.;

    // maximum number of steps, needed for dtraj (saving geodesics) 
    int maxsteplength = -1;

    // minimum time, used to track the longest geodesic
    double t0 = 0.;

    // time when first geodesic first goes inside of R<100 ball.
    // the difference (t0 - tgeoi) is the minimum amount of time
    // required to be within the dumps. i.e., if the integration
    // must go beyond this range, it is not required to load new
    // dumps here. similar for tgeof, as the diagram below shows
    //
    //  rmax_geo      r=100.  r=100.                    r=rcam
    //
    //      |-----------|--------|----------------------> [cam]
    //      t0        tgeof    tgeoi                       t=0
    //
    double tgeoi = -1.e100;
    double tgeof = 0.;

    int nsteps[nx][ny];

    // first pass through geodesics to find lengths
#pragma omp parallel for schedule(dynamic,2) collapse(2) reduction(max:tgeoi) reduction(max:maxsteplength) reduction(min:t0) reduction(min:tgeof) shared(nsteps)
    for (int i=0; i<nx; ++i) {
      for (int j=0; j<ny; ++j) {
        if (j==0) fprintf(stderr, "%d ", i);

        int nstep = 0;
        double dl;
        double X[NDIM], Xhalf[NDIM], Kcon[NDIM], Kconhalf[NDIM];
        init_XK(i,j, Xcam, fovx,fovy, X, Kcon, rotcam, xoff, yoff);

        for (int k=0; k<NDIM; ++k) Kcon[k] *= freq;

        double tgeoitmp = 1.;
        double tgeoftmp = 1.;

        MULOOP Xhalf[mu] = X[mu];
        while (!stop_backward_integration(X, Xhalf, Kcon, Xcam)) {
          dl = stepsize(X, Kcon);
          push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
          nstep++;

          if (nstep > MAXNSTEP - 2) {
            fprintf(stderr, "MAXNSTEP exceeded on j=%d i=%d\n", j,i);
            break;
          }

          if (tgeoitmp > 0. && X[1] < log(100.)) tgeoitmp = X[0];
          if (tgeoftmp > 0. && X[1] > log(100.) && Kcon[1] < 0.) tgeoftmp = X[0]; // Kcon is for forward integration
        }

        nsteps[i][j] = nstep--;

        if (nstep > maxsteplength) maxsteplength = nstep;
        if (X[0] < t0) t0 = X[0];
        if (tgeoitmp < 0. && tgeoitmp > tgeoi) tgeoi = tgeoitmp;
        if (tgeoftmp < 0. && tgeoftmp < tgeof) tgeof = tgeoftmp;
        else if (tgeoftmp > 0. && X[0] < tgeof) tgeof = X[0];
      }
    }

    fprintf(stderr, "\n");

    fprintf(stderr, "%g %g %g\n", t0, tgeoi, tgeof);

    // now create data structures for geodesics and light rays
    maxsteplength += 2;
    fprintf(stderr, "geodesic size = %g GB\n", 1. * sizeof(struct of_traj) * nx*ny * maxsteplength / 1024/1024/1024);
    struct of_traj *dtraj = malloc(sizeof(*dtraj) * nx*ny * maxsteplength);

    // TODO: if DTd is different, we'll have to calcluate this a different way
    int nconcurrentimgs = 2. + 1. * fabs(t0) / params.img_cadence;
    fprintf(stderr, "images size = %g GB\n", 1. * sizeof(struct of_image) * nx*ny * nconcurrentimgs / 1024/1024/1024);
    struct of_image *dimage = malloc(sizeof(*dimage) * nx*ny * nconcurrentimgs);

    // now populate the geodesic data structures
#pragma omp parallel for schedule(dynamic,2) collapse(2)
    for (int i=0; i<nx; ++i) {
      for (int j=0; j<ny; ++j) {
        if (j==0) fprintf(stderr, "%d ", i);

        int nstep = 0;
        double dl;
        double X[NDIM], Xhalf[NDIM], Kcon[NDIM], Kconhalf[NDIM];
        init_XK(i,j, Xcam, fovx,fovy, X, Kcon, rotcam, xoff, yoff);
        for (int k=0; k<NDIM; ++k) Kcon[k] *= freq;

        MULOOP Xhalf[mu] = X[mu];
        while (!stop_backward_integration(X, Xhalf, Kcon, Xcam)) {
          dl = stepsize(X, Kcon);
          push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
          nstep++;

          int stepidx = imgindex(nstep,i,j);

          dtraj[stepidx].dl = dl * L_unit * HPL / (ME * CL * CL);
          for (int k=0; k<NDIM; ++k) {
            dtraj[stepidx].X[k] = X[k];
            dtraj[stepidx].Kcon[k] = Kcon[k];
            dtraj[stepidx].Xhalf[k] = Xhalf[k];
            dtraj[stepidx].Kconhalf[k] = Kconhalf[k];
          }

          if (nstep > MAXNSTEP - 2) {
            fprintf(stderr, "MAXNSTEP exceeded on j=%d i=%d\n", j,i);
            break;
          }
        }

      }
    }

    fprintf(stderr, "\n\nnow beginning radiative transfer calculation ...\n");

    // initialize state
    int nimg = 0, nopenimgs = 0;
    double last_img_target = tA - tgeof;
    double *target_times = malloc(sizeof(*target_times) * nconcurrentimgs);
    int *valid_images = malloc(sizeof(*valid_images) * nconcurrentimgs);
    for (int i=0; i<nconcurrentimgs; ++i) valid_images[i] = 0;

    fprintf(stderr, "first image will be produced for t = %g\n", last_img_target);

    // optionally start from restart file
    if( access("restart.h5", F_OK) != -1 ) {

      fprintf(stderr, "attempting to load restart file...\n");
      double ttA, ttB;
      read_restart("restart.h5", &ttA, &ttB, &last_img_target, &nopenimgs, &nimg,
            nconcurrentimgs, nconcurrentimgs*nx*ny,
            target_times, valid_images, dimage);
      update_data_until(&tA, &tB, ttA);
    }

    // now do radiative transfer
    while (1) {

      while ( last_img_target + t0 < tB ) {
        // add a new image if the final time for that image is before the final valid dump
        target_times[nimg] = last_img_target;
        if ( last_img_target + tgeoi < tf ) {
          valid_images[nimg] = 1;
          nopenimgs++;
          for (int i=0; i<nx; ++i) {
            for (int j=0; j<ny; ++j) {
              int pxidx = imgindex(nimg,i,j);
              dimage[pxidx].nstep = nsteps[i][j];
              dimage[pxidx].intensity = 0.;
              dimage[pxidx].tau = 0.;
              dimage[pxidx].tauF = 0.;
            }
          }
          if (++nimg == nconcurrentimgs) nimg = 0;
        }
        last_img_target += params.img_cadence;
      }

      for (int k=0; k<nconcurrentimgs; ++k) {  // images

        if (valid_images[k] == 0) continue;

        int do_output = 1;

#pragma omp parallel for schedule(dynamic,2) collapse(2)
        for (int i=0; i<nx; ++i) {
          for (int j=0; j<ny; ++j) {
            if (j==0) fprintf(stderr, "%d ", i);

            double ji,ki, jf,kf;
            double Xi[NDIM],Xhalf[NDIM],Xf[NDIM];
            double Kconi[NDIM],Kconhalf[NDIM],Kconf[NDIM];

            int pxidx = imgindex(k,i,j);

            while (dimage[pxidx].nstep > 1) {

              int stepidx = imgindex(dimage[pxidx].nstep,i,j);
              int pstepidx = imgindex(dimage[pxidx].nstep-1,i,j);

              for (int l=0; l<NDIM; ++l) {
                Xi[l]       = dtraj[stepidx].X[l];
                Xhalf[l]    = dtraj[stepidx].Xhalf[l];
                Xf[l]       = dtraj[pstepidx].X[l];
                Kconi[l]    = dtraj[stepidx].Kcon[l];
                Kconhalf[l] = dtraj[stepidx].Kconhalf[l];
                Kconf[l]    = dtraj[pstepidx].Kcon[l];
              }

              // adjust time. constant factor is to avoid floating point precision issues
              Xi[0] += target_times[k] + 1.e-5;
              Xhalf[0] += target_times[k] + 1.e-5;
              Xf[0] += target_times[k] + 1.e-5;

              // only integrate for points within the bounds
              if (Xi[0] < tA) {
                // this should only fire when we first start the program
                Xf[0] += tA - Xi[0];
                Xhalf[0] += tA - Xi[0];
                Xi[0] = tA;
              }
              if (Xi[0] >= tB) {
                if (Xf[0] >= tf) {
                  Xi[0] += tf - Xf[0];
                  Xhalf[0] += tf - Xf[0];
                  Xf[0] = tf;
                } else {
                  break;
                }
              }

              get_jkinv(Xi, Kconi, &ji, &ki);
              get_jkinv(Xf, Kconf, &jf, &kf);

              dimage[pxidx].intensity = 
                approximate_solve(dimage[pxidx].intensity, ji,ki,jf,kf, dtraj[stepidx].dl, &(dimage[pxidx].tau));
              
              // polarized transport
              if (! only_unpolarized) {
                evolve_N(Xi, Kconi, Xhalf, Kconhalf, Xf, Kconf, dtraj[stepidx].dl, dimage[pxidx].N_coord, &(dimage[pxidx].tauF));
                if (isnan(creal(dimage[pxidx].N_coord[0][0]))) {
                  exit(-2);
                }
              }

              dimage[pxidx].nstep -= 1;
            }

            if (dimage[pxidx].nstep != 1) {
              do_output = 0;
            }

          }
        }

        if (do_output) {

          // image, imageS, taus
          for (int i=0; i<nx; ++i) {
            for (int j=0; j<ny; ++j) {

              int pxidx = imgindex(k,i,j);
              int pstepidx = imgindex(dimage[pxidx].nstep,i,j);

              image[i*ny+j] = dimage[pxidx].intensity * pow(freqcgs, 3.);
              taus[i*ny+j] = dimage[pxidx].tau;

              if (! only_unpolarized) {
                double Stokes_I, Stokes_Q, Stokes_U, Stokes_V;
                project_N(dtraj[pstepidx].X, dtraj[pstepidx].Kcon, dimage[pxidx].N_coord, &Stokes_I, &Stokes_Q, &Stokes_U, &Stokes_V, rotcam);
                imageS[(i*ny+j)*NIMG+0] = Stokes_I * pow(freqcgs, 3.);
                imageS[(i*ny+j)*NIMG+1] = Stokes_Q * pow(freqcgs, 3.);
                imageS[(i*ny+j)*NIMG+2] = Stokes_U * pow(freqcgs, 3.);
                imageS[(i*ny+j)*NIMG+3] = Stokes_V * pow(freqcgs, 3.);
                imageS[(i*ny+j)*NIMG+4] = dimage[pxidx].tauF;
                if (params.qu_conv == 0) {
                  imageS[(i*ny+j)*NIMG+1] *= -1;
                  imageS[(i*ny+j)*NIMG+2] *= -1;
                }
              } else {
                imageS[(i*ny+j)*NIMG+0] = 0.;
                imageS[(i*ny+j)*NIMG+1] = 0.;
                imageS[(i*ny+j)*NIMG+2] = 0.;
                imageS[(i*ny+j)*NIMG+3] = 0.;
                imageS[(i*ny+j)*NIMG+4] = 0.;
              }
            }
          }

          //fprintf(stderr, "saving image %d at t = %g\n", k, target_times[k]);
          char dfname[256];
          snprintf(dfname, 255, params.outf, target_times[k]);
          dump(image, imageS, taus, dfname, scale, Xcam, fovx, fovy, &params);
          valid_images[k] = 0;
          nopenimgs--;
        }

      }

      if (nopenimgs < 1) break;
      update_data(&tA, &tB);

      // write a restart file however frequency as desired
      if (params.restart_int > 0. && tA > next_restart_after) {
        while (tA > next_restart_after) next_restart_after += params.restart_int;
        char rfname[256];
        snprintf(rfname, 200, "restarts/restart_%05.1f.h5", tA);
        write_restart(rfname, tA, tB, last_img_target, nopenimgs, nimg, 
              nconcurrentimgs, nconcurrentimgs*nx*ny,
              target_times, valid_images, dimage);
      }

    }

  // FAST LIGHT
  } else {

#pragma omp parallel for schedule(dynamic,1) collapse(2) shared(image,imageS)
    for (int i=0; i<nx; ++i) {
      for (int j=0; j<ny; ++j) {
        if (j==0) fprintf(stderr, "%d ", i);

        double X[NDIM],Xhalf[NDIM],Xi[NDIM],Xf[NDIM],Kcon[NDIM],Kconhalf[NDIM],Kconi[NDIM],Kconf[NDIM];
        double dl, ji,ki, jf,kf;
        double complex N_coord[NDIM][NDIM];
        double Intensity, Stokes_I, Stokes_Q, Stokes_U, Stokes_V, tauF, Tau;

        // Setup geodesic
        init_XK(i,j, Xcam, fovx,fovy, X, Kcon, rotcam, xoff, yoff);
        MULOOP Kcon[mu] *= freq;
        MULOOP Xhalf[mu] = X[mu];
        int nstep = 0;
        struct of_traj *traj = calloc(MAXNSTEP, sizeof(struct of_traj));

        // Integrate backwards
        while (!stop_backward_integration(X, Xhalf, Kcon, Xcam)) {
          /* This stepsize function can be troublesome inside of R = 2M,
             and should be used cautiously in this region. */
          dl = stepsize(X, Kcon);

          /* move photon one step backwards, the procecure updates X
             and Kcon full step and returns also values in the middle */
          push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
          nstep++;

          traj[nstep].dl = dl * L_unit * HPL / (ME * CL * CL);

          MULOOP {
            traj[nstep].X[mu] = X[mu];
            traj[nstep].Kcon[mu] = Kcon[mu];
            traj[nstep].Xhalf[mu] = Xhalf[mu];
            traj[nstep].Kconhalf[mu] = Kconhalf[mu];
          }

          if (nstep > MAXNSTEP - 2) {
            fprintf(stderr, "MAXNSTEP exceeded on j=%d i=%d\n", j,i);
            break;
          }
        }

        nstep--; // don't record final step because it violated "stop" condition
        //fprintf(stderr, "Geodesic i: %d j: %d nsteps: %d\n", i, j, nstep);

        // integrate forwards along trajectory, including radiative transfer equation
        // initialize X, K
        MULOOP { Xi[mu] = traj[nstep].X[mu];
                 Kconi[mu] = traj[nstep].Kcon[mu]; }
        // Initialize emision variables
        init_N(Xi, Kconi, N_coord);
        Intensity = 0.0;
        Tau = 0.;
        tauF = 0.;
        get_jkinv(traj[nstep].X, traj[nstep].Kcon, &ji, &ki);
        
        // Option to save a variable along the geodesic
        int nstep_save = nstep+1;

        while (nstep > 1) {
          // initialize X,K
          MULOOP {
            Xi[mu]       = traj[nstep].X[mu];
            Kconi[mu]    = traj[nstep].Kcon[mu];
            Xhalf[mu]    = traj[nstep].Xhalf[mu];
            Kconhalf[mu] = traj[nstep].Kconhalf[mu];
            Xf[mu]       = traj[nstep - 1].X[mu];
            Kconf[mu]    = traj[nstep - 1].Kcon[mu];
          }

          // solve total intensity equation alone
#if THIN_DISK
          if (thindisk_region(Xi, Xf)) {
            get_model_i(Xf, Kconf, &Intensity);
          }
#else
          get_jkinv(Xf, Kconf, &jf, &kf);
          Intensity = approximate_solve(Intensity, ji, ki, jf, kf, traj[nstep].dl, &Tau);
          // swap start and finish
          ji = jf;
          ki = kf;
#endif

          // solve polarized transport
          if (! only_unpolarized ) {
            evolve_N(Xi, Kconi, Xhalf, Kconhalf, Xf, Kconf, traj[nstep].dl, N_coord, &tauF);
            if (isnan(creal(N_coord[0][0]))) {
              exit(-3);
            }
          }

          nstep--;
        }

        if (params.trace) {
          int stride = params.trace_stride;
          if (params.trace_i < 0 || params.trace_j < 0) { // If no single point is specified
            if (i % stride == 0 && j % stride == 0) { // Save every stride pixels
#pragma omp critical
              dump_var_along(i/stride, j/stride, nstep_save, traj, nx/stride, ny/stride, scale, Xcam, fovx, fovy, &params);
            }
          } else {
            if (i == params.trace_i && j == params.trace_j) { // Save just the one
#pragma omp critical
              dump_var_along(0, 0, nstep_save, traj, 1, 1, scale, Xcam, fovx, fovy, &params);
            }
          }
        }

        free(traj);

        // deposit intensity and Stokes parameter in pixels
        image[i*ny+j] = Intensity * pow(freqcgs, 3);
        taus[i*ny+j] = Tau;
        if (! only_unpolarized) {
          project_N(Xf, Kconf, N_coord, &Stokes_I, &Stokes_Q, &Stokes_U, &Stokes_V, rotcam);
          imageS[(i*ny+j)*NIMG+0] = Stokes_I * pow(freqcgs, 3);
          imageS[(i*ny+j)*NIMG+1] = Stokes_Q * pow(freqcgs, 3);
          imageS[(i*ny+j)*NIMG+2] = Stokes_U * pow(freqcgs, 3);
          imageS[(i*ny+j)*NIMG+3] = Stokes_V * pow(freqcgs, 3);
          imageS[(i*ny+j)*NIMG+4] = tauF;
          if (params.qu_conv == 0) {
            imageS[(i*ny+j)*NIMG+1] *= -1;
            imageS[(i*ny+j)*NIMG+2] *= -1;
          }
          if (isnan(imageS[(i*ny+j)*NIMG+0])) exit(-1);
        } else {
          imageS[(i*ny+j)*NIMG+0] = 0.;
          imageS[(i*ny+j)*NIMG+1] = 0.;
          imageS[(i*ny+j)*NIMG+2] = 0.;
          imageS[(i*ny+j)*NIMG+3] = 0.;
          imageS[(i*ny+j)*NIMG+4] = 0.;
        }
      }
    }

    // finish up
    fprintf(stderr, "\nscale = %e\n", scale);

    double Ftot = 0.;
    double Ftot_unpol = 0.;
    double Imax = 0.0;
    double Iavg = 0.0;
    int imax = 0;
    int jmax = 0;
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        Ftot_unpol += image[i*ny+j]*scale;
        Ftot += imageS[(i*ny+j)*NIMG+0] * scale;
        Iavg += imageS[(i*ny+j)*NIMG+0];
        if (imageS[(i*ny+j)*NIMG+0] > Imax) {
          imax = i;
          jmax = j;
          Imax = imageS[(i*ny+j)*NIMG+0];
        }
      }
    }

    fprintf(stderr, "imax=%d jmax=%d Imax=%g Iavg=%g\n", imax, jmax, Imax, Iavg/(nx*ny));
    fprintf(stderr, "freq: %g Ftot: %g (%g unpol) scale=%g\n", freqcgs, Ftot, Ftot_unpol, scale);
    fprintf(stderr, "nuLnu = %g\n", 4.*M_PI*Ftot * Dsource * Dsource * JY * freqcgs);

    // don't dump if we've been asked to quench output. useful for batch jobs
    // like when fitting light curve fluxes
    if (!quench_output) {
      // dump result. if specified, also output ppm image
      dump(image, imageS, taus, params.outf, scale, Xcam, fovx, fovy, &params);
      if (params.add_ppm) {
        // TODO respect filename from params?
        make_ppm(image, freq, nx, ny, "ipole_lfnu.ppm");
      }
    }
  }

  time = omp_get_wtime() - time;
  printf("Total wallclock time: %g s\n\n", time);

  return 0;
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
    double X[NDIM], double Kcon[NDIM], double rotcam, double xoff, double yoff)
{
  double Econ[NDIM][NDIM];
  double Ecov[NDIM][NDIM];
  double Kcon_tetrad[NDIM];

  make_camera_tetrad(Xcam, Econ, Ecov);

  // Construct outgoing wavevectors
  // xoff: allow arbitrary offset for e.g. ML training imgs
  // +0.5: project geodesics from px centers
  // -0.01: Prevent nasty artifact at 0 spin/phicam
  double dxoff = (xoff+i+0.5-0.01)/nx - 0.5;
  double dyoff = (yoff+j+0.5)/ny - 0.5;
  Kcon_tetrad[0] = 0.;
#if THIN_DISK
  // TODO correct this more generally and take as parameter
  Kcon_tetrad[1] = (dxoff*cos(rotcam) - dyoff*sin(rotcam)) * fovx  - (0.659986/40*fovx);
#else
  Kcon_tetrad[1] = (dxoff*cos(rotcam) - dyoff*sin(rotcam)) * fovx;
#endif
  Kcon_tetrad[2] = (dxoff*sin(rotcam) + dyoff*cos(rotcam)) * fovy;
  Kcon_tetrad[3] = 1.;

  /* normalize */
  null_normalize(Kcon_tetrad, 1.);

  /* translate into coordinate frame */
  tetrad_to_coordinate(Econ, Kcon_tetrad, Kcon);

  /* set position */
  for (int mu = 0; mu < NDIM; mu++) X[mu] = Xcam[mu];
}

/* must be a stable, approximate solution to radiative transfer
 * that runs between points w/ initial intensity I, emissivity
 * ji, opacity ki, and ends with emissivity jf, opacity kf.
 *
 * return final intensity
 */
double approximate_solve(double Ii, double ji, double ki, double jf,
    double kf, double dl, double *tau)
{
  double efac, If, javg, kavg, dtau;

  javg = (ji + jf) / 2.;
  kavg = (ki + kf) / 2.;

  dtau = dl * kavg;
  *tau += dtau;

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

