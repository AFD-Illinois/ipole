#include "decs.h"

#include "model.h"
#include "model_geodesics.h"
#include "model_radiation.h"
#include "model_tetrads.h"

#include "radiation.h"
#include "coordinates.h"
#include "debug_tools.h"
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

// Print a backtrace of a certain pixel, useful for debugging
// Set to -1 to disable
#define DIAG_PX_I -1
#define DIAG_PX_J -1

// When running in check mode, we downsample the image by
// some skip amount. Set that skip here.
#define CHECK_DOWNSAMPLE_SKIP 10

// Some useful blocks of code to re-use
// Note the difference between "int nx,ny" in get_pixel and "long int nx,ny" in save_pixel
// explained in parameter parsing.
void get_pixel(size_t i, size_t j, int nx, int ny, double Xcam[NDIM], Params params,
               double fovx, double fovy, double freq, int only_intensity, double scale,
               double *Intensity, double *Is, double *Qs, double *Us, double *Vs,
               double *Tau, double *tauF);
void save_pixel(double *image, double *imageS, double *taus, size_t i, size_t j, 
                size_t nx, size_t ny, int only_unpol,
                double Intensity, double Is, double Qs, double Us, double Vs,
                double freqcgs, double Tau, double tauF);
void print_image_stats(double *image, double *imageS, size_t nx, size_t ny, Params params, double scale);
void save_pixelTransfer(double *image, double *imageS, double *taus,
                        size_t iold, size_t jold, size_t inew, size_t jnew, size_t nx, size_t ny, int only_intensity); //nearest neighbor saving
void lininterp4(double *image, double *imageS, double *taus, size_t i1, size_t j1,
                size_t i2,size_t j2, size_t i3, size_t j3, size_t i4, size_t j4, size_t inew, size_t jnew,
                size_t nx, size_t ny, int only_intensity); //linear interpolation for floater case 
void lininterp2(double *image, double *imageS, double *taus, size_t i1, size_t j1,
                size_t i2,size_t j2,size_t inew, size_t jnew, size_t nx, size_t ny, int only_intensity); //linear interpolation for same row or column case


// TODO use this more.  Are long=only or upscaled ops faster?
static inline size_t imgindex(size_t n, size_t i, size_t j, size_t nx, size_t ny) {return (n*nx + i)*ny + j;}

// global variables. TODO scope into main
static double tf = 0.;

Params params = { 0 };

int main(int argc, char *argv[]) 
{
  // motd
  fprintf(stderr, "%s. githash: %s\n", VERSION_STRING, xstr(VERSION));
  fprintf(stderr, "notes: %s\n\n", xstr(NOTES));

  // initialization
  double time = omp_get_wtime();

  double tA, tB; // for slow light
  double Xcam[NDIM];
  double freq, scale;
  double DX, DY, fovx, fovy;

#pragma omp parallel
  if (omp_get_thread_num() == 0) {
    fprintf(stderr, "nthreads = %d\n", omp_get_num_threads());
  }

  // load values from parameter file. handle all actual
  // model parameter comprehension in the model/* files
  load_par_from_argv(argc, argv, &params);

  // now that we've loaded all parameters, tell our model about
  // them and use init_model to load the first dump
  init_model(&tA, &tB);

  // Adaptive resolution option
  // nx, ny are the resolution at maximum refinement level
  // nx_min, ny_min are at coarsest level
  // TODO check for obvious BS here
  if (params.nx_min < 0) {
    params.nx_min = params.nx;
    params.ny_min = params.ny;
  }

  int refine_level = log2((params.nx - params.nx % 2)/(params.nx_min - params.nx_min % 2))+1;
  // INTERNAL SIZE.  If nx or ny is even, compute an extra row/column to use the 2^N+1 scheme
  size_t nx, ny, nxmin, nymin;
  if (refine_level > 1 && params.nx % 2 == 0) {
    nx = params.nx + 1;
    nxmin = params.nx_min + 1;
  } else {
    nx = params.nx;
    nxmin = params.nx_min;
  }
  if (refine_level > 1 && params.ny % 2 == 0) {
    ny = params.ny + 1;
    nymin = params.ny_min + 1;
  } else {
    ny = params.ny;
    nymin = params.ny_min;
  }

  // normalize frequency to electron rest-mass energy
  double freqcgs = params.freqcgs;
  freq = params.freqcgs * HPL / (ME * CL * CL);
  if (freq == 0) {
    fprintf(stderr, "Frequency cannot be zero. Quitting!\n");
    exit(1);
  }

  // Initialize the camera
  params.rotcam *= M_PI/180.;

  // translate to geodesic coordinates
  native_coord(params.rcam, params.thetacam, params.phicam, Xcam);
  fprintf(stderr, "Xcam[] = %e %e %e %e\n", Xcam[0], Xcam[1], Xcam[2], Xcam[3]);

  params.dsource *= PC;
  double Dsource = params.dsource; // Shorthand

  // set DX/DY using fov_dsource if possible, otherwise DX, otherwise old default
  double fov_to_d = Dsource / L_unit / MUAS_PER_RAD;

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
  fovx = DX / params.rcam;
  fovy = DY / params.rcam;

  scale = (DX * L_unit / params.nx) * (DY * L_unit / params.ny) / (Dsource * Dsource) / JY;
  fprintf(stderr,"intensity [cgs] to flux per pixel [Jy] conversion: %g\n",scale);
  fprintf(stderr,"Dsource: %g [cm]\n",Dsource);
  fprintf(stderr,"Dsource: %g [kpc]\n",Dsource/(1.e3*PC));
  fprintf(stderr,"FOVx, FOVy: %g %g [GM/c^2]\n",DX,DY);
  fprintf(stderr,"FOVx, FOVy: %g %g [rad]\n",DX*L_unit/Dsource,DY*L_unit/Dsource);
  fprintf(stderr,"FOVx, FOVy: %g %g [muas]\n",DX*L_unit/Dsource * MUAS_PER_RAD ,DY*L_unit/Dsource * MUAS_PER_RAD);
  if (refine_level > 1) {
    fprintf(stderr,"Resolution: %dx%d, refined up to %dx%d (%d levels)\n",
        params.nx_min, params.ny_min, params.nx, params.ny, refine_level);
    fprintf(stderr,"Refinement when relative error is >%g%% or absolute error is >%g%% of estimated total flux, for pixel brightness > %f%% of average\n",
            params.refine_rel*100, params.refine_abs*100, params.refine_cut*100);
  } else {
    fprintf(stderr,"Resolution: %dx%d\n", params.nx, params.ny);
  }

  double *taus = calloc(nx*ny, sizeof(*taus));
  double *image = calloc(nx*ny, sizeof(*image));
  double *imageS = NULL;
  if (taus == NULL || image == NULL) {
    fprintf(stderr, "Could not allocate image memory!");
    exit(-1);
  }
  if(!params.only_unpolarized) {
    imageS = calloc(nx*ny*NIMG, sizeof(*imageS));
    if (imageS == NULL) {
      fprintf(stderr, "Could not allocate image memory!");
      exit(-1);
    }
  }

  if (params.perform_check) {
    fprintf(stderr, "Running in check mode...\n");
  }

  double initialization_time = omp_get_wtime() - time;

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

    // TODO disallows large slow-light images
    int nsteps[nx][ny];

    // first pass through geodesics to find lengths
#pragma omp parallel for schedule(dynamic,2) collapse(2) reduction(max:tgeoi) reduction(max:maxsteplength) reduction(min:t0) reduction(min:tgeof) shared(nsteps)
    for (size_t i=0; i<nx; ++i) {
      for (size_t j=0; j<ny; ++j) {
        if (j==0) fprintf(stderr, "%ld ", i);

        int nstep = 0;
        double dl;
        double X[NDIM], Xhalf[NDIM], Kcon[NDIM], Kconhalf[NDIM];
        init_XK(i,j, params.nx, params.ny, Xcam, params, fovx, fovy, X, Kcon);

        MULOOP Kcon[mu] *= freq;

        double tgeoitmp = 1.;
        double tgeoftmp = 1.;

        MULOOP Xhalf[mu] = X[mu];
        while (!stop_backward_integration(X, Xhalf, Kcon)) {
          dl = stepsize(X, Kcon, params.eps);
          push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
          nstep++;

          if (nstep > params.maxnstep - 2) {
            fprintf(stderr, "maxnstep exceeded on j=%ld i=%ld\n", j,i);
            break;
          }

          if (tgeoitmp > 0. && X[1] < log(100.)) tgeoitmp = X[0];
          if (tgeoftmp > 0. && X[1] > log(100.) && Kcon[1] < 0.) tgeoftmp = X[0]; // Kcon is for forward integration
        }

        nsteps[i][j] = --nstep;

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
    struct of_traj *dtraj = calloc(nx*ny * maxsteplength, sizeof(*dtraj));

    // TODO: if DTd is different, we'll have to calcluate this a different way
    int nconcurrentimgs = 2. + 1. * fabs(t0) / params.img_cadence;
    fprintf(stderr, "images size = %g GB\n", 1. * sizeof(struct of_image) * nx*ny * nconcurrentimgs / 1024/1024/1024);
    struct of_image *dimage = calloc(nx*ny * nconcurrentimgs, sizeof(*dimage));

    // now populate the geodesic data structures
#pragma omp parallel for schedule(dynamic,2) collapse(2)
    for (size_t i=0; i<nx; ++i) {
      for (size_t j=0; j<ny; ++j) {
        if (j==0) fprintf(stderr, "%ld ", i);

        int nstep = 0;
        double dl;
        double X[NDIM], Xhalf[NDIM], Kcon[NDIM], Kconhalf[NDIM];
        init_XK(i,j, params.nx, params.ny, Xcam, params, fovx,fovy, X, Kcon);
        for (int k=0; k<NDIM; ++k) Kcon[k] *= freq;

        MULOOP Xhalf[mu] = X[mu];
        while (!stop_backward_integration(X, Xhalf, Kcon)) {
          dl = stepsize(X, Kcon, params.eps);
          push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
          nstep++;

          int stepidx = imgindex(nstep,i,j,nx,ny);

          dtraj[stepidx].dl = dl * L_unit * HPL / (ME * CL * CL);
          for (int k=0; k<NDIM; ++k) {
            dtraj[stepidx].X[k] = X[k];
            dtraj[stepidx].Kcon[k] = Kcon[k];
            dtraj[stepidx].Xhalf[k] = Xhalf[k];
            dtraj[stepidx].Kconhalf[k] = Kconhalf[k];
          }

          if (nstep > params.maxnstep - 2) {
            fprintf(stderr, "maxnstep exceeded on j=%ld i=%ld\n", j,i);
            break;
          }
        }

      }
    }

    fprintf(stderr, "\n\nnow beginning radiative transfer calculation ...\n");

    // initialize state
    int nimg = 0, nopenimgs = 0;
    double last_img_target = tA - tgeof;
    double *target_times = calloc(nconcurrentimgs, sizeof(*target_times));
    int *valid_images = calloc(nconcurrentimgs, sizeof(*valid_images));
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
              int pxidx = imgindex(nimg,i,j,nx,ny);
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
        for (size_t i=0; i<nx; ++i) {
          for (size_t j=0; j<ny; ++j) {
            if (j==0) fprintf(stderr, "%ld ", i);

            double ji,ki, jf,kf;
            double Xi[NDIM],Xhalf[NDIM],Xf[NDIM];
            double Kconi[NDIM],Kconhalf[NDIM],Kconf[NDIM];

            int pxidx = imgindex(k,i,j,nx,ny);

            while (dimage[pxidx].nstep > 1) {

              int stepidx = imgindex(dimage[pxidx].nstep,i,j,nx,ny);
              int pstepidx = imgindex(dimage[pxidx].nstep-1,i,j,nx,ny);

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

              get_jkinv(Xi, Kconi, &ji, &ki, &params);
              get_jkinv(Xf, Kconf, &jf, &kf, &params);

              dimage[pxidx].intensity = 
                approximate_solve(dimage[pxidx].intensity, ji,ki,jf,kf, dtraj[stepidx].dl, &(dimage[pxidx].tau));
              
              // polarized transport
              if (! params.only_unpolarized) {
                evolve_N(Xi, Kconi, Xhalf, Kconhalf, Xf, Kconf, dtraj[stepidx].dl, dimage[pxidx].N_coord, &(dimage[pxidx].tauF), 0, &params);
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

              int pxidx = imgindex(k,i,j,nx,ny);
              int pstepidx = imgindex(dimage[pxidx].nstep,i,j,nx,ny);

              image[i*ny+j] = dimage[pxidx].intensity * pow(freqcgs, 3.);
              taus[i*ny+j] = dimage[pxidx].tau;

              if (! params.only_unpolarized) {
                double Stokes_I, Stokes_Q, Stokes_U, Stokes_V;
                project_N(dtraj[pstepidx].X, dtraj[pstepidx].Kcon, dimage[pxidx].N_coord, &Stokes_I, &Stokes_Q, &Stokes_U, &Stokes_V, params.rotcam);
                imageS[(i*ny+j)*NIMG+0] = Stokes_I * pow(freqcgs, 3.);
                imageS[(i*ny+j)*NIMG+1] = Stokes_Q * pow(freqcgs, 3.);
                imageS[(i*ny+j)*NIMG+2] = Stokes_U * pow(freqcgs, 3.);
                imageS[(i*ny+j)*NIMG+3] = Stokes_V * pow(freqcgs, 3.);
                imageS[(i*ny+j)*NIMG+4] = dimage[pxidx].tauF;
                if (params.qu_conv == 0) {
                  imageS[(i*ny+j)*NIMG+1] *= -1;
                  imageS[(i*ny+j)*NIMG+2] *= -1;
                }
              }
            }
          }

          //fprintf(stderr, "saving image %d at t = %g\n", k, target_times[k]);
          char dfname[256] = {0};
          snprintf(dfname, 255, params.outf, target_times[k]);
          dump(image, imageS, taus, dfname, scale, Xcam, fovx, fovy, nx, ny, &params);
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
    // BASE IMAGE at n_min
    // Allocate it, or use the existing allocation for just 1 level
    size_t initialspacingx = (nx - 1) / (nxmin - 1);
    size_t initialspacingy = (ny - 1) / (nymin - 1);

#if DEBUG
    fprintf(stderr, "Image dimensions: %d %d, memory dimensions %ld %ld, minimum %ld %ld\n",
           params.nx, params.ny, nx, ny, nxmin, nymin);

    fprintf(stderr, "Intial spacing: %ld\n", initialspacingx);
#endif

    int *interp_flag = calloc(nx * ny, sizeof(*interp_flag));
    double *prelimarray = NULL;
    if (interp_flag == NULL) {
      fprintf(stderr, "Could not allocate adaptive memory!\n");
      exit(-1);
    }
    if (refine_level > 1) {
      prelimarray = calloc(nxmin * nymin, sizeof(*prelimarray));
      if (prelimarray == NULL) {
        fprintf(stderr, "Could not allocate adaptive memory!\n");
        exit(-1);
      }
    }

    // Get the "base image" -- if not adaptively refining, this is just the whole image
    // Note the average of the interpolated image != the average of calculated pixels, though they're close.
#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (size_t i = 0; i < nx; i += initialspacingx) {
      for (size_t j = 0; j < ny; j += initialspacingy) {

        if (params.perform_check) {
          if ( i % CHECK_DOWNSAMPLE_SKIP != 0 || j % CHECK_DOWNSAMPLE_SKIP != 0 ) {
            continue;
          }
        }

        if (j==0) fprintf(stderr, "%ld ", i);
        size_t thislocation = i / initialspacingy * nymin + j / initialspacingx;
        double Intensity = 0;
        double Is = 0, Qs = 0, Us = 0, Vs = 0;
        double Tau = 0, tauF = 0;

        get_pixel(i, j, params.nx, params.ny, Xcam, params, fovx, fovy, freq,
                   params.only_unpolarized, scale, &Intensity, &Is, &Qs, &Us, &Vs,
                   &Tau, &tauF);

        if (params.perform_check) {

          save_pixel(image, imageS, taus, i/CHECK_DOWNSAMPLE_SKIP, j/CHECK_DOWNSAMPLE_SKIP, nx, ny, 
                     params.only_unpolarized, Intensity, Is, Qs, Us, Vs, freqcgs, Tau, tauF);
      
          int offset = nx % CHECK_DOWNSAMPLE_SKIP;
          if (offset == 0) offset = CHECK_DOWNSAMPLE_SKIP;

          int i2 = nx - i - offset;

          Intensity = 0;
          Is = 0; Qs = 0; Us = 0; Vs = 0;
          Tau = 0; tauF = 0;

          get_pixel(i2, j, params.nx, params.ny, Xcam, params, fovx, fovy, freq,
                     params.only_unpolarized, scale, &Intensity, &Is, &Qs, &Us, &Vs,
                     &Tau, &tauF);

          save_pixel(image, imageS, taus, (int)ceil(nx/CHECK_DOWNSAMPLE_SKIP) + i2/CHECK_DOWNSAMPLE_SKIP, j/CHECK_DOWNSAMPLE_SKIP, nx, ny, 
                     params.only_unpolarized, Intensity, Is, Qs, Us, Vs, freqcgs, Tau, tauF);

        } else {

          save_pixel(image, imageS, taus, i, j, nx, ny, params.only_unpolarized, Intensity, Is, Qs, Us,
                      Vs, freqcgs, Tau, tauF);

        }

        interp_flag[i*ny+j] = 0;

        if (refine_level > 1) {
          // computes a total interpolated flux from first pass
          // adds the total number of pixels adjacent to this one
          // assumes nx=ny and nx_min=ny_min
          
          if (i % (nx - 1) != 0 && j % (ny - 1) != 0) {
            // middle case
            prelimarray[thislocation] = Intensity * initialspacingx * initialspacingx;
          } else if ((i == 0 && j % (ny - 1) != 0) || (i % (nx - 1) != 0 && j == 0)) {
            // bottom or left vertical edge
            prelimarray[thislocation] = Intensity * (initialspacingx / 2 + 1) * initialspacingx;
          } else if (i % (nx - 1) != 0 || j % (ny - 1) != 0) {
            // top or right vertical edge
            prelimarray[thislocation] = Intensity * (initialspacingx / 2) * initialspacingx;
          } else {
            // corner case
            if (i == 0 && j == 0) {
              // bottom left corner
              prelimarray[thislocation] = Intensity * (initialspacingx / 2 + 1) * (initialspacingx / 2 + 1);
            } else if (i == 0 || j == 0) {
              // bottom right or top left
              prelimarray[thislocation] = Intensity * (initialspacingx / 2 + 1) * initialspacingx / 2;
            } else {
              // top right
              prelimarray[thislocation] = Intensity * (initialspacingx * initialspacingx) / 4;
            }
          }
        }
      }
    }

    // compute estimated flux total and intensity average
    double Iavg = 0;
    if (refine_level > 1) {
      double interp_tot = 0;
      for (size_t i = 0; i < nxmin * nymin; i++) {
        interp_tot += prelimarray[i];
      }

      // Average intensity per pixel
      Iavg = interp_tot * pow(freqcgs, 3) / (nx * ny);

      // Print calculated total intensity for debug
#if DEBUG
      fprintf(stderr, "\nInitial flux guess: %g", Iavg * (params.nx * params.ny) * scale);
#endif
      fprintf(stderr, "\n\n"); // TODO even for non-refined?
    }

    for (int refined_level = 1; refined_level < refine_level; refined_level++) {
      size_t newspacingx = initialspacingx / pow(2, refined_level);
      size_t newspacingy = initialspacingy / pow(2, refined_level);
      fprintf(stderr, "Refining level %d of %d, spacing %ld,%ld\n", refined_level+1, refine_level, newspacingx, newspacingy);

#pragma omp parallel for schedule(dynamic,1) collapse(2)
      for (size_t i = 0; i < nx; i += newspacingx) {
        for (size_t j = 0; j < ny; j += newspacingy) {

          if (j == 0) fprintf(stderr, "%ld ", i);

          double Intensity = 0;
          double Is = 0, Qs = 0, Us = 0, Vs = 0;
          double Tau = 0, tauF = 0;

          size_t previousspacingx = newspacingx * 2;
          size_t previousspacingy = newspacingy * 2;

          double I1, I2, I3, I4, err_abs, err_rel;
          if (i % previousspacingx == 0 && j % previousspacingy == 0) {
            // pixel has already been ray-traced
            continue;
          } else if (i % previousspacingx == 0 && j % previousspacingy != 0) {
            // pixel lies on pre-existing column

            I1 = image[i*ny+j-newspacingy]; // below
            I2 = image[i*ny+j+newspacingy]; // above
            err_abs = (I2 - I1) / 2 / Iavg;
            err_rel = (I2 - I1) / 2 / I1;

            if ((fabs(err_abs) > params.refine_abs && // could be changed to || if wanted
                fabs(err_rel) > params.refine_rel)
                && fabs(I1) > params.refine_cut) {
              // ray trace (tolerances exceeded)

              get_pixel(i, j, params.nx, params.ny, Xcam, params, fovx, fovy, freq,
                         params.only_unpolarized, scale, &Intensity, &Is, &Qs, &Us,
                         &Vs, &Tau, &tauF);

              save_pixel(image, imageS, taus, i, j, nx, ny, params.only_unpolarized, Intensity, Is,
                          Qs, Us, Vs, freqcgs, Tau, tauF);

              interp_flag[i*ny+j] = 0;
                  
              } else {
              // interpolate
              interp_flag[i*ny+j] = 1;
              if (params.nearest_neighbor) {
                // nearest interpolation
                save_pixelTransfer(image, imageS, taus, i, j - newspacingy, i,
                                    j, nx, ny, params.only_unpolarized);
                // fills in with the nearest neighbor (choosing one side)
              } else {
                // linear
                lininterp2(image, imageS, taus, i, j - newspacingy, i,
                            j + newspacingy, i, j, nx, ny, params.only_unpolarized);

              }
            }
          } else if (i % previousspacingx != 0 && j % previousspacingy == 0) {
            // pixel lies on pre-existing row
            I1 = image[(i-newspacingx)*ny+j]; //left
            I2 = image[(i+newspacingx)*ny+j]; //right
            err_abs = (I2 - I1) / 2 / Iavg;
            err_rel = (I2 - I1) / 2 / I1;

            if ((fabs(err_abs) > params.refine_abs && //could be changed back to || if wanted
                fabs(err_rel) > params.refine_rel)
                && fabs(I1) > params.refine_cut) {
              // ray trace (tolerances exceeded)

              get_pixel(i, j, params.nx, params.ny, Xcam, params, fovx, fovy, freq,
                         params.only_unpolarized, scale, &Intensity, &Is, &Qs, &Us,
                         &Vs, &Tau, &tauF);

              save_pixel(image, imageS, taus, i, j, nx, ny, params.only_unpolarized, Intensity, Is,
                          Qs, Us, Vs, freqcgs, Tau, tauF);

              interp_flag[i*ny+j] = 0;

            } else {
              // interpolate
              interp_flag[i*ny+j] = 1;
              if (params.nearest_neighbor) {
                // nearest interpolation
                save_pixelTransfer(image, imageS, taus, i - newspacingx, j, i,
                                    j, nx, ny, params.only_unpolarized);
                // fills in with the nearest neighbor (choosing one side)
              } else {
                // linear
                lininterp2(image, imageS, taus, i - newspacingx, j,
                            i + newspacingx, j, i, j, nx, ny, params.only_unpolarized);
              }
            }
          } else {
            // pixel lies equidistant from four corners
            I1 = image[(i-newspacingx)*ny+j-newspacingy]; // bottom left
            I2 = image[(i+newspacingx)*ny+j-newspacingy]; // bottom right
            I3 = image[(i-newspacingx)*ny+j+newspacingy]; // upper left
            I4 = image[(i+newspacingx)*ny+(j+newspacingy)]; // upper right
            (void)I4; // silence unused warning

            // Refinement criterion thanks to Zack Gelles: absolute & relative error of
            // central corner under nearest-neighbor, estimated by Taylor expanding at lower-left pixel
            // Make sure absolute error is in Jy/muas^2

            double err_abs = ((I2 + I3) / 2 - I1) / Iavg;
            double err_rel = (I2 + I3) / (2 * I1) - 1.;

            if ((fabs (err_abs) > params.refine_abs && //could be changed to && if wanted
                fabs (err_rel) > params.refine_rel)
                && fabs (I1) > params.refine_cut) {

              // ray trace (tolerances exceeded)

              get_pixel(i, j, params.nx, params.ny, Xcam, params, fovx, fovy, freq,
                         params.only_unpolarized, scale, &Intensity, &Is, &Qs, &Us,
                         &Vs, &Tau, &tauF);

              save_pixel(image, imageS, taus, i, j, nx, ny, params.only_unpolarized, Intensity, Is,
                          Qs, Us, Vs, freqcgs, Tau, tauF);

              interp_flag[i*ny+j] = 0;

            } else {
              // interpolate
              interp_flag[i*ny+j] = 1;
              if (params.nearest_neighbor) {
                // nearest interpolation
                save_pixelTransfer(image, imageS, taus, i - newspacingx,
                                    j - newspacingy, i, j, nx, ny, params.only_unpolarized);
                // fills in with the nearest neighbor (choosing one side)
              } else {
                // linear
                lininterp4(image, imageS, taus, i - newspacingx,
                            j - newspacingy, i + newspacingx, j - newspacingy,
                            i - newspacingx, j + newspacingy, i + newspacingx,
                            j + newspacingy, i, j, nx, ny, params.only_unpolarized);
              }
            }
          }
        }
      }

      // Print how many pixels were interpolated
      size_t total_interpolated = 0;
#pragma omp parallel for collapse(2) reduction(+:total_interpolated)
      for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
          total_interpolated += interp_flag[i*ny+j];
        }
      }
      // Report interpolation stats vs the number of computed, relevant pixels, not the total
      size_t nx_level = params.nx/newspacingx;
      size_t ny_level = params.ny/newspacingy;
      fprintf(stderr, "\n%ld of %ld (%f%%) of computed pixels at %ldx%ld were interpolated\n\n",
              total_interpolated, nx_level * ny_level, ((double) total_interpolated) / (nx_level * ny_level) * 100, nx_level, ny_level);
    }

    // TODO print only for "real" pixels
    print_image_stats(image, imageS, nx, ny, params, scale);

    // don't dump if we've been asked to quench output. useful for batch jobs
    // like when fitting light curve fluxes
    if (!params.quench_output) {
      // dump result. if specified, also output ppm image
      if (params.perform_check) {
        dump_check(image, imageS, taus, &params, CHECK_DOWNSAMPLE_SKIP);
      } else {
        dump(image, imageS, taus, params.outf, scale, Xcam, fovx, fovy, nx, ny, &params);
      }
      if (params.add_ppm) {
        // TODO respect filename from params?
        make_ppm(image, freq, nx, ny, "ipole_lfnu.ppm");
      }
    }
  } // SLOW_LIGHT

  time = omp_get_wtime() - time;
  fprintf(stderr, "Total wallclock time: %g s (%g s)\n\n", time, initialization_time);

  return 0;
}

// TODO Move these?
void get_pixel(size_t i, size_t j, int nx, int ny, double Xcam[NDIM], Params params,
               double fovx, double fovy, double freq, int only_intensity, double scale,
               double *Intensity, double *Is, double *Qs, double *Us, double *Vs,
               double *Tau, double *tauF)
{
  double X[NDIM], Kcon[NDIM];
  double complex N_coord[NDIM][NDIM];

  // Integrate backward to find geodesic trajectory
  init_XK(i,j, params.nx, params.ny, Xcam, params, fovx, fovy, X, Kcon);
  struct of_traj *traj = calloc(params.maxnstep, sizeof(struct of_traj));
#if !INTEGRATOR_TEST
  MULOOP Kcon[mu] *= freq;
#endif

  int nstep = trace_geodesic(X, Kcon, traj, params.eps, params.maxnstep, Xcam, 0);
  if (nstep >= params.maxnstep-1) {
    // You almost certainly don't want to continue if this happens
    fprintf(stderr, "\nMaxNStep exceeded in pixel %ld %ld!\n", i, j);
    exit(-10);
  }

  // Integrate emission forward along trajectory
  int oddflag = integrate_emission(traj, nstep, Intensity, Tau, tauF, N_coord, &params);

  if (!only_intensity) {
    project_N(X, Kcon, N_coord, Is, Qs, Us, Vs, params.rotcam);
  }

  // Catch anything with bad tetrads, anything we manually specify, and signficantly negative pixels
  if (oddflag != 0 || (i == DIAG_PX_I && j == DIAG_PX_J) || *Is < -1.e-10) {
    fprintf(stderr, "\nOddity in pixel %ld %ld (flag %d):\n", i, j, oddflag);
    //print_vector("Starting X", X);
    print_vector("Starting Kcon", Kcon);
    fprintf(stderr, "nstep: %d\n", nstep);
    fprintf(stderr, "Final Stokes parameters: [%g %g %g %g]\n", *Is, *Qs, *Us, *Vs);
  }

  // Record values along the geodesic if requested
  // TODO this is most likely not compatible with adaptive mode
  if (params.trace) {
    size_t stride = params.trace_stride;
    if (params.trace_i < 0 || params.trace_j < 0) { // If no single point is specified
      if (i % stride == 0 && j % stride == 0) { // Save every stride pixels
#pragma omp critical
        dump_var_along(i/stride, j/stride, nstep, traj, params.nx/stride, params.ny/stride, scale, Xcam, fovx, fovy, &params);
      }
    } else {
      if (i == params.trace_i && j == params.trace_j) { // Save just the one
#pragma omp critical
        {
          dump_var_along(0, 0, nstep, traj, 1, 1, scale, Xcam, fovx, fovy, &params);
        }
      }
    }
  }

  free(traj);
}

void save_pixel(double *image, double *imageS, double *taus, size_t i, size_t j, size_t nx, size_t ny, int only_intensity,
                double Intensity, double Is, double Qs, double Us, double Vs,
                double freqcgs, double Tau, double tauF)
{
  // deposit the intensity and Stokes parameter in pixel
  image[i*ny+j] = Intensity * pow(freqcgs, 3);
  taus[i*ny+j] = Tau;

  if (!only_intensity) {
    imageS[(i*ny+j)*NIMG+0] = Is * pow(freqcgs, 3);
    if (params.qu_conv == 0) {
      imageS[(i*ny+j)*NIMG+1] = -Qs * pow(freqcgs, 3);
      imageS[(i*ny+j)*NIMG+2] = -Us * pow(freqcgs, 3);
    } else {
      imageS[(i*ny+j)*NIMG+1] = Qs * pow(freqcgs, 3);
      imageS[(i*ny+j)*NIMG+2] = Us * pow(freqcgs, 3);
    }
    imageS[(i*ny+j)*NIMG+3] = Vs * pow(freqcgs, 3);
    imageS[(i*ny+j)*NIMG+4] = tauF;

    if (isnan(imageS[(i*ny+j)*NIMG+0])) {
      fprintf(stderr, "NaN in image! Exiting.\n");
      exit(-1);
    }
  }
}

void save_pixelTransfer(double *image, double *imageS, double *taus, size_t iold, size_t jold,
                        size_t inew, size_t jnew, size_t nx, size_t ny, int only_intensity)
{
  // deposit the intensity and Stokes parameter in pixel
    double Intensity=image[iold*ny+jold];
    image[inew*ny+jnew]=Intensity;

    if(!only_intensity) {
      double Tau=taus[iold*ny+jold];
      double Is=imageS[(iold*ny+jold)*NIMG+0];
      double Qs=imageS[(iold*ny+jold)*NIMG+1];
      double Us=imageS[(iold*ny+jold)*NIMG+2];
      double Vs=imageS[(iold*ny+jold)*NIMG+3];
      double tauF=imageS[(iold*ny+jold)*NIMG+4];

      taus[inew*ny+jnew] = Tau;
      imageS[(inew*ny+jnew)*NIMG+0] = Is;
      imageS[(inew*ny+jnew)*NIMG+1] = Qs;
      imageS[(inew*ny+jnew)*NIMG+2] = Us;
      imageS[(inew*ny+jnew)*NIMG+3] = Vs;
      imageS[(inew*ny+jnew)*NIMG+4] = tauF;

      if (isnan(imageS[(iold*ny+jold)*NIMG+0])) {
        fprintf(stderr, "NaN in image! Exiting.\n");
        exit(-1);
      }
    }
}

void lininterp2(double *image, double *imageS, double *taus, size_t i1, size_t j1,
                size_t i2, size_t j2, size_t inew, size_t jnew, size_t nx, size_t ny, int only_intensity)
{
  // deposit the intensity and Stokes parameter in pixel
  double Intensity = .5 * (image[i1*ny+j1] + image[i2*ny+j2]);
  image[inew*ny+jnew] = Intensity;
  double Tau=.5*(taus[i1*ny+j1]+taus[i2*ny+j2]);
  taus[inew*ny+jnew] = Tau;

  if(!only_intensity) {
    double Is=.5*(imageS[(i1*ny+j1)*NIMG+0]+imageS[(i2*ny+j2)*NIMG+0]);
    double Qs=.5*(imageS[(i1*ny+j1)*NIMG+1]+imageS[(i2*ny+j2)*NIMG+1]);
    double Us=.5*(imageS[(i1*ny+j1)*NIMG+2]+imageS[(i2*ny+j2)*NIMG+2]);
    double Vs=.5*(imageS[(i1*ny+j1)*NIMG+3]+imageS[(i2*ny+j2)*NIMG+3]);
    double tauF=.5*(imageS[(i1*ny+j1)*NIMG+4]+imageS[(i2*ny+j2)*NIMG+4]);

    imageS[(inew*ny+jnew)*NIMG+0] = Is;
    imageS[(inew*ny+jnew)*NIMG+1] = Qs;
    imageS[(inew*ny+jnew)*NIMG+2] = Us;
    imageS[(inew*ny+jnew)*NIMG+3] = Vs;
    imageS[(inew*ny+jnew)*NIMG+4] = tauF;

    if (isnan(imageS[(i1*ny+j1)*NIMG+0])||isnan(imageS[(i2*ny+j2)*NIMG+0])) {
      fprintf(stderr, "NaN in image! Exiting.\n");
      exit(-1);
    }
  }
}


void lininterp4(double *image, double *imageS, double *taus, size_t i1, size_t j1,
                size_t i2,size_t j2, size_t i3, size_t j3, size_t i4, size_t j4, size_t inew, size_t jnew,
                size_t nx, size_t ny, int only_intensity)
{
  // deposit the intensity and Stokes parameter in pixel
  double Intensity = .25 * (image[i1*ny+j1] + image[i2*ny+j2] + image[i3*ny+j3] + image[i4*ny+j4]);
  image[inew*ny+jnew] = Intensity;
  double Tau=.25*(taus[i1*ny+j1]+taus[i2*ny+j2]+taus[i3*ny+j3]+taus[i4*ny+j4]);
  taus[inew*ny+jnew] = Tau;

  if(!only_intensity) {
    double Is=.25*(imageS[(i1*ny+j1)*NIMG+0]+imageS[(i2*ny+j2)*NIMG+0]+imageS[(i3*ny+j3)*NIMG+0]+imageS[(i4*ny+j4)*NIMG+0]);
    double Qs=.25*(imageS[(i1*ny+j1)*NIMG+1]+imageS[(i2*ny+j2)*NIMG+1]+imageS[(i3*ny+j3)*NIMG+1]+imageS[(i4*ny+j4)*NIMG+1]);
    double Us=.25*(imageS[(i1*ny+j1)*NIMG+2]+imageS[(i2*ny+j2)*NIMG+2]+imageS[(i3*ny+j3)*NIMG+2]+imageS[(i4*ny+j4)*NIMG+2]);
    double Vs=.25*(imageS[(i1*ny+j1)*NIMG+3]+imageS[(i2*ny+j2)*NIMG+3]+imageS[(i3*ny+j3)*NIMG+3]+imageS[(i4*ny+j4)*NIMG+3]);
    double tauF=.25*(imageS[(i1*ny+j1)*NIMG+4]+imageS[(i2*ny+j2)*NIMG+4]+imageS[(i3*ny+j3)*NIMG+4]+imageS[(i4*ny+j4)*NIMG+4]);

    imageS[(inew*ny+jnew)*NIMG+0] = Is;
    imageS[(inew*ny+jnew)*NIMG+1] = Qs;
    imageS[(inew*ny+jnew)*NIMG+2] = Us;
    imageS[(inew*ny+jnew)*NIMG+3] = Vs;
    imageS[(inew*ny+jnew)*NIMG+4] = tauF;

    if (isnan(imageS[(i1*ny+j1)*NIMG+0])||isnan(imageS[(i2*ny+j2)*NIMG+0])||isnan(imageS[(i3*ny+j3)*NIMG+0])||isnan(imageS[(i4*ny+j4)*NIMG+0])) {
      fprintf(stderr, "NaN in image! Exiting.\n");
      exit(-1);
    }
  }
}

void print_image_stats(double *image, double *imageS, size_t nx, size_t ny, Params params, double scale)
{
  double Ftot = 0.;
  double Ftot_unpol = 0.;
  double Imax = 0.0;
  double Iavg = 0.0;
  double Qtot = 0.;
  double Utot = 0.;
  double Vtot = 0.;
  size_t imax = 0;
  size_t jmax = 0;
  for (size_t i = 0; i < params.nx; i++) {
    for (size_t j = 0; j < params.ny; j++) {
      Ftot_unpol += image[i*ny+j]*scale;

      if (!params.only_unpolarized) {
        Ftot += imageS[(i*ny+j)*NIMG+0] * scale;
        Iavg += imageS[(i*ny+j)*NIMG+0];
        Qtot += imageS[(i*ny+j)*NIMG+1] * scale;
        Utot += imageS[(i*ny+j)*NIMG+2] * scale;
        Vtot += imageS[(i*ny+j)*NIMG+3] * scale;
        if (imageS[(i*ny+j)*NIMG+0] > Imax) {
          imax = i;
          jmax = j;
          Imax = imageS[(i*ny+j)*NIMG+0];
        }
      }
    }
  }

  // output normal flux quantities
  fprintf(stderr, "\nscale = %e\n", scale);
  fprintf(stderr, "imax=%ld jmax=%ld Imax=%g Iavg=%g\n", imax, jmax, Imax, Iavg/(params.nx*params.ny));
  fprintf(stderr, "freq: %g Ftot: %g Jy (%g Jy unpol xfer) scale=%g\n", params.freqcgs, Ftot, Ftot_unpol, scale);
  fprintf(stderr, "nuLnu = %g erg/s\n", 4.*M_PI*Ftot * params.dsource * params.dsource * JY * params.freqcgs);

  // output polarized transport information
  double LPfrac = 100.*sqrt(Qtot*Qtot+Utot*Utot)/Ftot;
  double CPfrac = 100.*Vtot/Ftot;
  fprintf(stderr, "I,Q,U,V [Jy]: %g %g %g %g\n", Ftot, Qtot, Utot, Vtot);
  fprintf(stderr, "LP,CP [%%]: %g %g\n", LPfrac, CPfrac);
}
