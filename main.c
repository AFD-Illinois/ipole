#include "decs.h"
#include "defs.h"
#include <omp.h>
#include "hdf5_utils.h"

#define MAXNSTEP 50000

// QU convention sets whether QU/EVPA is measured as
//   - 0  "East of North"  or  "observer convention" 
//   - 1  "North of West"  
#define QU_CONVENTION 0 

#define imgindex(n,i,j) (((n)*NX+(i))*NY+(j))

struct of_traj {
  double dl;
  double X[NDIM];
  double Kcon[NDIM];
  double Xhalf[NDIM];
  double Kconhalf[NDIM];
} traj[MAXNSTEP];
#pragma omp threadprivate(traj)

// used for slow light
struct of_image {
  int nstep;
  double intensity;
  double tau;
  double tauF;
  double complex N_coord[NDIM][NDIM];
};

// global variables
double thetacam, freqcgs;
char fnam[STRLEN];
int quench_output = 0;
int only_unpolarized = 0;

// can be moved to decs if so desired
void write_restart(const char *fname, double tA, double tB, double last_img_target,
    int nopenimgs, int nimg, int nconcurrentimgs, int s2,
    double *target_times, int *valid_imgs, struct of_image *dimages);
void read_restart(const char *fname, double *tA, double *tB, double *last_img_target,
    int *nopenimgs, int *nimg, int nconcurrentimgs, int s2,
    double *target_times, int *valid_imgs, struct of_image *dimages);

double tf = 0.;

Params params = { 0 };

int main(int argc, char *argv[]) 
{
  // motd
  fprintf(stderr, "%s. githash: %s\n", VERSION_STR, xstr(VERSION));
  fprintf(stderr, "notes: %s\n\n", xstr(NOTES));

  // initialization
  double time = omp_get_wtime();
  set_levi_civita();

  double tA, tB; // for slow light
  double phicam, rotcam, rcam, Xcam[NDIM];
  double xoff, yoff;
  double freq, Dsource, scale;
  double DX, DY, fovx, fovy;

  double *taus = malloc(sizeof(*taus) * NX*NY);
  double *imageS = malloc(sizeof(*imageS) * NX*NY*NIMG);
  double *image = malloc(sizeof(*image) * NX*NY);

  // "forward declarations"
  double root_find(double X[NDIM]);

#pragma omp parallel
  if (omp_get_thread_num() == 0) {
    fprintf(stderr, "nthreads = %d\n", omp_get_num_threads());
  }

  // load values from parameter file. handle all actual
  // model parameter comprehension in the model/* files
  load_par_from_argv(argc, argv, &params);

  // figure out if we should run in a custom mode
  for (int i=0; i<argc; ++i) {
    if ( strcmp(argv[i], "-quench") == 0 ) quench_output = 1;
    else if ( strcmp(argv[i], "-unpol") == 0 ) only_unpolarized = 1;
  }

  // now that we've loaded all parameters, tell our model about
  // them and use init_model to load the first dump
  parse_input(argc, argv, &params);
  init_model(&tA, &tB);
  R0 = 0.0;

  // normalize frequency to electron rest-mass energy
  freq = freqcgs * HPL / (ME * CL * CL);

  // set camera location
  //   rcam       [ GM/c^2 ]
  //   thetacam   [ degrees ]
  //   phicam     [ degrees ]
  rcam = 1000.;
  phicam = 0.0;
  rotcam = 0.0;
  xoff = 0.0;
  yoff = 0.0;
  if (params.loaded) {
    phicam = params.phicam;
    rotcam = params.rotcam*M_PI/180.;
    xoff = params.xoff;
    yoff = params.yoff;
  }
  // translate to geodesic coordinates
  double x[NDIM] = {0., rcam, thetacam/180.*M_PI, phicam/180.*M_PI};
  Xcam[0] = 0.0;
  Xcam[1] = log(rcam);
  Xcam[2] = root_find(x);
  Xcam[3] = phicam/180.*M_PI;
  fprintf(stderr, "Xcam[] = %e %e %e %e\n", Xcam[0], Xcam[1], Xcam[2], Xcam[3]);
  fprintf(stderr,
      "cam_th_cal=%g [deg] th_beg=%g th_len=%g a=%g R0=%g hslope=%g\n",
      (th_beg + th_len * Xcam[2]) * 180. / M_PI, th_beg, th_len, a,
      R0, hslope);

  // choose source distance here
  Dsource = DM87;
  //Dsource = DM87_gas;
  //Dsource = DSGRA;

  // set DX/DY using fov_muas to ensure full fov is given uas x uas
  double fov_muas = 160.;
  DX = fov_muas * Dsource / L_unit / 2.06265e11;

  // alternatively, set DX, DY here to ensure full fov is M x M
  // DX = 40.; // fov = 40 GM/c^2 in plane of the hole

  // make image square and set other fov values
  DY = DX;
  fovx = DX / rcam;
  fovy = DY / rcam;

  // Maximum radius of radiation interactions (GM/c^2)
  //rmax = 50.;

  scale = (DX * L_unit / NX) * (DY * L_unit / NY) / (Dsource * Dsource) / JY;
  fprintf(stderr,"L_unit = %e DX = %e NX = %i Dsource = %e JY = %e\n", L_unit, DX, NX, Dsource,JY);
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

    int nsteps[NX][NY];

    // first pass through geodesics to find lengths
#pragma omp parallel for schedule(dynamic,2) collapse(2) reduction(max:tgeoi) reduction(max:maxsteplength) reduction(min:t0) reduction(min:tgeof) shared(nsteps)
    for (int i=0; i<NX; ++i) {
      for (int j=0; j<NY; ++j) {
        if (j==0) fprintf(stderr, "%d ", i);

        int nstep = 0;
        double dl;
        double X[NDIM], Xhalf[NDIM], Kcon[NDIM], Kconhalf[NDIM];
        init_XK(i,j, Xcam, fovx,fovy, X, Kcon, rotcam, xoff, yoff);

        for (int k=0; k<NDIM; ++k) Kcon[k] *= freq;

        double tgeoitmp = 1.;
        double tgeoftmp = 1.;

        while (!stop_backward_integration(X, Kcon, Xcam)) {
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
    fprintf(stderr, "geodesic size = %g GB\n", 1. * sizeof(struct of_traj) * NX*NY * maxsteplength / 1024/1024/1024);
    struct of_traj *dtraj = malloc(sizeof(*dtraj) * NX*NY * maxsteplength);

    // TODO: if DTd is different, we'll have to calcluate this a different way
    int nconcurrentimgs = 2. + 1. * fabs(t0) / params.img_cadence;
    fprintf(stderr, "images size = %g GB\n", 1. * sizeof(struct of_image) * NX*NY * nconcurrentimgs / 1024/1024/1024);
    struct of_image *dimage = malloc(sizeof(*dimage) * NX*NY * nconcurrentimgs);

    // now populate the geodesic data structures
#pragma omp parallel for schedule(dynamic,2) collapse(2)
    for (int i=0; i<NX; ++i) {
      for (int j=0; j<NY; ++j) {
        if (j==0) fprintf(stderr, "%d ", i);

        int nstep = 0;
        double dl;
        double X[NDIM], Xhalf[NDIM], Kcon[NDIM], Kconhalf[NDIM];
        init_XK(i,j, Xcam, fovx,fovy, X, Kcon, rotcam, xoff, yoff);
        for (int k=0; k<NDIM; ++k) Kcon[k] *= freq;

        while (!stop_backward_integration(X, Kcon, Xcam)) {
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
            nconcurrentimgs, nconcurrentimgs*NX*NY,
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
          for (int i=0; i<NX; ++i) {
            for (int j=0; j<NY; ++j) {
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
        for (int i=0; i<NX; ++i) {
          for (int j=0; j<NY; ++j) {

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
          for (int i=0; i<NX; ++i) {
            for (int j=0; j<NY; ++j) {

              int pxidx = imgindex(k,i,j);
              int pstepidx = imgindex(dimage[pxidx].nstep,i,j);

              image[i*NY+j] = dimage[pxidx].intensity * pow(freqcgs, 3.);
              taus[i*NY+j] = dimage[pxidx].tau;

              if (! only_unpolarized) {
                double Stokes_I, Stokes_Q, Stokes_U, Stokes_V;
                project_N(dtraj[pstepidx].X, dtraj[pstepidx].Kcon, dimage[pxidx].N_coord, &Stokes_I, &Stokes_Q, &Stokes_U, &Stokes_V, rotcam);
                imageS[(i*NY+j)*NIMG+0] = Stokes_I * pow(freqcgs, 3.);
                imageS[(i*NY+j)*NIMG+1] = Stokes_Q * pow(freqcgs, 3.);
                imageS[(i*NY+j)*NIMG+2] = Stokes_U * pow(freqcgs, 3.);
                imageS[(i*NY+j)*NIMG+3] = Stokes_V * pow(freqcgs, 3.);
                imageS[(i*NY+j)*NIMG+4] = dimage[pxidx].tauF * pow(freqcgs, 3.);
                if (QU_CONVENTION == 0) {
                  imageS[(i*NY+j)*NIMG+1] *= -1;
                  imageS[(i*NY+j)*NIMG+2] *= -1;
                }
              } else {
                imageS[(i*NY+j)*NIMG+0] = 0.;
                imageS[(i*NY+j)*NIMG+1] = 0.;
                imageS[(i*NY+j)*NIMG+2] = 0.;
                imageS[(i*NY+j)*NIMG+3] = 0.;
                imageS[(i*NY+j)*NIMG+4] = 0.;
              }
            }
          }

          //fprintf(stderr, "saving image %d at t = %g\n", k, target_times[k]);
          char dfname[256];
          snprintf(dfname, 255, params.outf, target_times[k]);
          dump(image, imageS, taus, dfname, scale, Dsource, Xcam, DX, DY,
               fovx, fovy, rcam, thetacam, phicam, rotcam, xoff, yoff);
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
              nconcurrentimgs, nconcurrentimgs*NX*NY,
              target_times, valid_images, dimage);
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
        double Intensity, Stokes_I, Stokes_Q, Stokes_U, Stokes_V, tauF, Tau;

        init_XK(i,j, Xcam, fovx,fovy, X, Kcon, rotcam, xoff, yoff);
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
            break;
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
        Tau = 0.;
        tauF = 0.;

        get_jkinv(traj[nstep].X, traj[nstep].Kcon, &ji, &ki);
        
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
          get_jkinv(Xf, Kconf, &jf, &kf);
          Intensity = approximate_solve(Intensity, ji, ki, jf, kf, traj[nstep].dl, &Tau);

          // solve polarized transport
          if (! only_unpolarized ) {
            evolve_N(Xi, Kconi, Xhalf, Kconhalf, Xf, Kconf, traj[nstep].dl, N_coord, &tauF);
            if (isnan(creal(N_coord[0][0]))) {
              exit(-3);
            }
          }

          // swap start and finish
          ji = jf;
          ki = kf;

          nstep--;
        }

        // deposit intensity and Stokes parameter in pixels
        image[i*NY+j] = Intensity * pow(freqcgs, 3);
        taus[i*NY+j] = Tau;
        if (! only_unpolarized) {
          project_N(Xf, Kconf, N_coord, &Stokes_I, &Stokes_Q, &Stokes_U, &Stokes_V, rotcam);
          imageS[(i*NY+j)*NIMG+0] = Stokes_I * pow(freqcgs, 3);
          imageS[(i*NY+j)*NIMG+1] = Stokes_Q * pow(freqcgs, 3);
          imageS[(i*NY+j)*NIMG+2] = Stokes_U * pow(freqcgs, 3);
          imageS[(i*NY+j)*NIMG+3] = Stokes_V * pow(freqcgs, 3);
          imageS[(i*NY+j)*NIMG+4] = tauF;
          if (QU_CONVENTION == 0) {
            imageS[(i*NY+j)*NIMG+1] *= -1;
            imageS[(i*NY+j)*NIMG+2] *= -1;
          }
          if (isnan(imageS[(i*NY+j)*NIMG+0])) exit(-1);
        } else {
          imageS[(i*NY+j)*NIMG+0] = 0.;
          imageS[(i*NY+j)*NIMG+1] = 0.;
          imageS[(i*NY+j)*NIMG+2] = 0.;
          imageS[(i*NY+j)*NIMG+3] = 0.;
          imageS[(i*NY+j)*NIMG+4] = 0.;
        }

        // print progress
        if ( (nprogress % NY) == 0 ) {
          fprintf(stderr, "%d ", nprogress/NY);
        }

#pragma omp atomic
        nprogress += 1;

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
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NY; j++) {
        Ftot_unpol += image[i*NY+j]*scale;
        Ftot += imageS[(i*NY+j)*NIMG+0] * scale;
        Iavg += imageS[(i*NY+j)*NIMG+0];
        if (imageS[(i*NY+j)*NIMG+0] > Imax) {
          imax = i;
          jmax = j;
          Imax = imageS[(i*NY+j)*NIMG+0];
        }
      }
    }

    fprintf(stderr, "imax=%d jmax=%d Imax=%g Iavg=%g\n", imax, jmax, Imax, Iavg/(NX*NY));
    fprintf(stderr, "freq: %g Ftot: %g (%g unpol) scale=%g\n", freqcgs, Ftot, Ftot_unpol, scale);
    fprintf(stderr, "nuLnu = %g\n", 4.*M_PI*Ftot * Dsource * Dsource * JY * freqcgs);

    // don't dump if we've been asked to quench output. useful for batch jobs
    // like when fitting light curve fluxes
    if (!quench_output) {
      // dump result. if parameters have been loaded, don't also
      // output image
      if (params.loaded) {
        dump(image, imageS, taus, params.outf, scale, Dsource, Xcam, DX, DY, 
             fovx, fovy, rcam, thetacam, phicam, rotcam, xoff, yoff);
      } else {
        dump(image, imageS, taus, "ipole.dat", scale, Dsource, Xcam, DX, DY, 
             fovx, fovy, rcam, thetacam, phicam, rotcam, xoff, yoff);
        IMLOOP image[i*NY+j] = log(image[i*NY+j] + 1.e-50);
        make_ppm(image, freq, "ipole_lfnu.ppm");
      }
    }

  }

  time = omp_get_wtime() - time;
  printf("Total wallclock time: %g s\n\n", time);

  return 0;
}

void read_restart(const char *fname, double *tA, double *tB, double *last_img_target,
        int *nopenimgs, int *nimg, int nconcurrentimgs, int s2,
        double *target_times, int *valid_imgs, struct of_image *dimages) {

  hdf5_open(fname);

  hdf5_set_directory("/");
  hdf5_read_single_val(tA, "/tA", H5T_IEEE_F64LE);
  hdf5_read_single_val(tB, "/tB", H5T_IEEE_F64LE);
  hdf5_read_single_val(last_img_target, "/last_img_target", H5T_IEEE_F64LE);
  hdf5_read_single_val(nimg, "/nimg", H5T_STD_I32LE);
  hdf5_read_single_val(nopenimgs, "/nopenimgs", H5T_STD_I32LE);

  int *tint = malloc(s2 * sizeof(*tint));
  double *tdbl = malloc(s2 * sizeof(*tdbl));

  for (int i=0; i<nconcurrentimgs; ++i) {
    tdbl[i] = target_times[i];
    tint[i] = valid_imgs[i];
  }

  hsize_t dims[] = { nconcurrentimgs };
  hsize_t start[] = { 0 };

  hdf5_read_array(target_times, "/target_times", 1, dims, start, dims, dims, start, H5T_IEEE_F64LE);
  hdf5_read_array(valid_imgs, "/valid_images", 1, dims, start, dims, dims, start, H5T_STD_I32LE);

  dims[0] = s2;  
  
  hdf5_read_array(tint, "/dimg/nstep", 1, dims, start, dims, dims, start, H5T_STD_I32LE);
  for (int i=0; i<s2; ++i) dimages[i].nstep = tint[i];

  hdf5_read_array(tdbl, "/dimg/intensity", 1, dims, start, dims, dims, start, H5T_IEEE_F64LE);
  for (int i=0; i<s2; ++i) dimages[i].intensity = tdbl[i];

  hdf5_read_array(tdbl, "/dimg/tau", 1, dims, start, dims, dims, start, H5T_IEEE_F64LE);
  for (int i=0; i<s2; ++i) dimages[i].tau = tdbl[i];

  hdf5_read_array(tdbl, "/dimg/tauF", 1, dims, start, dims, dims, start, H5T_IEEE_F64LE);
  for (int i=0; i<s2; ++i) dimages[i].tauF = tdbl[i];

  for (int mu=0; mu<4; ++mu) {
    for (int nu=0; nu<4; ++nu) {
      char tgt[20];
      snprintf(tgt, 18, "/dimg/Nr%d%d", mu, nu);
      hdf5_read_array(tdbl, tgt, 1, dims, start, dims, dims, start, H5T_IEEE_F64LE);
      for (int i=0; i<s2; ++i) dimages[i].N_coord[mu][nu] = tdbl[i];
      snprintf(tgt, 18, "/dimg/Ni%d%d", mu, nu);
      hdf5_read_array(tdbl, tgt, 1, dims, start, dims, dims, start, H5T_IEEE_F64LE);
      for (int i=0; i<s2; ++i) dimages[i].N_coord[mu][nu] += tdbl[i] * _Complex_I;
    }
  }

  free(tint);
  free(tdbl);

  hdf5_close();

}

void write_restart(const char *fname, double tA, double tB, double last_img_target,
        int nopenimgs, int nimg, int nconcurrentimgs, int s2, 
        double *target_times, int *valid_imgs, struct of_image *dimages) {

  hid_t fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0) {
    fprintf(stderr, "! unable to create hdf5 restart file.\n");
    exit(-4);
  }

  h5io_add_attribute_str(fid, "/", "githash", xstr(VERSION));
  h5io_add_data_dbl(fid, "/tA", tA);
  h5io_add_data_dbl(fid, "/tB", tB);
  h5io_add_data_dbl(fid, "/last_img_target", last_img_target);
  h5io_add_data_int(fid, "/nimg", nimg);
  h5io_add_data_int(fid, "/nopenimgs", nopenimgs);
  h5io_add_data_dbl_1d(fid, "/target_times", nconcurrentimgs, target_times);
  h5io_add_data_int_1d(fid, "/valid_images", nconcurrentimgs, valid_imgs);
  h5io_add_group(fid, "/dimg");

  // save dimg struct
  int *tint = malloc(s2 * sizeof(*tint));
  double *tdbl = malloc(s2 * sizeof(*tdbl));

  for (int i=0; i<s2; ++i) tint[i] = dimages[i].nstep;
  h5io_add_data_int_1d(fid, "/dimg/nstep", s2, tint);
  
  for (int i=0; i<s2; ++i) tdbl[i] = dimages[i].intensity;
  h5io_add_data_dbl_1d(fid, "/dimg/intensity", s2, tdbl);
    
  for (int i=0; i<s2; ++i) tdbl[i] = dimages[i].tau;
  h5io_add_data_dbl_1d(fid, "/dimg/tau", s2, tdbl);

  for (int i=0; i<s2; ++i) tdbl[i] = dimages[i].tauF;
  h5io_add_data_dbl_1d(fid, "/dimg/tauF", s2, tdbl);

  for (int mu=0; mu<4; ++mu) {
    for (int nu=0; nu<4; ++nu) {
      char tgt[20];
      snprintf(tgt, 18, "/dimg/Nr%d%d", mu, nu);
      for (int i=0; i<s2; ++i) tdbl[i] = creal(dimages[i].N_coord[mu][nu]);
      h5io_add_data_dbl_1d(fid, tgt, s2, tdbl); 
      snprintf(tgt, 18, "/dimg/Ni%d%d", mu, nu);
      for (int i=0; i<s2; ++i) tdbl[i] = cimag(dimages[i].N_coord[mu][nu]);
      h5io_add_data_dbl_1d(fid, tgt, s2, tdbl); 
    }
  }

  free(tint);
  free(tdbl);
  
  H5Fclose(fid);
}


void dump(double image[], double imageS[], double taus[],
    const char *fname, double scale, double Dsource, double cam[NDIM], double DX, 
    double DY, double fovx, double fovy, double rcam, double thetacam, double phicam,
    double rotcam, double xoff, double yoff)
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

  if (QU_CONVENTION == 0) h5io_add_data_str(fid, "/header/evpa_0", "N");
  else h5io_add_data_str(fid, "/header/evpa_0", "W");

  h5io_add_data_str(fid, "/header/version", VERSION_STR);

  h5io_add_group(fid, "/header/camera");
  h5io_add_data_int(fid, "/header/camera/nx", NX);
  h5io_add_data_int(fid, "/header/camera/ny", NY);
  h5io_add_data_dbl(fid, "/header/camera/dx", DX);
  h5io_add_data_dbl(fid, "/header/camera/dy", DY);
  h5io_add_data_dbl(fid, "/header/camera/fovx", fovx);
  h5io_add_data_dbl(fid, "/header/camera/fovy", fovy);
  h5io_add_data_dbl(fid, "/header/camera/rcam", rcam);
  h5io_add_data_dbl(fid, "/header/camera/thetacam", thetacam);
  h5io_add_data_dbl(fid, "/header/camera/phicam", phicam);
  h5io_add_data_dbl(fid, "/header/camera/rotcam", rotcam);
  h5io_add_data_dbl(fid, "/header/camera/xoff", xoff);
  h5io_add_data_dbl(fid, "/header/camera/yoff", yoff);
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
      Ftot_unpol += image[i*NY+j] * scale;
      Ftot += imageS[(i*NY+j)*NIMG+0] * scale;
    }
  }

  // output stuff
  h5io_add_data_dbl(fid, "/Ftot", Ftot);
  h5io_add_data_dbl(fid, "/Ftot_unpol", Ftot_unpol);
  h5io_add_data_dbl(fid, "/nuLnu", 4. * M_PI * Ftot * Dsource * Dsource * JY * freqcgs);
  h5io_add_data_dbl(fid, "/nuLnu_unpol", 4. * M_PI * Ftot_unpol * Dsource * Dsource * JY * freqcgs);

  h5io_add_data_dbl_2ds(fid, "/unpol", NX, NY, image);
  h5io_add_data_dbl_2ds(fid, "/tau", NX, NY, taus);
  h5io_add_data_dbl_3ds(fid, "/pol", NX, NY, 5, imageS);

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
    double X[NDIM], double Kcon[NDIM], double rotcam, double xoff, double yoff)
{
  double Econ[NDIM][NDIM];
  double Ecov[NDIM][NDIM];
  double Kcon_tetrad[NDIM];

  make_camera_tetrad(Xcam, Econ, Ecov);

  /* construct *outgoing* wavevectors */
  double dxoff = (xoff+i-0.01)/(double)NX - 0.5;
  double dyoff = (yoff+j)/(double)NY - 0.5;
  Kcon_tetrad[0] = 0.;
  Kcon_tetrad[1] = (dxoff*cos(rotcam) - dyoff*sin(rotcam)) * fovx;
  Kcon_tetrad[2] = (dxoff*sin(rotcam) + dyoff*cos(rotcam)) * fovy;
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

