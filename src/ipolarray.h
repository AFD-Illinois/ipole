/*
 * ipolarray.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef IPOLARRAY_H
#define IPOLARRAY_H

#include "decs.h"
#include "par.h"

#include <complex.h>

// Top-level functions for solving emission
int integrate_emission(struct of_traj *traj, int nsteps,
                    double *Intensity, double *Tau, double *tauF,
                    double complex N_coord[NDIM][NDIM], Params *params);

// Needed for slow light.  TODO extend above to use instead
int evolve_N(double Xi[NDIM],double Kconi[NDIM],
    double Xhalf[NDIM], double Kconhalf[NDIM],
    double Xf[NDIM],double Kconf[NDIM],
    double dlam,
    double complex N_coord[NDIM][NDIM],
    double *tauF, int ZERO_EMISSION, Params *params);
double approximate_solve (double Ii, double ji, double ki, double jf, double kf,
                   double dl, double *tau);

void project_N(double X[NDIM],double Kcon[NDIM],
    double complex Ncon[NDIM][NDIM],
    double *Stokes_I, double *Stokes_Q,double *Stokes_U,double *Stokes_V, double rotcam);

#endif /* IPOLARRAY_H */
