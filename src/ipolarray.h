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
                    double complex N_coord[NDIM][NDIM], Params *params,
                    int print);

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
void parallel_transport_vector(double Xi[NDIM], double Xm[NDIM], double Xf[NDIM],
    double Ki[NDIM], double Km[NDIM], double Kf[NDIM],
    double Ni[NDIM],
    double Nm[NDIM],
    double Nf[NDIM], double dl);

/* tensor tools */
void complex_lower(double complex N[NDIM][NDIM], double gcov[NDIM][NDIM],
    int low1, int low2, double complex Nl[NDIM][NDIM]);
void stokes_to_tensor(double fI, double fQ, double fU, double fV,
    double complex f_tetrad[NDIM][NDIM]);
void tensor_to_stokes(double complex f_tetrad[NDIM][NDIM], double *fI,
    double *fQ, double *fU, double *fV);
void any_tensor_to_stokes(double complex f_any[NDIM][NDIM], double gcov[NDIM][NDIM],
    double *fI, double *fQ, double *fU, double *fV);
void complex_coord_to_tetrad_rank2(double complex T_coord[NDIM][NDIM],
    double Ecov[NDIM][NDIM],
    double complex T_tetrad[NDIM][NDIM]);
void complex_tetrad_to_coord_rank2(double complex T_tetrad[NDIM][NDIM],
    double Econ[NDIM][NDIM],
    double complex T_coord[NDIM][NDIM]);

// Pointer for keeping/writing a histogram of zone contents
extern double *visible_emission_histogram;

#endif /* IPOLARRAY_H */
