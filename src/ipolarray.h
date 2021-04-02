/*
 * ipolarray.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef IPOLARRAY_H
#define IPOLARRAY_H

#include "decs.h"

#include <complex.h>

void init_N(double Xi[NDIM],double Kconi[NDIM],double complex Ncon[NDIM][NDIM]);
void evolve_N(double Xi[NDIM],double Kconi[NDIM],
    double Xf[NDIM],double Kconf[NDIM],
    double Xhalf[NDIM],double Kconhalf[NDIM],
    double dlam,
    double complex N_coord[NDIM][NDIM],
    double *tauF);
void project_N(double X[NDIM],double Kcon[NDIM],
    double complex Ncon[NDIM][NDIM],
    double *Stokes_I, double *Stokes_Q,double *Stokes_U,double *Stokes_V, double rotcam);

#endif /* IPOLARRAY_H */
