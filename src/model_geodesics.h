/*
 * model_geodesics.h
 *
 *  Created on: Sep 11, 2019
 *      Author: bprather
 */

#ifndef SRC_MODEL_GEODESICS_H_
#define SRC_MODEL_GEODESICS_H_

#include "decs.h"
#include "par.h"

int trace_geodesic(double Xi[NDIM], double Kconi[NDIM], struct of_traj *traj, double eps, int step_max, double Xcam[NDIM], int print);
void init_XK(long int i, long int j, int nx, int ny, double Xcam[NDIM],
             Params params, double fovx, double fovy,
             double X[NDIM], double Kcon[NDIM]);

// Internal utilities still used for slow light
int stop_backward_integration(double X[NDIM], double Xhalf[NDIM], double Kcon[NDIM]);
double stepsize(double X[NDIM], double K[NDIM], double eps);

#endif /* SRC_MODEL_GEODESICS_H_ */
