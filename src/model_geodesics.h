/*
 * model_geodesics.h
 *
 *  Created on: Sep 11, 2019
 *      Author: bprather
 */

#ifndef SRC_MODEL_GEODESICS_H_
#define SRC_MODEL_GEODESICS_H_

#include "decs.h"

int stop_backward_integration(double X[NDIM], double Xhalf[NDIM], double Kcon[NDIM], double Xcam[NDIM]);
double stepsize(double X[NDIM], double K[NDIM]);

#endif /* SRC_MODEL_GEODESICS_H_ */
