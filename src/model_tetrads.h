/*
 * model_tetrads.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef MODEL_TETRADS_H
#define MODEL_TETRADS_H

#include "decs.h"

// define CAMERA_CENTER_ZAMO to use force the middle camera pixel to have k_phi = 0
#define CAMERA_CENTER_ZAMO

void make_camera_tetrad(double X[NDIM], double Econ[NDIM][NDIM],
                    double Ecov[NDIM][NDIM]);
void make_plasma_tetrad(double Ucon[NDIM], double Kcon[NDIM], double Bcon[NDIM],
                    double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM],
                    double Ecov[NDIM][NDIM]);

#endif /* MODEL_TETRADS_H */
