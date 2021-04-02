/*
 * model_tetrads.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef MODEL_TETRADS_H
#define MODEL_TETRADS_H

#include "decs.h"

int make_camera_tetrad(double X[NDIM], double Econ[NDIM][NDIM],
                    double Ecov[NDIM][NDIM]);
int make_camera_tetrad_old(double X[NDIM], double Econ[NDIM][NDIM],
                        double Ecov[NDIM][NDIM]);
int make_plasma_tetrad(double Ucon[NDIM], double Kcon[NDIM], double Bcon[NDIM],
                    double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM],
                    double Ecov[NDIM][NDIM]);

#endif /* MODEL_TETRADS_H */
