/*
 * tetrads.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef TETRADS_H
#define TETRADS_H

#include "decs.h"

void coordinate_to_tetrad(double Ecov[NDIM][NDIM], double K[NDIM],
                      double K_tetrad[NDIM]);
void tetrad_to_coordinate(double Ecov[NDIM][NDIM], double K_tetrad[NDIM],
                      double K[NDIM]);
void set_Econ_from_trial(double Econ[4], int defdir, double trial[4]);
int check_handedness(double Econ[NDIM][NDIM], double Gcov[NDIM][NDIM], double *dot);
void project_out(double vcona[NDIM], double vconb[NDIM], double Gcov[4][4]);

#endif /* TETRADS_H */
