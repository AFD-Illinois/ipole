/*
 * radiation.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef RADIATION_H
#define RADIATION_H

#include "decs.h"

/* radiation */
double Bnu_inv(double nu, double Thetae);
double jnu_inv(double nu, double Thetae, double Ne, double B, double theta);
double get_fluid_nu(double Kcon[NDIM], double ucov[NDIM]);
double get_bk_angle(double X[NDIM], double Kcon[NDIM], double ucov[NDIM],
              double bcon[NDIM], double bcov[NDIM]);

#endif /* RADIATION_H */
