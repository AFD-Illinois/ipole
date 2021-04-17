/*
 * geometry.h
 *
 *  Created on: Sep 9, 2019
 *      Author: bprather
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "decs.h"

int gcon_func(double gcov[][NDIM], double gcon[][NDIM]);
double gdet_func(double gcov[][NDIM]);
void get_connection(double *X, double lconn[][NDIM][NDIM]);

void flip_index(double *ucon, double Gcov[NDIM][NDIM], double *ucov);
// Old names aliased
inline void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov) {flip_index(ucon, Gcov, ucov);};
inline void raise(double *ucov, double Gcon[NDIM][NDIM], double *ucon) {flip_index(ucov, Gcon, ucon);};

void null_normalize(double Kcon[NDIM], double fnorm);
void normalize(double *vcon, double gcov[][NDIM]);
void normalize_to(double vcon[NDIM], double Gcov[NDIM][NDIM], double target);

int invert_matrix(double Am[][NDIM], double Aminv[][NDIM]);
double theta_func(double X[NDIM]);
int levi_civita(int i, int j, int k, int l);

#endif /* GEOMETRY_H */
