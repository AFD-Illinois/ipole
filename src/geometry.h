/*
 * geometry.h
 *
 *  Created on: Sep 9, 2019
 *      Author: bprather
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "decs.h"

int gcon_func(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]);
double gdet_func(double gcov[NDIM][NDIM]);
void get_connection(double X[NDIM], double lconn[NDIM][NDIM][NDIM]);

void flip_index(double ucon[NDIM], double Gcov[NDIM][NDIM], double ucov[NDIM]);
// Old names aliased
inline void lower(double ucon[NDIM], double Gcov[NDIM][NDIM], double ucov[NDIM]) {flip_index(ucon, Gcov, ucov);};
inline void raise(double ucov[NDIM], double Gcon[NDIM][NDIM], double ucon[NDIM]) {flip_index(ucov, Gcon, ucon);};

void null_normalize(double Kcon[NDIM], double fnorm);
void normalize(double vcon[NDIM], double gcov[NDIM][NDIM]);
void normalize_to(double vcon[NDIM], double Gcov[NDIM][NDIM], double target);

int invert_matrix(double Am[NDIM][NDIM], double Aminv[NDIM][NDIM]);
double theta_func(double X[NDIM]);
int levi_civita(int i, int j, int k, int l);

#endif /* GEOMETRY_H */
