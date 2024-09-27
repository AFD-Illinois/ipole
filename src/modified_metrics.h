#ifndef MODIFIED_METRICS_H
#define MODIFIED_METRICS_H

#include "decs.h"
extern double a;
extern double Rh;
extern double L_unit;

void gcov_EdGB_ks(double r, double th, double gcov[NDIM][NDIM]);
void gcov_DCS_ks(double r, double th, double gcov[NDIM][NDIM]);
void matrix_multiply(double A[NDIM][NDIM], double B[NDIM][NDIM], double result[NDIM][NDIM]);
void matrix_exponential(double A[NDIM][NDIM], double expA[NDIM][NDIM]);
double get_EdGB_Event_Horizon();
double get_dCS_Event_Horizon();
double event_horizon();
#endif