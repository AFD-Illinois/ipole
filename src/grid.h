#ifndef GRID_H
#define GRID_H

#include "decs.h"

extern int N1, N2, N3;

double gdet_zone(int i, int j, int k);

void ijktoX(int i, int j, int k, double X[NDIM]);
void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]);

int X_in_domain(double X[NDIM]);

double interp_scalar(double X[NDIM], double ***var);
double interp_scalar_time(double X[NDIM], double ***varA, double ***varB, double tfac);

#endif // GRID_H
