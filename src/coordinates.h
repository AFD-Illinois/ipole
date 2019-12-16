#ifndef COORDINATES_H
#define COORDINATES_H

#include "decs.h"

extern int METRIC_eKS;
extern int METRIC_MKS, METRIC_FMKS, METRIC_MKS3;
extern double a, hslope;
extern double poly_norm, poly_xt, poly_alpha, mks_smooth; // mmks
extern double mks3R0, mks3H0, mks3MY1, mks3MY2, mks3MP0; // mks3
extern double startx[NDIM], stopx[NDIM], dx[NDIM];
extern double R0, Rin, Rout, Rh;

void bl_coord(double *X, double *r, double *th);
void bl_to_ks(double X[NDIM], double ucon_bl[NDIM], double ucon_ks[NDIM]);
void ks_to_bl(double X[NDIM], double ucon_ks[NDIM], double ucon_bl[NDIM]);
void gcov_func(double *X, double gcov[][NDIM]);
void gcov_ks(double r, double th, double gcov[NDIM][NDIM]);
void gcov_bl(double r, double th, double gcov[NDIM][NDIM]);

// Internal
void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM]);
void set_dXdx(double X[NDIM], double dXdx[NDIM][NDIM]);
void vec_to_ks(double X[NDIM], double v_nat[NDIM], double v_ks[NDIM]);
void vec_from_ks(double X[NDIM], double v_ks[NDIM], double v_nat[NDIM]);

// Translation to native coords
void native_coord(double r, double th, double phi, double X[NDIM]);
double root_find(double X[NDIM]);

#endif // COORDINATES_H
