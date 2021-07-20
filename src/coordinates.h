#ifndef COORDINATES_H
#define COORDINATES_H

#include "decs.h"

// Poor man's enum of coordinate systems
// Modified Kerr-Schild coordinates.  Gammie '03.
#define METRIC_MKS 0
// New-style BHAC MKS.  Just like MKS but with X2 in [0,pi] and hslope->(1-hslope)
#define METRIC_BHACMKS 1
// "Funky" MKS coordinates as defined on IL wiki, see
// https://github.com/AFD-Illinois/docs/wiki/Coordinates
#define METRIC_FMKS 2
// MKS3 coordinates from koral-light
#define METRIC_MKS3 3
// Spherical coordinates in Minkowski space
#define METRIC_MINKOWSKI 4
// Exponential spherical coordinates in Minkowski space
#define METRIC_EMINKOWSKI 5
// eKS exponential radial coordinate; KS otherwise. note not the same 
// as eKS_internal, which has X2 in [0, 1]
#define METRIC_EKS 6

// Coordinate parameters.  See 
extern int use_eKS_internal;
extern int metric;
extern double a, hslope; // mks
extern double poly_norm, poly_xt, poly_alpha, mks_smooth; // fmks
extern double mks3R0, mks3H0, mks3MY1, mks3MY2, mks3MP0; // mks3
extern double startx[NDIM], stopx[NDIM], dx[NDIM];
extern double cstartx[NDIM], cstopx[NDIM];
extern double R0, Rin, Rout, Rh;

void bl_coord(double *X, double *r, double *th);
void bl_to_ks(double X[NDIM], double ucon_bl[NDIM], double ucon_ks[NDIM]);
void ks_to_bl(double X[NDIM], double ucon_ks[NDIM], double ucon_bl[NDIM]);
void gcov_func(double *X, double gcov[][NDIM]);
// TODO privatize these, why are they needed in models?
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
