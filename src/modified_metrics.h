#ifndef MODIFIED_METRICS_H
#define MODIFIED_METRICS_H

#include "decs.h"
extern double a;
extern double Rh;
extern double L_unit;

void gcov_EdGB_ks(double r, double th, double gcov[NDIM][NDIM]);
void gcov_EdGB_bl(double r, double th, double gcov[NDIM][NDIM]);
void EdGB_bl_to_ks(double r, double th, double edgb_T[NDIM][NDIM]);
void gcov_DCS_ks(double r, double th, double gcov[NDIM][NDIM]);
void gcov_dCS_bl(double r, double th, double gcov[NDIM][NDIM]);
void dCS_bl_to_ks(double r, double th, double dcs_T[NDIM][NDIM]);
double get_EdGB_Event_Horizon();
double get_dCS_Event_Horizon();
double event_horizon();
#endif