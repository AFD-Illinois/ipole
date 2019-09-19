#ifndef MODEL_H
#define MODEL_H

#include "model_params.h"
#include "decs.h"
#include "par.h"
#include "hdf5_utils.h"

extern double M_unit;
extern double L_unit;
extern double T_unit;
extern double RHO_unit;
extern double U_unit;
extern double B_unit;
extern double Te_unit;

void try_set_model_parameter(const char *word, const char *value);
void init_model(double *tA, double *tB);

int radiating_region(double X[NDIM]);

double get_model_thetae(double X[NDIM]);
double get_model_b(double X[NDIM]);
double get_model_ne(double X[NDIM]);
void get_model_fourv(double X[NDIM], double Ucon[NDIM], double Ucov[NDIM],
                     double Bcon[NDIM], double Bcov[NDIM]);
void get_model_bcov(double X[NDIM], double Bcov[NDIM]);
void get_model_bcon(double X[NDIM], double Bcon[NDIM]);
void get_model_ucov(double X[NDIM], double Ucov[NDIM]);
void get_model_ucon(double X[NDIM], double Ucon[NDIM]);

void update_data(double *tA, double *tB);
void update_data_until(double *tA, double *tB, double tgt);

void output_hdf5();

// Optional function to be able to trace primitives along a geodesic
void get_model_primitives(double X[NDIM], double *p);

#endif // MODEL_H
