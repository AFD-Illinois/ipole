#ifndef MODEL_H
#define MODEL_H

#include "decs.h"
#include "par.h"
#include "hdf5_utils.h"

#define NX (160)
#define NY (160)

#define SLOW_LIGHT (0)

// Model-specific definitions and globals
#define KRHO 0
#define UU   1
#define U1   2
#define U2   3
#define U3   4
#define B1   5
#define B2   6
#define B3   7
#define KEL  8
#define KTOT 9

extern double M_unit;
extern double L_unit;
extern double T_unit;
extern double RHO_unit;
extern double U_unit;
extern double B_unit;
extern double Te_unit;

extern double freqcgs;
extern double thetacam;
extern double gam;
extern double DTd;
extern double t0;
extern double trat_j;
extern double theta_j;
extern double trat_d;
extern int counterjet;
extern double rmax_geo;

void set_units();

void parse_input(int argc, char *argv[], Params *params);
void init_storage(void);

double stepsize(double X[NDIM], double K[NDIM]);
void init_model(double *tA, double *tB);
void init_physical_quantities(int n);
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

// Private functions (move into .c?)
void load_iharm_data(int n, char *, int dumpidx, int verbose);
double get_dump_t(char *fnam, int dumpidx);
void init_iharm_grid(char *fnam, int dumpidx);

#endif // MODEL_H
