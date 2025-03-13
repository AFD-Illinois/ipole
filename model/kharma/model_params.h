#ifndef MODEL_PARAMS_H
#define MODEL_PARAMS_H

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
// These two will never be used simultaneously,
// and never with KEL.
#define TFLK 8  // temperature of fluid in Kelvin
#define THF  8  // fluid temperature in me c^2

extern double DTd;
extern double sigma_cut;

// Model-specific function declarations

double get_dump_time(char *fnam, int dumpidx);
void set_units();
void load_kharma_data(int n, char *, int dumpidx, int verbose);
void init_physical_quantities(int n, double rescale_factor);
void init_storage(void);

#endif // MODEL_PARAMS_H
