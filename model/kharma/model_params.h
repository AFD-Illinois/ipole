#ifndef MODEL_PARAMS_H
#define MODEL_PARAMS_H

#include <stddef.h>

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
int read_info_attribute(const char *filename, const char *attr_name, int type, void *value, size_t value_size);
char *get_parfile(const char *fname);
int get_parameter_value(const char *input_text, const char *block, const char *key, int type, void *value, size_t value_size);

void init_storage_prims(void);
void init_storage_derived_quantities(void);

double get_dump_time(char *fnam, int dumpidx);
int read_parameters_and_alloate_memory(char *fnam, int dumpidx);
void load_kharma_data(int n, char *fnam, int dumpidx, int verbose);

void set_units();
void populate_boundary_conditions(int n);
double get_code_dMact(int i, int n);
void init_physical_quantities(int n, double rescale_factor);

#endif // MODEL_PARAMS_H
