#ifndef MODEL_PARAMS_H
#define MODEL_PARAMS_H

#define SLOW_LIGHT (0)
#define EMHD_CONDUCTION (0)
#define EMHD_VISCOSITY  (0)

// Model-specific definitions and globals
#define KRHO 0
#define UU   1
#define U1   2
#define U2   3
#define U3   4
#if (EMHD_CONDUCTION)
#if (EMHD_VISCOSITY)
#define QTILDE 5
#define DPTILDE 6
#define B1   7
#define B2   8
#define B3   9
#define KEL  10
#define KTOT 11
#else
#define QTILDE 5
#define B1   6
#define B2   7
#define B3   8
#define KEL  9
#define KTOT 10
#endif
#elif (EMHD_VISCOSITY)
#define DPTILDE 5
#define B1   6
#define B2   7
#define B3   8
#define KEL  9
#define KTOT 10
#else
#define B1   5
#define B2   6
#define B3   7
#define KEL  8
#define KTOT 9
#endif

// These two will never be used simultaneously,
// and never with KEL.
#define TFLK 8  // temperature of fluid in Kelvin
#define THF  8  // fluid temperature in me c^2

extern double DTd;
extern double sigma_cut;

#endif // MODEL_PARAMS_H
