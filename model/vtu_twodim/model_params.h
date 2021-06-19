#ifndef MODEL_PARAMS_H
#define MODEL_PARAMS_H

#define SLOW_LIGHT (0)

// used in model_geodesic to know when to stop integrating
extern double rmax_geo;

// used in various places to determine how to scale from
// geodesics to cgs units
extern double L_unit;

#endif // MODEL_PARAMS_H
