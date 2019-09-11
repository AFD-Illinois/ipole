#ifndef PAR_H
#define PAR_H

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#define TYPE_INT (1)
#define TYPE_DBL (2)
#define TYPE_STR (3)

// feel free to change any part of this structure
typedef struct params_t {
  double tp_over_te;    // tp_over_te (defaults to 3.)
  double trat_small;    // when beta << 1 (defaults to 1.)
  double trat_large;    // when beta >> 1 (defaults to 10.)
  double thetacam;      // in degrees from the pole
  double phicam;        // in degrees
  double rotcam;        // in degrees
  double freqcgs;       // ... in cgs
  double MBH;           // in Msun
  double M_unit;        // in cgs
  int counterjet;       // 0,1,2
  int add_ppm;         // Whether to additionally make a ppm image of I
  const char dump[256];
  const char outf[256];
  char loaded;

  // ML parameters
  double xoff, yoff;    // in pixels

  // slow light
  int dump_min;
  int dump_max;
  int dump_skip;
  double img_cadence;
  double restart_int;
} Params;

// modify this to set default values
void load_par_from_argv(int, char *[], Params *);

// modify this to actually add the reading subroutines
void try_set_parameter(const char *, const char *, Params *);

// only modify if you add/modify types
void load_par(const char *, Params *);
void set_by_word_val (const char *word, const char *value, const char *key, void *val, int type);

#endif // PAR_H
