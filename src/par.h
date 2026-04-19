#ifndef PAR_H
#define PAR_H

#include "decs.h"

#define TYPE_INT (1)
#define TYPE_DBL (2)
#define TYPE_STR (3)

#define TYPE_INT_ARRAY (4)
#define TYPE_DBL_ARRAY (5)

// feel free to change any part of this structure
typedef struct params_t {
  double rcam;          // Camera radius in r_g
  double thetacam;      // in degrees from the pole
  double phicam;        // in degrees
  double rotcam;        // in degrees
  double dx, dy;        // FOV in-plane in r_g
  double fovx_dsource, fovy_dsource; // FOV (from Earth) in muas
  int nx, ny;           // image dimensions in px
  double dsource;       // in pc
  double freqcgs;       // ... in cgs
  int old_centering;    // 0 uses k_phi=0 "ZAMO" new centering, 1 uses k^phi=0 old centering

  // Geodesic accuracy
  double eps;
  int maxnstep;

  int add_ppm;          // Whether to additionally make a ppm image of I
  int qu_conv;          // Convention for Stokes Q,U.  0 (default) -> East of North (observer).  1 -> North of West
  int quench_output;    // Quench output, i.e. "quench" argument
  int only_unpolarized; // Unpolarized transport only
  int perform_check;    // Output check & recheck images

  // Which e- energy distributions/emissivities to use
  int emission_type;
  // Whether to apply I > 0 "floor" when integrating forward the Stokes parameters
  // Probably harmless!
  int stokes_floors;

  int isolate_counterjet;

  const char dump[STRLEN];
  const char outf[STRLEN];

  // Subrings
  int target_nturns;
  int subring_dtheta;

  // Adaptive tracing
  int nx_min, ny_min;   // dimensions of lowest resolution image
  double refine_abs, refine_rel; // Refinement tolerances
  double refine_cut;    // minimum intensity at which to bother refining
  int nearest_neighbor; // use nearest-neighbor instead of matching-order interpolation

  // ML parameters
  double xoff, yoff;    // in pixels

  // slow light
  double img_cadence;
  double restart_int;

  // Save out variables along a geodesic
  int trace;
  int trace_stride;
  int trace_i, trace_j;
  const char trace_outf[STRLEN];

  // Histogram observer-visible emission from each zone
  int histo, histo_polar;
  const char histo_outf[STRLEN];
} Params;

// modify this to set default values
void load_par_from_argv(int, char *[], Params *);

// modify this to actually add the reading subroutines
void try_set_parameter(const char *, const char *, Params *);

// only modify if you add/modify types
void load_par(const char *, Params *);
void set_by_word_val(const char *word, const char *value, const char *key, void *val, int type);

#endif // PAR_H
