
#include "par.h"
#include "decs.h"
#include "model.h"
#include "model_radiation.h"

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

/*
 * ## program loading order
 * 
 * Knowing this order may be useful for adding new variables that should be
 * loaded / overwritten before / during / after other parameters are set.
 * 
 * 1. ```par.c:load_par(...)``` may be called if triggered by the command line
 * arguments ```-par path/to/par/file```. This function populates the
 * ```Params``` struct with various user-defined variables (often similar to/a
 * superset of values you would pass at runtime) such as inclination angle,
 * M_unit, and so on. If this function is called, the ```params``` data
 * structure is marked as having been loaded.
 * 2. ```model.c:parse_input(...)``` is called. This function is responsible for
 * populating model-defined variables such as M_unit, freqcgs, and so on. These
 * variables tend to live within the model.c/h files (as static globals) and
 * thus should not be referenced from without. If the params data structure is
 * marked as loaded from above, then this function ignores command line input
 * and sets all model variables according to the params structure. (TODO in the
 * future, it might be nice to allow command line arguments to override
 * parameter file values.)
 * 3. ```model.c:init_model(...)``` is called. This function is responsible for
 * setting up the initial data structures associated with the fluid dump file(s)
 * and grid geometry. By the time this function has returned, all initial dump
 * files should be loaded into memory. Before returning, this function may call
 * ```load_*_data(...)``` below.
 * 4. ```model.c:load_*_data(...)``` may be called on initial load or during the
 * ```update_data``` stage, if ipole is running in the slow-light configuration.
 */

// sets default values for elements of params and then calls other
// loading functions
void load_par_from_argv(int argc, char *argv[], Params *params) {

  char *word, *value, *saveptr;

  // set default values here
  params->add_ppm = 0;
  params->qu_conv = 0;
  params->quench_output = 0;
  params->only_unpolarized = 0;
  params->perform_check = 0;
  params->old_centering = 0;

  params->emission_type = 4;
  params->stokes_floors = 0;

  params->isolate_counterjet = 0;

  params->rcam = 1000.;
  params->thetacam = 90.;
  params->phicam = 0.;
  params->rotcam = 0.;
  params->nx = 160;
  params->ny = 160;

  params->eps = 0.01;
  params->maxnstep = 50000;

  params->dsource = DM87_PC; // or DSGRA_PC

  params->restart_int = -1.;

  params->target_nturns = -1;
  params->subring_dtheta = 0;

  params->nx_min = -1;
  params->ny_min = -1;
  // Decent guesses for e.g. Zack Gelles's sample image.
  params->refine_rel = 1.e-1;
  params->refine_abs = 1.e-5;
  // Disable these unless requested
  params->refine_cut = 0.0;
  params->nearest_neighbor = 0;

  params->xoff = 0.0;
  params->yoff = 0.0;

  params->trace = 0;
  params->trace_stride = 1;
  params->trace_i = -1;
  params->trace_j = -1;

  // I'm not sure there's still any advantage to "const" if we do this,
  // but hey, no warnings
  sscanf("trace.h5", "%s", (char *) (void *) params->trace_outf);

  params->histo = 0;
  params->histo_polar = 0;
  sscanf("histo.h5", "%s", (char *) (void *) params->histo_outf);

  // process each command line argument
  for (int i=0; i<argc; ++i) {
    if ( strcmp(argv[i], "-quench") == 0 ) params->quench_output = 1;
    else if ( strcmp(argv[i], "-unpol") == 0 ) params->only_unpolarized = 1;
    else if ( strcmp(argv[i], "-check") == 0 ) params->perform_check = 1;

    // read parameter from command line
    if ( strlen(argv[i])>2 && argv[i][0]=='-' && argv[i][1]=='-' ) {
      word = strtok_r(argv[i]+2, "=", &saveptr);
      value = strtok_r(NULL, "=", &saveptr);
      if (value != NULL) {
        try_set_parameter(word, value, params);
      } else if (i == argc-1 || (argv[i+1][0]=='-' && argv[i+1][1]=='-')) {
        try_set_parameter(word, "1", params);
      } else {
        try_set_parameter(word, argv[i+1], params);
      }
    }
    
    // read parameter file
    if ( strcmp(argv[i],"-par") == 0 ) {
      if ( i == argc-1 ) {
        fprintf(stderr, "! parameter filename must follow -par argument\n");
        exit(1);
      }
      load_par(argv[i+1], params);
    }

  }

  if (params->histo != 0 || params->histo_polar != 0) {
    fprintf(stderr, "histogramming emission is nuanced. please either use the trace feature or uncomment this exit() at your own risk.\n");
    exit(42);
  }
}

// sets parameter file entry if any key==word
void try_set_parameter(const char *word, const char *value, Params *params) {
  // Required parameters
  set_by_word_val(word, value, "freqcgs", &(params->freqcgs), TYPE_DBL);

  // Optional (defaulted) parameters
  set_by_word_val(word, value, "add_ppm", &(params->add_ppm), TYPE_INT);
  set_by_word_val(word, value, "qu_conv", &(params->qu_conv), TYPE_INT);
  set_by_word_val(word, value, "quench_output", &(params->quench_output), TYPE_INT);
  set_by_word_val(word, value, "only_unpolarized", &(params->only_unpolarized), TYPE_INT);
  set_by_word_val(word, value, "emission_type", &(params->emission_type), TYPE_INT);
  set_by_word_val(word, value, "stokes_floors", &(params->stokes_floors), TYPE_INT);
  set_by_word_val(word, value, "counterjet", &(params->isolate_counterjet), TYPE_INT);
  set_by_word_val(word, value, "old_centering", &(params->old_centering), TYPE_INT);

  set_by_word_val(word, value, "rcam", &(params->rcam), TYPE_DBL);
  set_by_word_val(word, value, "thetacam", &(params->thetacam), TYPE_DBL);
  set_by_word_val(word, value, "phicam", &(params->phicam), TYPE_DBL);
  set_by_word_val(word, value, "rotcam", &(params->rotcam), TYPE_DBL);
  set_by_word_val(word, value, "dsource", &(params->dsource), TYPE_DBL);

  // There are many ways to specify a FOV
  set_by_word_val(word, value, "fovx", &(params->fovx_dsource), TYPE_DBL);
  set_by_word_val(word, value, "fovy", &(params->fovy_dsource), TYPE_DBL);
  set_by_word_val(word, value, "fovx_dsource", &(params->fovx_dsource), TYPE_DBL);
  set_by_word_val(word, value, "fovy_dsource", &(params->fovy_dsource), TYPE_DBL);
  set_by_word_val(word, value, "fov", &(params->fovx_dsource), TYPE_DBL);
  set_by_word_val(word, value, "fov", &(params->fovy_dsource), TYPE_DBL);
  // Even in the plane
  set_by_word_val(word, value, "dx", &(params->dx), TYPE_DBL);
  set_by_word_val(word, value, "dy", &(params->dy), TYPE_DBL);

  set_by_word_val(word, value, "nx", &(params->nx), TYPE_INT);
  set_by_word_val(word, value, "ny", &(params->ny), TYPE_INT);
  set_by_word_val(word, value, "nx_min", &(params->nx_min), TYPE_INT);
  set_by_word_val(word, value, "ny_min", &(params->ny_min), TYPE_INT);
  set_by_word_val(word, value, "refine_abs", &(params->refine_abs), TYPE_DBL);
  set_by_word_val(word, value, "refine_rel", &(params->refine_rel), TYPE_DBL);
  set_by_word_val(word, value, "refine_cut", &(params->refine_cut), TYPE_DBL);
  set_by_word_val(word, value, "use_nearest_neighbor", &(params->nearest_neighbor), TYPE_INT);

  set_by_word_val(word, value, "target_nturns", &(params->target_nturns), TYPE_INT);
  set_by_word_val(word, value, "subring_dtheta", &(params->subring_dtheta), TYPE_INT);

  set_by_word_val(word, value, "eps", &(params->eps), TYPE_DBL);
  set_by_word_val(word, value, "maxnstep", &(params->maxnstep), TYPE_INT);

  set_by_word_val(word, value, "xoff", &(params->xoff), TYPE_DBL);
  set_by_word_val(word, value, "yoff", &(params->yoff), TYPE_DBL);

  set_by_word_val(word, value, "outfile", (void *)(params->outf), TYPE_STR);

  // for slow light
  set_by_word_val(word, value, "img_cadence", &(params->img_cadence), TYPE_DBL);
  set_by_word_val(word, value, "restart_int", &(params->restart_int), TYPE_DBL);

  // Save out variables along paths
  set_by_word_val(word, value, "trace", &(params->trace), TYPE_INT);
  set_by_word_val(word, value, "trace_stride", &(params->trace_stride), TYPE_INT);
  set_by_word_val(word, value, "trace_i", &(params->trace_i), TYPE_INT);
  set_by_word_val(word, value, "trace_j", &(params->trace_j), TYPE_INT);
  set_by_word_val(word, value, "trace_outf", (void *)(params->trace_outf), TYPE_STR);

  set_by_word_val(word, value, "histo", &(params->histo), TYPE_INT);
  set_by_word_val(word, value, "histo_polar", &(params->histo_polar), TYPE_INT);
  set_by_word_val(word, value, "histo_outf", (void *)(params->histo_outf), TYPE_STR);

  // Let models add/parse their own parameters we don't understand
  try_set_model_parameter(word, value);

  // Let radiation model load its own parameters
  try_set_radiation_parameter(word, value);
}

// sets default values for elements of params (if desired) and loads from par file 'fname'
void load_par(const char *fname, Params *params) {
 
  char line[256] = {0}, word[256] = {0}, value[256] = {0};
  FILE *fp = fopen(fname, "r");

  if (fp == NULL) {
    fprintf(stderr, "! unable to load parameter file '%s'. (%d: %s)\n", fname, errno, strerror(errno));
    exit(errno);
  }

  // modify parameters/types below
  while (fgets(line, 255, fp) != NULL) {
    if (line[0] == '#') continue;
    sscanf(line, "%s %s", word, value);
    try_set_parameter(word, value, params);
  }

  fclose(fp);

}

// loads value -> (type *)*val if word==key
void set_by_word_val(const char *word, const char *value, const char *key, void *val, int type) {

  if (strcmp(word, key) == 0) {
    switch (type) {
    case TYPE_INT:
      sscanf(value, "%d", (int *)val);
      break;
    case TYPE_DBL:
      sscanf(value, "%lf", (double *)val);
      break;
    case TYPE_STR:
      sscanf(value, "%s", (char *)val);
      break;
    default:
      fprintf(stderr, "! attempt to load unknown type '%d'.\n", type);
    }
  }

}


