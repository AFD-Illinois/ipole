
#include "par.h"
#include "decs.h"
#include "model.h"

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

// sets default values for elements of params and then calls other
// loading functions
void load_par_from_argv(int argc, char *argv[], Params *params) {

  char *word, *value, *saveptr;

  // set default values here
  params->add_ppm = 0;
  params->qu_conv = 0;

  params->rcam = 1000.;
  params->thetacam = 90.;
  params->phicam = 0.;
  params->rotcam = 0.;
  params->nx = 160;
  params->ny = 160;

  params->dsource = DM87_PC; // or DSGRA_PC

  params->restart_int = -1.;

  params->xoff = 0.0;
  params->yoff = 0.0;

  // process each command line argument
  for (int i=0; i<argc; ++i) {

    // read parameter from command line
    if ( strlen(argv[i])>2 && argv[i][0]=='-' && argv[i][1]=='-' ) {
      word = strtok_r(argv[i]+2, "=", &saveptr);
      value = strtok_r(NULL, "=", &saveptr);
      try_set_parameter(word, value, params);
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
}

// sets parameter file entry if any key==word
void try_set_parameter(const char *word, const char *value, Params *params) {
  // Required parameters
  set_by_word_val(word, value, "freqcgs", &(params->freqcgs), TYPE_DBL);

  // Optional (defaulted) parameters
  set_by_word_val(word, value, "add_ppm", &(params->add_ppm), TYPE_INT);
  set_by_word_val(word, value, "qu_conv", &(params->qu_conv), TYPE_INT);

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

  set_by_word_val(word, value, "xoff", &(params->xoff), TYPE_DBL);
  set_by_word_val(word, value, "yoff", &(params->yoff), TYPE_DBL);

  set_by_word_val(word, value, "outfile", (void *)(params->outf), TYPE_STR);

  // for slow light
  set_by_word_val(word, value, "img_cadence", &(params->img_cadence), TYPE_DBL);
  set_by_word_val(word, value, "restart_int", &(params->restart_int), TYPE_DBL);

  // Let models add/parse their own parameters we don't understand
  try_set_model_parameter(word, value);
}

// sets default values for elements of params (if desired) and loads from par file 'fname'
void load_par(const char *fname, Params *params) {
 
  char line[256], word[256], value[256];
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


