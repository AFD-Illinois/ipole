
#include "par.h"

// sets default values for elements of params and then calls other
// loading functions
void load_par_from_argv (int argc, char *argv[], Params *params) {

  char *word, *value, *saveptr;

  // set default values here
  params->counterjet = 0;
  params->tp_over_te = 3.;
  params->trat_small = 1.;
  params->trat_large = 40.;
  params->phicam = 0.;
  params->rotcam = 0.;
  params->dump_skip = 1;
  params->restart_int = -1.;

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

  // just in case
  params->loaded = 1;

}

// sets parameter file entry if any key==word
void try_set_parameter (const char *word, const char *value, Params *params) {

  set_by_word_val(word, value, "counterjet", &(params->counterjet), TYPE_INT);

  set_by_word_val(word, value, "thetacam", &(params->thetacam), TYPE_DBL);
  set_by_word_val(word, value, "phicam", &(params->phicam), TYPE_DBL);
  set_by_word_val(word, value, "rotcam", &(params->rotcam), TYPE_DBL);
  set_by_word_val(word, value, "freqcgs", &(params->freqcgs), TYPE_DBL);
  set_by_word_val(word, value, "MBH", &(params->MBH), TYPE_DBL);
  set_by_word_val(word, value, "M_unit", &(params->M_unit), TYPE_DBL);
  set_by_word_val(word, value, "tp_over_te", &(params->tp_over_te), TYPE_DBL);
  set_by_word_val(word, value, "trat_small", &(params->trat_small), TYPE_DBL);
  set_by_word_val(word, value, "trat_large", &(params->trat_large), TYPE_DBL);

  set_by_word_val(word, value, "dump", (void *)(params->dump), TYPE_STR);
  set_by_word_val(word, value, "outfile", (void *)(params->outf), TYPE_STR);

  // for slow light
  set_by_word_val(word, value, "dump_min", &(params->dump_min), TYPE_INT);
  set_by_word_val(word, value, "dump_max", &(params->dump_max), TYPE_INT);
  set_by_word_val(word, value, "dump_skip", &(params->dump_skip), TYPE_INT);
  set_by_word_val(word, value, "img_cadence", &(params->img_cadence), TYPE_DBL);
  set_by_word_val(word, value, "restart_int", &(params->restart_int), TYPE_DBL);

}

// sets default values for elements of params (if desired) and loads from par file 'fname'
void load_par (const char *fname, Params *params) {
 
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

  params->loaded = 1;

}

// loads value -> (type *)*val if word==key
void set_by_word_val (const char *word, const char *value, const char *key, void *val, int type) {

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


