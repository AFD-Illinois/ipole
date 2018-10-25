
#include "par.h"

// sets default values for elements of params (if desired) and loads
// from par file 'fname'
void load_par (const char *fname, Params *params) {
 
  char line[256];
  FILE *fp = fopen(fname, "r");

  if (fp == NULL) {
    fprintf(stderr, "! unable to load parameter file '%s'. (%d: %s)\n", fname, errno, strerror(errno));
    exit(errno);
  }

  // set default values here
  params->counterjet = 0;
  params->tp_over_te = 3.;
  params->trat_small = 1.;
  params->trat_large = 30.;
  params->phicam = 0.;

  // modify parameters/types below
  while (fgets(line, 255, fp) != NULL) {

    if (line[0] == '#') continue;

    read_param(line, "counterjet", &(params->counterjet), TYPE_INT);

    read_param(line, "thetacam", &(params->thetacam), TYPE_DBL);
    read_param(line, "phicam", &(params->phicam), TYPE_DBL);
    read_param(line, "freqcgs", &(params->freqcgs), TYPE_DBL);
    read_param(line, "MBH", &(params->MBH), TYPE_DBL);
    read_param(line, "M_unit", &(params->M_unit), TYPE_DBL);
    read_param(line, "tp_over_te", &(params->tp_over_te), TYPE_DBL);
    read_param(line, "trat_small", &(params->trat_small), TYPE_DBL);
    read_param(line, "trat_large", &(params->trat_large), TYPE_DBL);

    read_param(line, "dump", (void *)(params->dump), TYPE_STR);
    read_param(line, "outfile", (void *)(params->outf), TYPE_STR);
    
  }

  fclose(fp);

  params->loaded = 1;

}

// loads value -> (type *)*val if line corresponds to (key,value)
void read_param (const char *line, const char *key, void *val, int type) {

  char word[256], value[256];

  sscanf(line, "%s %s", word, value);

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



