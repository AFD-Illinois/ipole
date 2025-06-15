#include "model.h"

#include "decs.h"
#include "hdf5_utils.h"
#include "debug_tools.h"
#include "dict.h"

#include "coordinates.h"
#include "geometry.h"
#include "grid.h"
#include "model_radiation.h"  // Only for outputting emissivities
#include "par.h"
#include "utils.h"

#include "debug_tools.h"

#include <assert.h>
#include <string.h>
#include <ctype.h>


// Macros
#define NSUP (3) //how many files to load for slow light tracing
#define NVAR (10)
#define USE_FIXED_TPTE (0)
#define USE_MIXED_TPTE (1)
#define USE_GEODESIC_SIGMACUT (1)
/* ELECTRONS
*    0 : constant TP_OVER_TE
*    1 : use dump file model (kawazura?)
*    2 : use mixed TP_OVER_TE (beta model)
*    3 : use mixed TP_OVER_TE (beta model) with fluid temperature
* TODO the way this is selected is horrid.  Make it a parameter.
*/
#define ELECTRONS_TFLUID (3)


// Debug reader
// Save metric and physical quantities
#define DEBUG_READER (0)

// Units
double M_unit;
double L_unit;
double T_unit;
double RHO_unit;
double U_unit;
double B_unit;
double Te_unit;

// Thermodynamic parameters
static int RADIATION, ELECTRONS;
static double gam = 1.444444444444444, game = 1.333333333333333, gamp = 1.666666666666667;
static double Thetae_unit, Mdotedd;
// Molecular weights
static double Ne_factor = 1.;  // e.g., used for He with 2 protons+neutrons per 2 electrons
static double mu_i, mu_e, mu_tot;
// Temperature ratio
static double tp_over_te = 3.;
static double trat_small = 1.;
static double trat_large = 40.;
double beta_crit = 1.0;
static char electron_subgrid_model_string[STRLEN] = "UNKNOWN";
// enum to represent electron heating model
enum ElectronSubgridModel {
  ELECTRONS_UNKNOWN,
  ELECTRONS_CONSTANT,
  ELECTRONS_HOWES,
  ELECTRONS_KAWAZURA,
  ELECTRONS_WERNER,
  ELECTRONS_ROWAN,
  ELECTRONS_SHARMA
};
// Converts input model string to an enum value for electron subgrid model
int get_electron_subgrid_model(const char* electron_subgrid_model_string) {
  if (strcmp(electron_subgrid_model_string, "CONSTANT") == 0) {
    return ELECTRONS_CONSTANT;
  } else if (strcmp(electron_subgrid_model_string, "HOWES") == 0) {
    return ELECTRONS_HOWES;
  } else if (strcmp(electron_subgrid_model_string, "KAWAZURA") == 0) {
    return ELECTRONS_KAWAZURA;
  } else if (strcmp(electron_subgrid_model_string, "WERNER") == 0) {
    return ELECTRONS_WERNER;
  } else if (strcmp(electron_subgrid_model_string, "ROWAN") == 0) {
    return ELECTRONS_ROWAN;
  } else if (strcmp(electron_subgrid_model_string, "SHARMA") == 0) {
    return ELECTRONS_SHARMA;
  } else {
    return ELECTRONS_UNKNOWN;
  }
}

// Geodesic parameters
double sigma_cut = 1.;
double sigma_cut_high = -1.0;
// Ignore radiation interactions within one degree of polar axis
static double polar_cut = -1;
static double th_beg = 0.0174;

// File loading parameters
static int nloaded = 0;
static char fnam[STRLEN] = "dump.h5";
static int dumpskip = 1;
static int dumpmin, dumpmax, dumpidx;
double DTd;
static char kharma_format[STRLEN] = "new";

// Black hole parameters
static double MBH_solar = 6.2e9;
static double MBH; // Set from previous
static double Mdot_dump;
static double MdotEdd_dump;
static double Ladv_dump;

// Other parameters
static int reverse_field = 0;
double tf;
static hdf5_blob fluid_header = {0};


// Struct to store model-specific data
struct of_data {
  double t;
  double ****p; // NVAR, N1, N2, N3
  double ***ne; // N1, N2, N3
  double ***thetae;
  double ***b;
  double ***sigma;
  double ***beta;
};
static struct of_data dataA, dataB, dataC; // additional structs for slow light
static struct of_data *data[NSUP];


// Dictionary for model parameters
dict *model_params = NULL;


/**
 * @brief Attempts to set a model parameter based on the provided key and value.
 *
 * Various model parameters:
 *
 * - **Mass Parameters:**  
 *   - "MBH": Mass of the black hole (in solar masses).  
 *   - "M_unit": Mass unit.
 *
 * - **Dump-related Parameters:**  
 *   - "dump": Dump file name.  
 *   - "dump_min": Minimum dump index (for slow light mode).  
 *   - "dump_max": Maximum dump index (for slow light mode).  
 *   - "dump_skip": Dump skip value (for slow light mode).  
 *     After setting these, the dump index (@c dumpidx) is initialized to dump_min.
 *
 * - **Thermodynamic Parameters:**  
 *   - "tp_over_te": Ratio of proton to electron temperature.  
 *   - "trat_small": Small temperature ratio.  
 *   - "trat_large": Large temperature ratio.  
 *   - "beta_crit": Critical plasma beta.
 *
 * - **Geodesic Parameters:**  
 *   - "rmax_geo": Maximum radius for geodesic integration.  
 *   - "rmin_geo": Minimum radius for geodesic integration.  
 *   - "sigma_cut": Sigma cutoff value.  
 *   - "sigma_cut_high": Higher sigma cutoff value.  
 *   - "polar_cut_deg": Angular cut (in degrees) to remove the spine.
 *
 * - **Other Parameters:**  
 *   - "reverse_field": Flag to reverse the field orientation.
 *
 * The function utilizes @c set_by_word_val to set a model parameter based on the provided key and value.
 *
 * @param word  The key identifying which model parameter to set.
 * @param value The value to be assigned to the model parameter.
 */
void try_set_model_parameter(const char *word, const char *value)
{
  // Mass parameters
  set_by_word_val(word, value, "MBH", &MBH_solar, TYPE_DBL);
  set_by_word_val(word, value, "M_unit", &M_unit, TYPE_DBL);

  // Dump-related parameters
  set_by_word_val(word, value, "dump", (void *)fnam, TYPE_STR);
  // for slow light
  set_by_word_val(word, value, "dump_min", &dumpmin, TYPE_INT);
  set_by_word_val(word, value, "dump_max", &dumpmax, TYPE_INT);
  set_by_word_val(word, value, "dump_skip", &dumpskip, TYPE_INT);
  dumpidx = dumpmin;

  // Thermodynamic parameters
  set_by_word_val(word, value, "tp_over_te", &tp_over_te, TYPE_DBL);
  set_by_word_val(word, value, "trat_small", &trat_small, TYPE_DBL);
  set_by_word_val(word, value, "trat_large", &trat_large, TYPE_DBL);
  set_by_word_val(word, value, "beta_crit", &beta_crit, TYPE_DBL);

  // Geodesic parameters
  set_by_word_val(word, value, "rmax_geo", &rmax_geo, TYPE_DBL);
  set_by_word_val(word, value, "rmin_geo", &rmin_geo, TYPE_DBL);
  set_by_word_val(word, value, "sigma_cut", &sigma_cut, TYPE_DBL);
  set_by_word_val(word, value, "sigma_cut_high", &sigma_cut_high, TYPE_DBL);
  // allow cutting out the spine
  set_by_word_val(word, value, "polar_cut_deg", &polar_cut, TYPE_DBL);

  // Other parameters
  set_by_word_val(word, value, "reverse_field", &reverse_field, TYPE_INT);

  // Electron heating model
  set_by_word_val(word, value, "electron_subgrid_model", (void *)electron_subgrid_model_string, TYPE_STR);

  // KHARMA file format. Older (v5) files have different data ordering
  set_by_word_val(word, value, "kharma_format", (void *)kharma_format, TYPE_STR);
}


/**
 * @brief Opens the file, enters the "Info" group, and reads the specified attribute.
 *
 * This function opens the specified KHARMA dump file, navigates to the "Info" group,
 * and reads the attribute identified by @p attr_name. The attribute value is stored
 * in the memory pointed to by @p value. In addition to scalar integer, double, and string
 * attributes, it now supports integer and double array attributes.
 *
 * For scalar attributes (TYPE_INT, TYPE_DBL), the attribute is expected to contain
 * exactly one element; for array attributes (TYPE_INT_ARRAY, TYPE_DBL_ARRAY), the caller
 * must allocate a buffer large enough for all elements (the number of elements is not returned).
 *
 * @param filename   Name of the KHARMA dump file.
 * @param attr_name  Name of the attribute to read from the "Info" group.
 * @param type       Expected type of the attribute:
 *                   - TYPE_INT: scalar integer
 *                   - TYPE_DBL: scalar double
 *                   - TYPE_STR: string (fixed-length assumed; for variable-length, additional handling is needed)
 *                   - TYPE_INT_ARRAY: integer array
 *                   - TYPE_DBL_ARRAY: double array
 * @param value      Pointer to memory where the attribute value will be stored.
 *                   For TYPE_INT and TYPE_DBL, pass a pointer to an int or double.
 *                   For TYPE_INT_ARRAY and TYPE_DBL_ARRAY, pass a pointer to a pre-allocated buffer.
 *                   For TYPE_STR, pass a character buffer.
 * @param value_size For strings, the size of the provided buffer; ignored for numeric types.
 *
 * @return 0 on success, -1 on error.
 */
int read_info_attribute(const char *filename, const char *attr_name, int type, void *value, size_t value_size)
{
    hid_t file_id = -1, group_id = -1, attr_id = -1, attr_type_id = -1;
    H5T_class_t type_class;
    herr_t status;
    int ret = 0;

    /* Open the file in read-only mode */
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return -1;
    }

    /* Open the "Info" group */
    group_id = H5Gopen(file_id, "Info", H5P_DEFAULT);
    if (group_id < 0) {
        fprintf(stderr, "Error opening group 'Info'\n");
        H5Fclose(file_id);
        return -1;
    }

    /* Open the specified attribute */
    attr_id = H5Aopen(group_id, attr_name, H5P_DEFAULT);
    if (attr_id < 0) {
        fprintf(stderr, "Error opening attribute: %s\n", attr_name);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return -1;
    }

    /* Get the attribute's datatype and class */
    attr_type_id = H5Aget_type(attr_id);
    type_class = H5Tget_class(attr_type_id);

    if (type == TYPE_INT) {
        if (type_class != H5T_INTEGER) {
            fprintf(stderr, "Attribute %s is not an integer.\n", attr_name);
            ret = -1;
            goto cleanup;
        }
        /* Ensure the attribute is scalar */
        {
            hid_t space_id = H5Aget_space(attr_id);
            hsize_t npoints = H5Sget_simple_extent_npoints(space_id);
            if (npoints != 1) {
                fprintf(stderr, "Expected scalar integer but found array with %llu elements.\n",
                        (unsigned long long)npoints);
                ret = -1;
                H5Sclose(space_id);
                goto cleanup;
            }
            H5Sclose(space_id);
        }
        status = H5Aread(attr_id, H5T_NATIVE_INT, value);
        if (status < 0) {
            fprintf(stderr, "Error reading integer attribute %s\n", attr_name);
            ret = -1;
            goto cleanup;
        }
    } else if (type == TYPE_DBL) {
        if (type_class != H5T_FLOAT) {
            fprintf(stderr, "Attribute %s is not a float/double.\n", attr_name);
            ret = -1;
            goto cleanup;
        }
        /* Ensure the attribute is scalar */
        {
            hid_t space_id = H5Aget_space(attr_id);
            hsize_t npoints = H5Sget_simple_extent_npoints(space_id);
            if (npoints != 1) {
                fprintf(stderr, "Expected scalar double but found array with %llu elements.\n",
                        (unsigned long long)npoints);
                ret = -1;
                H5Sclose(space_id);
                goto cleanup;
            }
            H5Sclose(space_id);
        }
        status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, value);
        if (status < 0) {
            fprintf(stderr, "Error reading double attribute %s\n", attr_name);
            ret = -1;
            goto cleanup;
        }
    } else if (type == TYPE_STR) {
        if (type_class != H5T_STRING) {
            fprintf(stderr, "Attribute %s is not a string.\n", attr_name);
            ret = -1;
            goto cleanup;
        }
        /* For fixed-length strings we can read directly.
           (For variable-length strings, additional handling is required.) */
        status = H5Aread(attr_id, attr_type_id, value);
        if (status < 0) {
            fprintf(stderr, "Error reading string attribute %s\n", attr_name);
            ret = -1;
            goto cleanup;
        }
        /* Ensure null termination in case the string isn’t terminated */
        ((char *)value)[value_size - 1] = '\0';
    } else if (type == TYPE_INT_ARRAY) {
        if (type_class != H5T_INTEGER) {
            fprintf(stderr, "Attribute %s is not an integer array.\n", attr_name);
            ret = -1;
            goto cleanup;
        }
        {
            hid_t space_id = H5Aget_space(attr_id);
            /* Optionally, you could query npoints here if needed */
            status = H5Aread(attr_id, H5T_NATIVE_INT, value);
            H5Sclose(space_id);
            if (status < 0) {
                fprintf(stderr, "Error reading integer array attribute %s\n", attr_name);
                ret = -1;
                goto cleanup;
            }
        }
    } else if (type == TYPE_DBL_ARRAY) {
        if (type_class != H5T_FLOAT) {
            fprintf(stderr, "Attribute %s is not a double array.\n", attr_name);
            ret = -1;
            goto cleanup;
        }
        {
            hid_t space_id = H5Aget_space(attr_id);
            status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, value);
            H5Sclose(space_id);
            if (status < 0) {
                fprintf(stderr, "Error reading double array attribute %s\n", attr_name);
                ret = -1;
                goto cleanup;
            }
        }
    } else {
        fprintf(stderr, "Unsupported parameter type.\n");
        ret = -1;
        goto cleanup;
    }

cleanup:
    if (attr_type_id >= 0)
        H5Tclose(attr_type_id);
    if (attr_id >= 0)
        H5Aclose(attr_id);
    if (group_id >= 0)
        H5Gclose(group_id);
    if (file_id >= 0)
        H5Fclose(file_id);

    return ret;
}


/**
 * @brief Reads the entire "File" attribute from the "Input" group in an HDF5 file.
 *
 * Given the native KHARMA dump file in @p fname, this function opens the file, navigates to the
 * "Input" group, and reads the attribute named "File". The attribute may be stored as a
 * fixed-length or variable-length string.
 *
 * @param fname Name of KAHRMA dump file.
 * @return Pointer to a null-terminated string containing input par file.
 */
char *get_parfile(const char *fname)
{
    hid_t file_id, group_id, attr_id, attr_type;
    char *attr_data = NULL;
    herr_t status;

    /* Open the HDF5 file in read-only mode */
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Failed to open file %s\n", fname);
        return NULL;
    }

    /* Open the "Input" group */
    group_id = H5Gopen(file_id, "Input", H5P_DEFAULT);
    if (group_id < 0) {
        fprintf(stderr, "Failed to open group 'Input'\n");
        H5Fclose(file_id);
        return NULL;
    }

    /* Open the "File" attribute */
    attr_id = H5Aopen(group_id, "File", H5P_DEFAULT);
    if (attr_id < 0) {
        fprintf(stderr, "Failed to open attribute 'File'\n");
        H5Gclose(group_id);
        H5Fclose(file_id);
        return NULL;
    }

    /* Get the attribute's datatype */
    attr_type = H5Aget_type(attr_id);
    if (attr_type < 0) {
        fprintf(stderr, "Failed to get attribute type\n");
        H5Aclose(attr_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return NULL;
    }

    /* Check if the attribute is a variable-length string */
    if (H5Tis_variable_str(attr_type)) {
        /* For variable-length strings, HDF5 allocates the memory */
        char *vldata = NULL;
        status = H5Aread(attr_id, attr_type, &vldata);
        if (status < 0) {
            fprintf(stderr, "Failed to read variable-length attribute 'File'\n");
            H5Tclose(attr_type);
            H5Aclose(attr_id);
            H5Gclose(group_id);
            H5Fclose(file_id);
            return NULL;
        }
        /* Duplicate the string into our own allocated memory */
        attr_data = strdup(vldata);
        /* Free the memory allocated by HDF5 */
        H5free_memory(vldata);
    } else {
        /* For fixed-length strings, determine the size */
        size_t type_size = H5Tget_size(attr_type);
        if (type_size == 0) {
            fprintf(stderr, "Attribute size is 0\n");
            H5Tclose(attr_type);
            H5Aclose(attr_id);
            H5Gclose(group_id);
            H5Fclose(file_id);
            return NULL;
        }
        attr_data = (char *)malloc(type_size + 1);
        if (attr_data == NULL) {
            fprintf(stderr, "Failed to allocate memory for attribute string\n");
            H5Tclose(attr_type);
            H5Aclose(attr_id);
            H5Gclose(group_id);
            H5Fclose(file_id);
            return NULL;
        }
        status = H5Aread(attr_id, attr_type, attr_data);
        if (status < 0) {
            fprintf(stderr, "Failed to read attribute 'File'\n");
            free(attr_data);
            attr_data = NULL;
        } else {
            /* Ensure null termination */
            attr_data[type_size] = '\0';
        }
    }

    /* Clean up HDF5 objects */
    H5Tclose(attr_type);
    H5Aclose(attr_id);
    H5Gclose(group_id);
    H5Fclose(file_id);

    return attr_data;
}


/**
 * @brief Extract the value for a given key within a block (optional) from the KHARMA par file
 *
 * If a block name is provided (e.g. "parthenon/mesh"), the function first locates the
 * block (assumed to appear as "<parthenon/mesh>") and then searches for the key within that block.
 * If no block name is passed (i.e. block == NULL or block[0]=='\0'), the entire input is searched,
 * returning the first occurrence of the key.
 *
 * @param input_text The full text of the par file attribute.
 * @param block Optional block name (without angle brackets). For example, "parthenon/mesh".
 *              If provided, only key–value pairs within that block are considered.
 * @param key The key to search for (e.g. "nx1").
 * @param type The expected type (TYPE_INT, TYPE_DBL, or TYPE_STR).
 * @param value Pointer to where the converted value should be stored.
 *              For TYPE_INT and TYPE_DBL, pass a pointer to an int or double.
 *              For TYPE_STR, pass a character buffer.
 * @param value_size For strings, the size of the provided buffer; ignored for numeric types.
 *
 * @return 0 on success, or -1 if the key (or block) is not found or on conversion error.
 */
int get_parameter_value(const char *input_text, const char *block, const char *key, int type, void *value, size_t value_size)
{
  const char *search_region = input_text;
  const char *region_end = input_text + strlen(input_text);

  /* If a block is specified, narrow the search region to that block */
  if (block != NULL && block[0] != '\0') {
    char block_tag[256];
    snprintf(block_tag, sizeof(block_tag), "<%s>", block);
    const char *block_pos = strstr(input_text, block_tag);
    if (block_pos == NULL) {
      fprintf(stderr, "Block '%s' not found in input.\n", block);
      return -1;
    }
    search_region = block_pos;
    /* Define the end of the block as the start of the next block marker, if any */
    const char *next_block = strstr(block_pos + strlen(block_tag), "<");
    if (next_block != NULL) {
      region_end = next_block;
    }
  }

  const char *pos = search_region;
  size_t key_len = strlen(key);
  const char *key_loc = NULL;

  /* Loop to find an occurrence of the key that is not part of a longer word
  and lies within the designated search region. */
  while ((pos = strstr(pos, key)) != NULL && pos < region_end) {
    /* Check that the character before the key is either the start of the region or whitespace */
    if (pos != search_region && !isspace((unsigned char)*(pos - 1))) {
      pos++; // Not a valid match, continue search
      continue;
    }
    /* Check that the character after the key is either whitespace or '=' */
    char after = pos[key_len];
    if (after != '\0' && !isspace((unsigned char)after) && after != '=') {
      pos++; // Not a valid match, continue search
      continue;
    }
    /* Valid key found */
    key_loc = pos;
    break;
  }

  if (key_loc == NULL || key_loc >= region_end) {
    fprintf(stderr, "Key '%s' not found in %sinput.\n", key, (block && block[0] != '\0') ? block : "default ");
    return -1;
  }

  /* Locate the '=' sign after the key */
  const char *equals_sign = strchr(key_loc, '=');
  if (equals_sign == NULL || equals_sign >= region_end) {
    fprintf(stderr, "No '=' found for key '%s'.\n", key);
    return -1;
  }

  /* Move past '=' and skip any leading whitespace */
  equals_sign++;
  while (equals_sign < region_end && *equals_sign && isspace((unsigned char)*equals_sign)) {
    equals_sign++;
  }

  /* Extract and convert the value based on the expected type */
  if (type == TYPE_INT) {
    *(int *)value = atoi(equals_sign);
  } else if (type == TYPE_DBL) {
    *(double *)value = atof(equals_sign);
  } else if (type == TYPE_STR) {
    char *str_value = (char *)value;
    size_t i = 0;
    while (equals_sign < region_end && equals_sign[i] && !isspace((unsigned char)equals_sign[i]) 
    && equals_sign[i] != '#' && i < value_size - 1) {
      str_value[i] = equals_sign[i];
      i++;
    }
    str_value[i] = '\0';
  } else {
    fprintf(stderr, "Unsupported parameter type.\n");
    return -1;
  }

  return 0;
}


/**
 * @brief Get the time of a dump file.
 *
 * The time is read from the "Time" attribute in the "Info" group.
 * 
 * @param fnam    Filename of the dump file.
 * @param dumpidx Index of the dump file.
 *
 * @return The time of the dump file.
 */
double get_dump_time(char *fnam, int dumpidx)
{
  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);
  
  double t = -1.;
  read_info_attribute(fname, "Time", TYPE_DBL, &t, 0);

  return t;
}


/**
 * @brief Advance through dumps until the current time is closer to the target.
 *
 * This function iterates through data dumps until the time indicated by @p tA
 * is near the target time (@p tgt). It is primarily used when restarting from
 * a slowlight restart file, ensuring that the restart process begins from a dump
 * corresponding to a time close to @p tgt.
 *
 * @param tA Pointer to the time variable for dump A.
 * @param tB Pointer to the time variable for dump B.
 * @param tgt Target time value to advance towards.
 */
void update_data_until(double *tA, double *tB, double tgt)
{
  double tC = data[2]->t;

  while (tC < tgt) {
    dumpidx += dumpskip;
    tC = get_dump_time(fnam, dumpidx);
  }

  // reset dump index, just to be safe, then load on out ...
  dumpidx -= dumpskip;
  while (*tA < tgt) update_data(tA, tB);
}


/**
 * @brief Load the next expected dump into memory.
 *
 * This function uses an internal dump index variable to load the next "expected"
 * dump into memory (used for slowlight mode). After calling this function, it is
 * guaranteed that the data dumps are ordered in time:
 *
 *   data[0]->t < data[1]->t < data[2]->t
 *
 * The function uses internal pointers (dataA, dataB, dataC) to store the actual
 * data locations and then swaps which members live where to maintain the correct
 * temporal order.
 *
 * @param tA Pointer to the current time variable for dump A.
 * @param tB Pointer to the current time variable for dump B.
 */
void update_data(double *tA, double *tB)
{
  #if SLOW_LIGHT
  // Reorder dataA, dataB, dataC in data[]
  if (nloaded % 3 == 0) {
    data[0] = &dataB;
    data[1] = &dataC;
    data[2] = &dataA;
  } else if (nloaded % 3 == 1) {
    data[0] = &dataC;
    data[1] = &dataA;
    data[2] = &dataB;
  } else {
    data[0] = &dataA;
    data[1] = &dataB;
    data[2] = &dataC;
  }
  int nextdumpidx = dumpidx;
  dumpidx += dumpskip;
  if (nextdumpidx > dumpmax) {
    load_kharma_data(2, fnam, --nextdumpidx, 0);
    data[2]->t += 1.;
  } else {
    load_kharma_data(2, fnam, nextdumpidx, 0);
  }
  *tA = data[0]->t;
  *tB = data[1]->t;
  fprintf(stderr, "loaded data (dump %d) (%g < t < %g)\n", nextdumpidx, *tA, *tB);
  #else // FAST LIGHT
  if (nloaded % 3 == 0) {
    data[0] = &dataA;
    data[1] = &dataB;
    data[2] = &dataC;
  } else if (nloaded % 3 == 1) {
    data[0] = &dataB;
    data[1] = &dataC;
    data[2] = &dataA;
  } else if (nloaded % 3 == 2) {
    data[0] = &dataC;
    data[1] = &dataA;
    data[2] = &dataB;
  } else {
    printf("Fail! nloaded = %i nloaded mod 3 = %i\n", nloaded, nloaded % 3);
  }
  data[2]->t = data[1]->t + DTd;
  #endif 
}


/**
 * @brief Allocate memory for primitives.
 * 
 * Once the grid parameters have been read from the dump file, this function allocates
 * memory only for the primitive variables defined over the entire mesh. This is to 
 * minimize memory usage.
 */
void init_storage_prims(void)
{ 
  /* One ghost zone on each side of the domain */
  for (int n = 0; n < NSUP; n++) {
    data[n]->p = malloc_rank4(NVAR,N1+2,N2+2,N3+2);
  }
}


/**
 * @brief Allocate memory for derived variables in the data struct.
 * 
 * Once the primitives are defined over the entire mesh, the derived variables:
 * electron number density, electron temperature, magnetic field, magnetization, 
 * and plasma beta over the entire grid, are allocated memory.
 */
void init_storage_derived_quantities(void)
{ 
  /* One ghost zone on each side of the domain */
  for (int n = 0; n < NSUP; n++) {
    data[n]->ne     = malloc_rank3(N1+2,N2+2,N3+2);
    data[n]->thetae = malloc_rank3(N1+2,N2+2,N3+2);
    data[n]->b      = malloc_rank3(N1+2,N2+2,N3+2);
    data[n]->sigma  = malloc_rank3(N1+2,N2+2,N3+2);
    data[n]->beta   = malloc_rank3(N1+2,N2+2,N3+2);
  }
}


/**
 * @brief Read relevant fluid and geometry parameters and allocate memory for 'data' object.
 *
 * The relevant parameters (data stoted in header in iharm3d format) are saved to the model_params
 * dict object. The function then allocates memory for the 'data' object.
 * 
 * @param fnam    Filename of the dump file.
 * @param dumpidx Index of the dump file.
 *
 * @return 0 on success, -1 on error.
 */
int read_parameters_and_allocate_memory(char *fnam, int dumpidx)
{
  /* Get input parameter file from dump*/
  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);
  char *parfile = get_parfile(fname);

  /* Declaring necessary parameters */
  /* Thermodynamic parameters*/
  char has_electrons[20];
  char problem_id[20];
  /* Grid parameter*/
  char coordinate_system[256], base[256];
  int nx1_mb, nx2_mb, nx3_mb; // Meshblock size
  int nghost;
  double x1min, x1max, x2min, x2max, x3min, x3max; // Domain limits
  double dx1, dx2, dx3;
  /* Temporal parameters */
  double tlim;
  double dt, time, cfl, dt_min;
  int ncycle;
  /* Physical parameters */
  double rin, rmax;
  char bfield_type[256];

  char buffer[100];
  /* Read parameters from Info group */
  read_info_attribute(fname, "dt", TYPE_DBL, &dt, 0);
  dict_add(model_params, "dt", (snprintf(buffer, sizeof(buffer), "%.8g", dt), buffer));
  read_info_attribute(fname, "NCycle", TYPE_INT, &ncycle, 0);
  dict_add(model_params, "ncycle", (snprintf(buffer, sizeof(buffer), "%d", ncycle), buffer));
  read_info_attribute(fname, "Time", TYPE_DBL, &time, 0);
  dict_add(model_params, "time", (snprintf(buffer, sizeof(buffer), "%.8g", time), buffer));

  /* Read parameters from par file */
  get_parameter_value(parfile, "parthenon/job", "problem_id", TYPE_STR, problem_id, sizeof(problem_id));
  dict_add(model_params, "problem_id", problem_id);
  if (strcmp(problem_id, "torus") == 0) {
    get_parameter_value(parfile, "torus", "rin", TYPE_DBL, &rin, 0);
    dict_add(model_params, "rin", (snprintf(buffer, sizeof(buffer), "%.8g", rin), buffer));
    get_parameter_value(parfile, "torus", "rmax", TYPE_DBL, &rmax, 0);
    dict_add(model_params, "rmax", (snprintf(buffer, sizeof(buffer), "%.8g", rmax), buffer));
    get_parameter_value(parfile, "b_field", "type", TYPE_STR, bfield_type, sizeof(bfield_type));
    dict_add(model_params, "bfield_type", bfield_type);
  }
  get_parameter_value(parfile, "parthenon/time", "tlim", TYPE_DBL, &tlim, 0);
  dict_add(model_params, "tlim", (snprintf(buffer, sizeof(buffer), "%.8g", tlim), buffer));
  get_parameter_value(parfile, "parthenon/time", "dt_min", TYPE_DBL, &dt_min, 0);
  dict_add(model_params, "dt_min", (snprintf(buffer, sizeof(buffer), "%.8g", dt_min), buffer));
  get_parameter_value(parfile, "electrons", "on", TYPE_STR, has_electrons, sizeof(has_electrons));
  dict_add(model_params, "has_electrons", has_electrons);
  if (strncmp(has_electrons, "true", 19) == 0) {
    get_parameter_value(parfile, "electrons", "gamma_e", TYPE_DBL, &game, 0);
    get_parameter_value(parfile, "electrons", "gamma_p", TYPE_DBL, &gamp, 0);
    dict_add(model_params, "game", (snprintf(buffer, sizeof(buffer), "%.8g", game), buffer));
    dict_add(model_params, "gamp", (snprintf(buffer, sizeof(buffer), "%.8g", gamp), buffer));

    enum ElectronSubgridModel electron_subgrid_model = get_electron_subgrid_model(electron_subgrid_model_string);
    dict_add(model_params, "electron_subgrid_model", (snprintf(buffer, sizeof(buffer), "%d", electron_subgrid_model), buffer));
    ELECTRONS = 1;
  } else {
    ELECTRONS = 0;
  }
  Te_unit = Thetae_unit;

  get_parameter_value(parfile, "GRMHD", "cfl", TYPE_DBL, &cfl, 0);
  dict_add(model_params, "cfl", (snprintf(buffer, sizeof(buffer), "%.8g", cfl), buffer));
  get_parameter_value(parfile, "GRMHD", "gamma", TYPE_DBL, &gam, 0);
  dict_add(model_params, "gam", (snprintf(buffer, sizeof(buffer), "%.8g", gam), buffer));
  get_parameter_value(parfile, "coordinates", "base", TYPE_STR, coordinate_system, sizeof(coordinate_system));
  dict_add(model_params, "base", coordinate_system);
  get_parameter_value(parfile, "coordinates", "transform", TYPE_STR, coordinate_system, sizeof(coordinate_system));
  dict_add(model_params, "coordinate_system", coordinate_system);
  if (strncmp(coordinate_system, "fmks", 19) == 0) {
    metric = METRIC_FMKS;
    cstopx[2] = 1.0;
  } else if (strncmp(coordinate_system, "mks", 19) == 0) {
    metric = METRIC_MKS;
    cstopx[2] = 1.0;
  } else if (strncmp(coordinate_system, "eks", 19) == 0) {
    metric = METRIC_EKS;
    cstopx[2] = M_PI;
  } else {
    fprintf(stderr, "Unknown coordinate system: %s\n", coordinate_system);
    return -1;
  }
  get_parameter_value(parfile, "parthenon/mesh", "nx1", TYPE_INT, &N1, 0);
  dict_add(model_params, "nx1", (snprintf(buffer, sizeof(buffer), "%d", N1), buffer));
  get_parameter_value(parfile, "parthenon/mesh", "nx2", TYPE_INT, &N2, 0);
  dict_add(model_params, "nx2", (snprintf(buffer, sizeof(buffer), "%d", N2), buffer));
  get_parameter_value(parfile, "parthenon/mesh", "nx3", TYPE_INT, &N3, 0);
  dict_add(model_params, "nx3", (snprintf(buffer, sizeof(buffer), "%d", N3), buffer));
  get_parameter_value(parfile, "parthenon/mesh", "nghost", TYPE_INT, &nghost, 0);
  dict_add(model_params, "nghost", (snprintf(buffer, sizeof(buffer), "%d", nghost), buffer));
  get_parameter_value(parfile, "parthenon/meshblock", "nx1", TYPE_INT, &nx1_mb, 0);
  dict_add(model_params, "nx1_mb", (snprintf(buffer, sizeof(buffer), "%d", nx1_mb), buffer));
  get_parameter_value(parfile, "parthenon/meshblock", "nx2", TYPE_INT, &nx2_mb, 0);
  dict_add(model_params, "nx2_mb", (snprintf(buffer, sizeof(buffer), "%d", nx2_mb), buffer));
  get_parameter_value(parfile, "parthenon/meshblock", "nx3", TYPE_INT, &nx3_mb, 0);
  dict_add(model_params, "nx3_mb", (snprintf(buffer, sizeof(buffer), "%d", nx3_mb), buffer));
  get_parameter_value(parfile, "parthenon/mesh", "x1min", TYPE_DBL, &x1min, 0);
  dict_add(model_params, "x1min", (snprintf(buffer, sizeof(buffer), "%.8g", x1min), buffer));
  get_parameter_value(parfile, "parthenon/mesh", "x1max", TYPE_DBL, &x1max, 0);
  dict_add(model_params, "x1max", (snprintf(buffer, sizeof(buffer), "%.8g", x1max), buffer));
  get_parameter_value(parfile, "parthenon/mesh", "x2min", TYPE_DBL, &x2min, 0);
  dict_add(model_params, "x2min", (snprintf(buffer, sizeof(buffer), "%.8g", x2min), buffer));
  get_parameter_value(parfile, "parthenon/mesh", "x2max", TYPE_DBL, &x2max, 0);
  dict_add(model_params, "x2max", (snprintf(buffer, sizeof(buffer), "%.8g", x2max), buffer));
  get_parameter_value(parfile, "parthenon/mesh", "x3min", TYPE_DBL, &x3min, 0);
  dict_add(model_params, "x3min", (snprintf(buffer, sizeof(buffer), "%.8g", x3min), buffer));
  get_parameter_value(parfile, "parthenon/mesh", "x3max", TYPE_DBL, &x3max, 0);
  dict_add(model_params, "x3max", (snprintf(buffer, sizeof(buffer), "%.8g", x3max), buffer));
  /* Set startx and dx */
  // TODO: dx may change if we are using meshblocks with refinement (SMR/AMR)
  startx[1] = x1min;
  startx[2] = x2min;
  startx[3] = x3min;
  dx[1] = (x1max - x1min) / N1;
  dx[2] = (x2max - x2min) / N2;
  dx[3] = (x3max - x3min) / N3;

  /* Which electron model to use if electrons are present in dump */
  if (!USE_FIXED_TPTE && !USE_MIXED_TPTE) {
    if (ELECTRONS != 1) {
      fprintf(stderr, "No electron temperature model specified! Cannot continue\n");
      exit(-3);
    }
    ELECTRONS = 1;
    Thetae_unit = MP/ME;
  } else if (USE_FIXED_TPTE && !USE_MIXED_TPTE) {
    ELECTRONS = 0; // force TP_OVER_TE to overwrite bad electrons
    fprintf(stderr, "Using fixed tp_over_te ratio = %g\n", tp_over_te);
    /* Thetae_unit = MP/ME*(gam-1.)*1./(1. + tp_over_te);
    * see, e.g., Eq. 8 of the EHT GRRT formula list. 
    * this formula assumes game = 4./3 and gamp = 5./3
    */
    Thetae_unit = 2./3. * MP/ME / (2. + tp_over_te);
  } else if (USE_MIXED_TPTE && !USE_FIXED_TPTE) {
    ELECTRONS = 2;
    fprintf(stderr, "Using mixed tp_over_te with trat_small = %g, trat_large = %g, and beta_crit = %g\n", 
      trat_small, trat_large, beta_crit);
    /* Thetae_unit set per-zone below */
  } else {
    fprintf(stderr, "Unknown electron model %d! Cannot continue.\n", ELECTRONS);
    exit(-3);
  }

  /* Set electron temperature unit */
  Te_unit = Thetae_unit;
  
  /* Print sigma cut value */
  fprintf(stderr, "sigma_cut = %g\n", sigma_cut);

  /* Print coordinate system used */
  switch (metric) {
    case METRIC_FMKS:
      fprintf(stderr, "Using Funky Modified Kerr-Schild coordinates FMKS\n");
      break;
    case METRIC_MKS:
      fprintf(stderr, "Using odified Kerr-Schild coordinates MKS\n");
      break;
    case METRIC_EKS:
      fprintf(stderr, "Using Kerr-Schild coordinates with exponential radial coordiante\n");
      break;
  }

  /* Load metric-specific parameters */
  if (metric == METRIC_EKS) {
    get_parameter_value(parfile, "coordinates", "a", TYPE_DBL, &a, 0);
    dict_add(model_params, "a", (snprintf(buffer, sizeof(buffer), "%.8g", a), buffer));
    get_parameter_value(parfile, "coordinates", "r_in", TYPE_DBL, &Rin, 0);
    dict_add(model_params, "r_in", (snprintf(buffer, sizeof(buffer), "%.8g", Rin), buffer));
    get_parameter_value(parfile, "coordinates", "r_out", TYPE_DBL, &Rout, 0);
    dict_add(model_params, "r_out", (snprintf(buffer, sizeof(buffer), "%.8g", Rout), buffer));
    fprintf(stderr, "eKS parameters a: %f r_in: %f r_out: %f\n", a, Rin, Rout);
  } else if (metric == METRIC_MKS) {
    get_parameter_value(parfile, "coordinates", "a", TYPE_DBL, &a, 0);
    dict_add(model_params, "a", (snprintf(buffer, sizeof(buffer), "%.8g", a), buffer));
    get_parameter_value(parfile, "coordinates", "r_in", TYPE_DBL, &Rin, 0);
    dict_add(model_params, "r_in", (snprintf(buffer, sizeof(buffer), "%.8g", Rin), buffer));
    get_parameter_value(parfile, "coordinates", "r_out", TYPE_DBL, &Rout, 0);
    dict_add(model_params, "r_out", (snprintf(buffer, sizeof(buffer), "%.8g", Rout), buffer));
    get_parameter_value(parfile, "coordinates", "hslope", TYPE_DBL, &hslope, 0);
    dict_add(model_params, "hslope", (snprintf(buffer, sizeof(buffer), "%.8g", hslope), buffer));
    fprintf(stderr, "MKS parameters a: %f hslope: %f r_in: %f r_out: %f\n", a, hslope, Rin, Rout);
  } else if (metric == METRIC_FMKS) {
    get_parameter_value(parfile, "coordinates", "a", TYPE_DBL, &a, 0);
    dict_add(model_params, "a", (snprintf(buffer, sizeof(buffer), "%.8g", a), buffer));
    get_parameter_value(parfile, "coordinates", "r_in", TYPE_DBL, &Rin, 0);
    dict_add(model_params, "r_in", (snprintf(buffer, sizeof(buffer), "%.8g", Rin), buffer));
    get_parameter_value(parfile, "coordinates", "r_out", TYPE_DBL, &Rout, 0);
    dict_add(model_params, "r_out", (snprintf(buffer, sizeof(buffer), "%.8g", Rout), buffer));
    get_parameter_value(parfile, "coordinates", "hslope", TYPE_DBL, &hslope, 0);
    dict_add(model_params, "hslope", (snprintf(buffer, sizeof(buffer), "%.8g", hslope), buffer));
    get_parameter_value(parfile, "coordinates", "mks_smooth", TYPE_DBL, &mks_smooth, 0);
    dict_add(model_params, "mks_smooth", (snprintf(buffer, sizeof(buffer), "%.8g", mks_smooth), buffer));
    get_parameter_value(parfile, "coordinates", "poly_xt", TYPE_DBL, &poly_xt, 0);
    dict_add(model_params, "poly_xt", (snprintf(buffer, sizeof(buffer), "%.8g", poly_xt), buffer));
    get_parameter_value(parfile, "coordinates", "poly_alpha", TYPE_DBL, &poly_alpha, 0);
    dict_add(model_params, "poly_alpha", (snprintf(buffer, sizeof(buffer), "%.8g", poly_alpha), buffer));
    poly_norm = 0.5 * M_PI * 1. / (1. + 1. / (poly_alpha + 1.) * 1. / pow(poly_xt, poly_alpha));
    fprintf(stderr, "FMKS parameters a: %f hslope: %f Rin: %f Rout: %f mks_smooth: %f poly_xt: %f poly_alpha: %f poly_norm: %f\n", 
      a, hslope, Rin, Rout, mks_smooth, poly_xt, poly_alpha, poly_norm);
  }

  /* Don't emit beyond specified limit or coordinate limit */
  rmax_geo = fmin(rmax_geo, Rout);
  rmin_geo = fmax(rmin_geo, Rin);

  /* The rest of the grid bounds */
  stopx[0] = 1.;
  stopx[1] = startx[1] + (N1 * dx[1]);
  stopx[2] = startx[2] + (N2 * dx[2]);
  stopx[3] = startx[3] + (N3 * dx[3]);
  /* Start & stop for specifically the *coordinate system*, for step sizes & so on */
  cstartx[0] = 0;
  cstartx[1] = 0;
  cstartx[2] = 0;
  cstartx[3] = 0;
  cstopx[0]  = 0;
  cstopx[1]  = log(Rout);
  cstopx[3]  = 2 * M_PI;

  /* Print grid bounds */
  fprintf(stderr, "Native coordinate start: %g %g %g stop: %g %g %g\n",
    cstartx[1], cstartx[2], cstartx[3], cstopx[1], cstopx[2], cstopx[3]);
  fprintf(stderr, "Grid start: %g %g %g stop: %g %g %g\n",
    startx[1], startx[2], startx[3], stopx[1], stopx[2], stopx[3]);

  /* Allocate memory for the primitive member of the data struct */
  init_storage_prims();

  return 0;
}


/**
 * @brief Set units
 *
 * Once the relevant parameres are read---MBH and Munit---this function sets various units
 * which fixes the scale of the system.
 */
void set_units()
{
  MBH = MBH_solar * MSUN; // Convert to CGS

  L_unit   = GNEWT * MBH / (CL * CL);
  T_unit   = L_unit / CL;
  RHO_unit = M_unit / pow(L_unit, 3);
  U_unit   = RHO_unit * CL * CL;
  B_unit   = CL * sqrt(4.*M_PI*RHO_unit);

  Mdotedd  = 4. * M_PI * GNEWT * MBH * MP / CL / 0.1 / SIGMA_THOMSON;

  fprintf(stderr,"MBH: %g [Msun]\n",MBH/MSUN);
  fprintf(stderr,"L,T,M units: %g [cm] %g [s] %g [g]\n",L_unit,T_unit,M_unit) ;
  fprintf(stderr,"rho,u,B units: %g [g cm^-3] %g [g cm^-1 s^-2] %g [G] \n",RHO_unit,U_unit,B_unit) ;
}


/**
 * @brief Populate ghost zones based on boundary conditions
 *
 * This function sets primitives and magnetic field strength in the ghost zones based on the boundary conditions.
 * Extend along radial direction, reflect about poles, and periodic in azimuthal direction.
 * 
 * @param n       Index of the 'data' struct.
 */
void populate_boundary_conditions(int n)
{
  /* Radial -- just extend zones */
#pragma omp parallel for collapse(2)
  for (int j=1; j<N2+1; ++j) {
    for (int k=1; k<N3+1; ++k) {
      for (int l=0; l<NVAR; ++l) {
        data[n]->p[l][0][j][k] = data[n]->p[l][1][j][k];
        data[n]->p[l][N1+1][j][k] = data[n]->p[l][N1][j][k];
      }
      data[n]->b[0][j][k] = data[n]->b[1][j][k];
      data[n]->b[N1+1][j][k] = data[n]->b[N1][j][k];
    }
  }

  /* Elevation -- flip (this is a rotation by pi) */
#pragma omp parallel for collapse(2)
  for (int i=0; i<N1+2; ++i) {
    for (int k=1; k<N3+1; ++k) {
      if (N3%2 == 0) {
        int kflip = ( (k - 1) + (N3/2) ) % N3 + 1;
        for (int l=0; l<NVAR; ++l) {
          data[n]->p[l][i][0][k] = data[n]->p[l][i][1][kflip];
          data[n]->p[l][i][N2+1][k] = data[n]->p[l][i][N2][kflip];
        }
        data[n]->b[i][0][k] = data[n]->b[i][1][kflip];
        data[n]->b[i][N2+1][k] = data[n]->b[i][N2][kflip];
      } else {
        int kflip1 = ( k + (N3/2) ) % N3;
        int kflip2 = ( k + (N3/2) + 1 ) % N3;
        for (int l=0; l<NVAR; ++l) {
          data[n]->p[l][i][0][k]    = ( data[n]->p[l][i][1][kflip1]
                                      + data[n]->p[l][i][1][kflip2] ) / 2.;
          data[n]->p[l][i][N2+1][k] = ( data[n]->p[l][i][N2][kflip1]
                                      + data[n]->p[l][i][N2][kflip2] ) / 2.;
        }
        data[n]->b[i][0][k]    = ( data[n]->b[i][1][kflip1]
                                 + data[n]->b[i][1][kflip2] ) / 2.;
        data[n]->b[i][N2+1][k] = ( data[n]->b[i][N2][kflip1]
                                 + data[n]->b[i][N2][kflip2] ) / 2.;
      }
    }
  }

  /* Azimuth -- periodic */
#pragma omp parallel for collapse(2)
  for (int i=0; i<N1+2; ++i) {
    for (int j=0; j<N2+2; ++j) {
      for (int l=0; l<NVAR; ++l) {
        data[n]->p[l][i][j][0] = data[n]->p[l][i][j][N3];
        data[n]->p[l][i][j][N3+1] = data[n]->p[l][i][j][1];
      }
      data[n]->b[i][j][0] = data[n]->b[i][j][N3];
      data[n]->b[i][j][N3+1] = data[n]->b[i][j][1];
    }
  }
}


/**
 * @brief Compute horizon accretion rate at zone @p i index along X1
 *
 * This function computes the horizon accretion rate at zone @p i index along X1 by evaluating,
 * rho * u^r * sqrt(-g) * dx^2 * dx^3, integrated over the shell at r[i]
 * 
 * @param i Index of radial coordinate
 * @param n Index of the 'data' struct.
 */
double get_code_dMact(int i, int n)
{
  i += 1;
  double dMact = 0;
#pragma omp parallel for collapse(1) reduction(+:dMact)
  for (int j = 1; j < N2+1; ++j) {
    double X[NDIM] = { 0. };
    double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
    double g, r, th;

    // this assumes axisymmetry in the coordinates
    ijktoX(i-1, j-1, 0, X);
    gcov_func(X, gcov);
    gcon_func(gcov, gcon);
    g = gdet_zone(i-1, j-1, 0);
    bl_coord(X, &r, &th);

    for (int k = 1; k < N3+1; ++k) {
      ijktoX(i-1, j-1, k, X);
      double UdotU = 0;

      for(int l = 1; l < NDIM; l++)
        for(int m = 1; m < NDIM; m++)
          UdotU += gcov[l][m]*data[n]->p[U1+l-1][i][j][k]*data[n]->p[U1+m-1][i][j][k];
      double ufac = sqrt(-1./gcon[0][0]*(1 + fabs(UdotU)));

      double ucon[NDIM] = {0.};
      ucon[0] = -ufac * gcon[0][0];

      for(int l = 1; l < NDIM; l++)
        ucon[l] = data[n]->p[U1+l-1][i][j][k] - ufac*gcon[0][l];

      double ucov[NDIM] = { 0. };
      flip_index(ucon, gcov, ucov);

      dMact += g * data[n]->p[KRHO][i][j][k] * ucon[1];
    }

  }
  return dMact;
}


/**
 * @brief Initialize derived variables in 'data' struct
 * 
 * This function initializes the derived variables in the 'data' struct, such as electron number density,
 * electron temperature, magnetic field, magnetization, and plasma beta.
 * 
 * @param n              Index of the 'data' struct.
 * @param rescale_factor Rescale factor for magnetic field.
 */
void init_physical_quantities(int n, double rescale_factor)
{
#if DEBUG
  int ceilings = 0;
#endif

  rescale_factor = sqrt(rescale_factor);

  /* Cover everything, even ghost zones */
#pragma omp parallel for collapse(3)
  for (int i = 0; i < N1+2; i++) {
    for (int j = 0; j < N2+2; j++) {
      for (int k = 0; k < N3+2; k++) {
        data[n]->ne[i][j][k] = data[n]->p[KRHO][i][j][k] * RHO_unit/(MP+ME) * Ne_factor;

        data[n]->b[i][j][k] *= rescale_factor;

        double bsq = data[n]->b[i][j][k] / B_unit;
        bsq = bsq*bsq;

        double sigma_m = bsq/data[n]->p[KRHO][i][j][k];
        double beta_m = data[n]->p[UU][i][j][k]*(gam-1.)/0.5/bsq;
#if DEBUG
        if(isnan(sigma_m)) {
          sigma_m = 0;
          fprintf(stderr, "Setting zero sigma!\n");
        }
        if(isnan(beta_m)) {
          beta_m = INFINITY;
          fprintf(stderr, "Setting INF beta!\n");
        }
#endif
        if (ELECTRONS == 1) {
          data[n]->thetae[i][j][k] = data[n]->p[KEL][i][j][k] * 
                                     pow(data[n]->p[KRHO][i][j][k],game-1.)*Thetae_unit;
        } else if (ELECTRONS == 2) {
          double betasq = (beta_m * beta_m) / (beta_crit * beta_crit);
          double trat = trat_large * betasq / (1. + betasq) + trat_small / (1. + betasq);
          //Thetae_unit = (gam - 1.) * (MP / ME) / trat;
          // see, e.g., Eq. 8 of the EHT GRRT formula list
          double lcl_Thetae_u = (MP/ME) * (game-1.) * (gamp-1.) / ( (gamp-1.) + (game-1.) * trat );
          Thetae_unit = lcl_Thetae_u;
          data[n]->thetae[i][j][k] = lcl_Thetae_u * data[n]->p[UU][i][j][k]/data[n]->p[KRHO][i][j][k];
        } else if (ELECTRONS == 9) {
          /* Convert Kelvin -> Thetae */
          data[n]->thetae[i][j][k] = data[n]->p[TFLK][i][j][k] * KBOL / ME / CL / CL;
        } else if (ELECTRONS == ELECTRONS_TFLUID) {
          double beta = data[n]->p[UU][i][j][k] * (gam-1.) / 0.5 / bsq;
          double betasq = (beta * beta) / (beta_crit * beta_crit);
          double trat = trat_large * betasq/(1. + betasq) + trat_small /(1. + betasq);
          double dfactor = mu_tot / mu_e + mu_tot / mu_i * trat;
          data[n]->thetae[i][j][k] = data[n]->p[THF][i][j][k] / dfactor;
        } else {
          data[n]->thetae[i][j][k] = Thetae_unit*data[n]->p[UU][i][j][k]/data[n]->p[KRHO][i][j][k];
        }

        // Apply floor last in case the above is a very restrictive ceiling
        data[n]->thetae[i][j][k] = fmax(data[n]->thetae[i][j][k], 1.e-3);

        // Preserve sigma for cutting along geodesics, and for variable-kappa model
        data[n]->sigma[i][j][k] = fmax(sigma_m, SMALL);
        // Also record beta, for variable-kappa model
        data[n]->beta[i][j][k] = fmax(beta_m, SMALL);

        // Cut Ne (i.e. emission) based on sigma, if we're not doing so along each geodesic
        // Strongly magnetized = empty, no shiny spine
        if (sigma_m > sigma_cut && !USE_GEODESIC_SIGMACUT) {
          data[n]->b[i][j][k]     = 0.0;
          data[n]->ne[i][j][k]    = 0.0;
          data[n]->thetae[i][j][k]= 0.0;
        }
      }
    }
  }
}


/**
 * @brief Read fluid primitives, compute four-vectors, populate ghost zones, and initialize the 'data' struct.
 *
 * This function reads the primitives from the native KHARMA dump file, assembles the mesh from the meshblocks,
 * calculates the four-vectors in the physical zones, and initializes the derived quantities.
 * 
 * @param n       Index of the 'data' struct.
 * @param fname   Filename of the dump file.
 * @param data    Pointer to the 'data' struct.
 * @param verbose Verbosity level (print values)
 */
void load_kharma_data(int n, char *fnam, int dumpidx, int verbose)
{
  /* File name */
  char fname[256];
  snprintf(fname, 255, fnam, dumpidx);

  nloaded++;

  double dMact, Ladv;

  /* Check file exists and can be accessed */
  if ( hdf5_open(fname) < 0 ) {
    fprintf(stderr, "Unable to open file %s. Exiting!\n", fname);
    exit(-1);
  }

  /* Get number of meshblocks and meshblock size */
  int num_meshblocks = 1;
  int meshblock_size[NDIM-1] = {0};
  read_info_attribute(fname, "NumMeshBlocks", TYPE_INT, &num_meshblocks, 0);
  read_info_attribute(fname, "MeshBlockSize", TYPE_INT_ARRAY, meshblock_size, 0);
  int nx1_mb = meshblock_size[0];
  int nx2_mb = meshblock_size[1];
  int nx3_mb = meshblock_size[2];
  
  /* Get meshblock ordering */
  int **mb_order = malloc_rank2_int(num_meshblocks, NDIM-1);
  int rank = 2;
  hsize_t fdims_2[] = {num_meshblocks, NDIM-1};
  hsize_t fstart_2[] = {0, 0};
  hsize_t fcount_2[] = {num_meshblocks, NDIM-1};
  hsize_t mdims_2[] = {num_meshblocks, NDIM-1};
  hsize_t mstart_2[] = {0, 0};
  hdf5_read_array(mb_order[0], "Blocks/loc.lx123", rank, fdims_2, fstart_2, fcount_2, mdims_2, mstart_2, H5T_STD_I32LE);

  /* Read primitives into buffer */
  double *****primitives_buffer; // NMB, nx3_mb, nx2_mb, nx1_mb, NVAR
  /* Allocate memory */
  primitives_buffer = malloc_rank5(num_meshblocks, nx3_mb, nx2_mb, nx1_mb, NVAR);

  if (strcmp(kharma_format, "new") == 0) {
    fprintf(stderr, "Using latest KHARMA output format\n");
    /* Read scalar fields */
    int frank = 4;
    int mrank = 5;
    hsize_t fdims_4[4]  = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb};
    hsize_t fstart_4[4] = {0, 0, 0, 0};
    hsize_t fcount_4[4] = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb}; // Read the entire dataset
    /* In the memory buffer, we want to store the file data into the slice corresponding to "rho".
      So we set mstart such that the last (5th) dimension starts at KRHO, and mcount to read 1 element along that axis. */
    hsize_t mdims_5[5]  = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb, NVAR};
    hsize_t mstart_5[5] = {0, 0, 0, 0, KRHO};
    hsize_t mcount_5[5] = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb, 1};
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.rho", frank, fdims_4, fstart_4, fcount_4, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    mstart_5[4] = UU;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.u", frank, fdims_4, fstart_4, fcount_4, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    /* Read vector fields */
    frank = 5;
    hsize_t fdims_5[5]  = {num_meshblocks, NDIM-1, nx3_mb, nx2_mb, nx1_mb};
    hsize_t fstart_5[5] = {0, 0, 0, 0, 0};
    hsize_t fcount_5[5] = {num_meshblocks, 1, nx3_mb, nx2_mb, nx1_mb};
    mstart_5[4] = U1;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.uvec", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    fstart_5[1] = 1;
    mstart_5[4] = U2;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.uvec", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    fstart_5[1] = 2;
    mstart_5[4] = U3;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.uvec", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    fstart_5[1] = 0;
    mstart_5[4] = B1;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.B", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    fstart_5[1] = 1;
    mstart_5[4] = B2;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.B", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    fstart_5[1] = 2;
    mstart_5[4] = B3;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.B", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
  } else if (strcmp(kharma_format, "old") == 0) {
    fprintf(stderr, "Using old (v5) KHARMA output format\n");
    /* Read scalar fields */
    int frank = 5;
    int mrank = 5;
    hsize_t fdims_4[5]  = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb, 1};
    hsize_t fstart_4[5] = {0, 0, 0, 0, 0};
    hsize_t fcount_4[5] = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb, 1}; // Read the entire dataset
    /* In the memory buffer, we want to store the file data into the slice corresponding to "rho".
      So we set mstart such that the last (5th) dimension starts at KRHO, and mcount to read 1 element along that axis. */
    hsize_t mdims_5[5]  = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb, NVAR};
    hsize_t mstart_5[5] = {0, 0, 0, 0, KRHO};
    hsize_t mcount_5[5] = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb, 1};
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.rho", frank, fdims_4, fstart_4, fcount_4, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    mstart_5[4] = UU;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.u", frank, fdims_4, fstart_4, fcount_4, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    /* Read vector fields */
    frank = 5;
    hsize_t fdims_5[5]  = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb, NDIM-1};
    hsize_t fstart_5[5] = {0, 0, 0, 0, 0};
    hsize_t fcount_5[5] = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb, 1};
    mstart_5[4] = U1;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.uvec", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    fstart_5[4] = 1;
    mstart_5[4] = U2;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.uvec", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    fstart_5[4] = 2;
    mstart_5[4] = U3;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.uvec", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    fstart_5[4] = 0;
    mstart_5[4] = B1;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.B", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    fstart_5[4] = 1;
    mstart_5[4] = B2;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.B", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
    fstart_5[4] = 2;
    mstart_5[4] = B3;
    hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.B", frank, fdims_5, fstart_5, fcount_5, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
  }
  if (ELECTRONS == 1) {
    /* Read electron entropy */
    int frank = 4;
    int mrank = 5;
    hsize_t fdims_4[4]  = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb};
    hsize_t fstart_4[4] = {0, 0, 0, 0};
    hsize_t fcount_4[4] = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb}; // Read the entire dataset
    hsize_t mdims_5[5]  = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb, NVAR};
    hsize_t mstart_5[5] = {0, 0, 0, 0, KRHO};
    hsize_t mcount_5[5] = {num_meshblocks, nx3_mb, nx2_mb, nx1_mb, 1};
    mstart_5[4] = KEL;
    /* Get electron subgrid model enum */
    enum ElectronSubgridModel electron_subgrid_model = get_electron_subgrid_model(electron_subgrid_model_string);
    switch (electron_subgrid_model) {
      case ELECTRONS_CONSTANT:
        fprintf(stderr, "Using CONSTANT electron temperature model\n");
        hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.Kel_Constant", frank, fdims_4, fstart_4, fcount_4, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
        break;
      case ELECTRONS_HOWES:
        fprintf(stderr, "Using HOWES electron temperature model\n");
        hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.Kel_Howes", frank, fdims_4, fstart_4, fcount_4, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
        break;
      case ELECTRONS_KAWAZURA:
        fprintf(stderr, "Using KAWAZURA electron temperature model\n");
        hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.Kel_Kawazura", frank, fdims_4, fstart_4, fcount_4, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
        break;
      case ELECTRONS_WERNER:
        fprintf(stderr, "Using WERNER electron temperature model\n");
        hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.Kel_Werner", frank, fdims_4, fstart_4, fcount_4, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
        break;
      case ELECTRONS_ROWAN:
        fprintf(stderr, "Using ROWAN electron temperature model\n");
        hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.Kel_Rowan", frank, fdims_4, fstart_4, fcount_4, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
        break;
      case ELECTRONS_SHARMA:
        fprintf(stderr, "Using SHARMA electron temperature model\n");
        hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.Kel_Sharma", frank, fdims_4, fstart_4, fcount_4, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
        break;
      case ELECTRONS_UNKNOWN:
        fprintf(stderr, "Unknown electron subgrid model: %s\n", electron_subgrid_model_string);
        break;
    }

    /* Read fluid entropy */
    // mstart_5[4] = KTOT;
    // hdf5_read_array_multidim(primitives_buffer[0][0][0][0], "prims.Ktot", frank, fdims_4, fstart_4, fcount_4, mrank, mdims_5, mstart_5, mcount_5, H5T_IEEE_F64LE);
  }
  
  /* Assemble mesh */
#pragma omp parallel for collapse(4)
  for (int mb = 0; mb < num_meshblocks; mb++) {
    for (int kb = 0; kb < meshblock_size[2]; kb++) {
      for (int jb = 0; jb < meshblock_size[1]; jb++) {
        for (int ib = 0; ib < meshblock_size[0]; ib++) {
          int *mbd_loc = mb_order[mb];
          data[n]->p[KRHO][1+mbd_loc[0]*meshblock_size[0]+ib]
                          [1+mbd_loc[1]*meshblock_size[1]+jb]
                          [1+mbd_loc[2]*meshblock_size[2]+kb] = primitives_buffer[mb][kb][jb][ib][KRHO];
            data[n]->p[UU][1+mbd_loc[0]*meshblock_size[0]+ib]
                          [1+mbd_loc[1]*meshblock_size[1]+jb]
                          [1+mbd_loc[2]*meshblock_size[2]+kb] = primitives_buffer[mb][kb][jb][ib][UU];
            data[n]->p[U1][1+mbd_loc[0]*meshblock_size[0]+ib]
                          [1+mbd_loc[1]*meshblock_size[1]+jb]
                          [1+mbd_loc[2]*meshblock_size[2]+kb] = primitives_buffer[mb][kb][jb][ib][U1];
            data[n]->p[U2][1+mbd_loc[0]*meshblock_size[0]+ib]
                          [1+mbd_loc[1]*meshblock_size[1]+jb]
                          [1+mbd_loc[2]*meshblock_size[2]+kb] = primitives_buffer[mb][kb][jb][ib][U2];
            data[n]->p[U3][1+mbd_loc[0]*meshblock_size[0]+ib]
                          [1+mbd_loc[1]*meshblock_size[1]+jb]
                          [1+mbd_loc[2]*meshblock_size[2]+kb] = primitives_buffer[mb][kb][jb][ib][U3];
            data[n]->p[B1][1+mbd_loc[0]*meshblock_size[0]+ib]
                          [1+mbd_loc[1]*meshblock_size[1]+jb]
                          [1+mbd_loc[2]*meshblock_size[2]+kb] = primitives_buffer[mb][kb][jb][ib][B1];
            data[n]->p[B2][1+mbd_loc[0]*meshblock_size[0]+ib]
                          [1+mbd_loc[1]*meshblock_size[1]+jb]
                          [1+mbd_loc[2]*meshblock_size[2]+kb] = primitives_buffer[mb][kb][jb][ib][B2];
            data[n]->p[B3][1+mbd_loc[0]*meshblock_size[0]+ib]
                          [1+mbd_loc[1]*meshblock_size[1]+jb]
                          [1+mbd_loc[2]*meshblock_size[2]+kb] = primitives_buffer[mb][kb][jb][ib][B3];

            if (ELECTRONS == 1) {
              data[n]->p[KEL][1+mbd_loc[0]*meshblock_size[0]+ib]
                             [1+mbd_loc[1]*meshblock_size[1]+jb]
                             [1+mbd_loc[2]*meshblock_size[2]+kb] = primitives_buffer[mb][kb][jb][ib][KEL];
              // data[n]->p[KTOT][1+mbd_loc[0]*meshblock_size[0]+ib]
              //                 [1+mbd_loc[1]*meshblock_size[1]+jb]
              //                 [1+mbd_loc[2]*meshblock_size[2]+kb] = primitives_buffer[mb][kb][jb][ib][KTOT];
            }

        }
      }
    }
  }

  /* Free memory */
  free_rank5(primitives_buffer, num_meshblocks, nx3_mb, nx2_mb, nx1_mb);

  /* Close file */
  hdf5_close();

  /* Get dump time */
  double t = -1.;
  read_info_attribute(fname, "Time", TYPE_DBL, &t, 0);
  data[n]->t = t;

  /* Reversing B Field */
  if(reverse_field) {
    double multiplier = -1.0;
    for(int i = 0; i < N1 + 2; i++){
      for(int j = 0; j < N2 + 2; j++){
        for(int k = 0; k< N3 + 2; k++){ 
          data[n]->p[B1][i][j][k] = multiplier*data[n]->p[B1][i][j][k];
          data[n]->p[B2][i][j][k] = multiplier*data[n]->p[B2][i][j][k];
          data[n]->p[B3][i][j][k] = multiplier*data[n]->p[B3][i][j][k];
        }
      }
    }
  }

  /* Allocate memory for derived quantities */
  init_storage_derived_quantities();

  dMact = Ladv = 0.;

  /* Construct four-vectors over "real" zones */
  #pragma omp parallel for collapse(2) reduction(+:dMact,Ladv)
  for(int i = 1; i < N1+1; i++) {
    for(int j = 1; j < N2+1; j++) {

      double X[NDIM] = {0.};
      double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
      double g, r, th;

      // this assumes axisymmetry in the coordinates
      ijktoX(i-1, j-1, 0, X);
      gcov_func(X, gcov);
      gcon_func(gcov, gcon);
      g = gdet_zone(i-1, j-1, 0);

      bl_coord(X, &r, &th);

      for(int k = 1; k < N3+1; k++){

        ijktoX(i-1,j-1,k,X);
        double UdotU = 0.;

        /* The four-vector reconstruction should have gcov and gcon and gdet using the modified coordinates
        interpolating the four vectors to the zone center !!!! */
        for(int l = 1; l < NDIM; l++) 
          for(int m = 1; m < NDIM; m++) 
            UdotU += gcov[l][m]*data[n]->p[U1+l-1][i][j][k]*data[n]->p[U1+m-1][i][j][k];
        double ufac = sqrt(-1./gcon[0][0]*(1 + fabs(UdotU)));

        // TODO: move this above
        double ucon[NDIM] = { 0. };
        ucon[0] = -ufac * gcon[0][0];

        for(int l = 1; l < NDIM; l++) 
          ucon[l] = data[n]->p[U1+l-1][i][j][k] - ufac*gcon[0][l];

        double ucov[NDIM] = { 0. };
        flip_index(ucon, gcov, ucov);

        /* Reconstruct the magnetic field three vectors */
        double udotB = 0.;

        for (int l = 1; l < NDIM; l++) {
          udotB += ucov[l]*data[n]->p[B1+l-1][i][j][k];
        }

        double bcon[NDIM] = {0.};
        double bcov[NDIM] = {0.};

        bcon[0] = udotB;
        for (int l = 1; l < NDIM; l++) {
          bcon[l] = (data[n]->p[B1+l-1][i][j][k] + ucon[l]*udotB)/ucon[0];
        }
        flip_index(bcon, gcov, bcov);

        double bsq = 0.;
        for (int l=0; l<NDIM; ++l) bsq += bcon[l] * bcov[l];
        data[n]->b[i][j][k] = sqrt(bsq) * B_unit;


        if(i <= 21) { dMact += g * data[n]->p[KRHO][i][j][k] * ucon[1]; }
        if(i >= 21 && i < 41 && 0) Ladv += g * data[n]->p[UU][i][j][k] * ucon[1] * ucov[0];
        if(i <= 21) Ladv += g * data[n]->p[UU][i][j][k] * ucon[1] * ucov[0];

      }
    }
  }

  /* Check if 21st zone (i.e., 20th without ghost zones) is beyond r_eh otherwise recompute */
  double r_eh = 1. + sqrt(1. - a * a);
  int N2_by_2 = (int)(N2 / 2);

  double X[NDIM] = {0.};
  ijktoX(20, N2_by_2, 0, X);

  double r, th;
  bl_coord(X, &r, &th);

  if (r < r_eh) {
    for (int i = 1; i < N1+1; i++) {
      ijktoX(i, N2_by_2, 0, X);
      bl_coord(X, &r, &th);
      if (r >= r_eh) {
        fprintf(stderr, "r_eh is beyond regular zones. recomputing at %g...\n", r);
        dMact = get_code_dMact(i, n) * 21;
        break;
      }
    }
  }

  /* Copy primitives and four-vectors to ghost zones according to boundary conditions */
  populate_boundary_conditions(n);

  dMact *= dx[3]*dx[2] ;
  dMact /= 21. ;
  Ladv *= dx[3]*dx[2] ;
  Ladv /= 21. ;

  Mdot_dump = -dMact * (M_unit/T_unit);
  MdotEdd_dump = Mdotedd;
  Ladv_dump =  Ladv;

  double rescale_factor = 1.;

  if (verbose == 2) {
    fprintf(stderr,"dMact: %g [code]\n",dMact);
    fprintf(stderr,"Ladv: %g [code]\n",Ladv_dump);
    fprintf(stderr,"Mdot: %g [g/s] \n",Mdot_dump);
    fprintf(stderr,"Mdot: %g [MSUN/YR] \n",Mdot_dump/(MSUN / YEAR));
    fprintf(stderr,"Mdot: %g [Mdotedd]\n",Mdot_dump/MdotEdd_dump);
    fprintf(stderr,"Mdotedd: %g [g/s]\n",MdotEdd_dump);
    fprintf(stderr,"Mdotedd: %g [MSUN/YR]\n",MdotEdd_dump/(MSUN/YEAR));
  } else if (verbose == 1) {
    fprintf(stderr,"Mdot: %g [MSUN/YR] \n",Mdot_dump/(MSUN / YEAR));
    fprintf(stderr,"Mdot: %g [Mdotedd]\n",Mdot_dump/MdotEdd_dump);
  }

  // now construct useful scalar quantities (over full (+ghost) zones of data)
  init_physical_quantities(n, rescale_factor);

}


/**
 * @brief Initialize KHARMA model
 *
 * This function reads the fluid data from the dump file and stores it in the 'data' struct.
 * 
 * @param tA Address of time variable for dataA (only for slowlight calculation)
 * @param tB Address of time variable for dataB (only for slowlight calculation)
 */
void init_model(double *tA, double *tB)
{
  /* Set up initial ordering of data[] */
  data[0] = &dataA;
  data[1] = &dataB;
  data[2] = &dataC;

  /* Create parameter dictionary */
  model_params = dict_new();

  /* Read relevant parameters from dump file and allocate memory for data struct */
  fprintf(stderr, "Reading parameters, allocating memory...\n");
  read_parameters_and_allocate_memory(fnam, dumpmin);

  /* Set all dimensional quantities from loaded parameters */
  set_units();

  /* Read fluid data */
  fprintf(stderr, "Reading data...\n");
  load_kharma_data(0, fnam, dumpidx, 2);  // replaced dumpmin -> 2 because apparently that argument was just .. removed
  dumpidx += dumpskip;
  #if SLOW_LIGHT
  update_data(tA, tB);
  update_data(tA, tB);
  tf = get_dump_time(fnam, dumpmax) - 1.e-5;
  #else // FAST LIGHT
  data[2]->t = 10000.;
  #endif // SLOW_LIGHT

  /* Event horizon radius */
  Rh = 1 + sqrt(1. - a * a);

  /* Possibly cut around the pole */
  if (polar_cut >= 0) th_beg = 0.0174 * polar_cut;

  #if DEBUG_READER
    /* Set filename */
    char debug_fname[256];
    snprintf(debug_fname, sizeof(debug_fname), "debug_reader_kharma.h5");

    /* Create HDF5 file*/
    hdf5_create(debug_fname);

    /* Compute gcov, gcon */
    size_t total_elements = NDIM * NDIM * (N2 + 2) * (N1 + 2);
    double *gcov_global = malloc(total_elements * sizeof(double));
    double *gcon_global = malloc(total_elements * sizeof(double));
    // Use the arrays via indexing. For example, to access element [mu][nu][j][i]:
    #define IDX(mu, nu, j, i) (((mu) * NDIM * (N2+2) * (N1+2)) + ((nu) * (N2+2) * (N1+2)) + ((j) * (N1+2)) + (i))

#pragma omp parallel for collapse(2)
  for (int i = 0; i < N1+2; i++) {
    for (int j = 0; j < N2+2; j++) {

      double X[NDIM] = {0.};
      ijktoX(i, j, 0, X);
      double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
      gcov_func(X, gcov);
      gcon_func(gcov, gcon);

      for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
          gcov_global[IDX(mu, nu, j, i)] = gcov[mu][nu];
          gcon_global[IDX(mu, nu, j, i)] = gcon[mu][nu];
        }
      }
    }
  }

    /* Write gcov, gcon to file */
    hsize_t dims_grid[4] = { NDIM, NDIM, N2+2, N1+2 };
    hdf5_write_full_array(gcov_global, "gcov", 4, dims_grid, H5T_NATIVE_DOUBLE);
    hdf5_write_full_array(gcon_global, "gcon", 4, dims_grid, H5T_NATIVE_DOUBLE);

    // Free the allocated memory when done
    free(gcov_global);
    free(gcon_global);

    /* Write physical quantities */
    hsize_t dims_phys[3] = { N1+2, N2+2, N3+2 };
    hdf5_write_full_array(data[0]->ne[0][0], "ne", 3, dims_phys, H5T_NATIVE_DOUBLE);
    hdf5_write_full_array(data[0]->thetae[0][0], "thetae", 3, dims_phys, H5T_NATIVE_DOUBLE);
    hdf5_write_full_array(data[0]->b[0][0], "b", 3, dims_phys, H5T_NATIVE_DOUBLE);
    hdf5_write_full_array(data[0]->sigma[0][0], "sigma", 3, dims_phys, H5T_NATIVE_DOUBLE);
    hdf5_write_full_array(data[0]->beta[0][0], "beta", 3, dims_phys, H5T_NATIVE_DOUBLE);

    /* Close HDF5 file */
    hdf5_close();

  #endif // DEBUG_READER
}


/** 
 * @brief How far have we interpolated between data[nA]->t and data[nB]->t
 * 
 * In slowlight mode, we perform linear interpolation in time. This function tells
 * us how far we've progressed from data[nA]->t to data[nB]->t but "in reverse" as
 * tinterp == 1 -> we're exactly on nA and tinterp == 0 -> we're exactly on nB.
 * 
 * @param X  local coordinates
 * @param nA dataA index
 * @param nB dataB index
 */
double set_tinterp_ns(double X[NDIM], int *nA, int *nB)
{
  #if SLOW_LIGHT
  if (X[0] < data[1]->t) {
    *nA = 0; *nB = 1;
  } else {
    *nA = 1; *nB = 2;
  }
  double tinterp = 1. - ( X[0] - data[*nA]->t ) / ( data[*nB]->t - data[*nA]->t );
  if (tinterp < 0.) tinterp = 0.; //  In slow light, when we reset based on tB, sometimes we overshoot
  if (tinterp > 1.) tinterp = 1.; //  TODO, this should really only happen at r >> risco, but still...
  return tinterp;
  #else
  *nA = 0;
  *nB = 0;
  return 0.;
  #endif // SLOW_LIGHT
}


// Calculate Ucon,Ucov,Bcon,Bcov from primitives at location X using 
// interpolation (on the primitives). This has been all wrapped into
// a single function because some calculations require each other.
void get_model_fourv(double X[NDIM], double Kcon[NDIM],
                     double Ucon[NDIM], double Ucov[NDIM],
                     double Bcon[NDIM], double Bcov[NDIM])
{
  double gcov[NDIM][NDIM], gcon[NDIM][NDIM];

  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  // If we're outside of the logical domain, default to
  // normal observer velocity for Ucon/Ucov and default
  // Bcon/Bcov to zero.
  if ( X_in_domain(X) == 0 ) {

    Ucov[0] = -1./sqrt(-gcon[0][0]);
    Ucov[1] = 0.;
    Ucov[2] = 0.;
    Ucov[3] = 0.;
    Ucon[0] = 0.;
    Ucon[1] = 0.;
    Ucon[2] = 0.;
    Ucon[3] = 0.;

    for (int mu=0; mu<NDIM; ++mu) {
      Ucon[0] += Ucov[mu] * gcon[0][mu];
      Ucon[1] += Ucov[mu] * gcon[1][mu];
      Ucon[2] += Ucov[mu] * gcon[2][mu];
      Ucon[3] += Ucov[mu] * gcon[3][mu];
      Bcon[mu] = 0.;
      Bcov[mu] = 0.;
    }

    return;
  }

  // Set Ucon and get Ucov by lowering

  // interpolate primitive variables first
  double Vcon[NDIM];
  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);
  Vcon[1] = interp_scalar_time(X, data[nA]->p[U1], data[nB]->p[U1], tfac);
  Vcon[2] = interp_scalar_time(X, data[nA]->p[U2], data[nB]->p[U2], tfac);
  Vcon[3] = interp_scalar_time(X, data[nA]->p[U3], data[nB]->p[U3], tfac);

  // translate to four velocity
  double VdotV = 0.;
  for (int i = 1; i < NDIM; i++)
    for (int j = 1; j < NDIM; j++)
      VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
  double Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
  Ucon[0] = -Vfac * gcon[0][0];
  for (int i = 1; i < NDIM; i++)
    Ucon[i] = Vcon[i] - Vfac * gcon[0][i];

  // lower (needed for Bcon)
  lower(Ucon, gcov, Ucov);

  // Now set Bcon and get Bcov by lowering

  // interpolate primitive variables first
  double Bcon1 = interp_scalar_time(X, data[nA]->p[B1], data[nB]->p[B1], tfac);
  double Bcon2 = interp_scalar_time(X, data[nA]->p[B2], data[nB]->p[B2], tfac);
  double Bcon3 = interp_scalar_time(X, data[nA]->p[B3], data[nB]->p[B3], tfac);

  // get Bcon
  Bcon[0] = Bcon1*Ucov[1] + Bcon2*Ucov[2] + Bcon3*Ucov[3];
  Bcon[1] = (Bcon1 + Ucon[1] * Bcon[0]) / Ucon[0];
  Bcon[2] = (Bcon2 + Ucon[2] * Bcon[0]) / Ucon[0];
  Bcon[3] = (Bcon3 + Ucon[3] * Bcon[0]) / Ucon[0];

  // lower
  lower(Bcon, gcov, Bcov);
}


// Get the primitive variables interpolated to a point X,
// And fill them in the next 8 array slots after p
// Not used for transport but useful for plotting along a geodesic later
void get_model_primitives(double X[NDIM], double *p)
{
  if ( X_in_domain(X) == 0 ) return;

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  for (int np=0; np<8; np++) {
    p[np] = interp_scalar_time(X, data[nA]->p[np], data[nA]->p[np], tfac);
  }
}


double get_model_thetae(double X[NDIM])
{
  if ( X_in_domain(X) == 0 ) return 0.;
  
  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);
  double thetae = interp_scalar_time(X, data[nA]->thetae, data[nB]->thetae, tfac);

#if DEBUG
  if (thetae < 0. || isnan(thetae)) {
    printf("thetae negative or NaN!\n");
    printf("X[] = %g %g %g %g\n", X[0], X[1], X[2], X[3]);
    printf("t = %e %e %e\n", data[0]->t, data[1]->t, data[2]->t);
    double thetaeA = interp_scalar(X, data[nA]->thetae);
    double thetaeB = interp_scalar(X, data[nB]->thetae);
    printf("thetaeA, thetaeB = %e %e", thetaeA, thetaeB);
    printf("thetae, tfac = %e %e\n", thetae, tfac);
  }
#endif

  return thetae;
}


//b field strength in Gauss
double get_model_b(double X[NDIM])
{
  // TODO how *should* we handle exiting the domain consistently?
  if ( X_in_domain(X) == 0 ) return 0.;
  
  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  return interp_scalar_time(X, data[nA]->b, data[nB]->b, tfac);
}


double get_model_sigma(double X[NDIM]) 
{
  if ( X_in_domain(X) == 0 ) return 0.;

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  return interp_scalar_time(X, data[nA]->sigma, data[nB]->sigma, tfac);
}


double get_model_beta(double X[NDIM]) 
{
  if ( X_in_domain(X) == 0 ) return 0.;  // TODO inf?

  double betaA, betaB, tfac;
  int nA, nB;
  tfac = set_tinterp_ns(X, &nA, &nB);
  betaA = interp_scalar(X, data[nA]->beta);
  betaB = interp_scalar(X, data[nB]->beta);
  return tfac*betaA + (1. - tfac)*betaB;
}


double get_sigma_smoothfac(double sigma)
{
  double sigma_above = sigma_cut;
  if (sigma_cut_high > 0) sigma_above = sigma_cut_high;
  if (sigma < sigma_cut) return 1;
  if (sigma >= sigma_above) return 0;
  double dsig = sigma_above - sigma_cut;
  return cos(M_PI / 2. / dsig * (sigma - sigma_cut));
}


double get_model_ne(double X[NDIM])
{
  if ( X_in_domain(X) == 0 ) return 0.;

  double sigma_smoothfac = 1;

#if USE_GEODESIC_SIGMACUT
  double sigma = get_model_sigma(X);
  if (sigma > sigma_cut) return 0.;
  sigma_smoothfac = get_sigma_smoothfac(sigma);
#endif

  int nA, nB;
  double tfac = set_tinterp_ns(X, &nA, &nB);

  return interp_scalar_time(X, data[nA]->ne, data[nB]->ne, tfac) * sigma_smoothfac;
}


void output_hdf5()
{
  hdf5_set_directory("/");

  hdf5_write_single_val(&Mdot_dump, "Mdot", H5T_IEEE_F64LE);
  hdf5_write_single_val(&MdotEdd_dump, "MdotEdd", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Ladv_dump, "Ladv", H5T_IEEE_F64LE);

  // Write fluid header
  // Obtain values from the dictionary `model_params`
  hdf5_make_directory("fluid_header");
  hdf5_set_directory("/fluid_header/");
  // Get root level parameters
  double a = atof(dict_get(model_params, "a", NULL));
  const char *bfield_type = dict_get(model_params, "bfield_type", NULL);
  char bfield_type_str[STRLEN];
  if (bfield_type) {
      strncpy(bfield_type_str, bfield_type, STRLEN - 1);
      bfield_type_str[STRLEN - 1] = '\0';
  } else {
      bfield_type_str[0] = '\0';
  }
  const char *coord_base = dict_get(model_params, "base", NULL);
  char base[STRLEN], metric[STRLEN];
  if (coord_base) {
      strncpy(base, coord_base, STRLEN - 1);
      base[STRLEN - 1] = '\0';
  } else {
      base[0] = '\0';
  }
  double cfl = atof(dict_get(model_params, "cfl", NULL));
  double cour = atof(dict_get(model_params, "cfl", NULL));
  double dt = atof(dict_get(model_params, "dt", NULL));
  double dt_min = atof(dict_get(model_params, "dt_min", NULL));
  // double dx1 = atof(dict_get(model_params, "dx1", NULL));
  // double dx2 = atof(dict_get(model_params, "dx2", NULL));
  // double dx3 = atof(dict_get(model_params, "dx3", NULL));
  double gam = atof(dict_get(model_params, "gam", NULL));
  double gamma = atof(dict_get(model_params, "gam", NULL));
  const char *electrons_on = dict_get(model_params, "has_electrons", NULL);
  char has_electrons[STRLEN];
  if (electrons_on) {
      strncpy(has_electrons, electrons_on, STRLEN - 1);
      has_electrons[STRLEN - 1] = '\0';
  } else {
      has_electrons[0] = '\0';
  }
  if (strncmp(has_electrons, "true", 19) == 0) {
    double game = atof(dict_get(model_params, "game", NULL));
    double gamp = atof(dict_get(model_params, "gamp", NULL));
  }
  double hslope = atof(dict_get(model_params, "hslope", NULL));

  const char *coord_system = dict_get(model_params, "coordinate_system", NULL);
  if (coord_system) {
      strncpy(metric, coord_system, STRLEN - 1);
      metric[STRLEN - 1] = '\0';
  } else {
      metric[0] = '\0';
  }
  double mks_smooth = atof(dict_get(model_params, "mks_smooth", NULL));
  int ncycle = atoi(dict_get(model_params, "ncycle", NULL));
  int nghost = atoi(dict_get(model_params, "nghost", NULL));
  int nx1 = atoi(dict_get(model_params, "nx1", NULL));
  int nx2 = atoi(dict_get(model_params, "nx2", NULL));
  int nx3 = atoi(dict_get(model_params, "nx3", NULL));
  int nx1_mb = atoi(dict_get(model_params, "nx1_mb", NULL));
  int nx2_mb = atoi(dict_get(model_params, "nx2_mb", NULL));
  int nx3_mb = atoi(dict_get(model_params, "nx3_mb", NULL));
  double poly_alpha = atof(dict_get(model_params, "poly_alpha", NULL));
  double poly_xt = atof(dict_get(model_params, "poly_xt", NULL));
  const char *problem_id = dict_get(model_params, "problem_id", NULL);
  char problem[STRLEN];
  if (problem_id) {
      strncpy(problem, problem_id, STRLEN - 1);
      problem[STRLEN - 1] = '\0';
  } else {
      problem[0] = '\0';
  }
  double r_eh = 1. + sqrt(1. - a * a);
  double r_in = atof(dict_get(model_params, "r_in", NULL));
  double r_out = atof(dict_get(model_params, "r_out", NULL));
  double rhor = r_eh;
  double rin = atof(dict_get(model_params, "rin", NULL));
  double rmax = atof(dict_get(model_params, "rmax", NULL));
  double startx1 = atof(dict_get(model_params, "x1min", NULL));
  double startx2 = atof(dict_get(model_params, "x2min", NULL));
  double startx3 = atof(dict_get(model_params, "x3min", NULL));
  double time = atof(dict_get(model_params, "time", NULL));
  double tlim = atof(dict_get(model_params, "tlim", NULL));
  double x1min = atof(dict_get(model_params, "x1min", NULL));
  double x1max = atof(dict_get(model_params, "x1max", NULL));
  double x2min = atof(dict_get(model_params, "x2min", NULL));
  double x2max = atof(dict_get(model_params, "x2max", NULL));
  double x3min = atof(dict_get(model_params, "x3min", NULL));
  double x3max = atof(dict_get(model_params, "x3max", NULL));
  // Write base fluid header
  hdf5_write_single_val(&a, "a", H5T_IEEE_F64LE);

  hid_t str_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(str_type, strlen(base) + 1);
  H5Tset_strpad(str_type, H5T_STR_NULLTERM);
  hdf5_write_single_val(bfield_type_str, "bfield_type", str_type);
  hdf5_write_single_val(base, "base", str_type);

  hdf5_write_single_val(&cfl, "cfl", H5T_IEEE_F64LE);
  hdf5_write_single_val(&cour, "cour", H5T_IEEE_F64LE);
  hdf5_write_single_val(&dt, "dt", H5T_IEEE_F64LE);
  hdf5_write_single_val(&dt_min, "dt_min", H5T_IEEE_F64LE);
  // hdf5_write_single_val(&dx1, "dx1", H5T_IEEE_F64LE);
  // hdf5_write_single_val(&dx2, "dx2", H5T_IEEE_F64LE);
  // hdf5_write_single_val(&dx3, "dx3", H5T_IEEE_F64LE);
  hdf5_write_single_val(&gam, "gam", H5T_IEEE_F64LE);
  hdf5_write_single_val(&game, "game", H5T_IEEE_F64LE);
  hdf5_write_single_val(&gamp, "gamp", H5T_IEEE_F64LE);
  hdf5_write_single_val(&gamma, "gamma", H5T_IEEE_F64LE);
  hdf5_write_single_val(&hslope, "hslope", H5T_IEEE_F64LE);
  hdf5_write_single_val(&mks_smooth, "mks_smooth", H5T_IEEE_F64LE);
  hdf5_write_single_val(&ncycle, "ncycle", H5T_STD_I32LE);
  hdf5_write_single_val(&nghost, "nghost", H5T_STD_I32LE);
  hdf5_write_single_val(&nx1, "nx1", H5T_STD_I32LE);
  hdf5_write_single_val(&nx2, "nx2", H5T_STD_I32LE);
  hdf5_write_single_val(&nx3, "nx3", H5T_STD_I32LE);
  hdf5_write_single_val(&nx1_mb, "nx1_mb", H5T_STD_I32LE);
  hdf5_write_single_val(&nx2_mb, "nx2_mb", H5T_STD_I32LE);
  hdf5_write_single_val(&nx3_mb, "nx3_mb", H5T_STD_I32LE);
  hdf5_write_single_val(&poly_alpha, "poly_alpha", H5T_IEEE_F64LE);
  hdf5_write_single_val(&poly_xt, "poly_xt", H5T_IEEE_F64LE);

  H5Tset_size(str_type, strlen(problem) + 1);
  H5Tset_strpad(str_type, H5T_STR_NULLTERM);
  hdf5_write_single_val(problem, "problem", str_type);

  hdf5_write_single_val(&r_eh, "r_eh", H5T_IEEE_F64LE);
  hdf5_write_single_val(&r_in, "r_in", H5T_IEEE_F64LE);
  hdf5_write_single_val(&r_out, "r_out", H5T_IEEE_F64LE);
  hdf5_write_single_val(&rhor, "rhor", H5T_IEEE_F64LE);
  hdf5_write_single_val(&rin, "rin", H5T_IEEE_F64LE);
  hdf5_write_single_val(&rmax, "rmax", H5T_IEEE_F64LE);
  hdf5_write_single_val(&startx1, "startx1", H5T_IEEE_F64LE);
  hdf5_write_single_val(&startx2, "startx2", H5T_IEEE_F64LE);
  hdf5_write_single_val(&startx3, "startx3", H5T_IEEE_F64LE);
  hdf5_write_single_val(&time, "t", H5T_IEEE_F64LE);
  hdf5_write_single_val(&tlim, "tlim", H5T_IEEE_F64LE);
  hdf5_write_single_val(&x1min, "x1min", H5T_IEEE_F64LE);
  hdf5_write_single_val(&x1max, "x1max", H5T_IEEE_F64LE);
  hdf5_write_single_val(&x2min, "x2min", H5T_IEEE_F64LE);
  hdf5_write_single_val(&x2max, "x2max", H5T_IEEE_F64LE);
  hdf5_write_single_val(&x3min, "x3min", H5T_IEEE_F64LE);
  hdf5_write_single_val(&x3max, "x3max", H5T_IEEE_F64LE);

  // Get geom parameters in fluid header
  hdf5_make_directory("geom");
  hdf5_set_directory("/fluid_header/geom/");
  // Write geom parameters in fluid header
  // hdf5_write_single_val(&dx1, "dx1", H5T_IEEE_F64LE);
  // hdf5_write_single_val(&dx2, "dx2", H5T_IEEE_F64LE);
  // hdf5_write_single_val(&dx3, "dx3", H5T_IEEE_F64LE);
  hdf5_write_single_val(&startx1, "startx1", H5T_IEEE_F64LE);
  hdf5_write_single_val(&startx2, "startx2", H5T_IEEE_F64LE);
  hdf5_write_single_val(&startx3, "startx3", H5T_IEEE_F64LE);

  // Write FMKS parameters in geom
  if (strcmp(metric, "fmks") == 0) {
    hdf5_make_directory("fmks");
    hdf5_set_directory("/fluid_header/geom/fmks/");

    hdf5_write_single_val(&a, "a", H5T_IEEE_F64LE);
    hdf5_write_single_val(&hslope, "hslope", H5T_IEEE_F64LE);
    hdf5_write_single_val(&mks_smooth, "mks_smooth", H5T_IEEE_F64LE);
    hdf5_write_single_val(&poly_alpha, "poly_alpha", H5T_IEEE_F64LE);
    hdf5_write_single_val(&poly_xt, "poly_xt", H5T_IEEE_F64LE);
    hdf5_write_single_val(&r_eh, "r_eh", H5T_IEEE_F64LE);
    hdf5_write_single_val(&r_in, "r_in", H5T_IEEE_F64LE);
    hdf5_write_single_val(&r_out, "r_out", H5T_IEEE_F64LE);

    hdf5_set_directory("/fluid_header/geom/");
    hdf5_make_directory("mmks");
    hdf5_set_directory("/fluid_header/geom/mmks/");
    hdf5_write_single_val(&a, "a", H5T_IEEE_F64LE);
    hdf5_write_single_val(&hslope, "hslope", H5T_IEEE_F64LE);
    hdf5_write_single_val(&mks_smooth, "mks_smooth", H5T_IEEE_F64LE);
    hdf5_write_single_val(&poly_alpha, "poly_alpha", H5T_IEEE_F64LE);
    hdf5_write_single_val(&poly_xt, "poly_xt", H5T_IEEE_F64LE);
    hdf5_write_single_val(&r_eh, "r_eh", H5T_IEEE_F64LE);
    hdf5_write_single_val(&r_in, "r_in", H5T_IEEE_F64LE);
    hdf5_write_single_val(&r_out, "r_out", H5T_IEEE_F64LE);
  }

  hdf5_set_directory("/header/");
#if SLOW_LIGHT
  hdf5_write_single_val(&(data[1]->t), "t", H5T_IEEE_F64LE);
#else // FAST LIGHT
  hdf5_write_single_val(&(data[0]->t), "t", H5T_IEEE_F64LE);
#endif

  hdf5_write_single_val(&sigma_cut, "sigma_cut", H5T_IEEE_F64LE);
  hdf5_make_directory("electrons");
  hdf5_set_directory("/header/electrons/");
  if (ELECTRONS == 0) {
    hdf5_write_single_val(&tp_over_te, "tp_over_te", H5T_IEEE_F64LE);
  } else if (ELECTRONS == 2) {
    hdf5_write_single_val(&trat_small, "rlow", H5T_IEEE_F64LE);
    hdf5_write_single_val(&trat_large, "rhigh", H5T_IEEE_F64LE);
    hdf5_write_single_val(&beta_crit, "beta_crit", H5T_IEEE_F64LE);
  } else if (ELECTRONS == ELECTRONS_TFLUID) {
    hdf5_write_single_val(&mu_i, "mu_i", H5T_IEEE_F64LE);
    hdf5_write_single_val(&mu_e, "mu_e", H5T_IEEE_F64LE);
    hdf5_write_single_val(&mu_tot, "mu_tot", H5T_IEEE_F64LE);
  }
  hdf5_write_single_val(&ELECTRONS, "type", H5T_STD_I32LE);

  hdf5_set_directory("/header/");
  hdf5_write_single_val(&reverse_field,"field_config",H5T_STD_I32LE);
  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  hdf5_write_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&M_unit, "M_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&T_unit, "T_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Te_unit, "Thetae_unit", H5T_IEEE_F64LE);

  hdf5_set_directory("/");
}


int radiating_region(double X[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);
  return (r > rmin_geo && r < rmax_geo && th > th_beg && th < (M_PI-th_beg));
}


// In case we want to mess with emissivities directly
void get_model_jar(double X[NDIM], double Kcon[NDIM],
    double *jI, double *jQ, double *jU, double *jV,
    double *aI, double *aQ, double *aU, double *aV,
    double *rQ, double *rU, double *rV) {return;}


void get_model_jk(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv) {return;}
