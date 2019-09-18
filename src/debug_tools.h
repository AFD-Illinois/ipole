/*
 * debug_tools.h
 *
 *  Created on: Sep 16, 2019
 *      Author: bprather
 */

#ifndef SRC_DEBUG_TOOLS_H_
#define SRC_DEBUG_TOOLS_H_

#include "decs.h"

void print_matrix(char *name, double g[NDIM][NDIM]);
void print_matrix_c(char *name, double complex g[NDIM][NDIM]);
void print_vector(char *name, double v[NDIM]);
void check_ortho(double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM]);
void check_u(double Ucon[NDIM], double Ucov[NDIM]);

#endif /* SRC_DEBUG_TOOLS_H_ */
