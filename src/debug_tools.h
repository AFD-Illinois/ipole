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

void dump_at_X(double X[NDIM]);


void check_ortho(double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM]);
void check_u(double Ucon[NDIM], double Ucov[NDIM]);

/*
 Check that the coherency tensor N satisfies certain
 basic properties:
 k . N = N . k = 0
 hermitian
 evaluate the invariants: I, Q^2 + U^2, V^2
 */
void check_N(double complex N[NDIM][NDIM], double Kcon[NDIM],
    double gcov[NDIM][NDIM]);

#endif /* SRC_DEBUG_TOOLS_H_ */
