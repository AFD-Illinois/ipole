/*
 * utils.h
 *
 *  Created on: Sep 9, 2019
 *      Author: bprather
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

int **malloc_rank2_int(int n1, int n2);

double ***malloc_rank3(int n1, int n2, int n3);
float ***malloc_rank3_float(int n1, int n2, int n3);

double ****malloc_rank4(int n1, int n2, int n3, int n4);
float ****malloc_rank4_float(int n1, int n2, int n3, int n4);

double *****malloc_rank5(int n1, int n2, int n3, int n4, int n5);
float *****malloc_rank5_float(int n1, int n2, int n3, int n4, int n5);

void free_rank5(double *****A, int n1, int n2, int n3, int n4);
#endif /* SRC_UTILS_H_ */
