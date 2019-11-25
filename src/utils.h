/*
 * utils.h
 *
 *  Created on: Sep 9, 2019
 *      Author: bprather
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

double ***malloc_rank3(int n1, int n2, int n3);
float ***malloc_rank3_float(int n1, int n2, int n3);

double ****malloc_rank4(int n1, int n2, int n3, int n4);
float ****malloc_rank4_float(int n1, int n2, int n3, int n4);

#endif /* SRC_UTILS_H_ */
