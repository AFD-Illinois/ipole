/*
 * utils.h
 *
 *  Created on: Sep 9, 2019
 *      Author: bprather
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

double ***malloc_rank3(u_int64_t n1, u_int64_t n2, u_int64_t n3);
float ***malloc_rank3_float(u_int64_t n1, u_int64_t n2, u_int64_t n3);

double ****malloc_rank4(u_int64_t n1, u_int64_t n2, u_int64_t n3, u_int64_t n4);
float ****malloc_rank4_float(u_int64_t n1, u_int64_t n2, u_int64_t n3, u_int64_t n4);

#endif /* SRC_UTILS_H_ */
