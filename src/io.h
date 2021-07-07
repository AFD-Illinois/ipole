/*
 * io.h
 *
 *  Created on: Sep 10, 2019
 *      Author: bprather
 */

#ifndef IO_H
#define IO_H

#include "decs.h"
#include "par.h"

void write_restart(const char *fname, double tA, double tB, double last_img_target,
                   int nopenimgs, int nimg, int nconcurrentimgs, int s2,
                   double *target_times, int *valid_imgs, struct of_image *dimages);
void read_restart(const char *fname, double *tA, double *tB, double *last_img_target,
                  int *nopenimgs, int *nimg, int nconcurrentimgs, int s2,
                  double *target_times, int *valid_imgs, struct of_image *dimages);

void dump_check(double image[], double imageS[], double taus[], Params *params, int skip);
void dump(double image[], double imageS[], double taus[],
          const char *fname, double scale, double cam[NDIM],
          double fovx, double fovy, size_t nx, size_t ny, Params *params);
void dump_var_along(int i, int j, int nstep, struct of_traj *traj, int nx, int ny,
                    double scale, double cam[NDIM], double fovx, double fovy, Params *params);

#endif // IO_H
