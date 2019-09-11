/*
 * io.h
 *
 *  Created on: Sep 10, 2019
 *      Author: bprather
 */

#ifndef IO_H
#define IO_H

#include "decs.h"

void write_restart(const char *fname, double tA, double tB, double last_img_target,
                   int nopenimgs, int nimg, int nconcurrentimgs, int s2,
                   double *target_times, int *valid_imgs, struct of_image *dimages);
void read_restart(const char *fname, double *tA, double *tB, double *last_img_target,
                  int *nopenimgs, int *nimg, int nconcurrentimgs, int s2,
                  double *target_times, int *valid_imgs, struct of_image *dimages);

void dump (double image[], double imageS[], double taus[], const char *fname,
      double scale, double Dsource, double cam[NDIM], double DX, double DY,
      double fovx, double fovy, double rcam, double thetacam, double phicam,
      double rotcam, double xoff, double yoff);

#endif // IO_H
