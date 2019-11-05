/*
 * image.h
 *
 *  Created on: Sep 9, 2019
 *      Author: bprather
 */

#ifndef IMAGE_H
#define IMAGE_H

/* imaging */
void make_ppm(double p[], int nx, int ny, double freq, char filename[]);
void rainbow_palette(double data, double min, double max, int *pRed,
                     int *pGreen, int *pBlue);

#endif // IMAGE_H
