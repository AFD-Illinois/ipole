#include "image.h"

#include "decs.h"
#include <ctype.h>

int compare_doubles(const void *a, const void *b)
{
  const double *da = (const double*) a;
  const double *db = (const double*) b;

  return (*da > *db) - (*da < *db);
}

#define TOP_FRAC	(0.005)	/* Brightest and dimmest pixels cut out of color scale */
#define BOT_FRAC	(0.005)

void make_ppm(double p[], int nx, int ny, double freq, char filename[])
{
  double q[nx*ny];

  int i, j, k;
  double min, max;
  FILE *fp;

  k = 0;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      q[k] = p[i * ny + j];
      k++;
    }

  int npixels = nx * ny;
  qsort(q, npixels, sizeof(double), compare_doubles);
  if (q[0] < 0) { /* must be log scaling */
    min = q[(int) (npixels * BOT_FRAC)];
  } else {
    i = 0;
    while (i < npixels - 1 && q[i] == 0.)
      i++;
    min = q[i];
  }
  max = q[npixels - (int) (npixels * TOP_FRAC)];

  if ((fp = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "Failed to open ppm file %s\n", filename);
    exit(125);
  }

  /* write out header information */
  fprintf(fp, "P6\n#  min=%g  , max=%g \n#  frequency=%g \n%d %d\n%d\n", min,
          max, freq, nx, ny, 255);
  fflush(fp);

  int red, green, blue;
  for (j = ny - 1; j >= 0; j--)
    for (i = 0; i < nx; i++) {
      rainbow_palette(p[i * ny + j], min, max, &red, &green, &blue);
      fputc((char) red, fp);
      fputc((char) green, fp);
      fputc((char) blue, fp);
    }

  fclose(fp);

  return;
}

#undef TOP_FRAC
#undef BOT_FRAC

/* 
 rainbow color palette, based on john.pal

 input: integer 0-255
 output: red, green, blue integers

 author: Bryan M. Johnson
 */

void rainbow_palette(double data, double min, double max, int *pRed,
                     int *pGreen, int *pBlue)
{
  double a, b, c, d, e, f;
  double x, y;
  double max_min = max - min;

  if (max_min > 0.0) { // trust no one
    x = (data - min) / (max_min);

    /* ========== Red ============ */
    a = 4.0 * x - 1.52549019607844;
    b = 4.52941176470589 - 4.0 * x;
    y = a < b ? a : b;
    *pRed = (int) (255.0 * y);
    *pRed = *pRed > 0 ? *pRed : 0;
    *pRed = *pRed < 255 ? *pRed : 255;

    /* ========== Green ========== */
    a = 4.0 * x - 0.521568627450979;
    b = 2.52549019607844 - 4.0 * x;
    c = a < b ? a : b;
    d = 4.0 * x - 1.53725490196073;
    e = 3.52941176470581 - 4.0 * x;
    f = d < e ? d : e;
    y = c > f ? c : f;
    *pGreen = (int) (255.0 * y);
    *pGreen = *pGreen > 0 ? *pGreen : 0;
    *pGreen = *pGreen < 255 ? *pGreen : 255;

    /* ========== Blue =========== */
    a = 4.0 * x + 0.498039215686276;
    b = 2.50980392156862 - 4.0 * x;
    y = a < b ? a : b;
    *pBlue = (int) (255.0 * y);
    *pBlue = *pBlue > 0 ? *pBlue : 0;
    *pBlue = *pBlue < 255 ? *pBlue : 255;

  } else {
    *pRed = *pGreen = *pBlue = (data > max ? 255 : 0);
  }

  return;
}
