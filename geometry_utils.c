
#include "decs.h"

/* 

   find determinant of covariant metric gcov

   find gcon, as numerical inverse of gcov

   These routines are taken directly out of HARM.

   CFG 21 July 06
   MM 11 July 17
   
*/

/* assumes gcov has been set first; returns sqrt{|g|} */
double gdet_func(double gcov[][NDIM])
{

    int i,j;
    int permute[NDIM];
    double gcovtmp[NDIM][NDIM];
    double gdet;
    int LU_decompose(double A[][NDIM], int permute[]);

    for (i = 0; i < NDIM; i++) 
    for (j = 0; j < NDIM; j++) {
	gcovtmp[i][j] = gcov[i][j];
    }
    if (LU_decompose(gcovtmp, permute) != 0) {
	fprintf(stderr, "gdet_func(): singular matrix encountered! \n");
	exit(1);
    }
    gdet = 1.;

    for (i = 0; i < NDIM; i++)
	gdet *= gcovtmp[i][i];

    return (sqrt(fabs(gdet)));
}


/* invert gcov to get gcon */
void gcon_func(double gcov[][NDIM], double gcon[][NDIM])
{
    int invert_matrix(double Am[][NDIM], double Aminv[][NDIM]);

    invert_matrix(gcov, gcon);

    /* done! */
}
