
#include "decs.h"

/* 

   set up the metric (indicies both up and down), and the
   connection coefficients on a grid, based
   on a HARM simulation grid.  
   
   In principle the geometry could be sampled on a finer or
   coarser grid than the simulation grid. 

   These routines are taken directly out of HARM.

   They require the code be compiled against the 
   Gnu scientific library (GSL).

   CFG 21 July 06
   
*/

gsl_matrix *gsl_gcov, *gsl_gcon;
gsl_permutation *perm;

/* assumes gcov has been set first; returns determinant */
double gdet_func(double gcov[][NDIM])
{
  /*
	double d;
	int k, l, signum;

	if (gsl_gcov == NULL) {
		gsl_gcov = gsl_matrix_alloc(NDIM, NDIM);
		gsl_gcon = gsl_matrix_alloc(NDIM, NDIM);
		perm = gsl_permutation_alloc(NDIM);
	}

	DLOOP gsl_matrix_set(gsl_gcov, k, l, gcov[k][l]);

	gsl_linalg_LU_decomp(gsl_gcov, perm, &signum);

	d = gsl_linalg_LU_det(gsl_gcov, signum);

	return (sqrt(fabs(d)));
  */

  //instead of gsl stuff we use harm3d stuff
	
	int i;
  int permute[NDIM];
  double gcovtmp[NDIM][NDIM];
  double detg;
  int LU_decompose( double A[][NDIM], int permute[] );

  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  gcovtmp[0][i] = gcov[0][i]; }
  if( LU_decompose( gcovtmp,  permute ) != 0  ) {
    fprintf(stderr, "gdet_func(): singular matrix encountered! \n");
    exit(1);
  }
  detg = 1.;
  for(i=0  ;i<NDIM ;i++) detg *= gcovtmp[i][i];
  return( sqrt(fabs(detg)) );
  

}


/* invert gcov to get gcon */

void gcon_func(double gcov[][NDIM], double gcon[][NDIM])
{
  /*
	int k, l, signum;
	
	if (gsl_gcov  == NULL) {
		gsl_gcov = gsl_matrix_alloc(NDIM, NDIM);
		gsl_gcon = gsl_matrix_alloc(NDIM, NDIM);
		perm = gsl_permutation_alloc(NDIM);
	}
	

	DLOOP gsl_matrix_set(gsl_gcov, k, l, gcov[k][l]);

	gsl_linalg_LU_decomp(gsl_gcov, perm, &signum);

	gsl_linalg_LU_invert(gsl_gcov, perm, gsl_gcon);

	DLOOP gcon[k][l] = gsl_matrix_get(gsl_gcon, k, l);
  */

  //instead of gsl stuff we use harm3d stuff
	
  int invert_matrix( double Am[][NDIM], double Aminv[][NDIM] )  ;
  
  invert_matrix( gcov, gcon );
  

  /* done! */
}


