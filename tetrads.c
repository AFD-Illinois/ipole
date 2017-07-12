
/*

all functions related to creation and manipulation of tetrads

*/

#include "decs.h"

/* input and vectors are contravariant (index up) */
void coordinate_to_tetrad(double Ecov[NDIM][NDIM], double K[NDIM],
			  double K_tetrad[NDIM])
{
    int k;

    for (k = 0; k < 4; k++) {
	K_tetrad[k] = Ecov[k][0] * K[0] 
	            + Ecov[k][1] * K[1] 
	            + Ecov[k][2] * K[2] 
	            + Ecov[k][3] * K[3];
    }
}

/* input and vectors are contravariant (index up) */
void tetrad_to_coordinate(double Econ[NDIM][NDIM], double K_tetrad[NDIM],
			  double K[NDIM])
{
    int l;

    for (l = 0; l < 4; l++) {
	K[l] = Econ[0][l] * K_tetrad[0] 
	     + Econ[1][l] * K_tetrad[1] 
	     + Econ[2][l] * K_tetrad[2] 
	     + Econ[3][l] * K_tetrad[3];
    }

    return;
}

#define SMALL_VECTOR	1.e-30

/** tetrad making routines **/

/* 

   econ/ecov index key:

   Econ[k][l]
   k: index attached to tetrad basis
   index down
   l: index attached to coordinate basis 
   index up

   Ecov[k][l]
   k: index attached to tetrad basis
   index up
   l: index attached to coordinate basis 
   index down

*/

/* 

   make orthonormal basis for plasma frame.

   e^0 along U
   e^2 along b
   e^3 along spatial part of K

*/

void make_plasma_tetrad(double Ucon[NDIM], double Kcon[NDIM],
			double Bcon[NDIM], double Gcov[NDIM][NDIM],
			double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM])
{
    int k, l;
    void normalize(double *vcon, double Gcov[4][4]);
    void project_out(double *vcona, double *vconb, double Gcov[4][4]);
    double check_handedness(double Econ[NDIM][NDIM],
			    double Gcov[NDIM][NDIM]);

    /* start w/ time component parallel to U */
    void set_Econ_from_trial(double Econ[4], int defdir, double trial[4]);
    set_Econ_from_trial(Econ[0], 0, Ucon);
    normalize(Econ[0], Gcov);

    /*** done w/ basis vector 0 ***/

    /* now use the trial vector in basis vector 3 */
    /* cast a suspicious eye on the trial vector... */
    set_Econ_from_trial(Econ[3], 3, Kcon);

    /* project out econ0 */
    project_out(Econ[3], Econ[0], Gcov);
    normalize(Econ[3], Gcov);

    /*** done w/ basis vector 3 ***/

    /* repeat for x2 unit basis vector */
    set_Econ_from_trial(Econ[2], 2, Bcon);

    /* project out econ0,3 */
    project_out(Econ[2], Econ[0], Gcov);
    project_out(Econ[2], Econ[3], Gcov);
    normalize(Econ[2], Gcov);

    /*** done w/ basis vector 2 ***/

    /* whatever is left is econ1 */
    for (k = 0; k < 4; k++)	/* trial vector */
        Econ[1][k] = 1.;
    /* project out econ[0-2] */
    project_out(Econ[1], Econ[0], Gcov);
    project_out(Econ[1], Econ[2], Gcov);
    project_out(Econ[1], Econ[3], Gcov);
    normalize(Econ[1], Gcov);

    /* check handedness */
    double dot = check_handedness(Econ, Gcov);

    if (fabs(fabs(dot) - 1.) > 1.e-10) {
	fprintf(stderr, "that's odd: %g\n", fabs(dot) - 1.);
    }

    /* we expect dot = 1. for right-handed system.  
       If not, flip Econ[1] to make system right-handed. */
    if (dot < 0.) {
	for (k = 0; k < 4; k++)
	    Econ[1][k] *= -1.;
    }

    /*** done w/ basis vector 3 ***/

    /* now make covariant version */
    for (k = 0; k < 4; k++) {

	/* lower coordinate basis index */
	lower(Econ[k], Gcov, Ecov[k]);
    }

    /* then raise tetrad basis index */
    for (l = 0; l < 4; l++) {
	Ecov[0][l] *= -1.;
    }

    /* paranoia: check orthonormality */
    /*
       double sum ;
       int m ;
       fprintf(stderr,"ortho check [plasma]:\n") ;
       for(k=0;k<NDIM;k++)
       for(l=0;l<NDIM;l++) {
       sum = 0. ;
       for(m=0;m<NDIM;m++) {
       sum += Econ[k][m]*Ecov[l][m] ;
       }
       fprintf(stderr,"sum: %d %d %g\n",k,l,sum) ;
       }
       fprintf(stderr,"\n") ;
       for(k=0;k<NDIM;k++)
       for(l=0;l<NDIM;l++) {
       fprintf(stderr,"%d %d %g\n",k,l,Econ[k][l]) ;
       }
       fprintf(stderr,"done ortho check.\n") ;
       fprintf(stderr,"\n") ;
     */

    /* done! */

}

/*
	make orthonormal basis for camera frame.

	e^0 along Ucam
	e^1 outward (!) along radius vector
	e^2 toward north pole of coordinate system
		("y" for the image plane)
	e^3 in the remaining direction
		("x" for the image plane)

        this needs to be arranged this way so that
        final Stokes parameters are measured correctly.

        essentially an interface to make_plasma_tetrad

*/

void make_camera_tetrad(double X[NDIM], double Econ[NDIM][NDIM],
			double Ecov[NDIM][NDIM])
{
    double Gcov[NDIM][NDIM];
    double Ucam[NDIM];
    double Kcon[NDIM];
    double trial[NDIM];

    gcov_func(X, Gcov);

    /* might be better to use normal observer here;
       at present, camera has dx^i/dtau = 0 */
    Ucam[0] = 1.;
    Ucam[1] = 0.;
    Ucam[2] = 0.;
    Ucam[3] = 0.;

    Kcon[0] = 0.;
    Kcon[1] = 1.;
    Kcon[2] = 0.;
    Kcon[3] = 0.;

    trial[0] = 0.;
    trial[1] = 0.;
    trial[2] = 1.;
    trial[3] = 0.;

    make_plasma_tetrad(Ucam, Kcon, trial, Gcov, Econ, Ecov);

    /* done! */
}

/*
    Kronecker delta
*/
double delta(int i, int j)
{
    if (i == j)
	return (1.);
    else
	return (0.);
}

void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov)
{

    ucov[0] = Gcov[0][0] * ucon[0]
	+ Gcov[0][1] * ucon[1]
	+ Gcov[0][2] * ucon[2]
	+ Gcov[0][3] * ucon[3];
    ucov[1] = Gcov[1][0] * ucon[0]
	+ Gcov[1][1] * ucon[1]
	+ Gcov[1][2] * ucon[2]
	+ Gcov[1][3] * ucon[3];
    ucov[2] = Gcov[2][0] * ucon[0]
	+ Gcov[2][1] * ucon[1]
	+ Gcov[2][2] * ucon[2]
	+ Gcov[2][3] * ucon[3];
    ucov[3] = Gcov[3][0] * ucon[0]
	+ Gcov[3][1] * ucon[1]
	+ Gcov[3][2] * ucon[2]
	+ Gcov[3][3] * ucon[3];

    return;
}

/* 

    normalize input vector (and overwrite)
    so that |v . v| = 1

*/

void normalize(double *vcon, double Gcov[4][4])
{
    int k, l;
    double norm;

    norm = 0.;
    for (k = 0; k < 4; k++)
	for (l = 0; l < 4; l++)
	    norm += vcon[k] * vcon[l] * Gcov[k][l];

    norm = sqrt(fabs(norm));
    for (k = 0; k < 4; k++)
	vcon[k] /= norm;

    return;
}

/*

    project out vconb from vcona

    both arguments are index up (contravariant)

    covariant metric is third argument.

    overwrite the first argument on return

*/

void project_out(double *vcona, double *vconb, double Gcov[4][4])
{

    double adotb, vconb_sq;
    int k, l;

    vconb_sq = 0.;
    for (k = 0; k < 4; k++)
	for (l = 0; l < 4; l++)
	    vconb_sq += vconb[k] * vconb[l] * Gcov[k][l];

    adotb = 0.;
    for (k = 0; k < 4; k++)
	for (l = 0; l < 4; l++)
	    adotb += vcona[k] * vconb[l] * Gcov[k][l];

    for (k = 0; k < 4; k++)
	vcona[k] -= vconb[k] * adotb / vconb_sq;

    return;
}

/* 

   copy the trial vector into a tetrad basis vector,
   checking to see if it is null, and if it is null
   setting to some default value 

*/
void set_Econ_from_trial(double Econ[4], int defdir, double trial[4])
{
    double norm = 0.;
    int k;

    for (k = 0; k < 4; k++)
	norm += fabs(trial[k]);
    for (k = 0; k < 4; k++)	/* trial vector */
	if (norm <= SMALL_VECTOR)	/* bad trial vector; default to radial direction */
	    Econ[k] = delta(k, defdir);
	else
	    Econ[k] = trial[k];

    return;
}

/* 
    check the handedness of a tetrad basis.

    basis is assumed to be in form e^\mu_{(a)} = Econ[a][mu]

    levi_(ijkl) e0^i e1^j e2^k e3^l will be +1 if spatial
    	components are right-handed, -1 if left-handed.

    experience suggests that roundoff produces errors of
    	order 10^{-12} in the result.

*/

double check_handedness(double Econ[NDIM][NDIM], double Gcov[NDIM][NDIM])
{
    int i, j, k, l;
    static int firstc = 1;
    void set_levi_civita(double levi_civita[NDIM][NDIM][NDIM][NDIM]);
    static double levi_civita[NDIM][NDIM][NDIM][NDIM];

    if (firstc) {
	firstc = 0;
	set_levi_civita(levi_civita);
    }

    double g = gdet_func(Gcov);

    /* check handedness */
    double dot = 0.;
    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    for (l = 0; l < 4; l++)
		for (k = 0; k < 4; k++) {
		    dot += g * levi_civita[i][j][k][l] *
			Econ[0][i] * Econ[1][j] * Econ[2][k] * Econ[3][l];
		}

    return (dot);
}

/* the completely antisymmetric symbol; not a tensor
   in the coordinate basis */
void set_levi_civita(double levi_civita[NDIM][NDIM][NDIM][NDIM])
{
    int i, j, k, l, n, do_sort, n_perm, val, n_swap;
    int index[NDIM];

    for (i = 0; i < NDIM; i++)
	for (j = 0; j < NDIM; j++)
	    for (k = 0; k < NDIM; k++)
		for (l = 0; l < NDIM; l++) {
		    if (i == j || i == k || i == l || j == k || j == l
			|| k == l) {
			levi_civita[i][j][k][l] = 0;
		    } else {
			index[0] = i;
			index[1] = j;
			index[2] = k;
			index[3] = l;
			do_sort = 1;
			n_perm = 0;
			while (do_sort) {
			    n_swap = 0;
			    for (n = 0; n < NDIM - 1; n++) {
				if (index[n] > index[n + 1]) {
				    n_perm++;
				    n_swap++;
				    val = index[n];
				    index[n] = index[n + 1];
				    index[n + 1] = val;
				}
			    }
			    do_sort = n_swap;
			}
			levi_civita[i][j][k][l] = (n_perm % 2) ? -1 : 1;
		    }
		}

    /* Test levi-civita : */
    /*                              
       for(i=0;i<NDIM;i++)  
       for(j=0;j<NDIM;j++)  
       for(k=0;k<NDIM;k++)  
       for(l=0;l<NDIM;l++) {                         
       fprintf(stdout,"levi-civita[%d%d%d%d] = %d \n", 
       i,j,k,l,levi_civita[i][j][k][l] );                                                
       fflush(stdout); 
       }       
     */
}
