
/*

model-dependent functions related to creation and manipulation of tetrads

*/

#include "decs.h"

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

#define SMALL_VECTOR (1.e-30)
int make_plasma_tetrad(double Ucon[NDIM], double Kcon[NDIM],
			double Bcon[NDIM], double Gcov[NDIM][NDIM],
			double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM])
{
    int k, l;
    void normalize(double *vcon, double Gcov[4][4]);
    void project_out(double *vcona, double *vconb, double Gcov[4][4]);
    double check_handedness(double Econ[NDIM][NDIM],
			    double Gcov[NDIM][NDIM]);

    void set_Econ_from_trial(double Econ[4], int defdir, double trial[4]);

#if 1 

    // alternative method inspired by bhlight
   
    // get character of Bcon and Kcon
    double Unorm = 0., Knorm = 0.;
    double UX1ness = 0., UX2ness = 0., UX3ness = 0.;
    double KX1ness = 0., KX2ness = 0., KX3ness = 0.;
    for (int mu=0; mu<NDIM; ++mu) {
        for (int nu=0; nu<NDIM; ++nu) {
            Knorm += Kcon[mu]*Kcon[nu]*Gcov[mu][nu];
            Unorm += Ucon[mu]*Ucon[nu]*Gcov[mu][nu];
        }
        UX1ness += Gcov[1][mu] * Ucon[mu];
        UX2ness += Gcov[2][mu] * Ucon[mu];
        UX3ness += Gcov[3][mu] * Ucon[mu];
        KX1ness += Gcov[1][mu] * Kcon[mu];
        KX2ness += Gcov[2][mu] * Kcon[mu];
        KX3ness += Gcov[3][mu] * Kcon[mu];
    }
    UX1ness = fabs(UX1ness);
    UX2ness = fabs(UX2ness);
    UX3ness = fabs(UX3ness);
    KX1ness = fabs(KX1ness);
    KX2ness = fabs(KX2ness);
    KX3ness = fabs(KX3ness);

    // set time component parallel to U and Econ[3] || Kcon
    set_Econ_from_trial(Econ[0], 0, Ucon);
    set_Econ_from_trial(Econ[2], 2, Bcon);
    set_Econ_from_trial(Econ[3], 3, Kcon);
    for (int mu=0; mu<NDIM; ++mu) Econ[1][mu] = 1.;

    normalize(Econ[0], Gcov);

    project_out(Econ[3], Econ[0], Gcov);
    normalize(Econ[3], Gcov);

    project_out(Econ[2], Econ[0], Gcov);
    project_out(Econ[2], Econ[3], Gcov);
    normalize(Econ[2], Gcov);

    project_out(Econ[1], Econ[0], Gcov);
    project_out(Econ[1], Econ[2], Gcov);
    project_out(Econ[1], Econ[3], Gcov);
    normalize(Econ[1], Gcov);

    // check handedness
    double dot = check_handedness(Econ, Gcov);

    /*
    int utetrad = 0;
    int ktetrad = 0;

    // if the dot condition fails, our alignment is 
    // probably degenerate. still attempt to align 
    // according to U and Kcon, then do best guess
    // following coordinates
    if ( fabs(fabs(dot) - 1.) > 1.e-10 ) {
    
        for (int mu=0; mu<NDIM; ++mu) {
            Econ[0][mu] = Ucon[mu];
            Econ[1][mu] = delta(1,mu);
            Econ[2][mu] = delta(2,mu);
            Econ[3][mu] = delta(3,mu);
        }

        if ( KX1ness > KX2ness && KX1ness > KX3ness ) {
            ktetrad = 1;
        } else if ( KX2ness > KX1ness && KX2ness > KX3ness ) {
            ktetrad = 2;
        } else { 
            ktetrad = 3;
        }

        for (int mu=0; mu<NDIM; ++mu) {
            Econ[ktetrad][mu] = Kcon[mu];
        }

        normalize(Econ[0], Gcov);

        project_out(Econ[3], Econ[0], Gcov);
        normalize(Econ[3], Gcov);

        project_out(Econ[2], Econ[0], Gcov);
        project_out(Econ[2], Econ[3], Gcov);
        normalize(Econ[2], Gcov);

        project_out(Econ[1], Econ[0], Gcov);
        project_out(Econ[1], Econ[2], Gcov);
        project_out(Econ[1], Econ[3], Gcov);
        normalize(Econ[1], Gcov);

        //fprintf(stderr, "%g %g %g %g %g %g\n", UX1ness, UX2ness, UX3ness, KX1ness, KX2ness, KX3ness);
        //fprintf(stderr, "is this odd? %g\n", fabs(dot) - 1.);

        //exit(-1);
    }
     */

    if (fabs(fabs(dot) - 1.) > 1.e-2) {
      fprintf(stderr, "\n");
      fprintf(stderr, "Ucon[] = %e %e %e %e\n", Ucon[0], Ucon[1], Ucon[2], Ucon[3]);
      fprintf(stderr, "Kcon[] = %e %e %e %e\n", Kcon[0], Kcon[1], Kcon[2], Kcon[3]);
      fprintf(stderr, "Bcon[] = %e %e %e %e\n", Bcon[0], Bcon[1], Bcon[2], Bcon[3]);
      fprintf(stderr, "nesses %g %g %g | %g %g %g\n", UX1ness, UX2ness, UX3ness, KX1ness, KX2ness, KX3ness);
      fprintf(stderr, "that's odd: %g\n", fabs(dot) - 1.);
     for (int i=0; i<4; ++i) {
          fprintf(stderr, "{%g, %g, %g, %g},\n", Econ[i][0], Econ[i][1], Econ[i][2], Econ[i][3]);
      }
      fprintf(stderr, "g = %g\n", gdet_func(Gcov));

      double udotk = 0.;
      for (int i=0; i<4; ++i) 
        for (int j=0; j<4; ++j) 
          udotk += Ucon[i] * Kcon[j] * Gcov[i][j];

      fprintf(stderr, "Knorm,Bnorm  %g %g\n\n\n", Knorm, Unorm);
      fprintf(stderr, "u.k = %g & debugflag %d\n", udotk, 0);
      exit(-1);
      return 1;
    }

    // dot = 1. for right-handed system. flip Econ[1] if necessary
    if ( dot < 0. ) {
      for (k = 0; k < 4; k++)
      Econ[1][k] *= -1.;
    }

    // now make covariant version by lowering coordinate basis index 
    for (k = 0; k < 4; k++) {
      lower(Econ[k], Gcov, Ecov[k]);
    }

    // and raise tetrad basis index
    for (l = 0; l < 4; l++) {
      Ecov[0][l] *= -1.;
    }

    return 0;

#else 

    /* start w/ time component parallel to U */
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
      printf("Ucon[] = %e %e %e %e\n", Ucon[0], Ucon[1],Ucon[2],Ucon[3]);
      printf("Kcon[] = %e %e %e %e\n", Kcon[0], Kcon[1],Kcon[2],Kcon[3]);
      printf("Bcon[] = %e %e %e %e\n", Bcon[0], Bcon[1], Bcon[2], Bcon[3]);
      fprintf(stderr, "that's odd: %g\n", fabs(dot) - 1.);
      exit(-1);
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

    return 0;

#endif // how to construct tetrads

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


    /* could use normal observer here; at present, camera has dx^i/dtau = 0 */
    Ucam[0] = 1.;
    Ucam[1] = 0.;
    Ucam[2] = 0.;
    Ucam[3] = 0.;

    /* this puts a photon with zero angular momentum in the center of
       the field of view */
    trial[0] = 1.;
    trial[1] = 1.;
    trial[2] = 0.;
    trial[3] = 0.;
    gcov_func(X, Gcov);
    double Gcon[NDIM][NDIM];
    gcon_func(Gcov, Gcon);
    lower(trial, Gcon, Kcon);

    trial[0] = 0.;
    trial[1] = 0.;
    trial[2] = 1.;
    trial[3] = 0.;

    make_plasma_tetrad(Ucam, Kcon, trial, Gcov, Econ, Ecov);

    /* done! */
}

