#ifndef SIMCOORDS_H
#define SIMCOORDS_H

#include "decs.h"
#include <stdio.h>

//                                 //
//  !!!!!!!!!! WARNING !!!!!!!!!!  //
//                                 //
//    We assume phi = x3, so x3    //
//    must span 0 -> 2 Pi !        //
//                                 //

// supported coordinate systems:
//   KORAL_MKS3 (forward and reverse available, implements reverse)
//   KORAL_JETCOORDS (only forward available)

#define SIMCOORDS_KORAL_MKS3 3
#define SIMCOORDS_KORAL_JETCOORDS 4

extern int use_simcoords;
extern int simcoords;

// interface
void load_simcoord_info_from_file(const char *fname);
void initialize_simgrid(size_t n1, size_t n2, double x1i, double x1f, double x2i, double x2f);
int simcoord_to_eks(double gridcoord[NDIM], double eKS[NDIM]);
int simcoordijk_to_eks(int i, int j, int k, double eKS[NDIM]);
double simcoordijk_to_gdet(int i, int j, int k);
void eks_to_simcoord(double eKS[NDIM], double gridcoord[NDIM]);
void finalize_simgrid();

// not particularly neat, but the idea is to "namespace" different
// metric parameters in structs, and then define one global struct
// instance for each coordinate system.
typedef struct metric_params_koral_mks3_struct {
  double r0;
  double h0;
  double my1;
  double my2;
  double mp0;
} metric_params_koral_mks3;
extern metric_params_koral_mks3 mp_koral_mks3;

typedef struct metric_params_koral_jetcoords_struct {
  double mksr0;
  double rbrk;
  double rmin;
  double rphotomax;
  double fjet;
  double fdisk;
  double runi;
  double rcoll_jet;
  double rcoll_disk;
  double rdecoll_jet;
  double rdecoll_disk;
  double alpha_1;
  double alpha_2;
  double cylindrify;
  double rcyl;
  double ncyl;
} metric_params_koral_jetcoords;
extern metric_params_koral_jetcoords mp_koral_jetcoords;

#endif // SIMCOORDS_H
