/*
 * utils.c
 *
 *  Created on: Sep 9, 2019
 *      Author: bprather
 */

#include <stdlib.h>
#include <stdio.h>

// malloc utilities
void *malloc_rank1(size_t n1, size_t size)
{
  void *A;

  if ((A = calloc(n1, size)) == NULL) {
    fprintf(stderr,"malloc failure in malloc_rank1\n");
    exit(123);
  }

  return A;
}


double ***malloc_rank3(int n1, int n2, int n3)
{

  double ***A;
  double *space;
  int i,j;

  space = malloc_rank1(n1*n2*n3, sizeof(double));
  A = malloc_rank1(n1, sizeof(double *));
  for(i = 0; i < n1; i++){
    A[i] = malloc_rank1(n2,sizeof(double *));
    for(j = 0; j < n2; j++){
      A[i][j] = &(space[n3*(j + n2*i)]);
    }
  }

  return A;
}

float ***malloc_rank3_float(int n1, int n2, int n3)
{

  float ***A;
  float *space;
  int i,j;

  space = malloc_rank1(n1*n2*n3, sizeof(float));
  A = malloc_rank1(n1, sizeof(float *));
  for(i = 0; i < n1; i++){
    A[i] = malloc_rank1(n2,sizeof(float *));
    for(j = 0; j < n2; j++){
      A[i][j] = &(space[n3*(j + n2*i)]);
    }
  }

  return A;
}

double ****malloc_rank4(int n1, int n2, int n3, int n4)
{

  double ****A;
  double *space;
  int i,j,k;

  space = malloc_rank1(n1*n2*n3*n4, sizeof(double));
  A = malloc_rank1(n1, sizeof(double *));
  for(i=0;i<n1;i++){
    A[i] = malloc_rank1(n2,sizeof(double *));
    for(j=0;j<n2;j++){
      A[i][j] = malloc_rank1(n3,sizeof(double *));
      for(k=0;k<n3;k++){
        A[i][j][k] = &(space[n4*(k + n3*(j + n2*i))]);
      }
    }
  }

  return A;
}

float ****malloc_rank4_float(int n1, int n2, int n3, int n4)
{

  float ****A;
  float *space;
  int i,j,k;

  space = malloc_rank1(n1*n2*n3*n4, sizeof(float));
  A = malloc_rank1(n1, sizeof(float *));
  for(i=0;i<n1;i++){
    A[i] = malloc_rank1(n2,sizeof(float *));
    for(j=0;j<n2;j++){
      A[i][j] = malloc_rank1(n3,sizeof(float *));
      for(k=0;k<n3;k++){
        A[i][j][k] = &(space[n4*(k + n3*(j + n2*i))]);
      }
    }
  }

  return A;
}
