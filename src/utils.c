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


int **malloc_rank2_int(int n1, int n2)
{
  int **A;
  int *space;
  int i;

  space = malloc_rank1(n1*n2, sizeof(int));
  A = malloc_rank1(n1, sizeof(int *));
  for(i = 0; i < n1; i++){
    A[i] = &(space[n2 * i]);
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


double *****malloc_rank5(int n1, int n2, int n3, int n4, int n5)
{

  double *****A;
  double *space;
  int i,j,k,l;

  space = malloc_rank1(n1*n2*n3*n4*n5, sizeof(double));
  A = malloc_rank1(n1, sizeof(double *));
  for(i=0;i<n1;i++){
    A[i] = malloc_rank1(n2,sizeof(double *));
    for(j=0;j<n2;j++){
      A[i][j] = malloc_rank1(n3,sizeof(double *));
      for(k=0;k<n3;k++){
        A[i][j][k] = malloc_rank1(n4,sizeof(double *));
        for (l=0;l<n4;l++){
          A[i][j][k][l] = &(space[n5*(l + n4*(k + n3*(j + n2*i)))]);
        }
      }
    }
  }

  return A;
}

float *****malloc_rank5_float(int n1, int n2, int n3, int n4, int n5)
{

  float *****A;
  float *space;
  int i,j,k,l;

  space = malloc_rank1(n1*n2*n3*n4*n5, sizeof(float));
  A = malloc_rank1(n1, sizeof(float *));
  for(i=0;i<n1;i++){
    A[i] = malloc_rank1(n2,sizeof(float *));
    for(j=0;j<n2;j++){
      A[i][j] = malloc_rank1(n3,sizeof(float *));
      for(k=0;k<n3;k++){
        A[i][j][k] = malloc_rank1(n4,sizeof(float *));
        for (l=0;l<n4;l++){
          A[i][j][k][l] = &(space[n5*(l + n4*(k + n3*(j + n2*i)))]);
        }
      }
    }
  }

  return A;
}


void free_rank5(double *****A, int n1, int n2, int n3, int n4) {
  if (A == NULL)
      return;
  // First free the contiguous data block.
  // Note: A[0][0][0][0] points to the block allocated by malloc_rank1 in malloc_rank5.
  free(A[0][0][0][0]);

  // Now free the pointer arrays allocated at the 4th level:
  for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
          for (int k = 0; k < n3; k++) {
              free(A[i][j][k]);  // This frees the array allocated for each A[i][j][k]
          }
      }
  }
  // Free the pointer arrays at the 3rd level:
  for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
          free(A[i][j]);  // This frees each array allocated for A[i][j]
      }
  }
  // Free the pointer arrays at the 2nd level:
  for (int i = 0; i < n1; i++) {
      free(A[i]);  // This frees each array allocated for A[i]
  }
  // Finally free the top-level pointer array:
  free(A);
}