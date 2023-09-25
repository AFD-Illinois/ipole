CC=h5cc
GSL_DIR=/sw/spack/deltacpu-2022-03/apps/gsl/2.7-gcc-11.2.0-h5vpe53/
HDF5_DIR = /sw/spack/deltacpu-2022-03/apps/hdf5/1.13.1-gcc-11.2.0-dschtbv/

CFLAGS = -march=native -mtune=native -std=gnu11 -O3 -flto -fopenmp -funroll-loops