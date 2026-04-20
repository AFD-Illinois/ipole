_delta_env := $(shell bash machines/delta.sh >/dev/null 2>&1)
CC=h5pcc
GSL_DIR=$(GSL_HOME)
#assumes you have loaded cray-hdf5 module which automatically sets HDF5_DIR
HDF5_DIR=$(HDF5_ROOT)

CFLAGS = -march=native -mtune=native -std=gnu11 -O3 -flto -fopenmp -funroll-loops
# $(info Using Delta settings: remember to prepend GSL lib location $(GSL_HOME)/lib to LD_LIBRARY_PATH for runtime)