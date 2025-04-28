# For oneapi. Load with ". /opt/intel/oneapi/setvars.sh"
#CC=/opt/intel/oneapi/compiler/2021.1.2/linux/bin/intel64/icc
CFLAGS=-xHost -Ofast -fstrict-aliasing -Wall -Werror -ipo -qopenmp
#HDF5_DIR=/usr/lib/hpc/gnu7/hdf5/1.10.5
