#!/bin/bash

# Load gsl module
module load cray-hdf5-parallel gsl

# Update the LD_LIBRARY_PATH to include the GSL lib directory
export LD_LIBRARY_PATH=$GSL_HOME/lib:$LD_LIBRARY_PATH

# Set the HDF5 directory
# export HDF5_DIR=/opt/cray/pe/hdf5-parallel/1.14.3.5/gnu/12.2

# Update the PATH to include the HDF5 bin directory
export PATH=$HDF5_DIR/bin:$PATH

# Update the LD_LIBRARY_PATH to include the HDF5 lib directory
export LD_LIBRARY_PATH=$HDF5_DIR/lib:$LD_LIBRARY_PATH

# Update the CPATH to include the HDF5 include directory
export CPATH=$HDF5_DIR/include:$CPATH

# Confirm the environment variables and module have been set
echo "GSL module loaded"
echo "HDF5_DIR is set to $HDF5_DIR"
echo "PATH is updated to include $HDF5_DIR/bin"
echo "LD_LIBRARY_PATH is updated to include $HDF5_DIR/lib and the GSL library path"
echo "CPATH is updated to include $HDF5_DIR/include"