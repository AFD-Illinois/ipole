# Before building, make sure you have got the required modules,
# and have set the appropriate paths
#
# 1. Load gsl.
# `module load gsl`
#
# 2. Update LD_LIBRARY_PATH so that ipole knows where to look for shared libraries.
# `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/spack/deltas11-2023-03/apps/linux-rhel8-zen3/gcc-11.4.0/gsl-2.7.1-ytg74v2/lib`
#
# 3. Download and compile your own HDF5 (Delta's HDF5 is perpetually broken).
# `mkdir -p ~/software/hdf5; cd ~/software/hdf5`
# `wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.2/src/hdf5-1.12.2.tar.gz`
# `tar -xvzf hdf5-1.12.2.tar.gz`
# `cd hdf5-1.12.2`
# `module load gcc`
# `module load cmake`
# `mkdir build; cd build`
# `cmake -DCMAKE_INSTALL_PREFIX=~/software/hdf5/install -DBUILD_SHARED_LIBS=ON ..`
# `make -j$(nproc)`
# `make install`
# 
# 4. Set environment variables for custom HDF5
# `export HDF5_DIR=$HOME/software/hdf5/install`
# `export PATH=$HDF5_DIR/bin:$PATH`
# `export LD_LIBRARY_PATH=$HDF5_DIR/lib:$LD_LIBRARY_PATH`
# `export CPATH=$HDF5_DIR/include:$CPATH`

CC=h5cc
GSL_DIR=/sw/spack/deltas11-2023-03/apps/linux-rhel8-zen3/gcc-11.4.0/gsl-2.7.1-ytg74v2
HDF5_DIR=/sw/spack/deltas11-2023-03/apps/linux-rhel8-zen3/gcc-11.4.0/hdf5-1.14.2-h7iecxa

CFLAGS = -march=native -mtune=native -std=gnu11 -O3 -flto -fopenmp -funroll-loops
