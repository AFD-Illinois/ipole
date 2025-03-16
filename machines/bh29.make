# For AOCC. Load with ". /opt/AMD/aocc-compiler-3.0.0/setenv_AOCC.sh"
CC=/opt/AMD/aocc-compiler-3.0.0/bin/clang
CFLAGS=-std=gnu99 -O3 -march=native -mtune=native -flto -fopenmp -funroll-loops
