GSL_DIR = /opt/sns/gsl/gcc-4.8.5/2.5/

CFLAGS_CUSTOM += -I$(GSL_DIR)/include
SYSTEM_LIBDIR += -L$(GSL_DIR)/lib
MATH_LIB += -lgsl -lgslcblas -lm