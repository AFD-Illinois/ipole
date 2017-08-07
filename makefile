#
# h5cc compiles for linking with HDF5 library
#
CC = h5cc -DH5_USE_16_API 
CFLAGS =  -fopenmp -I/usr/include -Wall
LDFLAGS = -lm -lgsl -lgslcblas 

SRCIPO = \
geodesics.c geodesics_gsl.c \
image.c \
main.c radiation.c tetrads.c ipolarray.c \
model_tetrads.c model_radiation.c model_geometry.c model_geodesics.c \
model_harm3d.c \
geometry.c

OBJIPO = \
geodesics.o geodesics_gsl.o \
image.o \
main.o radiation.o tetrads.o ipolarray.o \
model_tetrads.o model_radiation.o model_geometry.o model_geodesics.o \
model_harm3d.o \
geometry.o


ipole: $(OBJIPO) makefile 
	$(CC) $(CFLAGS) -o ipole $(OBJIPO) $(LDFLAGS)

$(OBJIPO) : makefile decs.h defs.h constants.h

clean:
	rm *.o 
cleanup:
	rm ipole*.ppm ipole.dat



