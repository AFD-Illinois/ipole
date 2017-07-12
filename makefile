#
# h5cc compiles for linking with HDF5 library
#
CC = h5cc -DH5_USE_16_API 
CFLAGS =  -fopenmp -I/usr/include -Wall
LDFLAGS = -lm -lgsl -lgslcblas 

SRCIPO = \
geodesics.c geodesics_gsl.c geometry_utils.c \
harm3d_model.c harm3d_utils.c init_harm3d_data.c  \
image.c jnu_mixed.c\
main.c radiation.c tetrads.c ipolarray.c lu.c

OBJIPO = \
geodesics.o geodesics_gsl.o geometry_utils.o \
harm3d_model.o harm3d_utils.o init_harm3d_data.o  \
image.o jnu_mixed.o\
main.o radiation.o tetrads.o ipolarray.o lu.o


ipole: $(OBJIPO) makefile 
	$(CC) $(CFLAGS) -o ipole $(OBJIPO) $(LDFLAGS)

$(OBJIPO) : makefile decs.h defs.h constants.h

clean:
	rm *.o 
cleanup:
	rm ipole*.ppm ipole.dat



