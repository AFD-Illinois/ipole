# h5cc (which uses icc) when hdf is used (i.e. for HARM3D data)
CC=/usr/bin/h5cc -DH5_USE_16_API
#CC = gcc
CFLAGS =  -fopenmp -I/usr/include -w
#CFLAGS =   -I/usr/include 
#-g -ggdb
LDFLAGS = -lm -lgsl -lgslcblas 
#
SRCMIB = \
geodesics.c geodesics_gsl.c geometry_utils.c \
harm3d_model.c harm3d_utils.c init_harm3d_data.c  \
image.c jnu_mixed.c\
main.c radiation.c tetrads_new.c brem.c ipolarray3.c lu.c
#harm_model.c harm_utils.c init_harm_data.c  \

OBJMIB = \
geodesics.o geodesics_gsl.o geometry_utils.o \
harm3d_model.o harm3d_utils.o init_harm3d_data.o  \
image.o jnu_mixed.o\
main.o radiation.o tetrads_new.o brem.o ipolarray3.o lu.o
#harm_model.o harm_utils.o init_harm_data.o  \


mibothros: $(OBJMIB) makefile 
	$(CC) $(CFLAGS) -o mibothros $(OBJMIB) $(LDFLAGS)
$(OBJMIB) : makefile decs.h defs.h constants.h
clean:
	rm *.o *.il
cleanup:
	rm mibothros*.ppm mibothros.dat


