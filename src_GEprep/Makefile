all : install

install : GEprep transform

clean : 
	rm *.o GEprep *.mod transform

OBJS = main.o error.o output_time.o inpinp.o store_element.o infil.o \
       distrib_CANO.o distrib_MULT.o distrib_CLMD.o distrib_ALSD.o \
       lstsquare.o minver.o lstsquare_w.o lstsquare_final.o

FC = ifort
ifeq ($(FC), ifort)
    LFLAGS = -O2 -ipo -no-prec-div -no-prec-sqrt -xHost -qopenmp
    #LFLAGS = -CB -traceback -pg -O2 -ipo -no-prec-div -xHost -openmp
else ifeq ($(FC), gfortran)
    LFLAGS = -O2 -march=native -fopenmp
else
    $(error Unsupported compiler: $(FC))
endif

CFLAGS = -c $(LFLAGS)

LIB = -lm

GEprep : COMDAT.f90 $(OBJS)
	$(FC) $(LFLAGS) -o GEprep COMDAT.f90 $(OBJS) $(LIB)

main.o : COMDAT.f90 main.f90
	$(FC) $(CFLAGS) COMDAT.f90 main.f90 $(LIB)

error.o : error.f90 COMDAT.f90
	$(FC) $(CFLAGS) COMDAT.f90 error.f90 $(LIB)

output_time.o : output_time.f90
	$(FC) $(CFLAGS) output_time.f90 $(LIB)

inpinp.o : inpinp.f90 COMDAT.f90
	$(FC) $(CFLAGS) COMDAT.f90 inpinp.f90 $(LIB)

store_element.o : store_element.f90 COMDAT.f90
	$(FC) $(CFLAGS) COMDAT.f90 store_element.f90 $(LIB)

infil.o : infil.f90
	$(FC) $(CFLAGS) infil.f90 $(LIB)

distrib_CANO.o : distrib_CANO.f90 COMDAT.f90
	$(FC) $(CFLAGS) COMDAT.f90 distrib_CANO.f90 $(LIB)

distrib_MULT.o : distrib_MULT.f90 COMDAT.f90
	$(FC) $(CFLAGS) COMDAT.f90 distrib_MULT.f90 $(LIB)

distrib_CLMD.o : distrib_CLMD.f90 COMDAT.f90
	$(FC) $(CFLAGS) COMDAT.f90 distrib_CLMD.f90 $(LIB)

distrib_ALSD.o : distrib_ALSD.f90 COMDAT.f90
	$(FC) $(CFLAGS) COMDAT.f90 distrib_ALSD.f90 $(LIB)

lstsquare.o : lstsquare.f90
	$(FC) $(CFLAGS) lstsquare.f90 $(LIB)

minver.o : minver.f90
	$(FC) $(CFLAGS) minver.f90 $(LIB)

lstsquare_w.o : lstsquare_w.f90
	$(FC) $(CFLAGS) lstsquare_w.f90 $(LIB)

lstsquare_final.o : lstsquare_final.f90
	$(FC) $(CFLAGS) lstsquare_final.f90 $(LIB)

transform.o : transform.f90
	$(FC) $(CFLAGS) transform.f90 $(LIB)

transform : transform.o lstsquare.o minver.o
	$(FC) $(LFLAGS) -o transform transform.o lstsquare.o minver.o $(LIB)
