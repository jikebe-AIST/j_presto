all : install

install : PCAaxis

clean : 
	rm *.o PCAaxis

OBJS = main.o error.o output_time.o store_element.o mkaxis.o eigen2.o

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

PCAaxis : COMDAT.f90 $(OBJS)
	$(FC) $(LFLAGS) -o PCAaxis COMDAT.f90 $(OBJS) $(LIB)

main.o : COMDAT.f90 main.f90
	$(FC) $(CFLAGS) COMDAT.f90 main.f90 $(LIB)

error.o : error.f90
	$(FC) $(CFLAGS) error.f90 $(LIB)

output_time.o : output_time.f90
	$(FC) $(CFLAGS) output_time.f90 $(LIB)

store_element.o : store_element.f90 COMDAT.f90
	$(FC) $(CFLAGS) COMDAT.f90 store_element.f90 $(LIB)

mkaxis.o : mkaxis.f90 COMDAT.f90
	$(FC) $(CFLAGS) COMDAT.f90 mkaxis.f90 $(LIB)

eigen2.o : eigen2.f
	$(FC) $(CFLAGS) eigen2.f $(LIB)
