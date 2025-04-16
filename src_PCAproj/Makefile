all : install

install : PCAproj

clean : 
	rm *.o *.mod PCAproj

OBJS = main.o error.o output_time.o store_element.o project.o

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


PCAproj : COMDAT.f90 $(OBJS)
	$(FC) $(LFLAGS) -o PCAproj COMDAT.f90 $(OBJS) $(LIB)

main.o : COMDAT.f90 main.f90
	$(FC) $(CFLAGS) COMDAT.f90 main.f90 $(LIB)

error.o : error.f90
	$(FC) $(CFLAGS) error.f90 $(LIB)

output_time.o : output_time.f90
	$(FC) $(CFLAGS) output_time.f90 $(LIB)

store_element.o : store_element.f90 COMDAT.f90
	$(FC) $(CFLAGS) COMDAT.f90 store_element.f90 $(LIB)

project.o : project.f90 COMDAT.f90
	$(FC) $(CFLAGS) COMDAT.f90 project.f90 $(LIB)
